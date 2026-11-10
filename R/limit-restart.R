# limit-restart.R -- a phase parameter running to its limit (#683, #729)
#
# Design: inst/dev/LIMIT-RESTART-DESIGN.md. A multiphase fit can walk a
# shape parameter toward a limit of the BNT parameterisation and return the
# last iterate with converged = TRUE. hazard() now restarts from a fixed grid
# and returns the best interior fit, or refuses. It never switches the model
# and never returns the drifted iterate. SAS/C has no such check: it chooses
# a limiting form only from values the analyst fixed (setg1.c:620-636).

# Boundary mechanisms that ARE a shape at a limit, with the parameter and the
# limit each one names. Anything else on $boundary (setup holds, the g3
# corner, CoE) is handled where it is raised and does not trigger.
.hzr_limit_mechanisms <- list(
  phase_beyond_data   = c("t_half", "Inf (the rise lies past the data)"),
  unbounded_phase     = c("t_half", "0 (below the first observed time)"),
  phase_discontinuity = c("nu", "0 (a step at t_half)")
)

# Internal shape name -> the hzr_phase() argument a user would fix.
.hzr_limit_user_name <- c(log_t_half = "t_half", nu = "nu", m = "m",
                          log_tau = "tau", gamma = "gamma", alpha = "alpha",
                          eta = "eta")

#' Every shape parameter's full name, mapped to its phase
#'
#' Built from the phases' own naming, never by suffix: a covariate column
#' called `nu` is `early.nu` too, and would be misread (classify by the
#' fit's structure, not by a string's shape).
#' @param phases Named list of `hzr_phase` objects.
#' @return Named character vector: names are `<phase>.<shape>`, values the
#'   phase name.
#' @keywords internal
#' @noRd
.hzr_limit_shape_names <- function(phases) {
  out <- character(0)
  for (nm in names(phases)) {
    shapes <- .hzr_phase_theta_names(phases[[nm]], nm)[-1L]  # drop log_mu
    out <- c(out, stats::setNames(rep(nm, length(shapes)), shapes))
  }
  out
}

#' Does a fitted multiphase object show a shape parameter at a limit?
#'
#' Reads only what the object carries: `$fit$boundary` records whose
#' mechanism is a limit, and a `$fit$weak` direction that names a shape.
#' A parameter the user fixed never triggers: holding it there was their
#' choice, as FIXNU is in SAS.
#'
#' @param object A `hazard` object, or a list with `$spec$phases` and `$fit`.
#' @return `NULL`, or a list of `list(phase, parameter, limit, source)`.
#' @keywords internal
#' @noRd
.hzr_limit_trigger <- function(object) {
  phases <- object$spec$phases
  fs <- object$fit
  if (!length(phases)) return(NULL)
  shape_phase <- .hzr_limit_shape_names(phases)
  theta_names <- names(fs$theta)
  fixed <- if (length(fs$fixed_mask) == length(theta_names)) {
    theta_names[as.logical(fs$fixed_mask)]
  } else {
    character(0)
  }
  internal_of <- function(phase, user) {
    key <- names(.hzr_limit_user_name)[.hzr_limit_user_name == user]
    paste0(phase, ".", key)
  }
  hits <- list()
  add <- function(phase, parameter, limit, source) {
    if (internal_of(phase, parameter) %in% fixed) return(invisible())
    id <- paste0(phase, ".", parameter)
    if (id %in% names(hits)) return(invisible())
    hits[[id]] <<- list(phase = phase, parameter = parameter, limit = limit,
                        source = source)
  }
  if (is.list(fs$boundary)) {
    for (r in fs$boundary) {
      m <- .hzr_limit_mechanisms[[r$mechanism %||% ""]]
      if (!is.null(m)) add(r$phase, m[[1L]], m[[2L]], r$mechanism)
    }
  }
  # A phase recorded as vanished has unidentified shapes because its mu went
  # to 0: a different limit, already warned (#565), and not a shape at a
  # limit. Its weak reading does not trigger (amended after Task 1: #444's
  # t_half = 3 fit vanished from every grid start).
  vanished <- if (is.list(fs$boundary)) {
    unlist(lapply(fs$boundary, function(r) {
      if (identical(r$mechanism, "coe_phase_vanished")) r$phase
    }))
  }
  if (is.list(fs$weak)) {
    for (p in intersect(fs$weak$params, names(shape_phase))) {
      if (shape_phase[[p]] %in% vanished) next
      short <- substring(p, nchar(shape_phase[[p]]) + 2L)
      add(shape_phase[[p]], .hzr_limit_user_name[[short]],
          "a flat direction (not identified)", "weak")
    }
  }
  if (!length(hits)) NULL else unname(hits)
}

# The grid's nu values. Spread across the sigmoid's usable range on the
# natural scale, away from both limits (0 and large).
.hzr_limit_grid_nu <- c(0.5, 1, 2)

#' Restart starts for a triggered fit
#'
#' Covariate coefficients are zeroed, as SAS's RETRY_T does
#' (hazrd3.c:343-361). For each triggered G1 phase, t_half is set to the
#' 25/50/75% quantiles of the event times and crossed with nu in
#' `.hzr_limit_grid_nu`; for a triggered g3 phase, tau takes the same
#' quantiles. Every other entry keeps `start`, and a fixed entry never moves.
#' Deterministic: no RNG, so a restarted fit is reproducible from its call.
#'
#' @param start Named full theta: the user's start.
#' @param trigger Output of `.hzr_limit_trigger()`.
#' @param phases The fit's phases.
#' @param event_times Times of the rows with an event.
#' @param fixed Names in `start` that must not move.
#' @return A list of named theta vectors, without duplicates.
#' @keywords internal
#' @noRd
.hzr_limit_grid <- function(start, trigger, phases, event_times, fixed) {
  et <- event_times[is.finite(event_times) & event_times > 0]
  if (!length(et)) return(list())
  q <- unique(unname(stats::quantile(et, c(0.25, 0.5, 0.75))))
  shape_phase <- .hzr_limit_shape_names(phases)
  base <- start
  is_cov <- !names(base) %in% c(names(shape_phase), paste0(names(phases), ".log_mu"))
  base[is_cov & !names(base) %in% fixed] <- 0
  set_free <- function(s, key, value) {
    if (key %in% names(s) && !key %in% fixed) s[[key]] <- value
    s
  }
  out <- list()
  for (ph in unique(vapply(trigger, `[[`, "", "phase"))) {
    if (identical(phases[[ph]]$type, "g3")) {
      for (tq in q) out[[length(out) + 1L]] <- set_free(base, paste0(ph, ".log_tau"), log(tq))
    } else {
      for (tq in q) for (nu in .hzr_limit_grid_nu) {
        s <- set_free(base, paste0(ph, ".log_t_half"), log(tq))
        out[[length(out) + 1L]] <- set_free(s, paste0(ph, ".nu"), nu)
      }
    }
  }
  out[!duplicated(lapply(out, unname))]
}

# TRUE while hazard() is evaluating a call on the limit layer's behalf, so
# the inner call fits and returns without checking again. Package-private;
# not a control element, so no user setting can turn the check off.
.hzr_limit_state <- new.env(parent = emptyenv())
.hzr_limit_state$inner <- FALSE

#' Evaluate a hazard() call as the limit layer's inner call
#'
#' Its warnings are held back (a discarded fit's warnings must not reach the
#' user as if they described the result) and its error is returned, not
#' raised.
#' @param cl A `hazard()` call. @param env Where to evaluate it.
#' @return `list(value, warnings)`; `value` is a fit or an error condition.
#' @keywords internal
#' @noRd
.hzr_limit_capture <- function(cl, env) {
  old <- .hzr_limit_state$inner
  .hzr_limit_state$inner <- TRUE
  on.exit(.hzr_limit_state$inner <- old, add = TRUE)
  caught <- list()
  value <- tryCatch(
    withCallingHandlers(eval(cl, env), warning = function(w) {
      caught[[length(caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }),
    error = function(e) e)
  list(value = value, warnings = caught)
}

.hzr_limit_describe <- function(trigger) {
  paste(vapply(trigger, function(r) {
    paste0(r$phase, ".", r$parameter, " ran to ", r$limit)
  }, character(1)), collapse = "; ")
}

#' The refusal condition
#' @keywords internal
#' @noRd
.hzr_limit_refusal <- function(trigger, best_objective, n_tried,
                               disabled = FALSE) {
  first <- trigger[[1L]]
  msg <- paste0(
    "hazard(): ", .hzr_limit_describe(trigger), ". ",
    if (disabled) {
      "No restart was tried (control$limit_restart = FALSE)"
    } else if (n_tried > 0L) {
      paste0("None of ", n_tried, " restarts reached a fit without a ",
             "parameter at a limit")
    } else {
      "No restart start could be built (no event times, or every shape fixed)"
    },
    "; the best objective reached was ", format(best_objective, digits = 8),
    ". The drifted estimate is not returned. Fix the parameter at a value ",
    "you choose, as SAS's FIX", toupper(first$parameter), " does: ",
    "hzr_phase(..., ", first$parameter, " = <value>, fixed = \"",
    first$parameter, "\"), or supply other starting values (`theta`).")
  structure(
    class = c("hzr_limit_refusal", "error", "condition"),
    list(message = msg, call = NULL, phase = first$phase,
         parameter = first$parameter, limit = first$limit,
         best_objective = best_objective, trigger = trigger))
}

.hzr_limit_usable <- function(v) {
  inherits(v, "hazard") && !isFALSE(v$fit$converged) &&
    isTRUE(is.finite(v$fit$objective))
}

#' hazard()'s limit layer: fit, check, restart or refuse
#'
#' @param cl `match.call()` of the user's `hazard()` call.
#' @param env Its `parent.frame()`.
#' @param control The call's evaluated `control`.
#' @return A `hazard` object; or an error is raised.
#' @keywords internal
#' @noRd
.hzr_limit_fit <- function(cl, env, control) {
  restart <- control$limit_restart %||% TRUE
  if (!(isTRUE(restart) || isFALSE(restart))) {
    stop("control$limit_restart must be TRUE or FALSE.", call. = FALSE)
  }
  first <- .hzr_limit_capture(cl, env)
  emit <- function(r) for (w in r$warnings) warning(w)
  if (inherits(first$value, "condition")) {
    emit(first)
    stop(first$value)
  }
  trigger <- .hzr_limit_trigger(first$value)
  if (is.null(trigger)) {
    emit(first)
    return(first$value)
  }
  # On a refusal the drifted fit's own warnings (weak, boundary) are raised
  # first: they describe the drift the refusal names.
  if (!restart) {
    emit(first)
    stop(.hzr_limit_refusal(trigger, first$value$fit$objective, 0L,
                            disabled = TRUE))
  }
  fit0 <- first$value
  prep <- .hzr_evaluate_prepare(fit0)
  # The user's start: their theta if they gave one, else the phase specs'
  # default, assembled as the optimizer assembles it.
  # cl$theta is an unevaluated call, so test its VALUE: a symbol that
  # evaluates to NULL means no theta was given.
  theta_given <- if (is.null(cl$theta)) NULL else eval(cl$theta, env)
  start <- if (!is.null(theta_given)) {
    stats::setNames(theta_given, prep$names)
  } else {
    stats::setNames(unlist(lapply(names(prep$phases), function(nm) {
      .hzr_phase_start(prep$phases[[nm]],
                       n_covariates = prep$covariate_counts[[nm]])
    }), use.names = FALSE), prep$names)
  }
  fixed <- prep$names[as.logical(fit0$fit$fixed_mask %||%
                                   logical(length(prep$names)))]
  grid <- .hzr_limit_grid(start, trigger, prep$phases,
                          event_times = fit0$data$time[fit0$data$status == 1],
                          fixed = fixed)
  ctrl <- control
  ctrl$n_starts <- 1L   # the grid IS the multi-start; 9 x n_starts is waste
  runs <- lapply(grid, function(s) {
    cl_i <- cl
    cl_i$theta <- s
    cl_i$control <- ctrl
    .hzr_limit_capture(cl_i, env)
  })
  rows <- data.frame(
    start = seq_along(runs),
    objective = vapply(runs, function(r) {
      if (.hzr_limit_usable(r$value)) r$value$fit$objective else NA_real_
    }, numeric(1)),
    triggered = vapply(runs, function(r) {
      .hzr_limit_usable(r$value) && !is.null(.hzr_limit_trigger(r$value))
    }, logical(1)),
    error = vapply(runs, function(r) {
      if (inherits(r$value, "condition")) conditionMessage(r$value) else NA_character_
    }, character(1)),
    stringsAsFactors = FALSE)
  # The winner is the best objective OVERALL, the original included. A clean
  # fit that is worse than a fit at a limit is not a rescue: there the data
  # support the limit, and the fit refuses (amended after Task 1: on #448's
  # avc fixture the best clean start was -208.31, against -21.66 at the step
  # and -78.15 returned today).
  best_any <- suppressWarnings(max(c(fit0$fit$objective, rows$objective),
                                   na.rm = TRUE))
  ok <- which(is.finite(rows$objective) & !rows$triggered)
  win <- if (length(ok)) ok[which.max(rows$objective[ok])] else NA_integer_
  # Relative tolerance: grid starts that reach the same optimum differ in
  # the last digits, and a tie must not read as "a limit fit was better".
  if (is.na(win) ||
        best_any - rows$objective[[win]] > 1e-6 * max(1, abs(best_any))) {
    emit(first)
    stop(.hzr_limit_refusal(trigger, best_any, length(runs)))
  }
  out <- runs[[win]]$value
  emit(runs[[win]])
  # The user's call, not the grid start's, so print() shows what they wrote
  # and hzr_bootstrap() refits from their start (through this layer again).
  out$call <- fit0$call
  out$call_env <- fit0$call_env
  out$fit$limit <- list(
    trigger = trigger,
    original = list(objective = fit0$fit$objective, theta = fit0$fit$theta,
                    weak = fit0$fit$weak, boundary = fit0$fit$boundary),
    grid = rows, winner = win)
  warning(structure(
    class = c("hzr_limit_restarted", "warning", "condition"),
    list(message = paste0(
      .hzr_limit_describe(trigger), " (objective ",
      format(fit0$fit$objective, digits = 8), "). Restarted from ",
      length(runs), " starts; returned start ", win, " (objective ",
      format(out$fit$objective, digits = 8), "). See fit$fit$limit."),
      call = NULL)))
  out
}
