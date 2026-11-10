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
