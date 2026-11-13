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
#' mechanism is a limit, a `$fit$weak` direction that names a shape, and a
#' `$fit$weak` of NA from an unusable Hessian (#734).
#' A parameter the user fixed never triggers: holding it there was their
#' choice, as FIXNU is in SAS.
#'
#' @param object A `hazard` object, or a list with `$spec$phases` and `$fit`.
#' @param unusable Whether an NA `$weak` from an unusable Hessian triggers.
#'   `FALSE` reads it as before #734; `.hzr_limit_fit()` judges a hard or
#'   flat-direction trigger's grid starts that way.
#' @return `NULL`, or a list of `list(phase, parameter, limit, source)`.
#' @keywords internal
#' @noRd
.hzr_limit_trigger <- function(object, unusable = TRUE) {
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
  # An NA reading on a fitted fit whose Hessian was unusable (#734): the
  # check could not look, and every measured drift into the m < 0 joint
  # limit ended there. It names no parameter, so each free G1 t_half is
  # named; the grid varies t_half and nu in that phase. Only on a fit nothing
  # else triggers on, so a hard trigger's grid, refusal and rescue are what
  # they were.
  if (unusable && !length(hits) && .hzr_limit_weak_unusable(object)) {
    for (ph in names(phases)) {
      if (!phases[[ph]]$type %in% c("cdf", "hazard") || ph %in% vanished) next
      if (!paste0(ph, ".log_t_half") %in% theta_names) next
      add(ph, "t_half", "no usable Hessian (identification not checked)",
          "weak_na")
    }
  }
  if (!length(hits)) NULL else unname(hits)
}

#' Is `$weak` NA because the Hessian at a fitted optimum was unusable?
#'
#' `$weak` is NA for any reason the check could not look, and only some of
#' those are about the fit (#734). Not one of them: an unfitted or imported
#' object ("fitting" on the record), a fit that did not converge or has no
#' objective (no optimum to take a Hessian at), and numDeriv not installed
#' (a missing Suggests says nothing about the fit, and the estimates must not
#' depend on it). The numDeriv reasons are the strings R/optimizer.R and
#' .hzr_optim_multiphase() record; `$weak` itself is unchanged.
#' @param object A `hazard` object.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
#' @noRd
.hzr_limit_weak_unusable <- function(object) {
  fs <- object$fit
  if (!.hzr_is_na_scalar(fs$weak)) return(FALSE)
  if ("fitting" %in% object$degraded) return(FALSE)
  if (isFALSE(fs$converged) || !isTRUE(is.finite(fs$objective))) return(FALSE)
  causes <- object$degraded_causes[c("standard_errors", "weak_direction_check")]
  !any(grepl("numDeriv not installed", causes, fixed = TRUE))
}

#' Does a trigger hold a record that is not weak-only?
#'
#' A flat direction (`"weak"`) and an unusable Hessian (`"weak_na"`) try the
#' grid and keep the original unless it does better; anything else restarts
#' or refuses.
#' @param trigger Output of `.hzr_limit_trigger()`.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
#' @noRd
.hzr_limit_is_hard <- function(trigger) {
  any(!vapply(trigger, `[[`, "", "source") %in% c("weak", "weak_na", "span"))
}

# The span check's threshold: a move of t_half across the data's whole
# log-span that costs less than this much log-likelihood leaves t_half
# unidentified (#734). It is half the median of a chi-squared on 1 df. A
# move that GAINS more than this shows the fit is not a maximum, and fires
# too, so the grid looks for the better fit (#739 review).
.hzr_limit_span_tau <- 0.5 * stats::qchisq(0.5, 1)

# The limit string a span record carries.
.hzr_limit_span_limit <- paste0(
  "a flat span (t_half moves across the data's whole range at no cost; ",
  "with m < 0 the phase is a power law in t)")

# The limit string when a move RAISES the log-likelihood by more than the
# threshold: the fit is not a maximum along t_half, which is a different
# conclusion from a flat span and must not read as one (#740 review).
.hzr_limit_span_gain <- paste0(
  "a span move that raises the log-likelihood (moving t_half across the ",
  "data's range does better, so the fit is not a maximum)")

# The limit string when no move could be computed and none was measured
# flat: the span check did not run to a result, so identification is not
# known (#739 review). Weak-only, as the flat span is.
.hzr_limit_span_unchecked <- paste0(
  "an unchecked span (moving t_half across the data's range could not be ",
  "computed, so whether it is identified is not known)")

#' Does moving a G1 phase's t_half across the data's span cost nothing?
#'
#' The m < 0 joint limit (#734): nu and t_half run to infinity together at a
#' fixed ((1 - 2^m)^-nu - 1) / t_half, and over the data the phase becomes a
#' power law in t. G at the last time can be far above sqrt(eps), so
#' `phase_beyond_data` does not fire, and the Hessian can be usable, so
#' neither does the weak check. This asks the likelihood directly: for each
#' free t_half of a cdf or hazard phase, move log t_half by L = log(t_max /
#' t_min) of the observed times, up and down; re-fit the phase's own other
#' free shape and scale parameters, first to its cumulative hazard over the
#' observed times (walking the move in 10 steps so a curved ridge is
#' followed), then on the likelihood; and fire when either move changes the
#' log-likelihood by less than `.hzr_limit_span_tau`, or raises it by more
#' (the fit is then not a maximum). Both moves are tried. A move measured
#' flat records a flat span (`.hzr_limit_span_limit`); a move that gains
#' records `.hzr_limit_span_gain` (#740 review); a move that could not
#' be computed, with no move measured flat, records an unchecked span
#' (`.hzr_limit_span_unchecked`), never a flat one (#739 review).
#' Run only where it can matter, a t_half outside the observed times or an
#' m below 0: on a fit inside the data with m >= 0 it never fired in the
#' measurements, and it costs from 30 to 2000 likelihood evaluations.
#'
#' A weak-only trigger: the fit is never refused on it. Measured on
#' 2026-11-10 at 1ed2e69a: it fired on all 14 fits the earlier checks
#' missed, and once on 38 G1 parity fits, at 0.082.
#'
#' @param object A fitted multiphase `hazard` object.
#' @param prep `.hzr_evaluate_prepare(object)`, or NULL to build it.
#' @return `NULL`, or trigger records as `.hzr_limit_trigger()` returns.
#' @keywords internal
#' @noRd
.hzr_limit_span <- function(object, prep = NULL) {
  if (!.hzr_limit_usable(object)) return(NULL)
  th <- object$fit$theta
  nm <- names(th)
  phases <- object$spec$phases
  if (is.null(nm) || !length(phases)) return(NULL)
  fixed <- if (length(object$fit$fixed_mask) == length(th)) {
    nm[as.logical(object$fit$fixed_mask)]
  } else {
    character(0)
  }
  cand <- Filter(function(ph) {
    keys <- paste0(ph, c(".log_t_half", ".nu", ".m"))
    phases[[ph]]$type %in% c("cdf", "hazard") && all(keys %in% nm) &&
      !keys[[1L]] %in% fixed
  }, names(phases))
  if (!length(cand)) return(NULL)
  if (is.null(prep)) {
    prep <- tryCatch(suppressWarnings(.hzr_evaluate_prepare(object)),
                     error = function(e) NULL)
  }
  if (is.null(prep)) return(NULL)
  # Weight-0 rows do not enter the likelihood, so they do not widen the
  # observed range, as in .hzr_limit_event_times().
  w <- prep$weights
  keep <- if (length(w) == length(prep$time)) !is.na(w) & w > 0 else TRUE
  rows <- function(v) if (length(v) == length(prep$time)) v[keep] else v
  ut <- c(rows(prep$time), rows(prep$time_lower), rows(prep$time_upper))
  ut <- sort(unique(ut[is.finite(ut) & ut > 0]))
  if (length(ut) < 2L) return(NULL)
  span <- log(max(ut) / min(ut))
  ll <- function(theta) {
    v <- tryCatch(suppressWarnings(.hzr_logl_multiphase(
      unname(theta), prep$time, prep$status, time_lower = prep$time_lower,
      time_upper = prep$time_upper, x = prep$x, weights = prep$weights,
      phases = prep$phases, covariate_counts = prep$covariate_counts,
      x_list = prep$x_list,
      objective = object$spec$objective %||% "likelihood")),
      error = function(e) NA_real_)
    if (length(v) == 1L && is.finite(v)) v else NA_real_
  }
  ll0 <- ll(th)
  if (is.na(ll0)) return(NULL)
  hits <- list()
  for (ph in cand) {
    pk <- paste0(ph, c(".log_mu", ".log_t_half", ".nu", ".m"))
    lth <- th[[pk[[2L]]]]
    if (lth >= log(min(ut)) && lth <= log(max(ut)) && th[[pk[[4L]]]] >= 0) {
      next
    }
    flat <- FALSE
    gain <- FALSE
    unchecked <- FALSE
    for (dir in c(1, -1)) {
      d <- .hzr_limit_span_move(th, pk, phases[[ph]]$type, lth + dir * span,
                                setdiff(pk, c(fixed, pk[[2L]])), ut, ll)
      # A move that cannot be computed measured nothing: try the other
      # direction before saying anything about the span.
      if (is.na(d)) {
        unchecked <- TRUE
      } else if (d - ll0 > .hzr_limit_span_tau) {
        # Tested first: a gain also passes the flat test below.
        gain <- TRUE
        break
      } else if (ll0 - d < .hzr_limit_span_tau) {
        flat <- TRUE
        break
      }
    }
    if (gain || flat || unchecked) {
      hits[[length(hits) + 1L]] <- list(
        phase = ph, parameter = "t_half",
        limit = if (gain) {
          .hzr_limit_span_gain
        } else if (flat) {
          .hzr_limit_span_limit
        } else {
          .hzr_limit_span_unchecked
        },
        source = "span")
    }
  }
  if (!length(hits)) NULL else hits
}

#' One span move: the best log-likelihood with log t_half held at `to`
#'
#' @param th The fit's named theta. @param pk The phase's four G1 names.
#' @param type `"cdf"` or `"hazard"`. @param to The moved log t_half.
#' @param free The phase's names to re-fit (its free ones but t_half).
#' @param ut The observed times. @param ll The log-likelihood, NA if none.
#' @return The log-likelihood after the move, or NA.
#' @keywords internal
#' @noRd
.hzr_limit_span_move <- function(th, pk, type, to, free, ut, ll) {
  if (!is.finite(to) || to > log(.Machine$double.xmax)) return(NA_real_)
  # The phase's log cumulative hazard over the observed times.
  log_h <- function(t) {
    tt <- exp(t[[pk[[2L]]]])
    if (!is.finite(tt)) return(rep(NA_real_, length(ut)))
    d <- tryCatch(hzr_decompos(ut, tt, t[[pk[[3L]]]], t[[pk[[4L]]]]),
                  error = function(e) NULL)
    if (is.null(d)) return(rep(NA_real_, length(ut)))
    t[[pk[[1L]]]] + if (identical(type, "cdf")) log(d$G) else log(-d$log_surv)
  }
  target <- log_h(th)
  ok <- is.finite(target)
  cur <- th
  if (length(free) && any(ok)) {
    match_cost <- function(v) {
      t2 <- cur
      t2[free] <- v
      # Both shapes negative is outside the G1 family.
      if (t2[[pk[[4L]]]] < 0 && t2[[pk[[3L]]]] < 0) return(1e100)
      h <- suppressWarnings(log_h(t2))[ok]
      if (any(!is.finite(h))) 1e100 else sum((h - target[ok])^2)
    }
    from <- th[[pk[[2L]]]]
    # Walk the move in steps, re-matching at each, so a curved ridge is
    # followed rather than jumped across.
    for (j in seq_len(10L)) {
      cur[[pk[[2L]]]] <- from + (to - from) * j / 10
      o <- tryCatch(stats::optim(cur[free], match_cost, method = "Nelder-Mead",
                                 control = list(maxit = 1000, reltol = 1e-12)),
                    error = function(e) NULL)
      if (!is.null(o)) {
        o2 <- tryCatch(stats::optim(o$par, match_cost, method = "BFGS",
                                    control = list(maxit = 300,
                                                   reltol = 1e-14)),
                       error = function(e) o)
        cur[free] <- if (o2$value <= o$value) o2$par else o$par
      }
    }
  }
  cur[[pk[[2L]]]] <- to
  matched <- ll(cur)
  polished <- NA_real_
  if (length(free)) {
    neg <- function(v) {
      t2 <- cur
      t2[free] <- v
      r <- ll(t2)
      if (is.na(r)) 1e100 else -r
    }
    o3 <- tryCatch(stats::optim(cur[free], neg, method = "BFGS",
                                control = list(maxit = 200, reltol = 1e-12)),
                   error = function(e) NULL)
    if (!is.null(o3) && o3$value < 1e100) polished <- -o3$value
  }
  best <- suppressWarnings(max(matched, polished, na.rm = TRUE))
  if (is.finite(best)) best else NA_real_
}

# The grid's nu values. Spread across the sigmoid's usable range on the
# natural scale, away from both limits (0 and large).
.hzr_limit_grid_nu <- c(0.5, 1, 2)

#' The event times the restart grid is placed on
#'
#' Every row with an event: an exact event at `time`, and a left- or
#' interval-censored one at its upper bound, which is SAS's TIME and where
#' `time_upper` holds it (an interval row's `time` is its lower bound). With
#' exact events only, a fit whose events are all censored to an interval
#' had no grid, and a hard trigger refused without trying a start. A row of
#' weight 0 is left out: the likelihood ignores it, so it must not move the
#' grid (#737 review).
#' @param data The rows the likelihood uses: `.hzr_evaluate_prepare()`'s
#'   output, which has dropped rows a phase design drops (NA covariates).
#' @return A numeric vector.
#' @keywords internal
#' @noRd
.hzr_limit_event_times <- function(data) {
  status <- data$status
  t <- data$time
  up <- data$time_upper
  bounded <- !is.na(status) & status %in% c(-1, 2)
  if (!is.null(up) && length(up) == length(t)) t[bounded] <- up[bounded]
  keep <- !is.na(status) & status %in% c(-1, 1, 2)
  w <- data$weights
  if (!is.null(w) && length(w) == length(t)) keep <- keep & !is.na(w) & w > 0
  t[keep]
}

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
# The grid start and control for the next inner call, read and cleared at its
# entry.
.hzr_limit_state$restart <- NULL
# The user's argument expressions, calling frame, call and its captured
# bindings, and a cache of the masked evaluations, for the next inner call;
# read and cleared at its entry (#683 review).
.hzr_limit_state$inputs <- NULL
# The outer call's cache while the layer runs, for the phase designs built
# where `inputs` does not reach (.hzr_phase_formula_design()).
.hzr_limit_state$cache <- NULL

#' Evaluate one masked argument of hazard(), once per hazard() call
#'
#' `time`, `status`, `time_lower`, `time_upper` and `weights` may name
#' columns of `data`, so hazard() evaluates the expression the user wrote
#' with `data` masking the calling frame. Under the limit layer each inner
#' fit would evaluate it again; `inputs$cache` holds the first evaluation
#' and its warnings, so the expression runs once and every inner fit reads
#' the same value, with the same warnings (#683 review). Every inner fit
#' sees the same `data` value, so the cached value is the one it would
#' compute. Outside the layer there is no cache and this is a plain eval().
#' @param name Argument name. @param exprs Named list of the expressions.
#' @param data The (masking) data. @param inputs `list(env, cache)`.
#' @return The argument's value.
#' @keywords internal
#' @noRd
.hzr_mask_eval <- function(name, exprs, data, inputs) {
  .hzr_limit_once(name, inputs, function() {
    eval(exprs[[name]], data, inputs$env)
  })
}

#' Compute a value once per hazard() call under the limit layer
#'
#' The cache behind .hzr_mask_eval(), for any step of hazard()'s body that
#' runs user code and gives every inner fit the same result: the formula's
#' parse runs the user's code in the formula, so each inner fit would run it
#' again (#683 review). The first call runs `fun()` and stores its value and
#' warnings under `key` in `inputs$cache`; a later call replays the warnings
#' and returns the stored value. Outside the layer there is no cache, and
#' this is `fun()`.
#' @param key Cache key. @param inputs `list(cache)`, or `NULL`.
#' @param fun A function of no arguments.
#' @return The value of `fun()`.
#' @keywords internal
#' @noRd
.hzr_limit_once <- function(key, inputs, fun) {
  cache <- inputs$cache
  if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
    hit <- cache[[key]]
    for (w in hit$warnings) warning(w)
    return(hit$value)
  }
  caught <- list()
  value <- withCallingHandlers(
    fun(),
    warning = function(w) caught[[length(caught) + 1L]] <<- w)
  if (!is.null(cache)) {
    assign(key, list(value = value, warnings = caught), envir = cache)
  }
  value
}

#' Build a phase formula's design once per hazard() call under the layer
#'
#' A phase formula runs the user's code each time its design is built: once
#' in each inner fit, once in the layer's own setup, and once more when a
#' grid start counts the phases' columns (#683 review). The design is built
#' deep in the fit, where the call's `inputs` list does not reach, so
#' .hzr_limit_fit() puts its cache in `.hzr_limit_state$cache` for the
#' duration of the call. An entry is reused only for the same phase,
#' formula and data, so a build on other rows (rows at time 0 dropped, say)
#' runs the formula as it would without the layer. Outside the layer this is
#' .hzr_formula_design().
#' @param nm Phase name. @param formula The phase formula.
#' @param data The data it is evaluated in.
#' @return As .hzr_formula_design().
#' @keywords internal
#' @noRd
.hzr_phase_formula_design <- function(nm, formula, data) {
  cache <- .hzr_limit_state$cache
  if (is.null(cache)) return(.hzr_formula_design(formula, data))
  key <- paste0("phase_design:", nm)
  if (exists(key, envir = cache, inherits = FALSE)) {
    hit <- cache[[key]]
    if (identical(hit$formula, formula) && identical(hit$data, data)) {
      for (w in hit$warnings) warning(w)
      return(hit$value)
    }
    return(.hzr_formula_design(formula, data))
  }
  caught <- list()
  value <- withCallingHandlers(
    .hzr_formula_design(formula, data),
    warning = function(w) caught[[length(caught) + 1L]] <<- w)
  assign(key, list(formula = formula, data = data, value = value,
                   warnings = caught), envir = cache)
  value
}

#' Evaluate a hazard() call as the limit layer's inner call
#'
#' Its warnings are held back (a discarded fit's warnings must not reach the
#' user as if they described the result) and its error is returned, not
#' raised.
#' @param cl A `hazard()` call. @param env Where to evaluate it.
#' @param inputs `NULL` (evaluate `cl` as written), or the user's argument
#'   expressions and call that the layer captured (see hazard()).
#' @param restart `NULL`, or `list(theta, control)` to replace the call's own.
#' @return `list(value, warnings)`; `value` is a fit or an error condition.
#' @keywords internal
#' @noRd
.hzr_limit_capture <- function(cl, env, inputs = NULL, restart = NULL) {
  old <- .hzr_limit_state$inner
  .hzr_limit_state$inner <- TRUE
  .hzr_limit_state$restart <- restart
  .hzr_limit_state$inputs <- inputs
  on.exit({
    .hzr_limit_state$inner <- old
    .hzr_limit_state$restart <- NULL
    .hzr_limit_state$inputs <- NULL
  }, add = TRUE)
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
    paste0(r$phase, ".", r$parameter,
           # A flat direction is not a value the parameter ran to (#736).
           if (r$source %in% c("weak", "weak_na", "span")) ": " else " ran to ",
           r$limit)
  }, character(1)), collapse = "; ")
}

#' One short, stable reason for a refusal, for failure tallies
#' @keywords internal
#' @noRd
.hzr_limit_reason <- function(e) {
  paste0("limit: ", .hzr_limit_describe(e$trigger))
}

#' The always-printed limit line (a line shown only on bad news cannot be
#' told from one its author forgot to write; R/degraded.R)
#' @keywords internal
#' @noRd
.hzr_limit_line <- function(limit) {
  paste0("  Limit check: ", if (is.null(limit)) {
    "none"
  } else if (!is.list(limit)) {
    "not run"
  } else if (identical(limit$outcome, "kept")) {
    paste0("restart tried (", .hzr_limit_describe(limit$trigger),
           "); no start clear of a limit did better, original kept")
  } else if (identical(limit$outcome, "skipped")) {
    paste0("restart not tried (", .hzr_limit_describe(limit$trigger),
           "); control$limit_restart = FALSE")
  } else if (identical(limit$outcome, "refused")) {
    paste0("REFUSED (", .hzr_limit_describe(limit$trigger),
           "); this is the drifted fit, kept for inspection, not a result")
  } else {
    paste0("restarted (", .hzr_limit_describe(limit$trigger), "); start ",
           limit$winner, " of ", nrow(limit$grid), " kept")
  })
}

#' The refusal condition
#' @keywords internal
#' @noRd
.hzr_limit_refusal <- function(trigger, best_objective, n_tried,
                               disabled = FALSE, n_errored = 0L,
                               grid = NULL, fit = NULL,
                               start_objective = best_objective) {
  first <- trigger[[1L]]
  # The objective from the user's own start leads. A restart that also ran to
  # a limit can reach a degenerate likelihood (a +1007 spike on tied times),
  # so the best over all starts is named as that, not as a result (#736).
  higher <- is.finite(best_objective) && (!is.finite(start_objective) ||
    best_objective - start_objective > 1e-6 * max(1, abs(start_objective)))
  msg <- paste0(
    "hazard(): ", .hzr_limit_describe(trigger), ". ",
    if (disabled) {
      "No restart was tried (control$limit_restart = FALSE)"
    } else if (n_tried > 0L) {
      paste0("None of ", n_tried, " restarts reached a fit without a ",
             "parameter at a limit",
             if (n_errored > 0L) {
               paste0(" (", n_errored, " ended in an error; see the error ",
                      "column of the restart grid)")
             })
    } else {
      "No restart start could be built (no event times, or every shape fixed)"
    },
    "; the objective from your start was ",
    format(start_objective, digits = 8),
    if (higher) {
      paste0(" (a restart that also ran to a limit reached ",
             format(best_objective, digits = 8), ")")
    },
    ". The drifted estimate is not returned. Fix the parameter at a value ",
    "you choose, as SAS's FIX", toupper(gsub("_", "", first$parameter)), " does: ",
    "hzr_phase(..., ", first$parameter, " = <value>, fixed = \"",
    first$parameter, "\"), or supply other starting values (`theta`).")
  structure(
    class = c("hzr_limit_refusal", "error", "condition"),
    list(message = msg, call = NULL, phase = first$phase,
         parameter = first$parameter, limit = first$limit,
         best_objective = best_objective,
         start_objective = start_objective, trigger = trigger,
         # The restart grid (start, objective, triggered, error), so a caller
         # can read why each start failed; NULL when none was tried.
         grid = grid,
         # The drifted fit hazard() used to return, for inspection only: its
         # $fit$limit says it was refused, so print() never shows it as a
         # result.
         fit = fit))
}

.hzr_limit_usable <- function(v) {
  inherits(v, "hazard") && !isFALSE(v$fit$converged) &&
    isTRUE(is.finite(v$fit$objective))
}

#' hazard()'s limit layer: fit, check, restart or refuse
#'
#' @param cl The inner `hazard()` call, over the outer call's own arguments.
#' @param env The outer `hazard()` frame, whose arguments are already forced.
#' @param inputs The user's argument expressions, calling frame and call.
#' @param theta,control The call's evaluated `theta` and `control`.
#' @return A `hazard` object; or an error is raised.
#' @keywords internal
#' @noRd
.hzr_limit_fit <- function(cl, env, inputs, theta, control) {
  # Phase designs are built once per call (.hzr_phase_formula_design()).
  old_cache <- .hzr_limit_state$cache
  .hzr_limit_state$cache <- inputs$cache
  on.exit(.hzr_limit_state$cache <- old_cache, add = TRUE)
  first <- .hzr_limit_capture(cl, env, inputs)
  emit <- function(r) for (w in r$warnings) warning(w)
  if (inherits(first$value, "condition")) {
    emit(first)
    stop(first$value)
  }
  # Read after the first fit, which forces `control` where the body does
  # (#737 review).
  restart <- control$limit_restart %||% TRUE
  if (!(isTRUE(restart) || isFALSE(restart))) {
    stop("control$limit_restart must be TRUE or FALSE.", call. = FALSE)
  }
  trigger <- .hzr_limit_trigger(first$value)
  # The m < 0 joint limit reads clean to every check above, so the
  # likelihood is asked directly, on a fit nothing else flagged (#734).
  if (is.null(trigger)) trigger <- .hzr_limit_span(first$value)
  if (is.null(trigger)) {
    emit(first)
    return(first$value)
  }
  # A weak-only trigger (a flat direction, or an unusable Hessian (#734), and
  # nothing else) is an ordinary ridge, not a drift: it tries the grid, and
  # keeps the original unless the grid does better. A hard trigger restarts
  # or refuses.
  hard <- .hzr_limit_is_hard(trigger)
  if (!restart && !hard) {
    # Not refitted, but flagged: the trigger goes on the record, so the fit
    # does not read as examined and clean.
    out <- first$value
    out$fit$limit <- list(trigger = trigger, grid = NULL,
                          winner = NA_integer_, outcome = "skipped")
    emit(first)
    return(out)
  }
  # On a refusal the drifted fit's own warnings (weak, boundary) are raised
  # first: they describe the drift the refusal names.
  # The fit the user's own start reached, marked refused, rides on the
  # refusal so it can be inspected, never returned (#737).
  drifted <- function(grid) {
    d <- first$value
    d$fit$limit <- list(trigger = trigger, grid = grid,
                        winner = NA_integer_, outcome = "refused")
    d
  }
  if (!restart) {
    emit(first)
    stop(.hzr_limit_refusal(trigger, first$value$fit$objective, 0L,
                            disabled = TRUE, fit = drifted(NULL)))
  }
  fit0 <- first$value
  # This rebuilds the first fit's phases, so any warning it raises is one the
  # first fit raised, and emit(first) reports those (#683 review).
  prep <- suppressWarnings(.hzr_evaluate_prepare(fit0))
  # The user's start: their theta if they gave one, else the phase specs'
  # default, assembled as the optimizer assembles it.
  start <- if (!is.null(theta)) {
    stats::setNames(theta, prep$names)
  } else {
    stats::setNames(unlist(lapply(names(prep$phases), function(nm) {
      .hzr_phase_start(prep$phases[[nm]],
                       n_covariates = prep$covariate_counts[[nm]])
    }), use.names = FALSE), prep$names)
  }
  fixed <- prep$names[as.logical(fit0$fit$fixed_mask %||%
                                   logical(length(prep$names)))]
  grid <- .hzr_limit_grid(start, trigger, prep$phases,
                          event_times = .hzr_limit_event_times(prep),
                          fixed = fixed)
  ctrl <- control
  ctrl$n_starts <- 1L   # the grid IS the multi-start; 9 x n_starts is waste
  runs <- lapply(grid, function(s) {
    .hzr_limit_capture(cl, env, inputs,
                       restart = list(theta = s, control = ctrl))
  })
  # A grid start is clean when it does not trigger. Under an unusable-Hessian
  # trigger, a start whose Hessian is unusable too is in the same unchecked
  # state, not a rescue (#734: on the #551 seed-15 fixture every start drifted
  # to the same limit, one 0.0017 above the original). A span trigger judges
  # its starts the same way: an unchecked start "rescued" a tied-times fit
  # with a degenerate spike (+772 against -147.6). Every other trigger judges
  # its starts as before #734.
  strict <- all(vapply(trigger, `[[`, "", "source") %in% c("weak_na", "span"))
  rows <- data.frame(
    start = seq_along(runs),
    objective = vapply(runs, function(r) {
      if (.hzr_limit_usable(r$value)) r$value$fit$objective else NA_real_
    }, numeric(1)),
    triggered = vapply(runs, function(r) {
      .hzr_limit_usable(r$value) &&
        !is.null(.hzr_limit_trigger(r$value, unusable = strict))
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
  tol <- 1e-6 * max(1, abs(best_any))
  # Under a span trigger a start that does better can sit at the same power
  # law; it wins only if its own span check is clean. Starts are examined
  # best first, and only those that beat the original (#734).
  if (any(vapply(trigger, `[[`, "", "source") == "span")) {
    for (i in ok[order(rows$objective[ok], decreasing = TRUE)]) {
      if (rows$objective[[i]] - fit0$fit$objective <= tol) break
      if (is.null(.hzr_limit_span(runs[[i]]$value))) break
      rows$triggered[[i]] <- TRUE
    }
    ok <- which(is.finite(rows$objective) & !rows$triggered)
    win <- if (length(ok)) ok[which.max(rows$objective[ok])] else NA_integer_
  }
  if (!hard) {
    # Weak-only: a rescue only when the grid beats the original; otherwise
    # keep the original, with the tried grid on the record.
    if (is.na(win) || rows$objective[[win]] - fit0$fit$objective <= tol) {
      out <- fit0
      out$fit$limit <- list(
        trigger = trigger,
        original = list(objective = fit0$fit$objective, theta = fit0$fit$theta,
                        weak = fit0$fit$weak, boundary = fit0$fit$boundary),
        grid = rows, winner = NA_integer_, outcome = "kept")
      emit(first)
      return(out)
    }
  } else if (is.na(win) || best_any - rows$objective[[win]] > tol) {
    emit(first)
    stop(.hzr_limit_refusal(trigger, best_any, length(runs),
                            n_errored = sum(!is.na(rows$error)), grid = rows,
                            fit = drifted(rows),
                            start_objective = fit0$fit$objective))
  }
  out <- runs[[win]]$value
  emit(runs[[win]])
  # The user's call, not the grid start's, so print() shows what they wrote
  # and hzr_bootstrap() refits from their start (through this layer again).
  out$call <- fit0$call
  out$call_env <- fit0$call_env
  # Likewise the user's control: the grid forced n_starts = 1, and
  # hzr_bootstrap() and hzr_stepwise() refit from spec$control. The state the
  # winner's own fit applied (conserve_applied and its reason) stays the
  # winner's, since CoE may apply in one fit and not the other.
  applied <- c("conserve_applied", "conserve_disabled_reason")
  ctl <- fit0$spec$control
  win_ctl <- out$spec$control
  for (k in applied) ctl[k] <- list(win_ctl[[k]])
  out$spec$control <- ctl
  out$fit$limit <- list(
    trigger = trigger,
    original = list(objective = fit0$fit$objective, theta = fit0$fit$theta,
                    weak = fit0$fit$weak, boundary = fit0$fit$boundary),
    grid = rows, winner = win, outcome = "restarted")
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
