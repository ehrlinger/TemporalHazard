# Evaluate a model at supplied parameters
#
# Parity work needs the likelihood at SAS's converged estimates, not at R's
# own optimum, and the model's shape there. Before this there was no way to
# ask for it: `hazard(fit = FALSE)` returns an object holding starting
# values, and `predict()` on one used to fail with "argument of length 0"
# (#144). It still refuses, by design -- predicting from starting values as
# though they were estimates is the defect this package exists to avoid --
# and this is the supported way to evaluate a model at parameters you supply.

#' Evaluate a hazard model at parameters you supply
#'
#' Computes a model's log-likelihood, and optionally its hazard and
#' cumulative hazard, at parameters **you** supply rather than at parameters
#' fitted from the data. This is what a parity check needs: the likelihood at
#' another program's converged estimates, evaluated by this package's own
#' likelihood.
#'
#' The result is not a fit and does not pretend to be one. It carries no
#' standard errors, no convergence status and no covariance matrix, because
#' nothing was estimated: [print()] says so on its first line, and the
#' parameters are labelled as supplied. Where a fitted model is what you
#' want, use `hazard(..., fit = TRUE)`.
#'
#' A phase built with `hzr_phase(constraint = )` derives one of its shapes
#' from the others, and that rule is applied to `theta` here as the fit
#' applies it, so the derived entry you pass is replaced rather than
#' used as given. At a fitted model's own estimates this returns that fit's
#' objective, with one exception: under Conservation of Events the fit
#' re-solves the conserved scale after recording its objective, so a fit
#' whose likelihood is steep in that scale can report a value it is not at.
#' Then this function returns the likelihood at the estimates, and the two
#' differ.
#'
#' @param object A `hazard` object, fitted or built with `fit = FALSE`. Its
#'   data, distribution and phase specification are used; its own `theta` is
#'   not.
#' @param theta Numeric vector of parameters to evaluate at, on the internal
#'   scale the model uses (see [hzr_theta_names()]). Its length must match
#'   the model's parameter count. Names, when supplied, must match too.
#' @param times Optional numeric vector of times, for multiphase models
#'   only. When given, the hazard and cumulative hazard at
#'   those times are returned for a covariate-free ("baseline") subject:
#'   every covariate at 0, so the curve is the shape the supplied parameters
#'   describe, not a prediction for any row of the data. The other families
#'   have no internal shape function that takes parameters directly, and
#'   writing their hazards out here would duplicate `predict()`; `times` is
#'   refused for them rather than mirrored.
#'
#' @return An object of class `hzr_evaluation`: a list with `theta`, the
#'   parameters the likelihood was evaluated at -- the vector supplied, with
#'   any constrained entry replaced by the value its phase derives, which is
#'   warned about as [hazard()] warns; `logLik`, the log-likelihood there;
#'   `dist`; `n_obs` and `n_events`, the rows the likelihood scored and the
#'   exact events among them (a left- or interval-censored row counts in
#'   `n_obs`, not in `n_events`); and, when `times` was given, `curve`, a
#'   data frame of `time`, `hazard` and `cumulative_hazard`.
#'
#' @seealso [hazard()] to fit a model, [hzr_theta_names()] for the parameter
#'   order.
#'
#' @examples
#' data(avc, package = "TemporalHazard")
#' spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
#'                dist = "weibull", theta = c(1, 1), fit = FALSE)
#' # The likelihood at parameters from somewhere else, e.g. a SAS run.
#' hzr_evaluate(spec, theta = c(0.05, 0.9))
#' @export
hzr_evaluate <- function(object, theta, times = NULL) {
  if (!inherits(object, "hazard")) {
    stop("'object' must be a 'hazard' object, from hazard().", call. = FALSE)
  }
  if (is.null(object$data) || is.null(object$data$time)) {
    stop("'object' carries no data to evaluate against. Build it with ",
         "hazard(formula, data = ) or hazard(time = , status = ).",
         call. = FALSE)
  }
  if (!is.numeric(theta) || !length(theta) || anyNA(theta) ||
        any(!is.finite(theta))) {
    stop("'theta' must be a finite numeric vector of parameters to ",
         "evaluate at.", call. = FALSE)
  }
  if (!is.null(times) &&
        (!is.numeric(times) || !length(times) || any(!is.finite(times)) ||
           any(times <= 0))) {
    # Strictly positive: every phase's cumulative hazard at 0 is 0, and a
    # shape's hazard there need not be finite, so a curve including 0 reads
    # as a number where there is none (#144 review).
    stop("'times' must be a numeric vector of finite positive times, or ",
         "NULL. Every phase's cumulative hazard at time 0 is 0.",
         call. = FALSE)
  }

  prepared <- .hzr_evaluate_prepare(object)
  dist <- object$spec$dist
  if (!length(prepared$time)) {
    # A log-likelihood of 0 over no rows is the best value there is, and it
    # reads as a parity success (#144 review).
    stop(if (length(object$data$time)) {
           paste0("No rows are left to evaluate: all ",
                  length(object$data$time), " were dropped for a missing ",
                  "value in a phase's design.")
         } else {
           "This model carries no observations to evaluate."
         },
         " There is no likelihood to report.", call. = FALSE)
  }
  if (length(theta) != prepared$n_par) {
    stored <- length(object$fit$theta)
    stop(if (identical(dist, "multiphase") && stored &&
               stored < prepared$n_par) {
           paste0("This model needs ", prepared$n_par, " parameters, and ",
                  "the object's own stored theta has ", stored,
                  ": hazard(fit = FALSE) does not resolve phase designs, so ",
                  "a specification whose phases carry covariates stores ",
                  "fewer than a fit would use. Supply the full vector. ")
         } else {
           ""
         },
         if (identical(dist, "multiphase")) {
           # The same count and sentence as hazard(fit = TRUE) (#408).
           .hzr_theta_length_message(
             length(theta),
             .hzr_phase_theta_counts(
               prepared$phases, covariate_counts = prepared$covariate_counts
             )
           )
         } else {
           paste0(
             "'theta' has ", length(theta), " ",
             if (length(theta) == 1L) "entry" else "entries",
             ", but this ", dist, " model takes ", prepared$n_par,
             if (is.null(prepared$names)) "." else
               paste0(": ", paste(utils::head(prepared$names, 3L),
                                  collapse = ", "),
                      if (length(prepared$names) > 3L) ", ..." else "", ".")
           )
         },
         call. = FALSE)
  }
  if (!is.null(names(theta)) && !is.null(prepared$names) &&
        !identical(names(theta), prepared$names)) {
    stop("'theta' is named, and its names are not the model's parameter ",
         "names (", paste(utils::head(prepared$names, 3L), collapse = ", "),
         if (length(prepared$names) > 3L) ", ..." else "",
         "). Supply it unnamed, or in that order.", call. = FALSE)
  }
  if (is.null(names(theta)) && !is.null(prepared$names)) {
    names(theta) <- prepared$names
  }

  if (identical(dist, "multiphase")) {
    # As hazard(fit = FALSE) does with a supplied theta: derive the
    # constrained entries and say which moved, so the vector reported is the
    # vector scored (#144 review).
    constrained <- .hzr_constrain_supplied_theta(
      unname(theta), prepared$phases, prepared$covariate_counts
    )
    names(constrained) <- names(theta)
    theta <- constrained
  }

  logl <- .hzr_logl_at(object, theta, prepared)
  curve <- if (is.null(times)) {
    NULL
  } else {
    .hzr_evaluate_curve(object, theta, times, prepared)
  }

  out <- list(
    theta = theta,
    logLik = logl,
    dist = dist,
    # The rows the likelihood actually scored, not the rows the object
    # carries: a phase design with an NA drops rows (#144 review).
    n_obs = length(prepared$time),
    n_events = sum(prepared$status == 1),
    curve = curve
  )
  structure(out, class = "hzr_evaluation")
}


#' What a model would be fitted with: designs, rows and parameter count
#'
#' `hazard()` expands the design for `time_windows` before fitting and
#' resolves per-phase designs at fit time, so an object built with
#' `fit = FALSE` does not carry either. Evaluating against what the object
#' happens to store scored a different model: on a `time_windows` fit the
#' log-likelihood came back -26390.82 where the fit reported -191.83, with
#' no warning (#144 review). This rebuilds what the fit would have used, by
#' the same functions the fit uses.
#'
#' @param object A `hazard` object.
#' @return A list with `time`, `status`, `time_lower`, `time_upper`, `x`,
#'   `weights` (row-aligned, window-expanded), `x_list` and
#'   `covariate_counts` for a multiphase model, `n_par`, the model's
#'   parameter count, and `names`, the parameter names where the model has
#'   them.
#' @keywords internal
#' @noRd
.hzr_evaluate_prepare <- function(object) {
  d <- object$data
  dist <- object$spec$dist
  x <- d$x
  if (!is.null(object$spec$time_windows) && !is.null(x)) {
    # As hazard() does before fitting.
    x <- .hzr_expand_time_varying_design(
      x = x, time = d$time, time_windows = object$spec$time_windows
    )
  }
  out <- list(time = d$time, status = d$status, time_lower = d$time_lower,
              time_upper = d$time_upper, x = x, weights = d$weights)
  if (identical(dist, "multiphase")) {
    built <- .hzr_multiphase_designs(
      d$time, d$status, time_lower = d$time_lower, time_upper = d$time_upper,
      x = x, weights = d$weights, phases = object$spec$phases,
      data = d$frame
    )
    out <- built[c("time", "status", "time_lower", "time_upper", "x",
                   "weights")]
    out$x_list <- built$x_list
    out$covariate_counts <- built$covariate_counts
    # The fit's own naming, not hzr_theta_names() over colnames(): a design
    # from the vector interface has no column names, and that route fell back
    # to the stored theta's length -- the vacuous check again (#144 review).
    # .hzr_phase_cov_names() synthesises x1, x2, ... for exactly that case.
    phases_v <- .hzr_validate_phases(object$spec$phases)
    out$names <- .hzr_theta_names_list(
      phases_v,
      cov_names = .hzr_phase_cov_names(phases_v, built$covariate_counts,
                                       built$x_list)
    )
    out$n_par <- length(out$names)
    out$phases <- phases_v
    return(out)
  }
  out$n_par <- .hzr_shape_parameter_count(dist, control = object$spec$control) +
    (if (is.null(x)) 0L else ncol(x))
  # Only if they describe THIS model: a stored theta of the wrong length
  # would otherwise be pasted onto a vector of another, and `names<-` errors
  # with a base-R message about attribute lengths (#144 review).
  stored_names <- names(object$fit$theta)
  out$names <- if (length(stored_names) == out$n_par) stored_names else NULL
  out
}


#' The log-likelihood of a model's data at supplied parameters
#'
#' One dispatcher for every distribution, so `hzr_evaluate()` scores with the
#' same likelihood the fit maximises, over the rows and designs
#' `.hzr_evaluate_prepare()` rebuilt.
#'
#' @param object A `hazard` object.
#' @param theta Parameters, on the internal scale.
#' @param prepared The result of `.hzr_evaluate_prepare()`.
#' @return A single log-likelihood.
#' @keywords internal
#' @noRd
.hzr_logl_at <- function(object, theta, prepared) {
  dist <- object$spec$dist
  args <- list(theta = unname(theta), time = prepared$time,
               status = prepared$status, time_lower = prepared$time_lower,
               time_upper = prepared$time_upper, weights = prepared$weights)
  if (identical(dist, "multiphase")) {
    return(.hzr_logl_multiphase(
      unname(theta), prepared$time, prepared$status,
      time_lower = prepared$time_lower, time_upper = prepared$time_upper,
      x = prepared$x, weights = prepared$weights,
      phases = .hzr_validate_phases(object$spec$phases),
      covariate_counts = prepared$covariate_counts, x_list = prepared$x_list,
      objective = object$spec$objective %||% "likelihood"
    ))
  }
  fn <- switch(
    dist,
    weibull = .hzr_logl_weibull,
    exponential = .hzr_logl_exponential,
    lognormal = .hzr_logl_lognormal,
    loglogistic = .hzr_logl_loglogistic,
    stop("hzr_evaluate() does not know how to evaluate dist = '", dist, "'.",
         call. = FALSE)
  )
  do.call(fn, c(args, list(x = prepared$x)))
}


#' Hazard and cumulative hazard at supplied parameters, for a baseline subject
#'
#' Every covariate is held at 0, so the curve is the shape the parameters
#' describe rather than a prediction for a row of the data. The covariate
#' COUNTS are the model's own: `theta` carries one slot per covariate per
#' phase, and splitting it as though the model had none read a later phase's
#' scale out of an earlier phase's coefficient -- 18 to 23 times wrong,
#' monotone and finite, with no warning (#144 review). Zero-valued designs
#' give a covariate-free subject without changing the split.
#'
#' @param object A `hazard` object.
#' @param theta Parameters, on the internal scale.
#' @param times Times to evaluate at.
#' @param prepared The result of `.hzr_evaluate_prepare()`.
#' @return A data frame of `time`, `hazard`, `cumulative_hazard`.
#' @keywords internal
#' @noRd
.hzr_evaluate_curve <- function(object, theta, times, prepared) {
  dist <- object$spec$dist
  if (!identical(dist, "multiphase")) {
    # Only the multiphase model has internal shape functions that take
    # parameters directly. Writing the other families' hazards out here
    # would duplicate what predict() computes inline, and two copies of the
    # same formula drift (#303 is that defect). So this is refused rather
    # than mirrored.
    stop("'times' is supported for dist = \"multiphase\" only; this model ",
         "is \"", dist, "\". Its log-likelihood at the supplied parameters ",
         "is returned either way. For a curve, fit the model and use ",
         "predict(), or evaluate the distribution directly.", call. = FALSE)
  }
  phases <- prepared$phases
  counts <- prepared$covariate_counts
  x_list <- stats::setNames(lapply(names(phases), function(nm) {
    k <- counts[[nm]] %||% 0L
    if (k > 0L) matrix(0, nrow = length(times), ncol = k) else NULL
  }), names(phases))
  data.frame(
    time = times,
    hazard = as.numeric(.hzr_multiphase_hazard(times, theta, phases, counts,
                                               x_list)),
    cumulative_hazard = as.numeric(.hzr_multiphase_cumhaz(times, theta,
                                                          phases, counts,
                                                          x_list))
  )
}


#' @export
print.hzr_evaluation <- function(x, ...) {
  cat("Model evaluated at SUPPLIED parameters -- not a fit\n")
  cat("  nothing was estimated here: no standard errors, no convergence,\n")
  cat("  no covariance. Use hazard(fit = TRUE) to fit.\n\n")
  cat("  distribution: ", x$dist, "\n", sep = "")
  cat("  observations: ", x$n_obs, " (", x$n_events, " events)\n", sep = "")
  cat("  logLik at the supplied parameters: ",
      format(x$logLik, digits = 8), "\n", sep = "")
  cat("\n  parameters supplied:\n")
  print(x$theta)
  if (!is.null(x$curve)) {
    cat("\n  hazard at ", nrow(x$curve),
        " time(s), for a covariate-free subject:\n", sep = "")
    print(utils::head(x$curve, 10L))
    if (nrow(x$curve) > 10L) {
      cat("  ... ", nrow(x$curve) - 10L, " more\n", sep = "")
    }
  }
  invisible(x)
}
