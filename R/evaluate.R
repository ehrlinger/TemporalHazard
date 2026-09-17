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
#' @param object A `hazard` object, fitted or built with `fit = FALSE`. Its
#'   data, distribution and phase specification are used; its own `theta` is
#'   not.
#' @param theta Numeric vector of parameters to evaluate at, on the internal
#'   scale the model uses (see [hzr_theta_names()]). Its length must match
#'   the model's parameter count. Names, when supplied, must match too.
#' @param times Optional numeric vector of times, for `dist =
#'   `"multiphase"` only. When given, the hazard and cumulative hazard at
#'   those times are returned for a covariate-free ("baseline") subject:
#'   every covariate at 0, so the curve is the shape the supplied parameters
#'   describe, not a prediction for any row of the data. The other families
#'   have no internal shape function that takes parameters directly, and
#'   writing their hazards out here would duplicate `predict()`; `times` is
#'   refused for them rather than mirrored.
#'
#' @return An object of class `hzr_evaluation`: a list with `theta` (the
#'   parameters supplied), `logLik` (the log-likelihood there), `dist`,
#'   `n_obs`, `n_events`, and, when `times` was given, `curve`, a data frame
#'   of `time`, `hazard` and `cumulative_hazard`.
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
           any(times < 0))) {
    stop("'times' must be a numeric vector of finite non-negative times, ",
         "or NULL.", call. = FALSE)
  }

  d <- object$data
  dist <- object$spec$dist
  expected <- object$fit$theta
  if (!is.null(expected) && length(theta) != length(expected)) {
    stop("'theta' has ", length(theta), " parameter",
         if (length(theta) == 1L) "" else "s", ", but this ", dist,
         " model has ", length(expected), ". See hzr_theta_names().",
         call. = FALSE)
  }
  known <- names(expected)
  if (is.null(known) && identical(dist, "multiphase")) {
    known <- tryCatch(
      hzr_theta_names(object$spec$phases,
                      covariates = lapply(object$fit$x_list, colnames)),
      error = function(e) NULL
    )
  }
  if (!is.null(names(theta)) && !is.null(known) &&
        !identical(names(theta), known)) {
    stop("'theta' is named, and its names are not the model's parameter ",
         "names (", paste(utils::head(known, 3L), collapse = ", "),
         if (length(known) > 3L) ", ..." else "",
         "). Supply it unnamed, or in that order.", call. = FALSE)
  }
  if (is.null(names(theta)) && !is.null(names(expected))) {
    names(theta) <- names(expected)
  }

  logl <- .hzr_logl_at(object, theta)

  out <- list(
    theta = theta,
    logLik = logl,
    dist = dist,
    n_obs = length(d$time),
    n_events = sum(d$status == 1),
    curve = if (is.null(times)) NULL else .hzr_evaluate_curve(object, theta,
                                                              times)
  )
  structure(out, class = "hzr_evaluation")
}


#' The log-likelihood of a model's data at supplied parameters
#'
#' One dispatcher for every distribution, so `hzr_evaluate()` scores with the
#' same likelihood the fit maximises. The multiphase designs come from
#' `.hzr_multiphase_designs()`, the function the optimizer itself uses, so an
#' unfitted object is evaluated against the design a fit would have built.
#'
#' @param object A `hazard` object.
#' @param theta Parameters, on the internal scale.
#' @return A single log-likelihood.
#' @keywords internal
#' @noRd
.hzr_logl_at <- function(object, theta) {
  d <- object$data
  dist <- object$spec$dist
  args <- list(theta = unname(theta), time = d$time, status = d$status,
               time_lower = d$time_lower, time_upper = d$time_upper,
               weights = d$weights)
  if (identical(dist, "multiphase")) {
    built <- .hzr_multiphase_designs(
      d$time, d$status, time_lower = d$time_lower, time_upper = d$time_upper,
      x = d$x, weights = d$weights, phases = object$spec$phases,
      data = d$frame
    )
    return(.hzr_logl_multiphase(
      unname(theta), built$time, built$status,
      time_lower = built$time_lower, time_upper = built$time_upper,
      x = built$x, weights = built$weights,
      phases = .hzr_validate_phases(object$spec$phases),
      covariate_counts = built$covariate_counts, x_list = built$x_list,
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
  do.call(fn, c(args, list(x = d$x)))
}


#' Hazard and cumulative hazard at supplied parameters, for a baseline subject
#'
#' Every covariate at 0, so the curve is the shape the parameters describe
#' rather than a prediction for a row of the data. That keeps the result from
#' reading as `predict()` output for a fit.
#'
#' @param object A `hazard` object.
#' @param theta Parameters, on the internal scale.
#' @param times Times to evaluate at.
#' @return A data frame of `time`, `hazard`, `cumulative_hazard`.
#' @keywords internal
#' @noRd
.hzr_evaluate_curve <- function(object, theta, times) {
  dist <- object$spec$dist
  if (identical(dist, "multiphase")) {
    phases <- .hzr_validate_phases(object$spec$phases)
    counts <- stats::setNames(integer(length(phases)), names(phases))
    x_list <- stats::setNames(vector("list", length(phases)), names(phases))
    haz <- .hzr_multiphase_hazard(times, unname(theta), phases, counts,
                                  x_list)
    cum <- .hzr_multiphase_cumhaz(times, unname(theta), phases, counts,
                                  x_list)
  } else {
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
  data.frame(time = times, hazard = as.numeric(haz),
             cumulative_hazard = as.numeric(cum))
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
