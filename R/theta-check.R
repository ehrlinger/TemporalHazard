# theta-check.R -- one definition of what a user-supplied theta must satisfy
# before the likelihood sees it, called from every entry point that takes
# one: hazard(), hzr_evaluate() and predict(). Separate copies of these rules
# drifted: hzr_evaluate() counted parameters as the score test does rather
# than as the likelihood does, and a Weibull scale or shape <= 0 reached the
# user as the likelihood's internal Inf sentinel.

#' Refuse a single-distribution theta the likelihood cannot use (#375, #383)
#'
#' Two rules, length first. **Length (#375):** one entry per parameter, the
#' distribution's shape parameters then one coefficient per design column,
#' counted as the likelihood counts them: `control$shape_param_count` is
#' ignored, as `wald.R`, the score test and the stepwise refit and trace
#' ignore it (#489). The only check was a lower bound against the
#' design, so a short theta fitted with a covariate silently dropped and a
#' long one returned its starting values unfitted. Checked for the four
#' distributions whose count is known; an unsupported one is left to the
#' refusal that already exists for it, rather than told it takes 0.
#'
#' **Positivity (#383):**
#' `.hzr_logl_weibull()` returns `Inf` for `mu <= 0` or `nu <= 0`: a sentinel
#' the optimizer's wrapper absorbs, but one that reached the user as
#' `logLik = Inf` from `hzr_evaluate()`, which reads as the best possible fit,
#' and as "non-finite value supplied by optim" from `hazard()`. The model is
#' defined only for positive `mu` and `nu`, since the likelihood takes their
#' logarithms, so a non-positive value is refused and named.
#'
#' **Representability (#566):** a fitted `mu` can be `Inf`, or a subnormal
#' number. The model depends on it only through `nu * log(mu)`, and with a
#' covariate far from zero that logarithm leaves the range a double holds at
#' a sound maximum. `Inf` passed the positivity rule, and `predict()` then
#' returned survival 0 and cumulative hazard `Inf` as if they were results; a
#' subnormal `mu` has lost most of its digits, and gave a cumulative hazard of
#' 0.4888 where the truth was 0.5243. A `mu` or `nu` that is not finite, or
#' is positive but below the smallest normal double, is refused and named,
#' ahead of positivity.
#'
#' @param theta The supplied parameter vector, natural scale.
#' @param dist The distribution.
#' @param n_coef Design columns, after any time-window expansion; `NULL`
#'   skips the length check.
#' @param windowed Whether the design was expanded by time windows, for the
#'   message.
#' @return `NULL`, invisibly; called for its error.
#' @noRd
.hzr_check_theta <- function(theta, dist, n_coef = NULL, windowed = FALSE) {
  if (!is.null(n_coef) &&
        dist %in% c("weibull", "exponential", "loglogistic", "lognormal")) {
    n_shape <- .hzr_shape_parameter_count(dist)
    if (length(theta) != n_shape + n_coef) {
      stop("'theta' has ", length(theta), " ",
           if (length(theta) == 1L) "entry" else "entries",
           ", but this ", dist, " model takes ", n_shape + n_coef, ": ",
           n_shape, " shape parameter", if (n_shape == 1L) "" else "s",
           ", then one coefficient per column of the design (", n_coef,
           " column", if (n_coef == 1L) "" else "s",
           if (windowed) ", one per covariate per time window",
           ").", call. = FALSE)
    }
  }
  if (!identical(dist, "weibull") || length(theta) < 2L) {
    return(invisible(NULL))
  }
  unrep <- c(if (.hzr_unrepresentable(theta[[1L]])) {
               paste0("scale mu = ", format(theta[[1L]]))
             },
             if (.hzr_unrepresentable(theta[[2L]])) {
               paste0("shape nu = ", format(theta[[2L]]))
             })
  if (length(unrep)) {
    stop("'theta' gives Weibull ", paste(unrep, collapse = " and "),
         ", which cannot be represented: the value is outside the range a ",
         "double holds at full precision. Where it came from a fit, a ",
         "covariate far from zero is the usual cause; centre or rescale the ",
         "covariates and refit.", call. = FALSE)
  }
  parts <- c(if (theta[[1L]] <= 0) {
               paste0("scale mu = ", format(theta[[1L]], digits = 4))
             },
             if (theta[[2L]] <= 0) {
               paste0("shape nu = ", format(theta[[2L]], digits = 4))
             })
  if (length(parts)) {
    both <- length(parts) == 2L
    stop("'theta' gives Weibull ", paste(parts, collapse = " and "),
         if (both) ", which must both be positive" else
           ", which must be positive",
         ": the likelihood takes ", if (both) "their logarithms" else
           "its logarithm", ".", call. = FALSE)
  }
  invisible(NULL)
}

#' Whether a positive parameter's value is outside what a double holds
#'
#' `TRUE` for a non-finite value, and for a positive one below
#' `.Machine$double.xmin`: subnormal numbers carry fewer than 53 bits, so
#' `log()` of one is not the logarithm the fit found (#566). Zero and
#' negative values are `FALSE` here; the positivity rule names those.
#' @noRd
.hzr_unrepresentable <- function(value) {
  !is.finite(value) || (value > 0 && value < .Machine$double.xmin)
}
