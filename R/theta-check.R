# theta-check.R -- one definition of what a user-supplied theta must satisfy
# before the likelihood sees it, called from every entry point that takes
# one: hazard(), hzr_evaluate() and predict(). Separate copies of these rules
# drifted: hzr_evaluate() counted parameters as the score test does rather
# than as the likelihood does, and a Weibull scale or shape <= 0 reached the
# user as the likelihood's internal Inf sentinel.

#' Refuse a Weibull theta whose scale or shape is not positive (#383)
#'
#' `.hzr_logl_weibull()` returns `Inf` for `mu <= 0` or `nu <= 0`: a sentinel
#' the optimizer's wrapper absorbs, but one that reached the user as
#' `logLik = Inf` from `hzr_evaluate()`, which reads as the best possible fit,
#' and as "non-finite value supplied by optim" from `hazard()`. The model is
#' defined only for positive `mu` and `nu`, since the likelihood takes their
#' logarithms, so a non-positive value is refused and named.
#'
#' @param theta The supplied parameter vector, natural scale.
#' @param dist The distribution; only `"weibull"` is checked.
#' @return `NULL`, invisibly; called for its error.
#' @noRd
.hzr_check_theta <- function(theta, dist) {
  if (!identical(dist, "weibull") || length(theta) < 2L) {
    return(invisible(NULL))
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
