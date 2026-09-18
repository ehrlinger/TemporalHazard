# wald.R -- Wald inference for stepwise covariate selection.
#
# The stepwise engine (hzr_stepwise, Phase 4b of the development plan)
# needs a lightweight Wald-test helper that works against any fitted
# hazard model:
#
#   - Scalar test for a single named coefficient  -> two-sided z-test
#   - Joint test for a set of named coefficients  -> chi-square Wald
#
# Both cases produce a p-value on the same scale so the stepwise loop
# can compare candidates uniformly.
#
# This helper is intentionally forgiving about boundary / fixed
# parameters: if the vcov matrix is absent or any requested coefficient
# has a non-finite variance, it returns p = NA rather than erroring,
# matching the summary.hazard() convention.

#' Wald test for one or more coefficients in a fitted hazard model
#'
#' Computes the Wald statistic and p-value for a subset of coefficients
#' from a fitted `hazard` object.  A single name yields a two-sided
#' z-test on the standardised coefficient; multiple names yield a joint
#' chi-square Wald test on the corresponding submatrix of the vcov.
#'
#' @param fit A fitted `hazard` object produced by [hazard()].
#' @param names Character vector of coefficient names to test.  For a
#'   multiphase fit, a subset of `names(coef(fit))`; for a
#'   single-distribution fit, the positional names (`mu`, `nu`, `beta1`,
#'   ...) whatever `theta` was named.
#'
#' @return A list with elements:
#' \describe{
#'   \item{stat}{Test statistic: `z` if `length(names) == 1L`,
#'     chi-square otherwise.}
#'   \item{df}{Degrees of freedom (1 for scalar, `length(names)` joint).}
#'   \item{p_value}{Two-sided p-value, or `NA_real_` if the test could
#'     not be computed (no vcov or non-finite SE for any coefficient in
#'     `names`).}
#'   \item{names}{Character vector of the coefficients tested.}
#'   \item{estimate}{Numeric vector of coefficient estimates for the
#'     requested `names`.}
#' }
#'
#' @keywords internal
#' @noRd
.hzr_wald_p <- function(fit, names) {
  if (!inherits(fit, "hazard")) {
    stop("`fit` must be a `hazard` object.", call. = FALSE)
  }
  if (!is.character(names) || length(names) == 0L || any(!nzchar(names))) {
    stop("`names` must be a non-empty character vector of coefficient names.",
         call. = FALSE)
  }

  theta <- fit$fit$theta
  if (is.null(theta)) {
    stop("`fit` has no fitted coefficients (did you call hazard(..., fit = TRUE)?).",
         call. = FALSE)
  }

  # Canonical coefficient names: multiphase stores them directly on
  # theta.  A single-distribution theta carries whatever names the user
  # gave its starting values (`c(mu =, nu =, age = 0)`), which the
  # positional `beta<k>` from .hzr_candidate_coef_name() never matches, and
  # a covariate named like a shape (`nu`) would match the shape (#304).  So
  # name those by position, as summary.hazard() does for an unnamed theta.
  # Position means something only when theta holds every coefficient:
  # hazard() fits a shape-only theta, and `beta1` would then name `mu`.
  if (fit$spec$dist != "multiphase") {
    # A windowed fit has a coefficient per covariate per window, and a
    # positional `beta<k>` names only the first window's.
    if (length(fit$spec$time_windows) > 0L) {
      stop("Wald tests of single covariates do not support `time_windows` ",
           "fits: each covariate has a coefficient per window.", call. = FALSE)
    }
    # The likelihood's own count: it ignores `control$shape_param_count`, so
    # honouring that here would let `beta1` name a shape parameter.
    n_shape <- .hzr_shape_parameter_count(fit$spec$dist)
    p <- if (is.null(fit$data$x)) 0L else ncol(fit$data$x)
    if (length(theta) != n_shape + p) {
      stop(
        "The fit's theta has ", length(theta), " value(s), but ",
        # Vowel-aware, not article-free: under hzr_bootstrap(scope = ) this
        # message becomes a replicate's failure-reason key, so every
        # consonant-initial dist keeps its key and only "a exponential" moves.
        if (grepl("^[aeiou]", fit$spec$dist)) "an " else "a ",
        fit$spec$dist, " fit with ", p, " covariate column(s) has ",
        n_shape + p, " coefficients, so they cannot be matched to the ",
        "covariates.  Refit with a starting value for each covariate.",
        call. = FALSE
      )
    }
  }
  coef_names <- if (fit$spec$dist == "multiphase") base::names(theta)
  if (is.null(coef_names) || !all(nzchar(coef_names))) {
    p <- if (is.null(fit$data$x)) 0L else ncol(fit$data$x)
    coef_names <- .hzr_parameter_names(
      theta = unname(theta),
      dist  = fit$spec$dist,
      p     = p
    )
  }

  missing_names <- setdiff(names, coef_names)
  if (length(missing_names) > 0L) {
    stop(
      "Unknown coefficient name(s): ",
      paste(sQuote(missing_names), collapse = ", "),
      ".  Available: ",
      paste(sQuote(coef_names), collapse = ", "),
      call. = FALSE
    )
  }

  idx <- match(names, coef_names)
  estimate <- setNames(theta[idx], names)

  # Access vcov directly on the fit rather than via vcov.hazard(), which
  # collapses to NA if any diagonal entry is NA (common for fixed shape
  # parameters). summary.hazard() handles partial NAs by masking, so
  # match that behaviour.
  vcov_mat <- fit$fit$vcov

  na_result <- function(stat = NA_real_, df = length(names)) {
    list(
      stat     = stat,
      df       = df,
      p_value  = NA_real_,
      names    = names,
      estimate = unname(estimate)
    )
  }

  if (is.null(vcov_mat) || !is.matrix(vcov_mat)) {
    return(na_result())
  }

  v_sub <- vcov_mat[idx, idx, drop = FALSE]
  diag_v <- diag(v_sub)
  if (any(!is.finite(diag_v)) || any(diag_v <= 0)) {
    return(na_result())
  }

  if (length(names) == 1L) {
    # Scalar z-test
    se <- sqrt(diag_v)
    z <- estimate / se
    p <- 2 * stats::pnorm(-abs(z))
    return(list(
      stat     = unname(z),
      df       = 1L,
      p_value  = unname(p),
      names    = names,
      estimate = unname(estimate)
    ))
  }

  # Joint chi-square Wald: W = beta' V^{-1} beta  ~  chi^2_{df}
  v_inv <- tryCatch(solve(v_sub), error = function(e) NULL)
  if (is.null(v_inv)) {
    return(na_result())
  }
  stat <- as.numeric(crossprod(estimate, v_inv %*% estimate))
  if (!is.finite(stat) || stat < 0) {
    return(na_result())
  }
  p <- stats::pchisq(stat, df = length(names), lower.tail = FALSE)

  list(
    stat     = stat,
    df       = length(names),
    p_value  = p,
    names    = names,
    estimate = unname(estimate)
  )
}
