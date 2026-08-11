#' @importFrom stats optim
#' @keywords internal
NULL

# optimizer.R -- Generic parametric hazard optimizer
#
# Consolidates the common optimization pattern shared across all distributions.
# Each distribution provides its own log-likelihood and gradient functions;
# this module handles the objective wrapping, optim() call, and post-fit
# Hessian/vcov computation.

#' Generic optimizer for parametric hazard likelihoods
#'
#' Maximises a log-likelihood by minimising its negation via \code{stats::optim}.
#' Handles constrained (L-BFGS-B) and unconstrained (BFGS) optimization,
#' post-fit Hessian computation, and vcov extraction.
#'
#' @param logl_fn Function(theta, time, status, time_lower, time_upper, x, ...)
#'   returning the scalar log-likelihood.  Must accept \code{return_gradient}
#'   argument (not used here directly, but the function signature must allow it).
#' @param gradient_fn Function(theta, time, status, time_lower, time_upper, x, ...)
#'   returning the score vector (gradient of the log-likelihood, NOT negated).
#' @param time Numeric vector of follow-up times.
#' @param status Numeric vector of event/censoring indicators.
#' @param time_lower Optional lower bounds for interval censoring.
#' @param time_upper Optional upper bounds for left/interval censoring.
#' @param x Optional design matrix.
#' @param theta_start Starting parameter vector.
#' @param control List of control options (maxit, reltol, abstol).
#' @param use_bounds Logical; if TRUE, use L-BFGS-B with lower bounds on
#'   constrained parameters.  Used for Weibull where mu/nu are on the natural
#'   scale.
#' @param lower_bounds Numeric vector of lower bounds for L-BFGS-B, the same
#'   length as \code{theta_start}.  Only used when \code{use_bounds = TRUE}.
#'   Defaults to \code{rep(1e-6, length(theta_start))} when NULL, but callers
#'   should supply distribution-specific bounds (e.g. \code{c(1e-6, 1e-6,
#'   rep(-Inf, p_cov))} for Weibull so that covariate betas are unconstrained).
#' @param hessian_fn Optional function(theta) returning the Hessian of the
#'   negative log-likelihood at \code{theta} (same scale as
#'   \code{numDeriv::hessian(objective)}).  When it returns a conformant square
#'   matrix (dimension equal to the number of parameters), that matrix is used
#'   for standard errors.  The numerical Hessian is used instead when
#'   \code{hessian_fn} is \code{NULL} (the default), when the function returns
#'   \code{NULL} (e.g. a censoring branch it does not cover analytically), or
#'   when it errors.  A non-NULL, non-conformant return raises a warning and
#'   also falls back to the numerical Hessian.
#'
#' @return List with par, value (log-likelihood), convergence, counts, message,
#'   hessian, vcov.
#' @noRd
.hzr_optim_generic <- function(
    logl_fn,
    gradient_fn,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    theta_start,
    weights = NULL,
    control = list(),
    use_bounds = FALSE,
    lower_bounds = NULL,
    hessian_fn = NULL) {

  control <- utils::modifyList(
    list(maxit = 1000, reltol = 1e-5, abstol = 1e-6),
    control
  )

  # Resolve lower bounds: default to 1e-6 on all params only when the caller
  # has not supplied explicit bounds.  Callers like .hzr_optim_weibull pass
  # c(1e-6, 1e-6, -Inf, ...) so that shape params are constrained but betas
  # are free.
  if (use_bounds && is.null(lower_bounds)) {
    lower_bounds <- rep(1e-6, length(theta_start))
  }

  # Negative log-likelihood for minimization
  objective <- function(theta) {
    if (!all(is.finite(theta))) return(1e10)
    # Only penalise infeasible shape params (those with finite lower bounds > -Inf).
    if (use_bounds && any(theta[is.finite(lower_bounds)] <= lower_bounds[is.finite(lower_bounds)])) return(1e10)

    ll <- -logl_fn(
      theta = theta, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x, weights = weights, return_gradient = FALSE
    )
    if (!is.finite(ll)) return(1e10)
    ll
  }

  # Negative score vector for minimization
  gradient <- function(theta) {
    if (!all(is.finite(theta))) return(rep(0, length(theta)))
    if (use_bounds && any(
      theta[is.finite(lower_bounds)] <= lower_bounds[is.finite(lower_bounds)]
    )) {
      return(rep(0, length(theta)))
    }

    grad <- tryCatch(
      gradient_fn(
        theta = theta, time = time, status = status,
        time_lower = time_lower, time_upper = time_upper,
        x = x, weights = weights
      ),
      error = function(e) rep(0, length(theta))
    )

    grad[!is.finite(grad)] <- 0
    -grad
  }

  # Dispatch to constrained or unconstrained optimizer
  if (use_bounds) {
    result <- optim(
      par = theta_start, fn = objective, gr = gradient,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = rep(Inf, length(theta_start)),
      control = list(
        maxit = control$maxit,
        factr = 1 / control$reltol,
        pgtol = control$abstol
      ),
      hessian = FALSE
    )
  } else {
    result <- optim(
      par = theta_start, fn = objective, gr = gradient,
      method = "BFGS",
      control = list(maxit = control$maxit, reltol = control$reltol),
      hessian = FALSE
    )
  }

  # Post-fit Hessian for standard errors.  Prefer the caller's analytic Hessian
  # (on the objective / negative-log-likelihood scale) when supplied; otherwise,
  # or when it declines by returning NULL, fall back to a numerical Hessian.
  hess_result <- NULL
  if (!is.null(hessian_fn)) {
    # A hook that ERRORS is a different thing from a hook that declines by
    # returning NULL, and swallowing the error made the two indistinguishable
    # downstream. Say which happened.
    hess_result <- tryCatch(
      hessian_fn(result$par),
      error = function(e) {
        warning("hessian_fn() errored, so the analytic Hessian was not used: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    # A non-NULL hook result must be a square matrix matching the parameter
    # dimension; a misbehaving hook should not silently produce a wrong vcov.
    p_dim <- length(result$par)
    if (!is.null(hess_result) &&
        !(is.matrix(hess_result) && all(dim(hess_result) == p_dim))) {
      warning("hessian_fn returned a non-conformant result; using numerical Hessian")
      hess_result <- NULL
    }
  }
  if (is.null(hess_result)) {
    # Reached when no analytic Hessian was supplied, or the hook declined by
    # returning NULL, or it returned something non-conformant, or it errored.
    # Each of those is reported at its own site above; this branch only knows
    # that it has no analytic Hessian, so the message says exactly that rather
    # than asserting a reason it cannot distinguish.
    #
    # numDeriv is the documented fallback, but it is a *Suggests*, so it is
    # legitimately absent on a machine installed without Suggests -- and then
    # there is no third option and no standard errors.
    #
    # That combination used to fail silently: rcond and pd came back NA and
    # vcov() returned a bare logical, with nothing anywhere naming numDeriv.
    # The user-visible symptom was diag(vcov(fit)) complaining about an
    # invalid 'nrow', which is unrecognisable from the cause.
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      warning("No analytic Hessian was obtained for this fit, and the ",
              "'numDeriv' fallback is not installed, so standard errors ",
              "cannot be computed. Install it with ",
              "install.packages(\"numDeriv\"); note numDeriv is a Suggests ",
              "dependency, so install.packages() and install_github() do ",
              "not pull it by default.", call. = FALSE)
    } else {
      hess_result <- tryCatch(
        numDeriv::hessian(objective, result$par),
        error = function(e) {
          warning("numDeriv::hessian() failed, so standard errors are ",
                  "unavailable: ", conditionMessage(e), call. = FALSE)
          NULL
        }
      )
    }
  }

  # Hardened inversion + conditioning diagnostics (Layer 1).
  inv <- if (is.matrix(hess_result)) {
    .hzr_safe_solve(hess_result)
  } else {
    # Reached only when no Hessian could be produced at all.
    # .hzr_safe_solve() warns on every degenerate path it handles, but this
    # branch bypasses it entirely, so it needs its own voice or the fit
    # returns NA diagnostics mutely.
    warning("No Hessian could be computed for this fit; standard errors, ",
            "rcond and pd are all NA.", call. = FALSE)
    list(vcov = NA, rcond = NA_real_, pd = NA)
  }

  list(
    par = result$par,
    value = -result$value,
    convergence = result$convergence,
    counts = result$counts,
    message = result$message,
    hessian = hess_result,
    vcov = inv$vcov,
    rcond = inv$rcond,
    pd = inv$pd
  )
}
