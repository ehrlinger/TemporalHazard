#' @importFrom stats optim
#' @keywords internal
NULL

# likelihood-exponential.R -- Exponential parametric hazard likelihood, gradient, and optimizer
#
# MODEL
# -----
# Constant-hazard (memoryless) proportional-hazards model:
#
#   h(t | x)  = lambda exp(eta)                    hazard (constant in t)
#   H(t | x)  = lambda t exp(eta)                  cumulative hazard
#   S(t | x)  = exp(-lambda t exp(eta))            survival
#
# where lambda > 0 (baseline rate), eta = x beta (linear predictor).
#
# THETA LAYOUT
# ------------
#   theta[1]   = log(lambda)   (unconstrained; lambda recovered via exp())
#   theta[2:p] = beta        (covariate coefficients; unrestricted)
#
# The exponential is a special case of Weibull with nu = 1.  It has no shape
# parameter so theta is one element shorter than Weibull for the same data.
# Unconstrained BFGS is used (no box bounds needed).
#
# FUNCTIONS
# ---------
#   .hzr_logl_exponential()     -- log-likelihood (optionally returning gradient)
#   .hzr_gradient_exponential() -- analytical score vector
#   .hzr_optim_exponential()    -- unconstrained BFGS wrapper

# Exponential Parametric Hazard Likelihood
#
# Evaluate the log-likelihood and its derivatives for exponential hazard models.
# The exponential is the simplest ALT model: constant hazard lambda, no shape parameter.

#' Log-likelihood for exponential hazard with covariates
#'
#' Computes the log-likelihood for right-censored data under the exponential
#' parametric hazard model with optional linear-predictor covariates.
#'
#' @param theta Vector of parameters:
#'   theta\[1\] = log(lambda) where lambda > 0 is the baseline hazard rate
#'   theta\[2:length\]: Covariate coefficients (linear on log-hazard scale)
#'
#' @param time Numeric vector of follow-up times (n)
#' @param status Numeric vector of event indicators: 1 = event, 0 = censored (n)
#' @param time_lower Optional numeric lower bound vector for interval-censored rows.
#'   Defaults to time if NULL.
#' @param time_upper Optional numeric upper bound vector for left/interval-censored rows.
#'   Defaults to time if NULL.
#' @param x Design matrix of covariates (n x p_coef); NULL for no covariates
#' @param return_gradient Logical; if TRUE, attach gradient vector as attribute
#'
#' @return Scalar log-likelihood value. If return_gradient = TRUE, gradient vector
#' is attached as attribute \code{"gradient"}.
#'
#' @details
#' The exponential hazard model is parameterized as:
#'
#' \deqn{h(t | x) = \lambda \exp(\eta)}
#'
#' where:
#' - \eqn{\lambda > 0} is baseline hazard rate
#' - \eqn{\eta = x \beta} is log-relative-hazard (covariates)
#'
#' The survival function is:
#'
#' \deqn{S(t | x) = \exp(-\lambda t \exp(\eta))}
#'
#' The log-likelihood for right-censored data with entry times \eqn{s_i} is:
#'
#' \deqn{\ell(\theta) = \sum_{i: \delta_i = 1} [\log\lambda + \eta_i]
#'   - \lambda \sum_i (t_i - s_i) \exp(\eta_i)}
#'
#' For status 0/1 rows, \eqn{s_i} is \code{time_lower} where
#' \code{time_lower < time} (counting-process entry, left truncation) and 0
#' otherwise, so each row contributes \eqn{H(t_i) - H(s_i)}. Status 2 rows
#' use \code{time_lower} as the interval lower bound instead.
#'
#' Reparameterization: theta\[1\] = log(lambda) avoids constrained optimization.
#'
#' Mixed censoring status coding:
#' - 1: exact event at time
#' - 0: right-censored at time
#' - -1: left-censored with upper bound time_upper \(or time\)
#' - 2: interval-censored in the interval \(time_lower, time_upper\)
#'
#' @noRd
.hzr_logl_exponential <- function(
    theta,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    weights = NULL,
    return_gradient = FALSE) {

  n <- length(time)

  # Default unit weights
  if (is.null(weights)) weights <- rep(1, n)

  # Extract parameters
  log_lambda <- theta[1]
  lambda <- exp(log_lambda)  # Always positive

  # Covariate coefficients (if any)
  if (!is.null(x)) {
    if (is.null(attr(x, "dimnames")[[2]])) {
      colnames(x) <- paste0("beta_", seq_len(ncol(x)))
    }
    beta <- theta[2:length(theta)]
    eta <- as.numeric(x %*% beta)  # n-vector
  } else {
    eta <- rep(0, n)
  }

  # Resolve censoring bounds once so each likelihood term can reference
  # a consistent lower/upper representation regardless of user input style.
  lower <- if (is.null(time_lower)) time else time_lower
  upper <- if (is.null(time_upper)) time else time_upper

  if (any(status %in% c(-1, 2)) && any(upper <= 0)) {
    return(Inf)
  }
  if (any(status == 2) && any(lower[status == 2] >= upper[status == 2])) {
    return(Inf)
  }

  # Cumulative hazards at event/lower/upper times; all terms can be
  # expressed using H(t) for the exponential model.
  cumhaz_event <- lambda * time * exp(eta)
  cumhaz_lower <- lambda * lower * exp(eta)
  cumhaz_upper <- lambda * upper * exp(eta)

  # Counting-process entry-time cumulative hazard H(start) for status 0/1
  # rows, as in .hzr_logl_weibull(): only genuine epoch rows (time_lower <
  # time) get a non-zero start.  Status 2 rows keep start = 0 because their
  # time_lower is the interval lower bound, used via cumhaz_lower.
  start_vec <- .hzr_entry_start_exponential(time, status, time_lower)
  cumhaz_start <- lambda * start_vec * exp(eta)

  idx_event <- status == 1
  idx_right <- status == 0
  idx_left <- status == -1
  idx_interval <- status == 2

  # Event: w * [log h(t) - (H(t) - H(start))]
  ll_event <- if (any(idx_event)) {
    sum(weights[idx_event] *
          ((log_lambda + eta[idx_event]) -
             (cumhaz_event[idx_event] - cumhaz_start[idx_event])))
  } else {
    0
  }

  # Right-censored: w * [-(H(t) - H(start))]
  ll_right <- if (any(idx_right)) {
    -sum(weights[idx_right] *
           (cumhaz_event[idx_right] - cumhaz_start[idx_right]))
  } else {
    0
  }

  # Left-censored: w * log(1 - exp(-H(u)))
  ll_left <- if (any(idx_left)) {
    sum(weights[idx_left] * hzr_log1mexp(cumhaz_upper[idx_left]))
  } else {
    0
  }

  # Interval-censored: w * [log(S(l)-S(u))]
  ll_interval <- if (any(idx_interval)) {
    delta_h <- cumhaz_upper[idx_interval] - cumhaz_lower[idx_interval]
    sum(weights[idx_interval] *
          (-cumhaz_lower[idx_interval] + hzr_log1mexp(delta_h)))
  } else {
    0
  }

  logl <- ll_event + ll_right + ll_left + ll_interval

  if (!is.finite(logl)) {
    return(Inf)
  }

  # If gradient requested, compute score vector
  if (return_gradient) {
    grad <- .hzr_gradient_exponential(
      theta = theta, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper, x = x,
      weights = weights,
      eta = eta, cumhaz = cumhaz_event,
      lambda = lambda, log_lambda = log_lambda
    )
    attr(logl, "gradient") <- grad
  }

  logl
}

#' Score vector (gradient of log-likelihood) for exponential
#'
#' Computes the score vector of the exponential log-likelihood w.r.t. all parameters.
#'
#' The log-likelihood is:
#'   \eqn{L = \sum(\delta_i * [\log(\lambda) + \eta_i]) - \sum(\lambda * (t_i - s_i) * \exp(\eta_i))}
#'
#' where \eqn{s_i} is the entry time (0 unless status 0/1 and time_lower < time).
#' With \eqn{H_i = \lambda (t_i - s_i) \exp(\eta_i)}, the derivatives are:
#' dL/d(log lambda) = sum(delta_i) - sum(H_i)
#' dL/dbeta_j    = t(X) %*% (delta - H)
#'
#' @noRd
.hzr_gradient_exponential <- function(
    theta,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    weights = NULL,
    eta = NULL,
    cumhaz = NULL,
    lambda = NULL,
    log_lambda = NULL) {

  n <- length(time)
  p <- length(theta)
  grad <- numeric(p)

  # Default unit weights -- keep gradient consistent with the LL when weights
  # are omitted.
  if (is.null(weights)) weights <- rep(1, n)

  # Numerical fallback for left/interval censoring contributions.
  # Right-censored/event-only rows still use the closed-form score below.
  if (any(status %in% c(-1, 2))) {
    return(.hzr_numeric_grad_exponential(theta, time, status, time_lower,
                                          time_upper, x, weights = weights))
  }

  # Recompute if not provided
  if (is.null(eta) || is.null(cumhaz) || is.null(lambda)) {
    log_lambda <- theta[1]
    lambda <- exp(log_lambda)

    if (!is.null(x)) {
      beta <- theta[2:length(theta)]
      eta <- as.numeric(x %*% beta)
    } else {
      eta <- rep(0, n)
    }

    cumhaz <- lambda * time * exp(eta)
  }

  # Entry-time H(start); zero except for genuine epoch rows (status 0/1 with
  # time_lower < time).  H is linear in lambda * exp(eta), so the entry term
  # enters every derivative as H(t) -> H(t) - H(start).
  cumhaz_start <- lambda * .hzr_entry_start_exponential(time, status,
                                                        time_lower) * exp(eta)

  # Weighted building blocks: every term below is the unweighted form with
  # `status` -> `w * status` and `cumhaz` -> `w * (H(t) - H(start))`.
  w_status <- weights * status
  w_cumhaz <- weights * (cumhaz - cumhaz_start)

  # ===== Gradient w.r.t. log(lambda) =====
  # dL/d(log lambda) = sum(w * delta) - sum(w * (H(t) - H(start)))
  grad[1] <- sum(w_status) - sum(w_cumhaz)

  # ===== Gradient w.r.t. beta (covariate coefficients) =====
  if (p > 1 && !is.null(x)) {
    # dL/dbeta = t(X) %*% (w * delta - w * H)
    residual <- w_status - w_cumhaz
    grad[2:p] <- as.numeric(crossprod(x, residual))
  }

  grad
}

#' Analytic Hessian of the exponential negative log-likelihood
#'
#' Returns the Hessian of the **objective** (negative log-likelihood) on the
#' internal parameter scale, matching \code{numDeriv::hessian(objective)} so it
#' can be used directly by \code{.hzr_optim_generic()}.  Coverage mirrors the
#' analytic gradient: closed form for event + right-censored rows only.  When any
#' row is left- or interval-censored (\code{status \%in\% c(-1, 2)}) it returns
#' \code{NULL} to request the numerical-Hessian fallback.
#'
#' Derivation: with \eqn{H_i = \lambda (t_i - s_i) e^{\eta_i}} (\eqn{s_i} the
#' entry time, 0 unless status 0/1 and time_lower < time) and design matrix
#' \eqn{\tilde X = [1 | X]}, the log-likelihood Hessian is
#' \eqn{-\tilde X^\top \mathrm{diag}(w H) \tilde X}; the objective Hessian is its
#' negation, \eqn{+\tilde X^\top \mathrm{diag}(w H) \tilde X}.
#'
#' @noRd
.hzr_hessian_exponential <- function(
    theta, time, status,
    time_lower = NULL, time_upper = NULL, x = NULL, weights = NULL) {

  n <- length(time)
  if (is.null(weights)) weights <- rep(1, n)

  # Mirror the analytic gradient's coverage: closed form is event + right
  # censored only.  Decline (NULL -> numDeriv fallback) on left/interval rows.
  if (any(status %in% c(-1, 2))) {
    return(NULL)
  }

  log_lambda <- theta[1]
  lambda <- exp(log_lambda)
  # Fail loud on inconsistent theta/x: a mismatch would otherwise yield a wrong
  # or non-conformant Hessian used for vcov without a clear error.
  p_cov <- length(theta) - 1L
  if (!is.null(x)) x <- as.matrix(x)
  n_x <- if (is.null(x)) 0L else NCOL(x)
  if (n_x != p_cov) {
    stop(".hzr_hessian_exponential(): ncol(x) (", n_x, ") must equal the number ",
         "of covariate parameters in theta (", p_cov, ").", call. = FALSE)
  }
  if (p_cov > 0L) {
    beta <- theta[2:length(theta)]
    eta <- as.numeric(x %*% beta)
  } else {
    eta <- rep(0, n)
  }

  # Weighted cumulative hazard per row, net of the entry-time H(start).
  start_vec <- .hzr_entry_start_exponential(time, status, time_lower)
  wH <- weights * lambda * (time - start_vec) * exp(eta)

  # Design matrix augmented with the log(lambda) "intercept" column.
  x_tilde <- if (is.null(x)) matrix(1, nrow = n, ncol = 1L) else cbind(1, x)

  # Objective (negative log-likelihood) Hessian: + X~' diag(wH) X~.
  hess <- crossprod(x_tilde, wH * x_tilde)
  dimnames(hess) <- list(names(theta), names(theta))
  hess
}

#' @noRd
.hzr_optim_exponential <- function(
    time, status, time_lower = NULL, time_upper = NULL,
    x = NULL, theta_start, weights = NULL, control = list()) {
  hessian_fn <- function(par) {
    .hzr_hessian_exponential(
      theta = par, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x, weights = weights
    )
  }
  .hzr_optim_generic(
    logl_fn = .hzr_logl_exponential,
    gradient_fn = .hzr_gradient_exponential,
    time = time, status = status,
    time_lower = time_lower, time_upper = time_upper,
    x = x, theta_start = theta_start, weights = weights,
    control = control, use_bounds = FALSE,
    hessian_fn = hessian_fn
  )
}

# Counting-process entry time per row: time_lower for status 0/1 rows with
# time_lower < time, 0 otherwise (status -1/2 rows, or time_lower NULL).
# Same rule as .hzr_logl_weibull().
.hzr_entry_start_exponential <- function(time, status, time_lower) {
  start_vec <- rep(0, length(time))
  if (!is.null(time_lower)) {
    epoch_idx <- status %in% c(0L, 1L) & time_lower < time
    start_vec[epoch_idx] <- time_lower[epoch_idx]
  }
  start_vec
}

.hzr_numeric_grad_exponential <- function(theta, time, status,
                                            time_lower = NULL, time_upper = NULL,
                                            x = NULL, weights = NULL) {
  # Reuse log-likelihood as scalar objective and differentiate numerically.
  objective <- function(par) {
    .hzr_logl_exponential(
      theta = par,
      time = time,
      status = status,
      time_lower = time_lower,
      time_upper = time_upper,
      x = x,
      weights = weights,
      return_gradient = FALSE
    )
  }

  if (requireNamespace("numDeriv", quietly = TRUE)) {
    return(as.numeric(numDeriv::grad(objective, theta)))
  }

  eps <- 1e-6
  grad <- numeric(length(theta))
  for (i in seq_along(theta)) {
    tp <- theta
    tm <- theta
    tp[i] <- tp[i] + eps
    tm[i] <- tm[i] - eps
    grad[i] <- (objective(tp) - objective(tm)) / (2 * eps)
  }
  grad
}
