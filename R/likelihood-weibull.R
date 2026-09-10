#' @importFrom stats optim setNames
#' @keywords internal
NULL

# likelihood-weibull.R -- Weibull parametric hazard likelihood, gradient, and optimizer
#
# MODEL
# -----
# Proportional-hazards (PH) parameterization:
#
#   H(t | x)  = (mu t)^nu exp(eta)              cumulative hazard
#   h(t | x)  = nu mu^nu t^(nu-1) exp(eta)      hazard = dH/dt = (nu/t) H
#   S(t | x)  = exp(-H(t | x))              survival
#
# where mu > 0 (scale), nu > 0 (shape), eta = x beta (linear predictor).
#
# THETA LAYOUT
# ------------
#   theta[1] = mu   (must be > 0; enforced by L-BFGS-B lower bound)
#   theta[2] = nu   (must be > 0; enforced by L-BFGS-B lower bound)
#   theta[3:p] = beta (covariate coefficients; unrestricted)
#
# Note: unlike the other distributions, Weibull optimises on the natural scale
# and uses L-BFGS-B with lower = 1e-6 instead of log-reparameterisation.
# This preserves direct interpretability of the mu/nu estimates.
#
# FUNCTIONS
# ---------
#   .hzr_logl_weibull()     -- log-likelihood (optionally returning gradient)
#   .hzr_gradient_weibull() -- analytical score vector
#   .hzr_optim_weibull()    -- L-BFGS-B wrapper with post-fit Hessian SE

# Weibull Parametric Hazard Likelihood
#
# Evaluate the log-likelihood and its derivatives for Weibull hazard models.
# This is the core mathematical engine for model fitting.

#' Log-likelihood for Weibull hazard with covariates
#'
#' Computes the negative log-likelihood for a sample under the Weibull
#' parametric hazard model with optional linear-predictor covariates.
#'
#' @param theta Vector of parameters:
#'   theta\[1:n_shape\]: Shape/scale parameters (distribution-specific)
#'   theta\[n_shape+1:length\]: Covariate coefficients (linear on log-hazard scale)
#'   where n_shape = number of shape parameters (e.g., 2 for Weibull scale/shape)
#'
#' @param time Numeric vector of follow-up times (n)
#' @param status Numeric vector of event indicators: 1 = event, 0 = censored (n)
#' @param time_lower Optional numeric lower bound vector for interval-censored rows.
#'   Defaults to time if NULL.
#' @param time_upper Optional numeric upper bound vector for left/interval-censored rows.
#'   Defaults to time if NULL.
#' @param x Design matrix of covariates (n x p_coef); NULL for no covariates
#' @param dist_name Character name of baseline distribution ("weibull", etc.)
#' @param return_gradient Logical; if TRUE, attach gradient vector as attribute
#' @param return_hessian Logical; if TRUE, attach Hessian matrix as attribute
#'
#' @return Scalar log-likelihood value. If return_gradient = TRUE, gradient vector
#' is attached as attribute \code{"gradient"}. If return_hessian = TRUE, Hessian
#' matrix is attached as \code{"hessian"}.
#'
#' @details
#' The Weibull hazard model is parameterized via its cumulative hazard
#' \eqn{H(t | x) = (\mu t)^{\nu} \exp(\eta)}, so the hazard is its exact
#' derivative:
#'
#' \deqn{h(t | x) = \nu \mu^{\nu} t^{\nu-1} \exp(\eta) = (\nu / t)\, H(t | x)}
#'
#' where:
#' - \eqn{\mu > 0} is scale (MU)
#' - \eqn{\nu > 0} is shape (NU)
#' - \eqn{\eta = x \beta} is log-relative-hazard (covariates)
#'
#' The log-likelihood for right-censored data is:
#'
#' \deqn{\ell(\theta) = \sum_{i: \delta_i = 1} \log h(t_i | x_i)
#'   - \sum_i S(t_i | x_i)}
#'
#' where S is the survival function (1 - CDF).
#'
#' Mixed censoring status coding:
#' - 1: exact event at time
#' - 0: right-censored at time
#' - -1: left-censored with upper bound time_upper \(or time\)
#' - 2: interval-censored in the interval \(time_lower, time_upper\)
#'
#' @noRd
.hzr_logl_weibull <- function(
    theta,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    weights = NULL,
    dist_name = "weibull",
    return_gradient = FALSE,
    return_hessian = FALSE) {

  # Shape parameter count for Weibull
  n_shape <- 2  # mu, nu (scale, shape)
  n <- length(time)

  # Default unit weights
  if (is.null(weights)) weights <- rep(1, n)

  # Extract parameters
  mu <- theta[1]    # scale > 0
  nu <- theta[2]    # shape > 0

  # Covariate coefficients (if any)
  if (!is.null(x)) {
    if (is.null(attr(x, "dimnames")[[2]])) {
      colnames(x) <- paste0("beta_", seq_len(ncol(x)))
    }
    beta <- theta[3:length(theta)]
    eta <- as.numeric(x %*% beta)  # n-vector
  } else {
    eta <- rep(0, n)
  }

  # Stability checks
  if (mu <= 0 || nu <= 0) {
    return(Inf)  # Infeasible: return large penalty
  }

  # Normalize censoring bounds to lower/upper vectors for unified formulas.
  lower <- if (is.null(time_lower)) time else time_lower
  upper <- if (is.null(time_upper)) time else time_upper

  # Exact events require strictly positive time because log(h(t)) is used.
  # Right-censored observations can be at time 0 (valid: log S(0) = 0).
  if (any(status == 1 & time <= 0)) {
    return(Inf)
  }
  if (any(status %in% c(-1, 2)) && any(upper <= 0)) {
    return(Inf)
  }
  if (any(status == 2) && any(lower[status == 2] >= upper[status == 2])) {
    return(Inf)
  }

  # Event-time hazard/cumulative-hazard (exact events only).
  # Hazard is the exact derivative of the cumulative hazard H = (mu t)^nu e^eta:
  #   h = dH/dt = nu * mu^nu * t^(nu-1) * e^eta = (nu / t) * H
  # (Form A, matching the C/SAS HAZARD reference where HF = MU * dG/dt.)
  haz_event <- nu * mu ^ nu * (time ^ (nu - 1)) * exp(eta)
  cumhaz_event <- (mu * time) ^ nu * exp(eta)

  # Lower/upper cumulative hazards for censoring contributions.
  cumhaz_lower <- (mu * lower) ^ nu * exp(eta)
  cumhaz_upper <- (mu * upper) ^ nu * exp(eta)

  # Counting-process entry-time cumulative hazard (H(start)) for
  # event/right-censored rows.  Only rows where time_lower < time and
  # status is 0 or 1 (genuine epoch entries) get a non-zero start.
  # Interval-censored rows (status == 2) always use start = 0 because
  # their contribution goes through cumhaz_lower/cumhaz_upper, not H(start).
  # This prevents time_lower being mistakenly interpreted as an epoch start
  # when it was supplied as the interval lower bound for status-2 rows.
  start_vec <- rep(0, n)
  if (!is.null(time_lower)) {
    epoch_idx <- status %in% c(0L, 1L) & time_lower < time
    start_vec[epoch_idx] <- time_lower[epoch_idx]
  }
  cumhaz_start <- (mu * start_vec) ^ nu * exp(eta)

  idx_event <- status == 1
  idx_right <- status == 0
  idx_left <- status == -1
  idx_interval <- status == 2

  # Exact events: w * [log h(stop) - (H(stop) - H(start))]
  ll_event <- if (any(idx_event)) {
    sum(weights[idx_event] *
          (log(haz_event[idx_event]) -
             (cumhaz_event[idx_event] - cumhaz_start[idx_event])))
  } else {
    0
  }

  # Right-censored: w * [-(H(stop) - H(start))]
  ll_right <- if (any(idx_right)) {
    -sum(weights[idx_right] *
           (cumhaz_event[idx_right] - cumhaz_start[idx_right]))
  } else {
    0
  }

  # Left-censored: w * [log(1 - exp(-H(u)))]
  ll_left <- if (any(idx_left)) {
    sum(weights[idx_left] * hzr_log1mexp(cumhaz_upper[idx_left]))
  } else {
    0
  }

  # Interval-censored: w * [-H(l) + log(1 - exp(-(H(u)-H(l))))]
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

  # If gradient requested, compute score vector.
  if (return_gradient) {
    grad <- .hzr_gradient_weibull(
      theta = theta, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper, x = x,
      weights = weights,
      eta = eta, cumhaz = cumhaz_event, haz = haz_event,
      mu = mu, nu = nu, n_shape = n_shape
    )
    attr(logl, "gradient") <- grad
  }

  # If Hessian requested (stub for M2)
  if (return_hessian) {
    attr(logl, "hessian") <- matrix(NA, length(theta), length(theta))
  }

  logl
}

#' Score vector (gradient of log-likelihood)
#'
#' Computes the score vector (gradient) of the Weibull log-likelihood w.r.t. all parameters.
#'
#' @param time_lower Optional lower bounds for interval-censored rows.
#' @param time_upper Optional upper bounds for left/interval-censored rows.
#'
#' The log-likelihood is:
#'   L = sum(delta_i * log h(t_i)) - sum(H(t_i))
#'
#' where h(t) = nu * mu^nu * t^(nu-1) * exp(eta) and H(t) = (mu*t)^nu * exp(eta)
#' (so log h = log(nu) + nu*log(mu) + (nu-1)*log(t) + eta). Derivatives are:
#'
#' dL/dmu  = (nu / mu) * sum(delta_i) - (nu / mu) * sum(H(t_i))
#' dL/dnu  = sum(delta_i / nu) + sum(delta_i) * log(mu) + sum(delta_i * log(t_i))
#'           - sum(log(mu * t_i) * H(t_i))
#' dL/dbeta_j = sum(delta_i * x_ij) - sum(H(t_i) * x_ij)  = t(X) %*% (delta - H)
#'
#' @noRd
.hzr_gradient_weibull <- function(
    theta,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    weights = NULL,
    eta = NULL,
    cumhaz = NULL,
    haz = NULL,
    mu = NULL,
    nu = NULL,
    n_shape = 2) {

  n <- length(time)
  p <- length(theta)
  grad <- numeric(p)

  # Default unit weights -- keep gradient consistent with the LL when weights
  # are omitted.
  if (is.null(weights)) weights <- rep(1, n)

  # Interval/left censoring currently uses a robust numerical gradient fallback.
  # The closed-form derivatives for these terms are deferred to a later pass.
  if (any(status %in% c(-1, 2))) {
    return(.hzr_numeric_grad_weibull(theta, time, status, time_lower,
                                       time_upper, x, weights = weights))
  }

  # Sanity check: if eta, cumhaz, haz not provided, compute them
  if (is.null(eta) || is.null(cumhaz) || is.null(haz)) {
    # Need to recompute from theta
    mu <- theta[1]
    nu <- theta[2]

    if (!is.null(x)) {
      beta <- theta[3:length(theta)]
      eta <- as.numeric(x %*% beta)
    } else {
      eta <- rep(0, n)
    }

    haz <- nu * mu ^ nu * (time ^ (nu - 1)) * exp(eta)
    cumhaz <- (mu * time) ^ nu * exp(eta)
  }

  # Counting-process entry-time cumulative hazard; H(start) = 0 except
  # for genuine epoch rows (status 0/1 with time_lower < time).
  start_vec <- rep(0, n)
  if (!is.null(time_lower)) {
    epoch_idx <- status %in% c(0L, 1L) & time_lower < time
    start_vec[epoch_idx] <- time_lower[epoch_idx]
  }
  cumhaz_start <- (mu * start_vec) ^ nu * exp(eta)

  # Weighted building blocks: each per-row contribution below is the
  # unweighted form scaled by `weights`, and the cumulative-hazard term is
  # the epoch difference `H(stop) - H(start)` (degenerates to `H(stop)` when
  # start = 0).
  w_status <- weights * status
  w_cumhaz_net <- weights * (cumhaz - cumhaz_start)

  # ===== Gradient w.r.t. mu (scale parameter) =====
  # log h = log(nu) + nu*log(mu) + (nu-1)*log(t) + eta, so the event term
  # contributes d log h/dmu = nu/mu.  dH(t)/dmu = (nu / mu) * H(t), so the
  # cumulative term subtracts -(nu/mu) * (H(stop) - H(start)).
  grad[1] <- (nu / mu) * (sum(w_status) - sum(w_cumhaz_net))

  # ===== Gradient w.r.t. nu (shape parameter) =====
  # dH(t)/dnu = log(mu*t) * H(t).  When start = 0, H(start) = 0 and the
  # log(mu*start) term is undefined; `ifelse` picks the well-defined 0
  # limit there.
  # log h = log(nu) + nu*log(mu) + (nu-1)*log(t) + eta, so the event term
  # contributes d log h/dnu = 1/nu + log(mu) + log(t); the log(mu) piece is
  # constant across rows.
  # Guard log(time): a legal right-censored time = 0 row has w_status = 0 and
  # cumhaz = 0, but an unguarded log(0) = -Inf turns those into 0 * -Inf = NaN,
  # poisoning the entire summed nu (shape) gradient component (start is already
  # guarded by the ifelse below).
  log_t <- log(pmax(time, .Machine$double.xmin))
  log_mu_start <- ifelse(start_vec > 0, log(mu * start_vec), 0)
  d_nu_start <- weights * log_mu_start * cumhaz_start
  grad[2] <- sum(w_status) / nu +
             sum(w_status) * log(mu) +
             sum(w_status * log_t) -
             (sum((log(mu) + log_t) * weights * cumhaz) - sum(d_nu_start))

  # ===== Gradient w.r.t. beta (covariate coefficients) =====
  if (n_shape < p && !is.null(x)) {
    # dL/dbeta = t(X) %*% (w * (delta - (H(stop) - H(start))))
    residual <- w_status - w_cumhaz_net
    grad[3:p] <- as.numeric(crossprod(x, residual))
  }

  grad
}

#' Analytic Hessian of the Weibull negative log-likelihood (internal scale)
#'
#' Returns the Hessian of the objective (negative log-likelihood) on the internal
#' \code{phi = (alpha, psi, beta)} parameterization that \code{.hzr_optim_weibull}
#' optimizes, where \code{g = exp(psi)} and \code{H(t) = exp(alpha + eta) t^g}.
#' This matches \code{numDeriv::hessian(objective)} so it can be used directly by
#' \code{.hzr_optim_generic()}; the downstream delta-method transform to the
#' natural \code{(mu, nu, beta)} scale is applied by the caller.  Coverage mirrors
#' the analytic gradient: event + right-censored rows only; returns \code{NULL}
#' for any left/interval-censored row (\code{status \%in\% c(-1, 2)}) to request
#' the numerical-Hessian fallback.
#'
#' @noRd
.hzr_hessian_weibull_internal <- function(
    theta, time, status,
    time_lower = NULL, time_upper = NULL, x = NULL, weights = NULL) {

  n <- length(time)
  if (is.null(weights)) weights <- rep(1, n)

  # Coverage contract: closed form is event + right censored only.
  if (any(status %in% c(-1, 2))) {
    return(NULL)
  }

  alpha <- theta[1]
  psi <- theta[2]
  g <- exp(psi)
  p <- length(theta)
  p_cov <- p - 2L
  # Fail loud on inconsistent theta/x: a silent has_cov = FALSE would zero the
  # covariate rows/cols and yield a conformant-but-wrong vcov.
  if (!is.null(x)) x <- as.matrix(x)
  n_x <- if (is.null(x)) 0L else NCOL(x)
  if (n_x != p_cov) {
    stop(".hzr_hessian_weibull_internal(): ncol(x) (", n_x, ") must equal the ",
         "number of covariate parameters in theta (", p_cov, ").", call. = FALSE)
  }
  has_cov <- p_cov > 0L
  beta <- if (has_cov) theta[3:p] else numeric(0)
  eta <- if (has_cov) as.numeric(x %*% beta) else rep(0, n)

  # Cumulative hazards and counting-process entry-time (epochs only).
  base <- exp(alpha + eta)
  cumhaz <- base * (time ^ g)
  start_vec <- rep(0, n)
  if (!is.null(time_lower)) {
    epoch_idx <- status %in% c(0L, 1L) & time_lower < time
    start_vec[epoch_idx] <- time_lower[epoch_idx]
  }
  cumhaz_start <- base * (start_vec ^ g)

  log_t <- log(pmax(time, .Machine$double.xmin))
  log_ts <- ifelse(start_vec > 0, log(start_vec), 0)

  w <- weights
  delta <- as.numeric(status == 1)
  net <- cumhaz - cumhaz_start                       # A_i

  # Weighted aggregate building blocks.
  wA   <- w * net
  wD   <- w * (cumhaz * log_t - cumhaz_start * log_ts)
  # psi,psi scalar.
  hpp <- -g * (sum(w * delta * log_t) - sum(w * cumhaz * log_t) +
                 sum(w * cumhaz_start * log_ts)) +
          g^2 * (sum(w * cumhaz * log_t^2) - sum(w * cumhaz_start * log_ts^2))

  # Design augmented with the alpha intercept column.
  x_tilde <- if (has_cov) cbind(1, x) else matrix(1, nrow = n, ncol = 1L)

  hq <- crossprod(x_tilde, wA * x_tilde)             # (alpha,beta) block
  v  <- g * as.numeric(crossprod(x_tilde, wD))        # psi cross (alpha,beta)

  # Assemble p x p with order (alpha = 1, psi = 2, beta = 3..p).
  idx_q <- if (has_cov) c(1L, 3:p) else 1L
  hess <- matrix(0, p, p)
  hess[idx_q, idx_q] <- hq
  hess[2L, idx_q] <- v
  hess[idx_q, 2L] <- v
  hess[2L, 2L] <- hpp
  dimnames(hess) <- list(names(theta), names(theta))
  hess
}

#' @noRd
.hzr_optim_weibull <- function(
    time, status, time_lower = NULL, time_upper = NULL,
    x = NULL, theta_start, weights = NULL, control = list()) {

  # -- Internal reparameterisation --------------------------------------------
  # The original hazard C code parameterises the cumulative hazard as
  #
  #     H(t | x) = exp(alpha + Xbeta) * t^gamma
  #
  # where alpha is a free intercept, gamma = exp(psi) is the shape, and beta enters
  # the same log-hazard-rate space as alpha.  This decouples scale from shape
  # and avoids the degenerate ridge that appears when covariates with large
  # means (e.g. age ~ 62) interact with the (log mu, log nu) parameterisation.
  #
  # We optimise over phi = (alpha, psi, beta_1, ...) -- all unconstrained -- then
  # back-transform to the user-facing (mu, nu, beta) scale:
  #
  #     nu = exp(psi)      mu = exp(alpha / nu)
  #
  # and apply the delta method for the variance-covariance matrix.
  # --------------------------------------------------------------------------

  n_shape <- 2L
  p <- length(theta_start)
  n <- length(time)

  # -- Convert user theta (mu, nu, beta) -> internal (alpha, psi, beta) ------
  mu_start    <- unname(theta_start[1])
  nu_start    <- unname(theta_start[2])
  alpha_start <- nu_start * log(mu_start)   # alpha = nu * log(mu)
  psi_start   <- log(nu_start)              # psi = log(nu)
  beta_start  <- if (p > n_shape) unname(theta_start[(n_shape + 1):p]) else numeric(0)
  phi_start   <- c(alpha_start, psi_start, beta_start)

  # -- Internal likelihood: H(t|x) = exp(alpha + Xbeta) * t^exp(psi) ----------------
  # `weights` is declared as a formal (with a default of NULL) so that when
  # .hzr_optim_generic forwards weights via ..., it matches this formal
  # rather than landing in ... -- which would collide with the explicit
  # `weights =` on the delegated call below.
  .outer_w <- if (is.null(weights)) rep(1, n) else weights

  logl_internal <- function(theta, time, status, time_lower, time_upper,
                            x, weights = NULL, ...) {
    w_use <- if (is.null(weights)) .outer_w else weights
    alpha <- theta[1]
    g     <- exp(theta[2])   # nu = shape
    beta  <- if (p > n_shape) theta[(n_shape + 1):p] else numeric(0)
    eta   <- if (!is.null(x) && length(beta) > 0) {
      as.numeric(x %*% beta)
    } else {
      rep(0, n)
    }

    # Left / interval censoring: delegate to full Weibull likelihood via
    # the natural-scale function, converting params on the fly.
    if (any(status %in% c(-1, 2))) {
      mu <- exp(alpha / g)
      return(.hzr_logl_weibull(c(mu, g, beta), time, status,
                               time_lower, time_upper, x,
                               weights = w_use, ...))
    }

    log_t <- log(pmax(time, .Machine$double.xmin))
    H     <- exp(alpha + eta) * (time ^ g)

    # Counting-process entry-time H(start); zero except for genuine epochs.
    start_vec <- rep(0, n)
    if (!is.null(time_lower)) {
      epoch_idx <- status %in% c(0L, 1L) & time_lower < time
      start_vec[epoch_idx] <- time_lower[epoch_idx]
    }
    H_start   <- exp(alpha + eta) * (start_vec ^ g)

    idx_event <- status == 1
    idx_right <- status == 0

    ll_event <- if (any(idx_event)) {
      sum(w_use[idx_event] * (alpha + eta[idx_event] + theta[2] +
          (g - 1) * log_t[idx_event] -
          (H[idx_event] - H_start[idx_event])))
    } else {
      0
    }

    ll_right <- if (any(idx_right)) {
      -sum(w_use[idx_right] * (H[idx_right] - H_start[idx_right]))
    } else {
      0
    }

    logl <- ll_event + ll_right
    if (!is.finite(logl)) return(Inf)
    logl
  }

  # -- Internal gradient -----------------------------------------------------
  # `weights` is captured from the enclosing optim frame via `.outer_w`
  # (same trick used by logl_internal above) so the analytic gradient stays
  # consistent with the weighted LL whether or not .hzr_optim_generic
  # forwards the argument via ...
  grad_internal <- function(theta, time, status, time_lower, time_upper, x,
                            weights = NULL, ...) {
    w_use <- if (is.null(weights)) .outer_w else weights
    alpha <- theta[1]
    psi   <- theta[2]
    g     <- exp(psi)
    beta  <- if (p > n_shape) theta[(n_shape + 1):p] else numeric(0)
    eta   <- if (!is.null(x) && length(beta) > 0) as.numeric(x %*% beta) else rep(0, n)

    # Left / interval censoring -> numerical gradient fallback (already weighted
    # because logl_internal propagates weights via the closure).
    if (any(status %in% c(-1, 2))) {
      obj <- function(th) {
        logl_internal(th, time, status, time_lower, time_upper, x,
                      weights = w_use)
      }
      if (requireNamespace("numDeriv", quietly = TRUE))
        return(as.numeric(numDeriv::grad(obj, theta)))
      eps <- 1e-6
      grd <- numeric(p)
      for (i in seq_along(theta)) {
        tp <- tm <- theta
        tp[i] <- tp[i] + eps
        tm[i] <- tm[i] - eps
        grd[i] <- (obj(tp) - obj(tm)) / (2 * eps)
      }
      return(grd)
    }

    log_t <- log(pmax(time, .Machine$double.xmin))
    H     <- exp(alpha + eta) * (time ^ g)

    # H(start) for genuine epoch rows; zero otherwise.
    # `log_t_start` is only used where start > 0 (guarded below).
    start_vec <- rep(0, n)
    if (!is.null(time_lower)) {
      epoch_idx <- status %in% c(0L, 1L) & time_lower < time
      start_vec[epoch_idx] <- time_lower[epoch_idx]
    }
    H_start   <- exp(alpha + eta) * (start_vec ^ g)
    log_t_start <- ifelse(start_vec > 0, log(pmax(start_vec, .Machine$double.xmin)), 0)

    grad  <- numeric(p)

    # Weighted building blocks.  Every term below is the unweighted form with
    # `status` -> `w * status` and `H` -> `w * (H(stop) - H(start))`.
    w_status   <- w_use * status
    w_H_net    <- w_use * (H - H_start)
    w_H        <- w_use * H
    w_H_start  <- w_use * H_start

    # dL/dalpha = sum(w * delta) - sum(w * (H(stop) - H(start)))
    grad[1] <- sum(w_status) - sum(w_H_net)

    # dL/dpsi = sum(w * delta) + g * sum(w * delta * log(stop))
    #          - g * [sum(w * H(stop) * log(stop)) - sum(w * H(start) * log(start))]
    grad[2] <- sum(w_status) +
               g * sum(w_status * log_t) -
               g * (sum(w_H * log_t) - sum(w_H_start * log_t_start))

    # dL/dbeta_j = X_j' (w * (delta - (H(stop) - H(start))))
    if (p > n_shape && !is.null(x)) {
      grad[(n_shape + 1):p] <- as.numeric(crossprod(x, w_status - w_H_net))
    }

    grad
  }

  # -- Analytic Hessian closure (internal scale) --------------------------------
  hessian_fn <- function(par) {
    .hzr_hessian_weibull_internal(
      theta = par, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x, weights = weights
    )
  }

  # -- Optimise (all unconstrained) ------------------------------------------
  result <- .hzr_optim_generic(
    logl_fn     = logl_internal,
    gradient_fn = grad_internal,
    time        = time, status = status,
    time_lower  = time_lower, time_upper = time_upper,
    x           = x,
    theta_start = phi_start,
    weights     = weights,
    control     = control,
    use_bounds  = FALSE,
    hessian_fn  = hessian_fn
  )

  # -- Back-transform to natural scale (mu, nu, beta) ----------------------------
  alpha_hat <- result$par[1]
  psi_hat   <- result$par[2]
  nu_hat    <- exp(psi_hat)
  mu_hat    <- exp(alpha_hat / nu_hat)
  beta_hat  <- if (p > n_shape) result$par[(n_shape + 1):p] else numeric(0)

  result$par <- setNames(c(mu_hat, nu_hat, beta_hat), names(theta_start))

  # Delta method: Cov(mu, nu, beta) = J * Cov(alpha, psi, beta) * J^T
  #
  #   dmu/dalpha = mu / nu             dmu/dpsi = -mu alpha / nu
  #   dnu/dalpha = 0                 dnu/dpsi = nu
  #   dbeta/dalpha = 0                 dbeta/dpsi = 0         dbeta/dbeta = I
  #
  # .hzr_safe_solve() masks a non-positive variance by setting its row and
  # column to NA. Transform the finite block and carry the mask through J: a
  # reported parameter is unavailable when it depends on any masked internal
  # one. Skipping the transform instead would leave the surviving entries on
  # the (alpha, psi) scale beside (mu, nu) estimates.
  if (is.matrix(result$vcov)) {
    J <- diag(p)
    J[1, 1] <- mu_hat / nu_hat
    J[1, 2] <- -mu_hat * alpha_hat / nu_hat
    J[2, 1] <- 0
    J[2, 2] <- nu_hat
    bad_int <- !is.finite(diag(result$vcov))
    v_int <- result$vcov
    v_int[bad_int, ] <- 0
    v_int[, bad_int] <- 0
    v_nat <- J %*% v_int %*% t(J)
    bad_nat <- as.vector((J != 0) %*% bad_int) > 0
    v_nat[bad_nat, ] <- NA_real_
    v_nat[, bad_nat] <- NA_real_
    result$vcov <- v_nat
  }

  result
}

.hzr_numeric_grad_weibull <- function(theta, time, status, time_lower = NULL,
                                       time_upper = NULL, x = NULL,
                                       weights = NULL) {
  objective <- function(par) {
    .hzr_logl_weibull(
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

  # Fallback central differences if numDeriv is unavailable.
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
