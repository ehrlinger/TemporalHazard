#' @keywords internal
NULL

# predict-cl.R -- Delta-method confidence limits for predict.hazard()
#
# PROBLEM
# -------
# predict.hazard() returns a point estimate eta(theta) for each row.  We want
# a standard error and a (1 - alpha) confidence interval around each one,
# matching SAS HAZARD's `hzp_calc_haz_CL.c` / `hzp_calc_srv_CL.c`.
#
# METHOD
# ------
# Delta method: if theta_hat ~ N(theta, V) then any smooth function eta(theta)
# is approximately Gaussian with
#
#   var(eta_hat)  ~  (d eta/d theta)^T  V  (d eta/d theta).
#
# For a batch of n prediction points, we assemble J (n x p) where
# J[i, ] = d eta_i / d theta, then
#
#   se_i  =  sqrt( (J V J^T)_{ii} )  =  sqrt( rowSums( (J V) * J ) ).
#
# SCALE CHOICE
# ------------
# Symmetric CLs on the natural scale can leave a [0, 1] survival or a
# [0, Inf) hazard, so we build them on a transformed scale and back-transform:
#
#   hazard / cumhaz   -> log scale:
#     se(log eta)  =  se(eta) / eta,
#     [eta * exp(-z se_log),  eta * exp(+z se_log)]
#
#   survival          -> log-log scale (log(-log S) = log H):
#     se(log(-log S))  =  se(H) / H,
#     S_lwr  =  exp(-H * exp(+z se_logH)),
#     S_upr  =  exp(-H * exp(-z se_logH))
#
#   linear_predictor  -> natural scale (symmetric):
#     [eta - z se, eta + z se]
#
# ENTRY POINT
# -----------
# `.hzr_predict_with_se()` takes a fitted `hazard` object, a type, and the
# materialised `x` / `time` for prediction.  It dispatches to the appropriate
# distribution-specific jacobian builder, runs the delta-method sandwich, and
# returns a list(fit, se.fit, lower, upper).

# ---------------------------------------------------------------------------
# Jacobian: Weibull
# ---------------------------------------------------------------------------
#
# Parameter layout (natural scale, matching .hzr_logl_weibull):
#   theta = [mu, nu, beta_1, ..., beta_p]
#
# H(t|x) = (mu t)^nu exp(eta),   h_rel = exp(eta),   eta = x beta.
#
# Analytical derivatives used below:
#   dH/dmu       = (nu / mu) H
#   dH/dnu       = log(mu t) H
#   dH/dbeta_j   = x_ij H
#
#   d exp(eta)/dbeta_j = x_ij exp(eta)    (other entries are 0)
#
#   eta itself is linear: d eta / d beta_j = x_ij, rest 0.

#' Jacobian of Weibull predictions with respect to theta
#'
#' @param type Prediction type (see `.hzr_predict_with_se`).
#' @param theta MLE parameter vector `c(mu, nu, beta_1, ...)`.
#' @param time Prediction times (may be NULL for `hazard` / `linear_predictor`).
#' @param x Design matrix (n x p_cov) or NULL.
#' @param p Length of theta.
#' @return Numeric n x p Jacobian.
#' @keywords internal
.hzr_predict_jacobian_weibull <- function(type, theta, time, x, p) {
  n_shape <- 2L
  mu <- theta[1]
  nu <- theta[2]
  beta <- if (p > n_shape) theta[(n_shape + 1L):p] else numeric(0)

  if (!is.null(x) && length(beta) > 0L) {
    eta <- as.numeric(x %*% beta)
  } else if (!is.null(x)) {
    eta <- rep(0, nrow(x))
  } else {
    # No covariates: treat all rows as a single reference
    eta <- if (!is.null(time)) rep(0, length(time)) else 0
  }
  n <- length(eta)
  J <- matrix(0, nrow = n, ncol = p)

  if (type %in% c("cumulative_hazard", "survival")) {
    # On the log scale, as predict() computes H: mu * time can overflow where
    # H is finite (#566).
    log_mu_t <- log(mu) + log(time)
    H <- exp(nu * log_mu_t + eta)
    J[, 1L] <- (nu / mu) * H
    J[, 2L] <- log_mu_t * H
    if (length(beta) > 0L && !is.null(x)) {
      J[, (n_shape + 1L):p] <- x * H
    }
  } else if (type == "hazard") {
    # Current API returns relative hazard = exp(eta); shape params drop out.
    h_rel <- exp(eta)
    if (length(beta) > 0L && !is.null(x)) {
      J[, (n_shape + 1L):p] <- x * h_rel
    }
  } else if (type == "linear_predictor") {
    if (length(beta) > 0L && !is.null(x)) {
      J[, (n_shape + 1L):p] <- x
    }
  }

  J
}

# ---------------------------------------------------------------------------
# Jacobian: multiphase
# ---------------------------------------------------------------------------
#
# H(t|x) = sum_j mu_j(x) Phi_j(t),   mu_j(x) = exp(log_mu_j + x_j beta_j).
#
# For each phase j, its sub-vector in theta contains (in order):
#   log_mu_j, (log_t_half_j, nu_j, m_j) or (log_tau_j, gamma_j, alpha_j, eta_j)
#   or () for a "constant" phase, then beta_j covariates.
#
# Derivatives:
#   dH/d log_mu_j    = mu_j Phi_j
#   dH/d shape_s_j   = mu_j dPhi_j/d shape_s     (from .hzr_phase_derivatives
#                                                 / .hzr_g3_phase_derivatives)
#   dH/d beta_jk     = x_ik mu_j Phi_j

#' Jacobian of multiphase cumulative-hazard predictions
#'
#' Only `cumulative_hazard` and `survival` are delegated here; `hazard`
#' and `linear_predictor` are rejected upstream for multiphase models.
#'
#' @param theta MLE parameter vector.
#' @param time Prediction times (length n).
#' @param phases Named list of `hzr_phase` objects.
#' @param covariate_counts Named integer vector.
#' @param x_list Named list of per-phase design matrices.
#' @param p Length of theta.
#' @param per_phase Logical; if `TRUE`, return a named list of per-phase
#'   `n x p` Jacobians (each with only that phase's columns nonzero) instead
#'   of their sum.
#' @return Numeric `n x p` Jacobian of `H(t|x)` (default), or a named list of
#'   per-phase `n x p` Jacobians when `per_phase = TRUE`.
#' @keywords internal
.hzr_predict_jacobian_multiphase <- function(theta, time, phases,
                                              covariate_counts, x_list, p,
                                              per_phase = FALSE) {
  n <- length(time)
  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  # One n x p matrix per phase; each phase fills only its own columns.
  J_list <- stats::setNames(
    lapply(seq_along(phases), function(i) matrix(0, nrow = n, ncol = p)),
    names(phases)
  )

  pos <- 1L
  for (nm in names(phases)) {
    ph <- phases[[nm]]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], ph)
    Jp <- J_list[[nm]]

    # mu_j(x_i)
    if (length(pars$beta) > 0L && !is.null(x_list[[nm]])) {
      eta_j <- pars$log_mu + as.numeric(x_list[[nm]] %*% pars$beta)
    } else {
      eta_j <- rep(pars$log_mu, n)
    }
    mu_j <- exp(eta_j)

    # Phi_j and shape derivatives
    if (ph$type == "constant") {
      Phi_j <- time
      Jp[, pos] <- mu_j * Phi_j
      pos <- pos + 1L
    } else if (ph$type == "g3") {
      tau_j <- exp(pars$log_tau)
      pd <- .hzr_g3_phase_derivatives(time, tau = tau_j,
                                        gamma = pars$gamma,
                                        alpha = pars$alpha,
                                        eta = pars$eta)
      Phi_j <- pd$Phi
      Jp[, pos] <- mu_j * Phi_j
      Jp[, pos + 1L] <- mu_j * pd$dPhi_dlog_tau
      Jp[, pos + 2L] <- mu_j * pd$dPhi_dgamma
      Jp[, pos + 3L] <- mu_j * pd$dPhi_dalpha
      Jp[, pos + 4L] <- mu_j * pd$dPhi_deta
      pos <- pos + 5L
    } else {
      t_half_j <- exp(pars$log_t_half)
      pd <- .hzr_phase_derivatives(time, t_half = t_half_j,
                                    nu = pars$nu, m = pars$m,
                                    type = ph$type)
      Phi_j <- pd$Phi
      Jp[, pos] <- mu_j * Phi_j
      Jp[, pos + 1L] <- mu_j * pd$dPhi_dlog_thalf
      Jp[, pos + 2L] <- mu_j * pd$dPhi_dnu
      Jp[, pos + 3L] <- mu_j * pd$dPhi_dm
      pos <- pos + 4L
    }

    # Covariate betas
    n_beta <- covariate_counts[[nm]]
    if (n_beta > 0L && !is.null(x_list[[nm]])) {
      x_phase <- x_list[[nm]]
      for (k in seq_len(n_beta)) {
        Jp[, pos] <- x_phase[, k] * mu_j * Phi_j
        pos <- pos + 1L
      }
    } else if (n_beta > 0L) {
      pos <- pos + n_beta
    }

    J_list[[nm]] <- Jp
  }

  if (per_phase) {
    return(J_list)
  }
  Reduce(`+`, J_list)
}

# ---------------------------------------------------------------------------
# Numeric jacobian fallback for exp / loglogistic / lognormal
# ---------------------------------------------------------------------------
#
# Rather than hand-code each distribution's closed form, we differentiate the
# prediction function numerically through the existing single-distribution
# cumhaz/survival formulas in predict.hazard().  The prediction function is
# built once per call and returns a length-n vector for any candidate theta.

#' Numeric Jacobian of a single-distribution prediction w.r.t. theta
#'
#' @param predict_fn Function(theta) -> numeric vector of length n.
#' @param theta MLE parameter vector.
#' @return Numeric n x p Jacobian.
#' @keywords internal
.hzr_predict_jacobian_numeric <- function(predict_fn, theta) {
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    return(numDeriv::jacobian(predict_fn, theta))
  }
  eps <- (.Machine$double.eps) ^ (1 / 3)
  p <- length(theta)
  fit <- predict_fn(theta)
  n <- length(fit)
  J <- matrix(0, nrow = n, ncol = p)
  for (j in seq_len(p)) {
    h <- eps * max(abs(theta[j]), 1)
    tp <- theta
    tm <- theta
    tp[j] <- tp[j] + h
    tm[j] <- tm[j] - h
    J[, j] <- (predict_fn(tp) - predict_fn(tm)) / (2 * h)
  }
  J
}

# ---------------------------------------------------------------------------
# Delta-method sandwich + scale transform
# ---------------------------------------------------------------------------

#' Per-row standard errors via the delta method
#'
#' @param J Jacobian (n x p).
#' @param vcov p x p variance-covariance matrix.
#' @return Numeric vector of length n.
#' @keywords internal
.hzr_predict_se_from_jacobian <- function(J, vcov) {
  JV <- J %*% vcov
  se <- sqrt(pmax(rowSums(JV * J), 0))
  se
}

#' Apply a scale transform + back-transform for CLs
#'
#' @param fit Point estimate (numeric vector length n).
#' @param se_nat SE on the natural scale (numeric vector length n).
#' @param level Confidence level.
#' @param scale One of "natural", "log", "loglog_survival".
#' @return data.frame with columns fit, se.fit, lower, upper.
#' @keywords internal
.hzr_predict_cl_from_se <- function(fit, se_nat, level, scale) {
  alpha <- 1 - level
  z <- stats::qnorm(1 - alpha / 2)

  lower <- rep(NA_real_, length(fit))
  upper <- rep(NA_real_, length(fit))

  if (scale == "natural") {
    lower <- fit - z * se_nat
    upper <- fit + z * se_nat
  } else if (scale == "log") {
    # eta > 0 required; any non-positive fit gets NA CLs.
    pos <- is.finite(fit) & fit > 0
    se_log <- ifelse(pos, se_nat / abs(fit), NA_real_)
    lower[pos] <- fit[pos] * exp(-z * se_log[pos])
    upper[pos] <- fit[pos] * exp(z * se_log[pos])
  } else if (scale == "loglog_survival") {
    # fit is S = exp(-H); se_nat is the SE of H (not S).  Build CLs on
    # log(-log S) = log H.
    #   S_lwr = exp(-H * exp(+z * se(log H)))
    #   S_upr = exp(-H * exp(-z * se(log H)))
    # The caller passes `fit = S` and `se_nat = se(H)`; we reconstruct
    # H = -log(S) to compute se(log H) = se(H) / H.
    H <- -log(pmin(pmax(fit, .Machine$double.xmin), 1))
    pos <- is.finite(H) & H > 0
    se_logH <- ifelse(pos, se_nat / H, NA_real_)
    lower[pos] <- exp(-H[pos] * exp(z * se_logH[pos]))
    upper[pos] <- exp(-H[pos] * exp(-z * se_logH[pos]))
  } else if (scale == "logit_survival") {
    # SAS HAZPRED's survival CL: build on the logit of cumulative incidence,
    # Z = logit(1 - S) = log(e^H - 1), with se(Z) = se(H) / (1 - S), then
    # back-transform S = 1 / (e^Z + 1).  `fit` is S = exp(-H); `se_nat` is
    # se(H).  (Reproduces hzp_calc_srv_CL.c.)
    H <- -log(pmin(pmax(fit, .Machine$double.xmin), 1))
    Fc <- 1 - fit
    pos <- is.finite(H) & H > 0 & Fc > 0
    Z <- ifelse(pos, log(expm1(H)), NA_real_)
    seZ <- ifelse(pos, se_nat / Fc, NA_real_)
    # S = 1 / (e^Z + 1) = plogis(-Z); use plogis for a numerically stable
    # logistic back-transform (avoids exp() overflow at large |Z +/- z*seZ|).
    lower[pos] <- stats::plogis(-(Z[pos] + z * seZ[pos]))
    upper[pos] <- stats::plogis(-(Z[pos] - z * seZ[pos]))
  } else {
    stop("Unknown CL scale: '", scale, "'.", call. = FALSE)
  }

  data.frame(fit = fit, se.fit = se_nat, lower = lower, upper = upper)
}

# ---------------------------------------------------------------------------
# Free-parameter vcov helper (shared by aggregate and decomposed se.fit paths)
# ---------------------------------------------------------------------------

#' Which variances overflowed or underflowed (#566)
#'
#' `TRUE` for a variance that is `NaN`, infinite, or below the smallest normal
#' double. `NA` is `FALSE`: that is how a fixed or masked parameter is marked.
#' @param d Numeric vector of variances, a vcov diagonal.
#' @return Logical vector, never `NA`.
#' @noRd
.hzr_variance_unrepresentable <- function(d) {
  # A negative variance is not an overflow or underflow; .hzr_free_vcov()
  # names it separately.
  is.nan(d) | is.infinite(d) |
    (!is.na(d) & d >= 0 & d < .Machine$double.xmin)
}

#' Which parameters a Weibull prediction depends on (#566)
#'
#' Read off the prediction's form, not the Jacobian's values. The relative
#' hazard and linear predictor never read mu or nu; a coefficient is unused
#' only where its covariate is 0 in every row.
#' @return Logical vector of length `p`.
#' @noRd
.hzr_weibull_used <- function(type, x, p) {
  used <- rep(TRUE, p)
  if (!type %in% c("cumulative_hazard", "survival")) used[1:2] <- FALSE
  if (p > 2L) {
    used[3:p] <- if (is.null(x)) FALSE else colSums(is.na(x) | x != 0) > 0
  }
  used
}

#' Which parameters each phase of a multiphase prediction depends on (#566)
#'
#' A phase depends on its own shape slots, and on each of its coefficients
#' whose covariate is not 0 in every row. The walk follows
#' .hzr_predict_jacobian_multiphase(). If it does not account for exactly `p`
#' parameters, every parameter counts as used, which withholds rather than
#' drops.
#' @return Named list of logical vectors of length `p`, one per phase.
#' @noRd
.hzr_multiphase_used <- function(phases, covariate_counts, x_list, p) {
  out <- list()
  pos <- 1L
  for (nm in names(phases)) {
    u <- rep(FALSE, p)
    n_shape <- switch(phases[[nm]]$type, constant = 1L, g3 = 5L, 4L)
    u[pos:(pos + n_shape - 1L)] <- TRUE
    pos <- pos + n_shape
    n_beta <- covariate_counts[[nm]]
    if (n_beta > 0L) {
      x <- x_list[[nm]]
      if (!is.null(x)) {
        u[pos:(pos + n_beta - 1L)] <- colSums(is.na(x) | x != 0) > 0
      }
      pos <- pos + n_beta
    }
    out[[nm]] <- u
  }
  if (pos - 1L != p) {
    out <- lapply(out, function(u) rep(TRUE, p))
  }
  out
}

#' Whether a vcov has the shape .hzr_free_vcov() can use
#' @noRd
.hzr_vcov_shape_ok <- function(vcov_mat, p) {
  !is.null(vcov_mat) && is.matrix(vcov_mat) &&
    nrow(vcov_mat) == p && ncol(vcov_mat) == p
}

#' Vcov submatrix for the delta-method sandwich
#'
#' Restricts the sandwich to the parameters the prediction uses. A parameter
#' held fixed (`fixed`, from the fit's `fixed_mask`; e.g. `fixed = "shapes"`)
#' carries an NA row and is dropped as known. (The CoE-conserved `log_mu` is
#' one when its full-information recompute failed; the caller warns.) An NA
#' variance on a parameter that is NOT fixed was masked by `.hzr_safe_solve()`
#' as non-positive: it withholds the standard error if the prediction uses
#' it (#586). Returns `NULL` (with a warning) when CLs cannot be computed.
#' Shared by the aggregate and decomposed se.fit paths.
#'
#' @param vcov_mat The fitted vcov (or NULL / wrong shape).
#' @param p Length of the parameter vector.
#' @param unused Indices of parameters the prediction does not depend on, by
#'   the prediction's form (`.hzr_weibull_used()`, `.hzr_multiphase_used()`),
#'   never read off numeric zeros in the Jacobian. They are dropped exactly,
#'   whatever their variance.
#' @param fixed Logical, `TRUE` for a parameter held fixed in the fit, or
#'   `NULL` when nothing is fixed (every single-distribution fit).
#' @param param_names Names for the warning, or `NULL`.
#' @return `list(vcov_use, free_idx)`, or `NULL` if unusable.
#' @keywords internal
.hzr_free_vcov <- function(vcov_mat, p, unused = integer(0), fixed = NULL,
                           param_names = NULL) {
  if (!.hzr_vcov_shape_ok(vcov_mat, p)) {
    warning("Variance-covariance matrix is unavailable; ",
            "standard errors and CLs will be NA.", call. = FALSE)
    return(NULL)
  }
  is_fixed <- if (length(fixed) == p) as.logical(fixed) %in% TRUE else
    rep(FALSE, p)
  # The sandwich runs over these: neither unused by the prediction nor fixed.
  free_idx <- setdiff(seq_len(p), union(unused, which(is_fixed)))
  # A fixed or masked parameter carries an NA variance and is dropped from the
  # sandwich below. A variance that overflowed or underflowed is not that. A
  # Weibull mu's variance carries mu^2: it is Inf for mu near exp(600), and
  # subnormal or 0 for mu near exp(-400), while its covariances are still
  # ordinary numbers. Dropping such a parameter computed the standard error
  # as if it were known exactly, and keeping a variance of 0 loses the
  # sandwich's positive term; both gave finite, wrong standard errors (#566).
  d <- diag(vcov_mat)[free_idx]
  # An NA here is not a fixed parameter: .hzr_safe_solve() masks a
  # non-positive variance with NA, the same mark. Dropping it computed the
  # standard error as if the parameter were known exactly (#586).
  masked <- free_idx[is.na(d) & !is.nan(d)]
  if (length(masked)) {
    nms <- if (length(param_names) == p) param_names[masked] else
      paste0("par", masked)
    warning("Variance-covariance matrix has no variance for an estimated ",
            "parameter (", paste(nms, collapse = ", "), "): it was masked ",
            "as non-positive when the Hessian was inverted; standard errors ",
            "and CLs will be NA.", call. = FALSE)
    return(NULL)
  }
  # A negative variance is a different fault: the covariance is not positive
  # definite (possible after the Weibull back-transform when the internal one
  # is indefinite). Name it; main returned an SE of 0 there.
  if (any(!is.na(d) & d < 0)) {
    warning("Variance-covariance matrix has a negative variance, so it is not ",
            "positive definite; standard errors and CLs will be NA.",
            call. = FALSE)
    return(NULL)
  }
  # A variance that overflowed or underflowed is not that either. A Weibull
  # mu's variance carries mu^2: it is Inf for mu near exp(600), and subnormal
  # or 0 for mu near exp(-400), while its covariances are still ordinary
  # numbers (#566).
  if (any(.hzr_variance_unrepresentable(d))) {
    warning("Variance-covariance matrix has a variance that cannot be ",
            "represented (it overflowed or underflowed, as for a parameter ",
            "far outside the usual range); standard errors and CLs will be ",
            "NA. Centre or rescale the covariates and refit.", call. = FALSE)
    return(NULL)
  }
  vcov_use <- vcov_mat[free_idx, free_idx, drop = FALSE]
  # Unused and fixed parameters are already out, so an infinite or NA
  # covariance here is one the prediction reads (#587 item 1).
  if (!all(is.finite(vcov_use))) {
    warning("Variance-covariance matrix has a covariance that is not finite ",
            "among the parameters this prediction uses; standard errors and ",
            "CLs will be NA.", call. = FALSE)
    return(NULL)
  }
  # A positive diagonal does not make a covariance: an indefinite one gave a
  # negative quadratic form, clamped to an SE of 0, or a positive one that is
  # no variance (#586). Screened on the correlation scale, so the tolerance
  # does not depend on the parameters' units.
  if (length(free_idx) > 0L) {
    s <- sqrt(diag(vcov_use))
    e <- eigen(vcov_use / outer(s, s), symmetric = TRUE,
               only.values = TRUE)$values
    if (min(e) < -length(e) * .Machine$double.eps * max(abs(e))) {
      warning("Variance-covariance matrix is not positive definite over the ",
              "parameters this prediction uses (the fit's Hessian was not at ",
              "a proper maximum); standard errors and CLs will be NA.",
              call. = FALSE)
      return(NULL)
    }
  }
  list(vcov_use = vcov_use, free_idx = free_idx)
}

#' Warn when a failed CoE recompute leaves a phase's variance out (#586)
#'
#' The conserved `log_mu` is then held fixed in `fixed_mask` and dropped from
#' the sandwich, and the object does not record which position it is.
#' @noRd
.hzr_warn_conserved_variance <- function(object) {
  if ("conserved_phase_variance" %in% object$degraded) {
    warning("Conservation of Events could not recompute the conserved ",
            "phase's variance for this fit, so these standard errors leave ",
            "it out and may be understated.", call. = FALSE)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

#' Compute point predictions + delta-method CLs
#'
#' Dispatched from predict.hazard() when `se.fit = TRUE`.
#'
#' DELTA-METHOD TARGET
#' -------------------
#' The quantity we actually differentiate depends on the prediction type and
#' the CL scale the type uses:
#'
#'   type = "cumulative_hazard" -> target = H, log-scale CL
#'   type = "survival"          -> target = H, log-log-survival CL
#'                                 (final `fit` is exp(-H); reported
#'                                 `se.fit` is S * se(H), the SE of S)
#'   type = "hazard"            -> target = exp(eta) (single-dist only),
#'                                 log-scale CL
#'   type = "linear_predictor"  -> target = eta, natural-scale CL
#'
#' The caller supplies `diff_fn(theta)` that returns the target vector of
#' length n.  For Weibull and multiphase we build J analytically; for
#' exp / loglogistic / lognormal we fall through to a numeric jacobian of
#' `diff_fn`.
#'
#' @param object Fitted `hazard` object.
#' @param type One of "hazard", "linear_predictor", "survival",
#'   "cumulative_hazard".
#' @param time Numeric prediction times or NULL (for type = "hazard" /
#'   "linear_predictor").
#' @param x Design matrix (single-distribution) or NULL.  Ignored for
#'   multiphase; use `x_list` instead.
#' @param x_list Named list of per-phase design matrices (multiphase only).
#' @param cov_counts Named integer covariate counts (multiphase only).
#' @param phases Named list of `hzr_phase` objects (multiphase only).
#' @param level Confidence level (default 0.95).
#' @param diff_fn Required function(theta) -> numeric vector of length n
#'   returning the delta-method target (H, exp(eta), or eta depending on
#'   `type`).  Used for the point estimate AND for the numeric jacobian
#'   fallback in exp / loglogistic / lognormal.
#' @param conf_type Survival CL transform: `"log-log"` (default, on
#'   `log(-log S)`) or `"logit"` (on `logit(1 - S)`, HAZPRED-compatible).
#'   Only consulted when `type == "survival"`.
#' @return data.frame with columns `fit`, `se.fit`, `lower`, `upper`.
#' @keywords internal
.hzr_predict_with_se <- function(object, type, time = NULL,
                                   x = NULL, x_list = NULL,
                                   cov_counts = NULL, phases = NULL,
                                   level = 0.95, diff_fn,
                                   conf_type = c("log-log", "logit")) {
  # Only survival CLs use conf_type; validate only then so an ignored value
  # (hazard / cumulative_hazard / linear_predictor) is not an error.
  if (type == "survival") conf_type <- match.arg(conf_type)

  theta <- object$fit$theta
  p <- length(theta)
  dist <- object$spec$dist
  target <- diff_fn(theta)

  # --- Build Jacobian of the delta-method target ---------------------------
  # Built before the vcov is screened, so the screen knows which parameters
  # this prediction depends on (#566).
  jacobian <- function() {
    if (dist == "weibull") {
      .hzr_predict_jacobian_weibull(type, theta, time, x, p)
    } else if (dist == "multiphase") {
      # The analytic multiphase Jacobian is for the cumulative hazard; the
      # instantaneous hazard has no analytic Jacobian here, so fall back to
      # a numeric Jacobian of its evaluator. (linear_predictor is rejected
      # upstream for multiphase.)
      if (type == "hazard") {
        .hzr_predict_jacobian_numeric(diff_fn, theta)
      } else {
        .hzr_predict_jacobian_multiphase(theta, time, phases,
                                          cov_counts, x_list, p)
      }
    } else {
      .hzr_predict_jacobian_numeric(diff_fn, theta)
    }
  }
  # An absent or wrong-sized vcov withholds every SE; say so before paying
  # for a Jacobian.
  if (.hzr_vcov_shape_ok(object$fit$vcov, p)) {
    J <- jacobian()
    # A parameter this prediction does not depend on, by its form: a Weibull
    # mu or nu for the relative hazard and linear predictor, or a coefficient
    # whose covariate is 0 in every row of newdata. Its variance, representable
    # or not, cannot change the standard error, so it is dropped exactly
    # rather than withholding it. Read off the form, not J's values: a column
    # can be zero numerically (exp(eta) underflowing, mu * time == 1) where
    # the prediction still depends on that parameter (#566).
    used <- if (dist == "weibull") {
      .hzr_weibull_used(type, x, p)
    } else if (dist == "multiphase") {
      Reduce(`|`, .hzr_multiphase_used(phases, cov_counts, x_list, p))
    } else {
      rep(TRUE, p)
    }
    fv <- .hzr_free_vcov(object$fit$vcov, p, unused = which(!used),
                         fixed = object$fit$fixed_mask,
                         param_names = names(theta))
  } else {
    fv <- .hzr_free_vcov(object$fit$vcov, p)
  }
  .hzr_warn_conserved_variance(object)
  if (is.null(fv)) {
    n <- length(target)
    fit <- if (type == "survival") exp(-target) else target
    na_vec <- rep(NA_real_, n)
    return(data.frame(fit = fit, se.fit = na_vec,
                      lower = na_vec, upper = na_vec))
  }
  vcov_use <- fv$vcov_use
  free_idx <- fv$free_idx

  # Restrict J to the free columns so the sandwich dimensions match.
  J <- J[, free_idx, drop = FALSE]

  # --- Sandwich -------------------------------------------------------------
  se_target <- .hzr_predict_se_from_jacobian(J, vcov_use)

  # --- Scale transform -----------------------------------------------------
  if (type == "survival") {
    fit <- exp(-target)
    surv_scale <- if (conf_type == "logit") "logit_survival" else "loglog_survival"
    res <- .hzr_predict_cl_from_se(fit, se_target, level, surv_scale)
    # The limits need se(H); the reported SE is of `fit` = S itself.
    # dS/dH = -S, so se(S) = S * se(H), as summary.survfit's std.err.
    res$se.fit <- fit * se_target
    res
  } else if (type == "cumulative_hazard" || type == "hazard") {
    .hzr_predict_cl_from_se(target, se_target, level, "log")
  } else if (type == "linear_predictor") {
    .hzr_predict_cl_from_se(target, se_target, level, "natural")
  } else {
    stop("Unknown prediction type '", type, "' for se.fit.", call. = FALSE)
  }
}

#' Per-phase + total cumulative-hazard predictions with delta-method CLs
#'
#' Long-format companion to [.hzr_predict_with_se()] for multiphase models under
#' `decompose = TRUE`. Computes log-scale CLs for the total cumulative hazard
#' and for each phase's additive contribution `H_j(t) = mu_j(x) Phi_j(t)`.
#' Per-phase CLs use only that phase's parameter block, so they do NOT sum to
#' the total CL (cross-phase covariance contributes only to the total).
#'
#' @param object Fitted multiphase `hazard` object.
#' @param time Numeric prediction times (length n).
#' @param x_list Named list of per-phase design matrices.
#' @param cov_counts Named integer covariate counts.
#' @param phases Named list of `hzr_phase` objects.
#' @param level Confidence level.
#' @return Long `data.frame(time, component, fit, se.fit, lower, upper)`;
#'   `component` is an ordered factor with levels `c("total", names(phases))`.
#' @keywords internal
.hzr_predict_with_se_decomposed <- function(object, time, x_list,
                                             cov_counts, phases,
                                             level = 0.95) {
  theta <- object$fit$theta
  p <- length(theta)
  components <- c("total", names(phases))

  # Point estimates: total + per-phase cumulative hazard.
  ph <- .hzr_multiphase_cumhaz(time, theta, phases, cov_counts, x_list,
                               per_phase = TRUE)
  fit_by_comp <- c(list(total = ph$total),
                   stats::setNames(lapply(names(phases), function(nm) ph[[nm]]),
                                   names(phases)))

  make_long <- function(cl_list) {
    rows <- lapply(components, function(cmp) {
      data.frame(time = time, component = cmp,
                 fit = cl_list[[cmp]]$fit, se.fit = cl_list[[cmp]]$se.fit,
                 lower = cl_list[[cmp]]$lower, upper = cl_list[[cmp]]$upper,
                 stringsAsFactors = FALSE)
    })
    res <- do.call(rbind, rows)
    res$component <- factor(res$component, levels = components, ordered = TRUE)
    rownames(res) <- NULL
    res
  }

  na_cl_of <- function(cmp) {
    n <- length(time)
    data.frame(fit = fit_by_comp[[cmp]], se.fit = rep(NA_real_, n),
               lower = rep(NA_real_, n), upper = rep(NA_real_, n))
  }

  vcov_mat <- object$fit$vcov
  if (!.hzr_vcov_shape_ok(vcov_mat, p)) {
    .hzr_free_vcov(vcov_mat, p)  # warns
    na_cl <- lapply(components, na_cl_of)
    names(na_cl) <- components
    return(make_long(na_cl))
  }

  # Per-phase Jacobians; the total is their sum.
  J_list <- .hzr_predict_jacobian_multiphase(theta, time, phases, cov_counts,
                                             x_list, p, per_phase = TRUE)
  J_by_comp <- c(list(total = Reduce(`+`, J_list)), J_list)
  used_list <- .hzr_multiphase_used(phases, cov_counts, x_list, p)
  used_by_comp <- c(list(total = Reduce(`|`, used_list)), used_list)

  # Screen the vcov per component, against the parameters that component
  # depends on (#566): a bad variance in one phase's block withholds that
  # phase's SE and the total's, not every phase's. Every warning of
  # .hzr_free_vcov() is followed by a NULL return, so catching it loses
  # nothing; the distinct messages are raised once below.
  cl_list <- lapply(components, function(cmp) {
    J <- J_by_comp[[cmp]]
    fv <- tryCatch(.hzr_free_vcov(vcov_mat, p,
                                  unused = which(!used_by_comp[[cmp]]),
                                  fixed = object$fit$fixed_mask,
                                  param_names = names(theta)),
                   warning = function(w) conditionMessage(w))
    if (!is.list(fv)) {
      return(list(cl = na_cl_of(cmp), msg = fv))
    }
    se <- .hzr_predict_se_from_jacobian(J[, fv$free_idx, drop = FALSE],
                                        fv$vcov_use)
    list(cl = .hzr_predict_cl_from_se(fit_by_comp[[cmp]], se, level, "log"),
         msg = character(0))
  })
  msgs <- unique(unlist(lapply(cl_list, `[[`, "msg")))
  for (m in msgs) warning(m, call. = FALSE)
  .hzr_warn_conserved_variance(object)
  cl_list <- lapply(cl_list, `[[`, "cl")
  names(cl_list) <- components
  make_long(cl_list)
}
