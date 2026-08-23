# likelihood-multiphase.R -- Multiphase additive hazard likelihood and optimizer
#
# REFERENCES
# ----------
# Blackstone EH, Naftel DC, Turner ME Jr. The decomposition of time-varying
# hazard into phases, each incorporating a separate stream of concomitant
# information. J Am Stat Assoc. 1986;81(395):615-624.
#
# Rajeswaran J, Blackstone EH, Ehrlinger J, Li L, Ishwaran H, Parides MK.
# Probability of atrial fibrillation after ablation: Using a parametric
# nonlinear temporal decomposition mixed effects model. Stat Methods Med Res.
# 2018;27(1):126-141.
#
# MODEL
# -----
# Additive cumulative hazard across J phases:
#
#   H(t | x) = sum_{j=1}^{J}  mu_j(x) * Phi_j(t; t_half_j, nu_j, m_j)
#
# where:
#   mu_j(x) = exp(alpha_j + x_j * beta_j)  -- phase-specific log-linear scale
#   Phi_j(t) depends on phase type:
#     "cdf"      -> G(t)             bounded [0, 1]   (early risk)
#     "hazard"   -> -log(1 - G(t))   monotone incr.   (late risk)
#     "constant" -> t                flat rate         (background)
#
# Total hazard:
#   h(t | x) = sum_j  mu_j(x) * phi_j(t)
#
# where phi_j = dPhi_j/dt.
#
# Survival:
#   S(t | x) = exp(-H(t | x))
#
# THETA LAYOUT (internal / estimation scale)
# ------------------------------------------
# For each phase j, the sub-vector is:
#   [log_mu_j, log_t_half_j, nu_j, m_j, beta_j_1, ..., beta_j_pj]
#   (constant phases omit log_t_half, nu, m)
#
# The full theta is the concatenation of all phase sub-vectors.
#
# FUNCTIONS
# ---------
#   .hzr_logl_multiphase()     -- log-likelihood
#   .hzr_optim_multiphase()    -- optimizer using .hzr_optim_generic()
#   .hzr_multiphase_cumhaz()   -- compute total H(t|x) from theta + phases
#   .hzr_multiphase_hazard()   -- compute total h(t|x) from theta + phases
#   .hzr_split_theta()         -- split theta into per-phase sub-vectors


# ============================================================================
# Theta vector manipulation
# ============================================================================

#' Split the full theta vector into per-phase sub-vectors
#'
#' @param theta Numeric vector -- full parameter vector (internal scale).
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector -- number of covariates per phase.
#' @return Named list of numeric vectors, one per phase.
#' @keywords internal
.hzr_split_theta <- function(theta, phases, covariate_counts) {
  result <- vector("list", length(phases))
  names(result) <- names(phases)
  pos <- 1L

  for (nm in names(phases)) {
    np <- .hzr_phase_n_params(phases[[nm]], n_covariates = covariate_counts[[nm]])
    result[[nm]] <- theta[pos:(pos + np - 1L)]
    pos <- pos + np
  }

  result
}


#' Extract parameters from a phase sub-vector (internal scale)
#'
#' @param theta_j Numeric sub-vector for one phase.
#' @param phase An `hzr_phase` object.
#' @return Named list with `log_mu`, `log_t_half`, `nu`, `m`, `beta`.
#' @keywords internal
.hzr_unpack_phase_theta <- function(theta_j, phase) {
  pos <- 1L
  log_mu <- theta_j[pos]
  pos <- pos + 1L

  if (phase$type == "g3") {
    log_tau <- theta_j[pos]
    pos <- pos + 1L
    gamma_  <- theta_j[pos]
    pos <- pos + 1L
    alpha_  <- theta_j[pos]
    pos <- pos + 1L
    eta_    <- theta_j[pos]
    pos <- pos + 1L
    log_t_half <- NA_real_
    nu <- NA_real_
    m <- NA_real_
  } else if (phase$type != "constant") {
    log_t_half <- theta_j[pos]
    pos <- pos + 1L
    nu <- theta_j[pos]
    pos <- pos + 1L
    m <- theta_j[pos]
    pos <- pos + 1L
    log_tau <- NA_real_
    gamma_ <- NA_real_
    alpha_ <- NA_real_
    eta_ <- NA_real_
  } else {
    log_t_half <- NA_real_
    nu <- NA_real_
    m <- NA_real_
    log_tau <- NA_real_
    gamma_ <- NA_real_
    alpha_ <- NA_real_
    eta_ <- NA_real_
  }

  beta <- if (pos <= length(theta_j)) theta_j[pos:length(theta_j)] else numeric(0)

  list(log_mu = log_mu, log_t_half = log_t_half, nu = nu, m = m,
       log_tau = log_tau, gamma = gamma_, alpha = alpha_, eta = eta_,
       beta = beta)
}


# ============================================================================
# Cumulative hazard and hazard computations
# ============================================================================

#' Compute total cumulative hazard H(t|x) for multiphase model
#'
#' @param time Numeric vector of times (n).
#' @param theta Full parameter vector (internal scale).
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector of covariate counts per phase.
#' @param x_list Named list of design matrices, one per phase. Each element is
#'   an n x p_j matrix or NULL.
#' @param per_phase Logical; if TRUE return a list with total and per-phase
#'   contributions.
#' @return If `per_phase = FALSE`: numeric vector of length n (total H(t|x)).
#'   If `per_phase = TRUE`: named list with `$total` and one element per phase.
#' @keywords internal
.hzr_multiphase_cumhaz <- function(time, theta, phases, covariate_counts,
                                    x_list, per_phase = FALSE) {
  n <- length(time)
  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  total <- rep(0, n)
  phase_contributions <- if (per_phase) vector("list", length(phases)) else NULL

  for (i in seq_along(phases)) {
    nm <- names(phases)[i]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])

    # mu_j(x) = exp(log_mu + x_j * beta_j)
    if (length(pars$beta) > 0 && !is.null(x_list[[nm]])) {
      eta_j <- pars$log_mu + as.numeric(x_list[[nm]] %*% pars$beta)
    } else {
      eta_j <- rep(pars$log_mu, n)
    }
    mu_j <- exp(eta_j)

    # Phi_j(t)
    if (phases[[nm]]$type == "constant") {
      phi_j <- time
    } else if (phases[[nm]]$type == "g3") {
      tau_j <- exp(pars$log_tau)
      d3 <- hzr_decompos_g3(time, tau = tau_j, gamma = pars$gamma,
                              alpha = pars$alpha, eta = pars$eta)
      phi_j <- d3$G3
    } else {
      t_half_j <- exp(pars$log_t_half)
      # Guard against invalid decomposition parameters.  Answer in the shape the
      # caller asked for: per_phase callers (the CoE adjustment, which runs
      # inside the objective on every evaluation) index the result by name, and
      # a bare vector there errors out of optim() and loses the whole start
      # rather than merely penalising the step.
      if (pars$m < 0 && pars$nu < 0) {
        inf_n <- rep(Inf, n)
        if (!per_phase) return(inf_n)
        infeasible <- rep(list(inf_n), length(phases))
        names(infeasible) <- names(phases)
        infeasible$total <- inf_n
        return(infeasible)
      }
      phi_j <- hzr_phase_cumhaz(time, t_half = t_half_j, nu = pars$nu,
                                  m = pars$m, type = phases[[nm]]$type)
    }

    contrib <- mu_j * phi_j
    total <- total + contrib

    if (per_phase) phase_contributions[[i]] <- contrib
  }

  if (per_phase) {
    names(phase_contributions) <- names(phases)
    phase_contributions$total <- total
    return(phase_contributions)
  }

  total
}


#' Compute total instantaneous hazard h(t|x) for multiphase model
#'
#' @inheritParams .hzr_multiphase_cumhaz
#' @return Numeric vector of length n.
#' @keywords internal
.hzr_multiphase_hazard <- function(time, theta, phases, covariate_counts,
                                    x_list) {
  n <- length(time)
  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  total <- rep(0, n)

  for (i in seq_along(phases)) {
    nm <- names(phases)[i]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])

    if (length(pars$beta) > 0 && !is.null(x_list[[nm]])) {
      eta_j <- pars$log_mu + as.numeric(x_list[[nm]] %*% pars$beta)
    } else {
      eta_j <- rep(pars$log_mu, n)
    }
    mu_j <- exp(eta_j)

    if (phases[[nm]]$type == "constant") {
      dphi_j <- rep(1, n)
    } else if (phases[[nm]]$type == "g3") {
      tau_j <- exp(pars$log_tau)
      d3 <- hzr_decompos_g3(time, tau = tau_j, gamma = pars$gamma,
                              alpha = pars$alpha, eta = pars$eta)
      dphi_j <- d3$g3
    } else {
      t_half_j <- exp(pars$log_t_half)
      if (pars$m < 0 && pars$nu < 0) return(rep(Inf, n))
      dphi_j <- hzr_phase_hazard(time, t_half = t_half_j, nu = pars$nu,
                                   m = pars$m, type = phases[[nm]]$type)
    }

    total <- total + mu_j * dphi_j
  }

  total
}


# ============================================================================
# Conservation of Events (CoE) -- Turner's theorem
# ============================================================================
#
# For the additive hazard model H(t|x) = Sigma_j mu_j(x) * Phi_j(t),
# the MLE constraint (score equation for mu) implies:
#
#   Sigma_i E(i) = Sigma_i H(t_i | x_i)
#
# i.e., total observed events equals total predicted cumulative hazard.
#
# This allows one phase's log_mu to be solved analytically at each
# optimizer iteration, reducing the optimization dimension by 1.
#
# Reference: Turner ME Jr., as implemented in C HAZARD setcoe.c / consrv.c.

#' Find the position of each phase's log_mu in the full theta vector
#'
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector of per-phase covariate counts.
#' @return Named integer vector: position of log_mu for each phase.
#' @keywords internal
.hzr_log_mu_positions <- function(phases, covariate_counts) {
  positions <- setNames(integer(length(phases)), names(phases))
  pos <- 1L
  for (nm in names(phases)) {
    positions[[nm]] <- pos
    np <- .hzr_phase_n_params(phases[[nm]],
                               n_covariates = covariate_counts[[nm]])
    pos <- pos + np
  }
  positions
}


#' Select the phase whose log_mu will be solved by conservation
#'
#' Chooses the phase contributing the largest share of total cumulative
#' hazard at the current theta, matching the C HAZARD SETCOE strategy.
#'
#' When one phase dominates by an extreme margin (>10x the median contribution)
#' it is typically due to poorly-scaled shape parameters at the starting theta
#' rather than a genuine signal that the phase should absorb all events.
#' Selecting such a phase as the fixmu phase traps the optimizer: CoE keeps
#' rescaling that phase's log_mu upward while the optimizer tries to drive it
#' to near-zero.  We therefore exclude extreme outliers and select among the
#' remaining phases.  If all phases are outliers (or only one phase exists)
#' the largest is used as a fallback.
#'
#' @param theta Full parameter vector (internal scale).
#' @param time Numeric vector of follow-up times.
#' @param status Numeric event indicator.
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector.
#' @param x_list Named list of per-phase design matrices.
#' @param weights Optional numeric vector of row weights (length n). Defaults to
#'   unit weights. Applied when summing per-phase cumhaz so that selection
#'   happens on the same scale as the (weighted) observed event count.
#' @param time_lower Optional numeric vector of counting-process entry (start)
#'   times. When supplied, phases are ranked by entry-time cumulative hazard
#'   `H(stop) - H(start)`, the scale on which events are conserved. `NULL` (the
#'   default) means no truncation, i.e. `H(start) = 0`.
#' @return Character: name of the phase to fix.
#' @keywords internal
.hzr_select_fixmu_phase <- function(theta, time, status,
                                     phases, covariate_counts, x_list,
                                     weights = NULL, time_lower = NULL) {
  decomp <- .hzr_multiphase_cumhaz(time, theta, phases,
                                     covariate_counts, x_list,
                                     per_phase = TRUE)
  # Left truncation (counting-process entry times): rank phases by the cumulative
  # hazard accrued over follow-up, H(stop) - H(start), the same scale on which
  # events are conserved. See `.hzr_conserve_events`.
  if (!is.null(time_lower)) {
    decomp_start <- .hzr_multiphase_cumhaz(time_lower, theta, phases,
                                            covariate_counts, x_list,
                                            per_phase = TRUE)
    for (nm in names(phases)) {
      decomp[[nm]] <- decomp[[nm]] - decomp_start[[nm]]
    }
  }
  if (is.null(weights)) weights <- rep(1, length(time))
  phase_sums <- vapply(names(phases), function(nm) {
    sum(weights * decomp[[nm]])
  }, numeric(1))

  # Exclude extreme outliers: a phase whose contribution exceeds 10x the
  # median is likely dominated by large shape-parameter values at the starting
  # theta rather than reflecting the true relative importance of the phase.
  # Fixing such a phase traps the optimizer when the true MLE has that phase
  # near zero.  Among the remaining (non-outlier) phases, pick the largest.
  if (length(phase_sums) > 1L) {
    finite_sums <- phase_sums[is.finite(phase_sums)]
    if (length(finite_sums) >= 1L) {
      med <- stats::median(finite_sums, na.rm = TRUE)
      if (is.finite(med) && med > 0) {
        non_outlier <- phase_sums[is.finite(phase_sums) & phase_sums <= 10 * med]
        if (length(non_outlier) >= 1L) {
          return(names(which.max(non_outlier)))
        }
      }
    }
  }

  names(which.max(phase_sums))
}


#' Apply the Conservation of Events adjustment to one phase's log_mu
#'
#' Given the current theta vector, analytically solve the fixmu phase's
#' log_mu so that total predicted events = total observed events.
#'
#' This is called BEFORE each likelihood evaluation inside the optimizer,
#' matching the C HAZARD `CONSRV` entry point.
#'
#' @param theta Full parameter vector (internal scale).
#' @param fixmu_phase Character: name of the phase whose log_mu is solved.
#' @param fixmu_pos Integer: position of that log_mu in theta.
#' @param time Numeric vector of follow-up times.
#' @param status Numeric event indicator.
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector.
#' @param x_list Named list of per-phase design matrices.
#' @param total_events Numeric: sum of observed events (precomputed, weighted
#'   when `weights` is supplied).
#' @param weights Optional numeric vector of row weights (length n). Defaults to
#'   unit weights. Applied when summing per-phase cumhaz so Turner's adjustment
#'   is computed on the same scale as `total_events`.
#' @param time_lower Optional numeric vector of counting-process entry (start)
#'   times. When supplied, conservation is enforced on the entry-time scale --
#'   `Sum E = Sum [H(stop) - H(start)]` -- by subtracting the entry-time
#'   cumulative hazard, matching the multiphase likelihood (and C HAZARD
#'   `setcoe` under `LCENSOR`/`STARTTME`). `NULL` (the default) means no
#'   truncation, i.e. `H(start) = 0`.
#' @return Updated theta vector with fixmu phase's log_mu adjusted.
#' @keywords internal
.hzr_conserve_events <- function(theta, fixmu_phase, fixmu_pos,
                                  time, status,
                                  phases, covariate_counts, x_list,
                                  total_events, weights = NULL,
                                  time_lower = NULL) {
  # Compute per-phase cumulative hazard contributions
  decomp <- .hzr_multiphase_cumhaz(time, theta, phases,
                                     covariate_counts, x_list,
                                     per_phase = TRUE)

  # Left truncation: subtract the per-phase entry-time cumulative hazard so the
  # conserved quantity is H(stop) - H(start), the same constraint the score
  # equation for mu satisfies. Without this, CoE conserves Sum H(stop) and the
  # fixmu phase absorbs the spurious Sum H(start), biasing its intercept.
  if (!is.null(time_lower)) {
    decomp_start <- .hzr_multiphase_cumhaz(time_lower, theta, phases,
                                            covariate_counts, x_list,
                                            per_phase = TRUE)
    decomp$total <- decomp$total - decomp_start$total
    decomp[[fixmu_phase]] <- decomp[[fixmu_phase]] - decomp_start[[fixmu_phase]]
  }

  if (is.null(weights)) weights <- rep(1, length(time))

  # Total predicted events (weighted sum of cumhaz across all observations)
  sumcz <- sum(weights * decomp$total)

  # Weighted contribution from the fixmu phase alone
  sumcj <- sum(weights * decomp[[fixmu_phase]])

  if (sumcj <= 0 || !is.finite(sumcj)) return(theta)

  # Discrepancy: how many events are unaccounted for
  devent <- total_events - sumcz

  # Events the fixmu phase should absorb
  jevent <- sumcj + devent

  if (jevent <= 0 || !is.finite(jevent)) return(theta)

  # Multiplicative adjustment in log scale
  lfactor <- log(jevent / sumcj)

  if (!is.finite(lfactor)) return(theta)

  theta[fixmu_pos] <- theta[fixmu_pos] + lfactor
  theta
}


# ============================================================================
# Log-likelihood
# ============================================================================

#' Log-likelihood for multiphase additive hazard model
#'
#' @param theta Full parameter vector (internal scale).
#' @param time Numeric vector of follow-up times (n).
#' @param status Numeric event indicator: 1 = event, 0 = right-censored,
#'   -1 = left-censored, 2 = interval-censored.
#' @param time_lower Optional lower bounds for interval censoring.
#' @param time_upper Optional upper bounds for left/interval censoring.
#' @param x Design matrix (unused directly; kept for interface compatibility).
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector of per-phase covariate counts.
#' @param x_list Named list of per-phase design matrices.
#' @param return_gradient Logical (ignored; gradient via separate function).
#' @param return_hessian Logical (ignored).
#' @param ... Ignored.
#' @return Scalar log-likelihood. Returns -Inf for infeasible parameters.
#' @keywords internal
.hzr_logl_multiphase <- function(theta, time, status,
                                  time_lower = NULL, time_upper = NULL,
                                  x = NULL, weights = NULL,
                                  phases, covariate_counts, x_list,
                                  return_gradient = FALSE,
                                  return_hessian = FALSE, ...) {

  n <- length(time)
  if (is.null(weights)) weights <- rep(1, n)

  # Feasibility: check parameter constraints

  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  for (nm in names(phases)) {
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])
    if (phases[[nm]]$type %in% c("cdf", "hazard")) {
      if (pars$m < 0 && pars$nu < 0) return(-Inf)
    }
    if (phases[[nm]]$type == "g3") {
      # G3 requires gamma > 0, alpha >= 0, eta > 0
      if (pars$gamma <= 0 || pars$eta <= 0 || pars$alpha < 0) return(-Inf)
    }
  }

  # Validate event times
  if (any(status == 1 & time <= 0)) return(-Inf)

  # Normalize censoring bounds
  lower <- if (is.null(time_lower)) time else time_lower
  upper <- if (is.null(time_upper)) time else time_upper

  # Cumulative hazard at event/censoring times
  cumhaz <- .hzr_multiphase_cumhaz(time, theta, phases, covariate_counts, x_list)
  if (any(!is.finite(cumhaz))) return(-Inf)

  # Counting-process entry-time cumulative hazard; H(start) = 0 when no
  # `time_lower` supplied (plain right-censored data).  For status in
  # {0, 1}, the contribution becomes H(stop) - H(start); for status == 2,
  # this value is unused (interval-censoring handled separately below).
  if (is.null(time_lower)) {
    cumhaz_start <- rep(0, n)
  } else {
    cumhaz_start <- .hzr_multiphase_cumhaz(time_lower, theta, phases,
                                             covariate_counts, x_list)
    if (any(!is.finite(cumhaz_start))) return(-Inf)
  }

  # Instantaneous hazard at event times (only needed for exact events)
  idx_event <- status == 1

  logl <- 0

  # Exact events: w * [log h(stop) - (H(stop) - H(start))]
  if (any(idx_event)) {
    haz <- .hzr_multiphase_hazard(time, theta, phases, covariate_counts, x_list)
    if (any(!is.finite(haz[idx_event])) || any(haz[idx_event] <= 0)) return(-Inf)
    logl <- logl + sum(weights[idx_event] *
                         (log(haz[idx_event]) -
                            (cumhaz[idx_event] - cumhaz_start[idx_event])))
  }

  # Right-censored: w * [-(H(stop) - H(start))]
  idx_right <- status == 0
  if (any(idx_right)) {
    logl <- logl - sum(weights[idx_right] *
                         (cumhaz[idx_right] - cumhaz_start[idx_right]))
  }

  # Left-censored: w * [log(1 - exp(-H(u)))]
  idx_left <- status == -1
  if (any(idx_left)) {
    cumhaz_upper <- .hzr_multiphase_cumhaz(upper, theta, phases,
                                             covariate_counts, x_list)
    logl <- logl + sum(weights[idx_left] *
                         hzr_log1mexp(cumhaz_upper[idx_left]))
  }

  # Interval-censored: w * [-H(l) + log(1 - exp(-(H(u) - H(l))))]
  idx_interval <- status == 2
  if (any(idx_interval)) {
    cumhaz_lower <- .hzr_multiphase_cumhaz(lower, theta, phases,
                                             covariate_counts, x_list)
    cumhaz_upper_iv <- .hzr_multiphase_cumhaz(upper, theta, phases,
                                                covariate_counts, x_list)
    delta_h <- cumhaz_upper_iv[idx_interval] - cumhaz_lower[idx_interval]
    logl <- logl + sum(weights[idx_interval] *
                         (-cumhaz_lower[idx_interval] +
                            hzr_log1mexp(delta_h)))
  }

  if (!is.finite(logl)) return(-Inf)

  logl
}


# ============================================================================
# Analytic gradient
# ============================================================================

#' Gradient of the multiphase log-likelihood
#'
#' Computes the score vector \eqn{d\ell / d\theta} using analytic chain-rule
#' formulas for `log_mu` and `beta` parameters, and central-difference
#' derivatives (via `.hzr_phase_derivatives()`) for shape parameters
#' (`log_t_half`, `nu`, `m`).
#'
#' @inheritParams .hzr_logl_multiphase
#' @return Numeric vector of length `length(theta)` -- the gradient.
#'   Returns a zero vector if any component is non-finite (guards optimizer).
#' @keywords internal
.hzr_gradient_multiphase <- function(theta, time, status,
                                      time_lower = NULL, time_upper = NULL,
                                      x = NULL, weights = NULL,
                                      phases, covariate_counts, x_list,
                                      ...) {
  n <- length(time)
  p <- length(theta)
  grad <- numeric(p)
  if (is.null(weights)) weights <- rep(1, n)

  # Feasibility check
  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  for (nm in names(phases)) {
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])
    if (phases[[nm]]$type %in% c("cdf", "hazard")) {
      if (pars$m < 0 && pars$nu < 0) return(grad)
    }
    if (phases[[nm]]$type == "g3") {
      if (pars$gamma <= 0 || pars$eta <= 0 || pars$alpha < 0) return(grad)
    }
  }

  # Normalize censoring bounds
  lower <- if (is.null(time_lower)) time else time_lower
  upper <- if (is.null(time_upper)) time else time_upper

  # Observation type masks
  idx_event    <- which(status == 1)
  idx_right    <- which(status == 0)
  idx_left     <- which(status == -1)

  idx_interval <- which(status == 2)

  has_left     <- length(idx_left) > 0
  has_interval <- length(idx_interval) > 0

  # -- Per-phase quantities --------------------------------------------------
  # We need, for each phase j and each observation i:
  #   mu_j(x_i), Phi_j(t_i), phi_j(t_i), and shape derivatives of Phi and phi
  #
  # Then the per-observation gradient contributions are assembled from:
  #   H(t) = sum_j mu_j * Phi_j    (cumulative hazard)
  #   h(t) = sum_j mu_j * phi_j    (instantaneous hazard)

  # Pre-compute total H(t) and h(t) at observation times (and censoring bounds)
  H_t <- rep(0, n)      # total cumulative hazard at time[i]
  h_t <- rep(0, n)      # total instantaneous hazard at time[i]

  # Storage for per-phase intermediate quantities
  n_phases <- length(phases)
  phase_mu         <- vector("list", n_phases)  # mu_j(x_i) for each i
  phase_Phi        <- vector("list", n_phases)  # Phi_j(t_i) for each i
  phase_phi        <- vector("list", n_phases)  # phi_j(t_i) for each i
  phase_deriv      <- vector("list", n_phases)  # derivative list at time
  phase_Phi_start  <- vector("list", n_phases)  # Phi_j(start_i), if needed
  phase_deriv_start <- vector("list", n_phases) # derivative list at start
  # The start time must be defined EXACTLY as the log-likelihood defines it:
  # `.hzr_logl_multiphase()` subtracts H(time_lower) from every status 0/1 row
  # whenever `time_lower` is supplied, with no further condition. An earlier
  # `time_lower < time` filter here excluded rows entering at their own event
  # or censoring time, so the derivative was taken of a different function from
  # the one being evaluated: the reported gradient pointed somewhere the
  # objective did not, and the optimizer followed it out of the region where
  # the fit means anything. The
  # `> 0` test is only a skip: H(0) = 0, so those rows contribute nothing
  # either way, and it must match `has_start` below or the term is weighted in
  # while its derivative is left at zero.
  need_start <- !is.null(time_lower) &&
                any(time_lower > 0 & status %in% c(0L, 1L))
  start_vec <- if (need_start) {
    sv <- rep(0, n)
    epoch_idx <- status %in% c(0L, 1L)
    sv[epoch_idx] <- time_lower[epoch_idx]
    sv
  } else {
    NULL
  }

  for (j in seq_along(phases)) {
    nm <- names(phases)[j]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])

    # mu_j(x_i) = exp(log_mu + x * beta)
    if (length(pars$beta) > 0 && !is.null(x_list[[nm]])) {
      eta_j <- pars$log_mu + as.numeric(x_list[[nm]] %*% pars$beta)
    } else {
      eta_j <- rep(pars$log_mu, n)
    }
    mu_j <- exp(eta_j)
    phase_mu[[j]] <- mu_j

    # Phi_j, phi_j, and shape derivatives at observation time
    if (phases[[nm]]$type == "constant") {
      phase_Phi[[j]] <- time
      phase_phi[[j]] <- rep(1, n)
      phase_deriv[[j]] <- NULL  # no shape params for constant
    } else if (phases[[nm]]$type == "g3") {
      tau_j <- exp(pars$log_tau)
      d3 <- hzr_decompos_g3(time, tau = tau_j, gamma = pars$gamma,
                              alpha = pars$alpha, eta = pars$eta)
      phase_Phi[[j]] <- d3$G3
      phase_phi[[j]] <- d3$g3
      # G3 shape derivatives via finite differences (for now)
      phase_deriv[[j]] <- .hzr_g3_phase_derivatives(
        time, tau = tau_j, gamma = pars$gamma,
        alpha = pars$alpha, eta = pars$eta)
    } else {
      t_half_j <- exp(pars$log_t_half)
      pd <- .hzr_phase_derivatives(time, t_half = t_half_j, nu = pars$nu,
                                    m = pars$m, type = phases[[nm]]$type)
      phase_Phi[[j]]   <- pd$Phi
      phase_phi[[j]]   <- pd$phi
      phase_deriv[[j]] <- pd
    }

    # Matching quantities at start_vec -- only Phi_j and its shape
    # derivatives matter for the H(start) gradient contribution.
    if (need_start) {
      if (phases[[nm]]$type == "constant") {
        phase_Phi_start[[j]] <- start_vec
        phase_deriv_start[[j]] <- NULL
      } else if (phases[[nm]]$type == "g3") {
        tau_j <- exp(pars$log_tau)
        d3s <- hzr_decompos_g3(start_vec, tau = tau_j, gamma = pars$gamma,
                                 alpha = pars$alpha, eta = pars$eta)
        phase_Phi_start[[j]] <- d3s$G3
        phase_deriv_start[[j]] <- .hzr_g3_phase_derivatives(
          start_vec, tau = tau_j, gamma = pars$gamma,
          alpha = pars$alpha, eta = pars$eta)
      } else {
        t_half_j <- exp(pars$log_t_half)
        pds <- .hzr_phase_derivatives(start_vec, t_half = t_half_j,
                                        nu = pars$nu, m = pars$m,
                                        type = phases[[nm]]$type)
        phase_Phi_start[[j]]   <- pds$Phi
        phase_deriv_start[[j]] <- pds
      }
    }

    H_t <- H_t + mu_j * phase_Phi[[j]]
    h_t <- h_t + mu_j * phase_phi[[j]]
  }

  # Guard: if total hazard or cumhaz is non-finite, return zero gradient
  if (any(!is.finite(H_t)) || any(!is.finite(h_t))) return(grad)

  # -- Per-observation weight vectors ----------------------------------------
  # dLogl/dH(t_i) depends on observation type:
  #   event:         d/dH [ log h(t) - H(t) ] = -1        (H part)
  #   right-cens:    d/dH [ -H(t) ]            = -1
  #   left-cens:     d/dH [ log(1-exp(-H(u))) ] = exp(-H) / (1 - exp(-H))
  #   interval-cens: more complex (handled separately)
  #
  # For events, there's also the d/d[theta] of log h(t):
  #   d(log h) / d(theta_j) = (1/h) * dh/d(theta_j)

  # Weight for the cumulative hazard part at `time[i]`: w_H_i = dLogl/dH(t_i),
  # scaled by Fisher weights so the score matches the weighted log-likelihood.
  w_H <- numeric(n)
  w_H[idx_event] <- -weights[idx_event]
  w_H[idx_right] <- -weights[idx_right]

  if (has_left) {
    H_upper <- .hzr_multiphase_cumhaz(upper, theta, phases, covariate_counts,
                                       x_list)
    # d/dH log(1 - exp(-H)) = exp(-H) / (1 - exp(-H))
    exp_neg_H_u <- exp(-H_upper[idx_left])
    one_minus   <- pmax(1 - exp_neg_H_u, .Machine$double.xmin)
    w_H[idx_left] <- weights[idx_left] * (exp_neg_H_u / one_minus)
    # Left-censored uses H(upper), not H(time), so we handle below
  }

  # H(start) contribution for counting-process rows (status in {0, 1}).
  # Per-row contribution is `+H(start_i)` inside the log-likelihood, so
  # the derivative w.r.t. theta picks up `+weights_i * dH(start)/dtheta`.
  has_start <- !is.null(time_lower) && any(time_lower > 0 & status %in% c(0L, 1L))
  w_H_start <- numeric(n)
  if (has_start) {
    w_H_start[idx_event] <- weights[idx_event]
    w_H_start[idx_right] <- weights[idx_right]
  }

  # Weight for the instantaneous hazard part (events only): w_i / h(t_i).
  # This multiplies dh/d(theta_j).
  inv_h <- numeric(n)
  if (length(idx_event) > 0) {
    h_event <- h_t[idx_event]
    h_event <- pmax(h_event, .Machine$double.xmin)
    inv_h[idx_event] <- weights[idx_event] / h_event
  }

  # -- Assemble gradient per phase -------------------------------------------
  pos <- 1L  # position in the full theta vector

  for (j in seq_along(phases)) {
    nm <- names(phases)[j]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])
    mu_j  <- phase_mu[[j]]
    Phi_j <- phase_Phi[[j]]
    phi_j <- phase_phi[[j]]
    Phi_j_start <- phase_Phi_start[[j]]  # NULL when need_start = FALSE

    # -- d(logl) / d(log_mu_j) ------------------------------------------
    # From H part: w_H_i * dH/d(log_mu) = w_H_i * mu_j * Phi_j
    # From h part (events): inv_h_i * dh/d(log_mu) = inv_h_i * mu_j * phi_j
    dlogl_dlog_mu <- sum(w_H * mu_j * Phi_j) +
                     sum(inv_h[idx_event] * mu_j[idx_event] * phi_j[idx_event])

    # Counting-process entry-time contribution: +H(start) inside the LL,
    # so its derivative picks up +w_H_start * dH(start)/d(log_mu).
    if (need_start) {
      dlogl_dlog_mu <- dlogl_dlog_mu +
        sum(w_H_start * mu_j * Phi_j_start)
    }

    # Left-censored: uses H(upper), need separate Phi_j evaluation there
    if (has_left) {
      if (phases[[nm]]$type == "constant") {
        Phi_j_upper <- upper[idx_left]
      } else if (phases[[nm]]$type == "g3") {
        tau_j <- exp(pars$log_tau)
        d3_upper <- hzr_decompos_g3(upper[idx_left], tau = tau_j,
                                      gamma = pars$gamma, alpha = pars$alpha,
                                      eta = pars$eta)
        Phi_j_upper <- d3_upper$G3
      } else {
        t_half_j <- exp(pars$log_t_half)
        Phi_j_upper <- hzr_phase_cumhaz(upper[idx_left], t_half = t_half_j,
                                          nu = pars$nu, m = pars$m,
                                          type = phases[[nm]]$type)
      }
      dlogl_dlog_mu <- dlogl_dlog_mu +
        sum(w_H[idx_left] * mu_j[idx_left] * Phi_j_upper)
      # Subtract the w_H * mu * Phi at time[i] that we already added above
      # (for left-censored, the H contribution is at upper, not time)
      dlogl_dlog_mu <- dlogl_dlog_mu -
        sum(w_H[idx_left] * mu_j[idx_left] * Phi_j[idx_left])
    }

    # Interval-censored: use numerical fallback for these observations
    # (handled at end of loop via finite-difference correction)

    grad[pos] <- dlogl_dlog_mu
    pos <- pos + 1L

    # -- Shape parameters (non-constant phases only) ---------------------
    if (phases[[nm]]$type %in% c("cdf", "hazard")) {
      pd <- phase_deriv[[j]]

      # Compute shape derivatives at upper[idx_left] for left-censoring
      pd_upper <- NULL
      if (has_left) {
        t_half_j <- exp(pars$log_t_half)
        pd_upper <- .hzr_phase_derivatives(upper[idx_left],
          t_half = t_half_j, nu = pars$nu, m = pars$m,
          type = phases[[nm]]$type)
      }

      pd_start <- if (need_start) phase_deriv_start[[j]] else NULL

      # For each shape param s in {log_t_half, nu, m}:
      shape_derivs <- list(
        list(dPhi = pd$dPhi_dlog_thalf, dphi = pd$dphi_dlog_thalf,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_dlog_thalf,
             dPhi_start = if (!is.null(pd_start)) pd_start$dPhi_dlog_thalf),
        list(dPhi = pd$dPhi_dnu,        dphi = pd$dphi_dnu,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_dnu,
             dPhi_start = if (!is.null(pd_start)) pd_start$dPhi_dnu),
        list(dPhi = pd$dPhi_dm,         dphi = pd$dphi_dm,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_dm,
             dPhi_start = if (!is.null(pd_start)) pd_start$dPhi_dm)
      )

      for (s in seq_along(shape_derivs)) {
        dPhi_ds <- shape_derivs[[s]]$dPhi
        dphi_ds <- shape_derivs[[s]]$dphi

        dlogl_ds <- sum(w_H * mu_j * dPhi_ds) +
                    sum(inv_h[idx_event] * mu_j[idx_event] * dphi_ds[idx_event])

        # Counting-process start contribution: +w_H_start * mu * dPhi(start)
        if (need_start && !is.null(shape_derivs[[s]]$dPhi_start)) {
          dlogl_ds <- dlogl_ds +
            sum(w_H_start * mu_j * shape_derivs[[s]]$dPhi_start)
        }

        # Left-censoring correction: swap dPhi at time for dPhi at upper
        if (has_left && !is.null(shape_derivs[[s]]$dPhi_upper)) {
          dlogl_ds <- dlogl_ds +
            sum(w_H[idx_left] * mu_j[idx_left] * shape_derivs[[s]]$dPhi_upper) -
            sum(w_H[idx_left] * mu_j[idx_left] * dPhi_ds[idx_left])
        }

        grad[pos] <- dlogl_ds
        pos <- pos + 1L
      }
    } else if (phases[[nm]]$type == "g3") {
      pd <- phase_deriv[[j]]

      # Compute shape derivatives at upper[idx_left] for left-censoring
      pd_upper <- NULL
      if (has_left) {
        tau_j <- exp(pars$log_tau)
        pd_upper <- .hzr_g3_phase_derivatives(upper[idx_left],
          tau = tau_j, gamma = pars$gamma,
          alpha = pars$alpha, eta = pars$eta)
      }

      pd_start <- if (need_start) phase_deriv_start[[j]] else NULL

      # For each shape param s in {log_tau, gamma, alpha, eta}:
      shape_derivs <- list(
        list(dPhi = pd$dPhi_dlog_tau, dphi = pd$dphi_dlog_tau,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_dlog_tau,
             dPhi_start = if (!is.null(pd_start)) pd_start$dPhi_dlog_tau),
        list(dPhi = pd$dPhi_dgamma,   dphi = pd$dphi_dgamma,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_dgamma,
             dPhi_start = if (!is.null(pd_start)) pd_start$dPhi_dgamma),
        list(dPhi = pd$dPhi_dalpha,   dphi = pd$dphi_dalpha,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_dalpha,
             dPhi_start = if (!is.null(pd_start)) pd_start$dPhi_dalpha),
        list(dPhi = pd$dPhi_deta,     dphi = pd$dphi_deta,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_deta,
             dPhi_start = if (!is.null(pd_start)) pd_start$dPhi_deta)
      )

      for (s in seq_along(shape_derivs)) {
        dPhi_ds <- shape_derivs[[s]]$dPhi
        dphi_ds <- shape_derivs[[s]]$dphi

        dlogl_ds <- sum(w_H * mu_j * dPhi_ds) +
                    sum(inv_h[idx_event] * mu_j[idx_event] * dphi_ds[idx_event])

        # Counting-process start contribution
        if (need_start && !is.null(shape_derivs[[s]]$dPhi_start)) {
          dlogl_ds <- dlogl_ds +
            sum(w_H_start * mu_j * shape_derivs[[s]]$dPhi_start)
        }

        # Left-censoring correction: swap dPhi at time for dPhi at upper
        if (has_left && !is.null(shape_derivs[[s]]$dPhi_upper)) {
          dlogl_ds <- dlogl_ds +
            sum(w_H[idx_left] * mu_j[idx_left] * shape_derivs[[s]]$dPhi_upper) -
            sum(w_H[idx_left] * mu_j[idx_left] * dPhi_ds[idx_left])
        }

        grad[pos] <- dlogl_ds
        pos <- pos + 1L
      }
    }

    # -- Covariate coefficients beta_jk ----------------------------------
    # dH/d(beta_jk) = x_ik * mu_j_i * Phi_j_i
    # dh/d(beta_jk) = x_ik * mu_j_i * phi_j_i
    if (length(pars$beta) > 0 && !is.null(x_list[[nm]])) {
      x_j <- x_list[[nm]]
      for (k in seq_along(pars$beta)) {
        x_k <- x_j[, k]
        dlogl_dbeta <- sum(w_H * x_k * mu_j * Phi_j) +
                       sum(inv_h[idx_event] * x_k[idx_event] *
                           mu_j[idx_event] * phi_j[idx_event])

        # Counting-process start contribution: +w_H_start * x_k * mu_j * Phi_start
        if (need_start) {
          dlogl_dbeta <- dlogl_dbeta +
            sum(w_H_start * x_k * mu_j * Phi_j_start)
        }

        if (has_left) {
          dlogl_dbeta <- dlogl_dbeta +
            sum(w_H[idx_left] * x_k[idx_left] *
                mu_j[idx_left] * Phi_j_upper) -
            sum(w_H[idx_left] * x_k[idx_left] *
                mu_j[idx_left] * Phi_j[idx_left])
        }

        grad[pos] <- dlogl_dbeta
        pos <- pos + 1L
      }
    } else {
      pos <- pos + covariate_counts[[nm]]
    }
  }

  # -- Interval-censored fallback: add correction via finite difference ---
  # Interval-censored observations are uncommon; their gradient contribution
  # is added via one-sided numerical difference on the weighted
  # log-likelihood of just those observations.
  if (has_interval) {
    w_iv <- weights[idx_interval]
    logl_iv <- function(th) {
      cumhaz_l <- .hzr_multiphase_cumhaz(lower, th, phases, covariate_counts,
                                          x_list)
      cumhaz_u <- .hzr_multiphase_cumhaz(upper, th, phases, covariate_counts,
                                          x_list)
      delta <- cumhaz_u[idx_interval] - cumhaz_l[idx_interval]
      sum(w_iv *
            (-cumhaz_l[idx_interval] + hzr_log1mexp(delta)))
    }
    eps_rel <- sqrt(.Machine$double.eps)
    ll0_iv <- logl_iv(theta)
    for (i in seq_len(p)) {
      h_i <- eps_rel * max(abs(theta[i]), 1)
      theta_p <- theta
      theta_p[i] <- theta_p[i] + h_i
      grad[i] <- grad[i] + (logl_iv(theta_p) - ll0_iv) / h_i
    }
  }

  # Safety: zero out non-finite entries
  grad[!is.finite(grad)] <- 0

  grad
}


# ============================================================================
# Optimizer
# ============================================================================

#' Deterministic multi-start perturbations
#'
#' Draws the offsets that shift the starting values for every multiphase start
#' after the first.  The draw is taken from an internally seeded stream and the
#' ambient `.Random.seed` is restored on exit, so a fit is reproducible without
#' `set.seed()` and never advances the caller's stream.  Both halves matter:
#' seeding alone would still let an intervening draw change the fit, and
#' restoring alone would leave the fit dependent on whatever the stream
#' happened to hold when it was called.
#'
#' @param n_starts Number of optimization starts.
#' @param n_par Length of the free-parameter vector being perturbed.
#' @param seed Single whole number selecting the ensemble of offsets.
#' @param sd Standard deviation of the offsets (default 0.5).
#' @return A list of `max(0, n_starts - 1)` numeric vectors of length `n_par`;
#'   empty for a single start, which is never perturbed.
#' @keywords internal
.hzr_start_perturbations <- function(n_starts, n_par, seed, sd = 0.5) {
  n_draw <- max(0L, as.integer(n_starts) - 1L)
  if (n_draw == 0L) return(list())

  # Restore the caller's stream exactly as we found it -- including the case
  # where the session had never drawn from it and `.Random.seed` did not
  # exist, where leaving ours behind would make the caller's next draw follow
  # from our seed instead of from R's own initialization.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
            add = TRUE)
  }

  set.seed(seed)
  lapply(seq_len(n_draw), function(i) stats::rnorm(n_par, sd = sd))
}


#' Fit a multiphase additive hazard model via maximum likelihood
#'
#' Assembles starting values from phase specifications, resolves per-phase
#' design matrices, and delegates to `.hzr_optim_generic()`.
#'
#' @param time Numeric vector of follow-up times.
#' @param status Numeric event indicator vector.
#' @param time_lower Optional lower bounds for interval censoring.
#' @param time_upper Optional upper bounds for left/interval censoring.
#' @param x Global design matrix (n x p) or NULL.
#' @param theta_start Starting parameter vector (full internal scale).
#'   If NULL, assembled automatically from phase specs.
#' @param control Named list of control options.
#' @param phases Named list of validated `hzr_phase` objects.
#' @param formula_global The global formula (used when phases have no
#'   phase-specific formula).
#' @param data Data frame containing covariates (needed for phase-specific
#'   formula evaluation).
#' @return List with par (internal scale), value, convergence, vcov, etc.
#' @keywords internal
.hzr_optim_multiphase <- function(time, status,
                                   time_lower = NULL, time_upper = NULL,
                                   x = NULL,
                                   theta_start = NULL,
                                   weights = NULL,
                                   control = list(),
                                   phases,
                                   formula_global = NULL,
                                   data = NULL) {

  if (is.null(weights)) weights <- rep(1, length(time))

  phases <- .hzr_validate_phases(phases)

  # --- Resolve per-phase design matrices and covariate counts ----------------
  # Build phase-specific design matrices with na.action = na.pass so that
  # rows with NA covariates are preserved (not silently dropped).  We then
  # align all observation-aligned inputs (time, status, weights, ...) by
  # removing any row that has NA in any phase's design matrix.  Without this
  # alignment, model.matrix()'s default na.omit shortens x_j relative to
  # time, and mu_j * phi_j emits a vector-recycling warning on every
  # likelihood evaluation.
  x_list <- vector("list", length(phases))
  names(x_list) <- names(phases)
  covariate_counts <- setNames(integer(length(phases)), names(phases))

  for (nm in names(phases)) {
    ph <- phases[[nm]]
    if (!is.null(ph$formula) && !is.null(data)) {
      # Phase-specific formula: build design matrix from data.  We route
      # through model.frame(na.action = na.pass) explicitly because
      # model.matrix() consults getOption("na.action") when it constructs
      # the frame internally, which defaults to na.omit and silently drops
      # rows with NA covariates.
      mf_j <- stats::model.frame(ph$formula, data = data,
                                   na.action = stats::na.pass)
      x_j <- stats::model.matrix(ph$formula, data = mf_j)[, -1L,
                                                           drop = FALSE]
      x_list[[nm]] <- x_j
      covariate_counts[[nm]] <- ncol(x_j)
    } else if (!is.null(x)) {
      # Inherit global design matrix
      x_list[[nm]] <- x
      covariate_counts[[nm]] <- ncol(x)
    } else {
      x_list[[nm]] <- NULL
      covariate_counts[[nm]] <- 0L
    }
  }

  # Align rows: drop any observation with an NA in any phase-specific design
  # matrix (or a row-count mismatch relative to time).
  n_time <- length(time)
  drop_rows <- logical(n_time)
  for (nm in names(phases)) {
    xm <- x_list[[nm]]
    if (is.null(xm) || ncol(xm) == 0L) next
    if (nrow(xm) != n_time) {
      stop("Design matrix for phase '", nm, "' has ", nrow(xm),
           " rows but length(time) is ", n_time,
           ". Check that 'data' matches 'time'/'status'.", call. = FALSE)
    }
    drop_rows <- drop_rows |
      apply(xm, 1L, function(row) any(is.na(row)))
  }
  if (any(drop_rows)) {
    keep <- !drop_rows
    time   <- time[keep]
    status <- status[keep]
    if (!is.null(time_lower)) time_lower <- time_lower[keep]
    if (!is.null(time_upper)) time_upper <- time_upper[keep]
    if (!is.null(weights))    weights    <- weights[keep]
    if (!is.null(x))          x          <- x[keep, , drop = FALSE]
    for (nm in names(phases)) {
      if (!is.null(x_list[[nm]])) {
        x_list[[nm]] <- x_list[[nm]][keep, , drop = FALSE]
      }
    }
  }

  # --- Assemble starting values if not provided ------------------------------
  if (is.null(theta_start)) {
    theta_start <- unlist(lapply(names(phases), function(nm) {
      .hzr_phase_start(phases[[nm]], n_covariates = covariate_counts[[nm]])
    }))
  }

  # Build parameter names for the full theta vector
  theta_names <- unlist(lapply(names(phases), function(nm) {
    cov_names <- if (covariate_counts[[nm]] > 0 && !is.null(x_list[[nm]])) {
      colnames(x_list[[nm]])
    } else {
      character(0)
    }
    # Ensure we have names even if colnames are NULL
    if (covariate_counts[[nm]] > 0 && length(cov_names) == 0) {
      cov_names <- paste0("x", seq_len(covariate_counts[[nm]]))
    }
    .hzr_phase_theta_names(phases[[nm]], nm, cov_names)
  }))

  names(theta_start) <- theta_names

  # --- Likelihood wrapper matching .hzr_optim_generic() signature -----------
  # NOTE: `weights` is an explicit formal so that .hzr_optim_generic can pass
  # it through without colliding with a closure capture via `...`
  # (which would trigger "formal argument 'weights' matched by multiple
  # actual arguments").
  logl_fn <- function(theta, time, status, time_lower, time_upper, x,
                      weights = NULL, ...) {
    .hzr_logl_multiphase(
      theta = theta, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x, weights = weights,
      phases = phases, covariate_counts = covariate_counts,
      x_list = x_list, ...
    )
  }

  # --- Analytic gradient (semi-analytic: chain-rule for mu/beta, central ----
  # differences for shape parameters via .hzr_phase_derivatives()).
  # Falls back to numerical gradient if analytic returns all zeros.
  #
  # IMPORTANT: Save a reference to the ORIGINAL (unwrapped) logl_fn for use

  # inside the gradient function's numerical fallback.  After this point,
  # logl_fn will be re-bound by CoE and fixed-mask wrappers, but the gradient
  # function always receives a FULL theta vector from its own wrapper chain,
  # so it must call the unwrapped version that also expects full theta.
  logl_fn_unwrapped <- logl_fn

  gradient_fn <- function(theta, time, status, time_lower, time_upper, x,
                          weights = NULL, ...) {
    grad <- .hzr_gradient_multiphase(
      theta = theta, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper, x = x,
      weights = weights,
      phases = phases, covariate_counts = covariate_counts, x_list = x_list
    )

    # Fallback: if gradient is all zero (e.g. at infeasible point), try
    # numerical gradient of the *weighted* LL to keep optimizer moving.
    if (all(grad == 0)) {
      eps_rel <- sqrt(.Machine$double.eps)
      p <- length(theta)
      ll0 <- logl_fn_unwrapped(theta, time, status, time_lower,
                                time_upper, x, weights = weights, ...)
      for (i in seq_len(p)) {
        h_i <- eps_rel * max(abs(theta[i]), 1)
        theta_plus <- theta
        theta_plus[i] <- theta_plus[i] + h_i
        ll_plus <- logl_fn_unwrapped(theta_plus, time, status,
                                      time_lower, time_upper, x,
                                      weights = weights, ...)
        grad[i] <- (ll_plus - ll0) / h_i
      }
    }

    grad
  }

  # --- Conservation of Events (CoE) setup -------------------------------------
  use_conserve <- if (!is.null(control$conserve)) control$conserve else TRUE
  control$conserve <- NULL  # remove before passing to optim

  fixmu_phase <- NULL
  fixmu_pos <- NULL
  total_events <- NULL
  log_mu_positions <- .hzr_log_mu_positions(phases, covariate_counts)

  # CoE is only valid for right-censored + exact-event data.  Interval and
  # left censoring require a different event-counting formulation.
  coe_supported_data <- all(status %in% c(0, 1))

  if (use_conserve && coe_supported_data && length(phases) >= 2L) {
    total_events <- sum(weights[status == 1])
    if (total_events > 0) {
      # Initial CoE scaling: adjust all log_mu proportionally (weighted sums
      # so the scaling is consistent with total_events).
      decomp_init <- .hzr_multiphase_cumhaz(
        time, theta_start, phases, covariate_counts, x_list,
        per_phase = TRUE
      )
      # Entry-time scale under left truncation: conserve H(stop) - H(start).
      init_total <- decomp_init$total
      if (!is.null(time_lower)) {
        decomp_init_start <- .hzr_multiphase_cumhaz(
          time_lower, theta_start, phases, covariate_counts, x_list,
          per_phase = TRUE
        )
        init_total <- init_total - decomp_init_start$total
      }
      sumcz_init <- sum(weights * init_total)
      if (sumcz_init > 0 && is.finite(sumcz_init)) {
        init_factor <- log(total_events / sumcz_init)
        if (is.finite(init_factor)) {
          for (nm in names(phases)) {
            theta_start[log_mu_positions[[nm]]] <-
              theta_start[log_mu_positions[[nm]]] + init_factor
          }
        }
      }

      # Select fixmu phase (largest weighted cumhaz contributor after scaling)
      fixmu_phase <- .hzr_select_fixmu_phase(
        theta_start, time, status, phases, covariate_counts, x_list,
        weights = weights, time_lower = time_lower
      )
      fixmu_pos <- log_mu_positions[[fixmu_phase]]

      # Apply precise CoE adjustment to the fixmu phase
      theta_start <- .hzr_conserve_events(
        theta_start, fixmu_phase, fixmu_pos,
        time, status, phases, covariate_counts, x_list, total_events,
        weights = weights, time_lower = time_lower
      )

      # Wrap logl and gradient: apply CoE before each evaluation
      logl_fn_pre_coe <- logl_fn
      logl_fn <- function(theta, time, status, time_lower,
                          time_upper, x, weights = NULL, ...) {
        theta <- .hzr_conserve_events(
          theta, fixmu_phase, fixmu_pos,
          time, status, phases, covariate_counts, x_list, total_events,
          weights = weights, time_lower = time_lower
        )
        logl_fn_pre_coe(theta, time, status, time_lower,
                        time_upper, x, weights = weights, ...)
      }

      gradient_fn_pre_coe <- gradient_fn
      gradient_fn <- function(theta, time, status, time_lower,
                              time_upper, x, weights = NULL, ...) {
        theta <- .hzr_conserve_events(
          theta, fixmu_phase, fixmu_pos,
          time, status, phases, covariate_counts, x_list, total_events,
          weights = weights, time_lower = time_lower
        )
        gradient_fn_pre_coe(theta, time, status, time_lower,
                            time_upper, x, weights = weights, ...)
      }
    } else {
      use_conserve <- FALSE
    }
  } else {
    use_conserve <- FALSE
  }

  # --- Fixed-parameter masking ------------------------------------------------
  # Build a logical mask: TRUE = free (optimized), FALSE = fixed (held constant)
  free_mask <- .hzr_phase_free_mask(phases, covariate_counts)

  # If CoE is active, also fix the fixmu phase's log_mu
  if (use_conserve && !is.null(fixmu_pos)) {
    free_mask[fixmu_pos] <- FALSE
  }

  any_fixed <- !all(free_mask)
  theta_fixed <- theta_start  # full vector; fixed entries stay at these values

  # When some parameters are fixed, wrap logl/gradient to operate on the
  # reduced (free-only) parameter vector.
  if (any_fixed) {
    free_idx <- which(free_mask)
    free_idx_eff <- free_idx

    # Expand reduced theta to full theta
    expand_theta <- function(theta_free) {
      theta_full <- theta_fixed
      theta_full[free_idx] <- theta_free
      theta_full
    }

    logl_fn_full <- logl_fn
    logl_fn <- function(theta, time, status, time_lower, time_upper, x,
                        weights = NULL, ...) {
      logl_fn_full(expand_theta(theta), time, status, time_lower,
                   time_upper, x, weights = weights, ...)
    }

    gradient_fn_full <- gradient_fn
    gradient_fn <- function(theta, time, status, time_lower, time_upper, x,
                            weights = NULL, ...) {
      grad_full <- gradient_fn_full(expand_theta(theta), time, status,
                                     time_lower, time_upper, x,
                                     weights = weights, ...)
      grad_full[free_idx]
    }

    theta_start_optim <- theta_start[free_idx]
  } else {
    free_idx_eff <- seq_along(theta_start)
    theta_start_optim <- theta_start
  }

  # --- Multi-start optimization -----------------------------------------------
  n_starts <- if (!is.null(control$n_starts)) control$n_starts else 5L
  control$n_starts <- NULL  # remove before passing to optim

  # Seed for the multi-start perturbations.  Fixed by default so that the
  # identical call returns the identical fit; vary it to probe a different set
  # of starts when a fit is suspected of sitting in a local optimum.
  #
  # The seed no longer decides whether there is a fit at all.  It once did:
  # roughly a quarter of seeds failed outright on the two-phase gradient
  # fixture, because the cumulative hazard's infeasible-shape guard threw
  # instead of penalising and took whole starts with it.  With that fixed, all
  # 50 seeds swept in test-multiphase-reproducibility.R converge.
  #
  # What the seed still decides is how good an objective the fit reaches, and
  # `starts` on the fit records which start the answer came from.  On real data
  # that means which optimum: the `avc` early+constant profile has an interior
  # maximum in m near 1 and falls away on both sides.
  #
  # On the two-phase gradient fixture it does NOT mean that, and the fixture
  # should not be read as an example of competing optima.  That likelihood has
  # no interior maximum in m.  Profiling it along m*nu, the objective climbs
  # past the -158.30 that a perturbed start reports (-158.29 at m = 50, -158.08
  # at m = 100) toward a finite supremum of -157.88 reached only as m -> Inf
  # with m*nu -> 0.753; the shape converges to 0.5*(t/t_half)^(1/k) truncated at
  # t_half*2^k.  So -158.30 is a point on a flat ridge where the optimizer hit
  # its tolerance, not a maximum -- its SE(m) is 42.8 on an m of 27.2 -- and
  # -159.15 is a near-boundary solution with nu driven to zero, whose vcov is
  # part non-finite.  Two seeds stopping at m = 27.1580 and 27.1501 look like
  # agreement on an optimum and are agreement on a stopping rule.
  #
  # Changing the default is a behavior change: re-check the three reference
  # problems (the gradient fixture, the `avc` early+constant fit, and the
  # synthetic two-phase fit in test-multiphase-reproducibility.R) before doing so.
  start_seed <- if (!is.null(control$start_seed)) control$start_seed else 3L
  control$start_seed <- NULL  # remove before passing to optim
  # Validated against what set.seed() can actually take, and rejected rather
  # than coerced.  Negative seeds are fine (set.seed(-1) is valid and
  # deterministic), but a non-whole number is not: set.seed() truncates it, so
  # 3.9 and 3 select the SAME ensemble, and a user sweeping 3.1/3.5/3.9 would
  # believe they had tried three sets of starts having tried one.  Out of
  # integer range, set.seed() fails with "supplied seed is not a valid integer"
  # and never names the argument that caused it.
  if (length(start_seed) != 1L || !is.numeric(start_seed) ||
      !is.finite(start_seed) || start_seed != trunc(start_seed) ||
      abs(start_seed) > .Machine$integer.max) {
    stop("'control$start_seed' must be a single whole number no larger in ",
         "magnitude than ", .Machine$integer.max, ".", call. = FALSE)
  }
  start_seed <- as.integer(start_seed)

  best_result <- NULL
  best_value <- -Inf
  best_start <- NA_integer_

  # Per-start outcomes.  The loop discards a start that errors, and without a
  # record of that an error raised inside the objective is indistinguishable
  # from a start that merely optimized poorly -- which is how a hard error in
  # the cumulative hazard was read as a convergence failure for as long as it
  # was.  Reported on the fit as `starts`.
  start_status  <- rep(NA_character_, n_starts)
  start_value   <- rep(NA_real_, n_starts)
  start_message <- rep(NA_character_, n_starts)

  # When shapes are fixed the free-parameter count is small (typically just the

  # mu intercepts).  BFGS can struggle on this reduced problem because the
  # likelihood surface is poorly scaled across phases (e.g. the late-phase mu
  # may be ~1e-8).  Strategy: run Nelder-Mead first (derivative-free, robust to
  # scaling) to find a good basin, then polish with BFGS for accurate SE.
  # Nelder-Mead requires >= 2 parameters; use Brent for 1D
  use_nelder_mead_warmup <- any_fixed && length(theta_start_optim) >= 2L &&
    length(theta_start_optim) <= 10L

  # Analytic Hessian closure -- covers event + right-censored rows.
  # Returns NULL for left/interval-censored data to trigger numDeriv fallback.
  # Defined OUTSIDE the n_starts loop: it is identical across all starts and
  # captures only stable (pre-loop) variables.
  hessian_fn_mp <- function(theta_free) {
    # theta_free is on the REDUCED (free-parameter) scale; expand to full.
    theta_full <- if (any_fixed) {
      th <- theta_fixed
      th[free_idx_eff] <- theta_free
      th
    } else {
      theta_free
    }
    H_full <- .hzr_hessian_multiphase(
      theta_full, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x, weights = weights,
      phases = phases, covariate_counts = covariate_counts, x_list = x_list
    )
    if (is.null(H_full)) return(NULL)
    # Restrict to the free-parameter block that the optimizer sees.
    # drop = FALSE keeps a single free parameter (e.g. a 2-phase fixed-shape
    # CoE fit, where CoE fixes one of the two log_mu) as a 1x1 matrix rather
    # than a scalar; otherwise .hzr_optim_generic rejects it as non-conformant
    # and needlessly discards the analytic Hessian for a numerical one.
    if (any_fixed) H_full[free_idx_eff, free_idx_eff, drop = FALSE] else H_full
  }

  # Drawn once, up front, so the optimization loop itself never touches the RNG.
  start_offsets <- .hzr_start_perturbations(
    n_starts, length(theta_start_optim), seed = start_seed
  )

  for (start_i in seq_len(n_starts)) {
    if (start_i == 1L) {
      theta_try <- theta_start_optim
    } else {
      # Deterministic perturbation of the starting values: the same data and
      # control give the same fit, and the caller's RNG stream is untouched.
      theta_try <- theta_start_optim + start_offsets[[start_i - 1L]]
    }

    # Nelder-Mead warm-up when shapes are fixed
    if (use_nelder_mead_warmup) {
      # Negative log-likelihood for minimization (matches .hzr_optim_generic)
      nm_objective <- function(theta) {
        if (!all(is.finite(theta))) return(1e10)
        ll <- -logl_fn(
          theta = theta, time = time, status = status,
          time_lower = time_lower, time_upper = time_upper,
          x = x, weights = weights, return_gradient = FALSE
        )
        if (!is.finite(ll)) return(1e10)
        ll
      }
      nm_result <- tryCatch(
        stats::optim(
          par = theta_try, fn = nm_objective,
          method = "Nelder-Mead",
          control = list(maxit = 5000L, reltol = 1e-10)
        ),
        error = function(e) NULL
      )
      if (!is.null(nm_result) && is.finite(nm_result$value)) {
        theta_try <- nm_result$par
      }
    }

    result <- tryCatch(
      .hzr_optim_generic(
        logl_fn     = logl_fn,
        gradient_fn = gradient_fn,
        time        = time, status = status,
        time_lower  = time_lower, time_upper = time_upper,
        x           = x,
        theta_start = theta_try,
        weights     = weights,
        control     = control,
        use_bounds  = FALSE,
        hessian_fn  = hessian_fn_mp
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      start_status[start_i]  <- "error"
      start_message[start_i] <- conditionMessage(result)
    } else if (!is.finite(result$value)) {
      start_status[start_i] <- "nonfinite"
    } else {
      start_status[start_i] <- "ok"
      start_value[start_i]  <- result$value
      if (result$value > best_value) {
        best_value  <- result$value
        best_result <- result
        best_start  <- start_i
      }
    }
  }

  starts <- data.frame(
    start     = seq_len(n_starts),
    status    = start_status,
    objective = start_value,
    best      = seq_len(n_starts) == best_start & !is.na(best_start),
    message   = start_message,
    stringsAsFactors = FALSE
  )

  if (is.null(best_result)) {
    # Name the errors.  "Failed to converge" is the right description of a start
    # that optimized badly and the wrong one for a start that threw, and the two
    # arrive here identically; saying which happened is the difference between a
    # numerical-tuning problem and a defect to go and find.
    errs <- unique(start_message[!is.na(start_message)])
    stop(
      "Multiphase optimization failed on all ", n_starts,
      if (n_starts == 1L) " start" else " starts", ": ",
      sum(start_status == "error", na.rm = TRUE), " errored, ",
      sum(start_status == "nonfinite", na.rm = TRUE),
      " returned a non-finite objective.",
      if (length(errs)) paste0(" Error(s) raised: ",
                               paste(errs, collapse = "; "), ".") else "",
      call. = FALSE
    )
  }

  # A start that throws is a defect, not a tuning problem, and the surviving
  # starts hide it: the fit comes back looking ordinary.  Measured on the
  # reference fixture, a healthy multi-start errors on none of its starts, so
  # this stays quiet unless something is actually wrong.  Losing start 1 to a
  # better perturbed optimum is routine and is reported in `starts`, not here.
  n_errored <- sum(start_status == "error", na.rm = TRUE)
  if (n_errored > 0L) {
    warning(
      n_errored, " of ", n_starts,
      if (n_starts == 1L) " optimization start" else " optimization starts",
      " errored and were discarded; the reported fit comes from start ",
      best_start, ". Error(s) raised: ",
      paste(unique(start_message[!is.na(start_message)]), collapse = "; "), ".",
      call. = FALSE
    )
  }

  # Expand optimized free params back to full theta vector
  if (any_fixed) {
    theta_full <- theta_fixed
    theta_full[free_idx] <- best_result$par
    best_result$par <- theta_full

    # Final CoE adjustment: solve fixmu log_mu at the optimum. Must pass
    # time_lower so the conserved log_mu in the returned coefficients is solved
    # on the same entry-time scale the objective was optimized on; otherwise a
    # left-truncated fit returns a conserved mu inconsistent with its objective.
    if (use_conserve && !is.null(fixmu_pos)) {
      best_result$par <- .hzr_conserve_events(
        best_result$par, fixmu_phase, fixmu_pos,
        time, status, phases, covariate_counts, x_list, total_events,
        weights = weights, time_lower = time_lower
      )
    }

    # Expand vcov to full dimension (NA for fixed params -- not estimated)
    if (!is.null(best_result$vcov) && is.matrix(best_result$vcov)) {
      p_full <- length(theta_fixed)
      vcov_full <- matrix(NA_real_, p_full, p_full)
      vcov_full[free_idx, free_idx] <- best_result$vcov
      best_result$vcov <- vcov_full
    }

    best_result$fixed_mask <- !free_mask

    # --- Full-information vcov for Conservation-of-Events fits ---------------
    # CoE removes the conserved phase's log_mu from the *search* (its score
    # equation IS the CoE constraint), so the optimizer's Hessian leaves that
    # parameter with no variance. But at the optimum the CoE solution is the
    # unconstrained MLE, so the reported *uncertainty* should be the full-
    # information vcov -- exactly what an all-mu-free (conserve = FALSE) fit
    # gives at the same point. Recompute the vcov from the unconstrained-
    # objective Hessian over the free set with the conserved log_mu restored,
    # so the conserved phase (and anything depending on it, e.g. se(H) and
    # survival CLs) gets a proper standard error.
    if (use_conserve && !is.null(fixmu_pos)) {
      recomputed <- FALSE
      free_unc <- free_mask
      free_unc[fixmu_pos] <- TRUE
      idx_unc <- which(free_unc)
      base_theta <- best_result$par
      neg_ll_unc <- function(th_free) {
        th <- base_theta
        th[idx_unc] <- th_free
        -logl_fn_pre_coe(th, time, status, time_lower, time_upper, x,
                         weights = weights, return_gradient = FALSE)
      }
      # Analytic Hessian attempt -- no numDeriv dependency
      H_unc <- tryCatch(
        .hzr_hessian_multiphase(
          base_theta, time = time, status = status,
          time_lower = time_lower, time_upper = time_upper,
          x = x, weights = weights,
          phases = phases, covariate_counts = covariate_counts,
          x_list = x_list
        ),
        error = function(e) NULL
      )
      if (is.matrix(H_unc)) {
        # Restrict to the unconstrained free set
        H_unc <- H_unc[idx_unc, idx_unc]
      }
      # numDeriv fallback only when analytic declines (left/interval rows)
      if (is.null(H_unc) && requireNamespace("numDeriv", quietly = TRUE)) {
        H_unc <- tryCatch(
          numDeriv::hessian(neg_ll_unc, base_theta[idx_unc]),
          error = function(e) NULL
        )
      }
      if (is.matrix(H_unc)) {
        inv_unc <- .hzr_safe_solve(H_unc)
        if (is.matrix(inv_unc$vcov)) {
          p_full <- length(base_theta)
          vcov_full <- matrix(NA_real_, p_full, p_full)
          vcov_full[idx_unc, idx_unc] <- inv_unc$vcov
          best_result$vcov  <- vcov_full
          best_result$rcond <- inv_unc$rcond
          best_result$pd    <- inv_unc$pd
          # The conserved log_mu now carries a variance; it is fixed only for
          # the search, not for inference.
          best_result$fixed_mask <- !free_unc
          recomputed <- TRUE
        }
      }
      if (!recomputed) {
        # Fall back to the reduced (search-only) vcov, in which the conserved
        # log_mu has no variance. Warn loudly so the understated uncertainty is
        # visible rather than silent.
        warning(
          "Could not compute the full-information variance for the ",
          "Conservation-of-Events-conserved phase 'log_mu' (numDeriv ",
          "unavailable or the Hessian was not invertible). Its standard error ",
          "stays NA and downstream standard errors / confidence limits for ",
          "that phase may be understated.", call. = FALSE)
      }
    }
  }

  # Restore parameter names
  names(best_result$par) <- theta_names

  # Store phase metadata for downstream use (predict, summary)
  best_result$phases <- phases
  best_result$covariate_counts <- covariate_counts
  best_result$x_list <- x_list

  # Which starts survived, and which one the reported fit came from.
  best_result$starts <- starts

  best_result
}


# ============================================================================
# G3 phase derivatives (central finite differences)
# ============================================================================

#' Finite-difference derivatives of G3 phase Phi and phi w.r.t. shape params
#'
#' Computes the G3 cumulative intensity, its time derivative, and their
#' partial derivatives with respect to log_tau, gamma, alpha, and eta using
#' central finite differences.
#'
#' @param time Numeric vector of positive times.
#' @param tau Positive scalar scale parameter.
#' @param gamma Positive scalar time exponent.
#' @param alpha Non-negative scalar shape parameter.
#' @param eta Positive scalar outer exponent.
#' @param h Relative step size for finite differences (default 1e-5).
#' @return Named list with Phi, phi, and 8 derivative vectors.
#' @keywords internal
.hzr_g3_phase_derivatives <- function(time, tau, gamma, alpha, eta,
                                       h = 1e-5) {

  # Base evaluation
  d0 <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                          alpha = alpha, eta = eta)
  Phi0 <- d0$G3
  phi0 <- d0$g3

  # Central differences for log_tau: tau * d/d(tau) = d/d(log(tau))
  eps_tau <- max(abs(tau) * h, 1e-10)
  if (tau - eps_tau > 0) {
    d_plus  <- hzr_decompos_g3(time, tau = tau + eps_tau, gamma = gamma,
                                 alpha = alpha, eta = eta)
    d_minus <- hzr_decompos_g3(time, tau = tau - eps_tau, gamma = gamma,
                                 alpha = alpha, eta = eta)
    dPhi_dtau <- (d_plus$G3 - d_minus$G3) / (2 * eps_tau)
    dphi_dtau <- (d_plus$g3 - d_minus$g3) / (2 * eps_tau)
  } else {
    d_plus  <- hzr_decompos_g3(time, tau = tau + eps_tau, gamma = gamma,
                                 alpha = alpha, eta = eta)
    dPhi_dtau <- (d_plus$G3 - Phi0) / eps_tau
    dphi_dtau <- (d_plus$g3 - phi0) / eps_tau
  }
  # Chain rule: d/d(log_tau) = tau * d/d(tau)
  dPhi_dlog_tau <- tau * dPhi_dtau
  dphi_dlog_tau <- tau * dphi_dtau

  # Central differences for gamma (must stay positive)
  eps_g <- max(abs(gamma) * h, 1e-10)
  if (gamma - eps_g > 0) {
    d_plus  <- hzr_decompos_g3(time, tau = tau, gamma = gamma + eps_g,
                                 alpha = alpha, eta = eta)
    d_minus <- hzr_decompos_g3(time, tau = tau, gamma = gamma - eps_g,
                                 alpha = alpha, eta = eta)
    dPhi_dgamma <- (d_plus$G3 - d_minus$G3) / (2 * eps_g)
    dphi_dgamma <- (d_plus$g3 - d_minus$g3) / (2 * eps_g)
  } else {
    d_plus  <- hzr_decompos_g3(time, tau = tau, gamma = gamma + eps_g,
                                 alpha = alpha, eta = eta)
    dPhi_dgamma <- (d_plus$G3 - Phi0) / eps_g
    dphi_dgamma <- (d_plus$g3 - phi0) / eps_g
  }

  # Central differences for alpha
  if (alpha > h) {
    eps_a <- max(abs(alpha) * h, 1e-10)
    d_plus  <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                                 alpha = alpha + eps_a, eta = eta)
    d_minus <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                                 alpha = alpha - eps_a, eta = eta)
    dPhi_dalpha <- (d_plus$G3 - d_minus$G3) / (2 * eps_a)
    dphi_dalpha <- (d_plus$g3 - d_minus$g3) / (2 * eps_a)
  } else {
    # alpha near 0: use forward difference
    eps_a <- h
    d_plus <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                                alpha = alpha + eps_a, eta = eta)
    dPhi_dalpha <- (d_plus$G3 - Phi0) / eps_a
    dphi_dalpha <- (d_plus$g3 - phi0) / eps_a
  }

  # Central differences for eta (must stay positive)
  eps_e <- max(abs(eta) * h, 1e-10)
  if (eta - eps_e > 0) {
    d_plus  <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                                 alpha = alpha, eta = eta + eps_e)
    d_minus <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                                 alpha = alpha, eta = eta - eps_e)
    dPhi_deta <- (d_plus$G3 - d_minus$G3) / (2 * eps_e)
    dphi_deta <- (d_plus$g3 - d_minus$g3) / (2 * eps_e)
  } else {
    d_plus  <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                                 alpha = alpha, eta = eta + eps_e)
    dPhi_deta <- (d_plus$G3 - Phi0) / eps_e
    dphi_deta <- (d_plus$g3 - phi0) / eps_e
  }

  list(
    Phi = Phi0,
    phi = phi0,
    dPhi_dlog_tau = dPhi_dlog_tau,
    dPhi_dgamma   = dPhi_dgamma,
    dPhi_dalpha   = dPhi_dalpha,
    dPhi_deta     = dPhi_deta,
    dphi_dlog_tau = dphi_dlog_tau,
    dphi_dgamma   = dphi_dgamma,
    dphi_dalpha   = dphi_dalpha,
    dphi_deta     = dphi_deta
  )
}


# ============================================================================
# Back-transformation helpers
# ============================================================================

#' Transform multiphase theta from internal to user-facing scale
#'
#' Converts log_mu -> mu, log_t_half -> t_half; nu, m, beta pass through.
#'
#' @param theta Full parameter vector (internal scale).
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector of per-phase covariate counts.
#' @return Named numeric vector on the reporting scale.
#' @keywords internal
.hzr_multiphase_theta_user <- function(theta, phases, covariate_counts) {
  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  values <- numeric(0)
  nms    <- character(0)

  for (nm in names(phases)) {
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])

    if (phases[[nm]]$type == "constant") {
      values <- c(values, exp(pars$log_mu))
      nms    <- c(nms, paste0(nm, ".mu"))
    } else if (phases[[nm]]$type == "g3") {
      values <- c(values, exp(pars$log_mu), exp(pars$log_tau),
                  pars$gamma, pars$alpha, pars$eta)
      nms    <- c(nms, paste0(nm, ".mu"), paste0(nm, ".tau"),
                  paste0(nm, ".gamma"), paste0(nm, ".alpha"), paste0(nm, ".eta"))
    } else {
      values <- c(values, exp(pars$log_mu), exp(pars$log_t_half),
                  pars$nu, pars$m)
      nms    <- c(nms, paste0(nm, ".mu"), paste0(nm, ".t_half"),
                  paste0(nm, ".nu"), paste0(nm, ".m"))
    }

    if (length(pars$beta) > 0) {
      values <- c(values, pars$beta)
      # Try to derive covariate names from the theta_split names
      n_skip <- switch(phases[[nm]]$type, constant = 1L, g3 = 5L, 4L)
      raw_names <- names(theta_split[[nm]])
      if (!is.null(raw_names) && length(raw_names) > n_skip) {
        beta_labels <- raw_names[(n_skip + 1L):length(raw_names)]
        # Strip phase prefix if present
        beta_labels <- sub(paste0("^", nm, "\\."), "", beta_labels)
      } else {
        beta_labels <- paste0("x", seq_along(pars$beta))
      }
      nms <- c(nms, paste0(nm, ".", beta_labels))
    }
  }

  names(values) <- nms
  values
}
