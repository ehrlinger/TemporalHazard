#' @keywords internal
NULL

# hessian-multiphase.R -- Analytic NLL Hessian for multiphase hazard models.
#
# Three helpers then the main function:
#   .hzr_phase_second_derivatives()     -- d2Phi/dss' and d2phi/dss' for CDF/hazard
#   .hzr_g3_phase_second_derivatives()  -- same for G3 (late) phases
#   .hzr_hessian_multiphase()           -- full p x p NLL Hessian

# Step size for 4-point second-difference formulas: h ~ eps^(1/4).
# Optimal for central second differences.
.hzr_h2 <- .Machine$double.eps^(1 / 4)

# ---------------------------------------------------------------------------
# CDF / hazard phase second derivatives
# ---------------------------------------------------------------------------

#' Second derivatives of CDF/hazard phase shape functions
#'
#' Computes the six unique mixed second partial derivatives of
#' \eqn{\Phi_j(t)} (cumulative shape) and \eqn{\phi_j(t)} (instantaneous
#' shape) with respect to the internal shape parameters
#' \eqn{\{log\_t\_half, \nu, m\}} via second-order central differences.
#'
#' This is the second-derivative analogue of \code{.hzr_phase_derivatives()}.
#'
#' @param time  Numeric vector of observation times.
#' @param t_half,nu,m  Phase shape parameters (natural scale).
#' @param type  Phase type: \code{"cdf"} or \code{"hazard"}.
#' @return A list with elements \code{d2Phi_dlog_thalf2},
#'   \code{d2Phi_dnu2}, \code{d2Phi_dm2},
#'   \code{d2Phi_dlog_thalf_dnu}, \code{d2Phi_dlog_thalf_dm},
#'   \code{d2Phi_dnu_dm}, and the matching \code{d2phi_*} entries.
#' @noRd
.hzr_phase_second_derivatives <- function(time, t_half, nu, m, type) {
  h <- .hzr_h2

  # Evaluate Phi and phi at given (t_half_p, nu_p, m_p)
  eval_at <- function(th, nup, mp) {
    d <- .hzr_phase_derivatives(time, t_half = th, nu = nup, m = mp, type = type)
    list(Phi = d$Phi, phi = d$phi)
  }

  d00 <- eval_at(t_half, nu, m)

  # hzr_decompos() rejects m < 0 && nu < 0 as mathematically undefined. A
  # fixed shape sitting exactly at m = 0 (nu < 0, Case 3L) or nu = 0 (m < 0,
  # Case 2L) is a valid, documented limiting case, but the central-difference
  # step below would probe the OTHER parameter's -h side, which crosses into
  # that undefined region. Detect it and fall back to a forward-only
  # (one-sided) difference in just that one direction; every other
  # direction/pair keeps the ordinary central-difference formula. (The two
  # conditions can't both hold: that would require m < 0 && nu < 0
  # simultaneously, which d00 above would already have rejected.)
  m_boundary  <- nu < 0 && abs(m)  < h
  nu_boundary <- m  < 0 && abs(nu) < h

  # Perturbed values on each internal scale:
  #   log_t_half: steps are multiplicative -> t_half * exp(+/-h)
  #   nu, m:      additive steps
  t_p  <- t_half * exp(h)
  t_m  <- t_half * exp(-h)
  nu_p <- nu + h
  nu_m <- if (nu_boundary) NA_real_ else nu - h
  m_p  <- m  + h
  m_m  <- if (m_boundary)  NA_real_ else m  - h

  # Six single-parameter perturbations for diagonal second differences
  # (skip the -h evaluation in a boundary direction -- it would error)
  d_tp  <- eval_at(t_p,    nu,   m)
  d_tm  <- eval_at(t_m,    nu,   m)
  d_np  <- eval_at(t_half, nu_p, m)
  d_nm  <- if (nu_boundary) NULL else eval_at(t_half, nu_m, m)
  d_mp  <- eval_at(t_half, nu,   m_p)
  d_mm  <- if (m_boundary)  NULL else eval_at(t_half, nu,   m_m)

  h2 <- h * h

  # Diagonal second differences. At a boundary, use the forward second
  # difference (f(x+2h) - 2f(x+h) + f(x)) / h^2 -- valid but O(h) rather
  # than O(h^2) locally, exactly at these two named limiting cases.
  d2Phi_tt <- (d_tp$Phi - 2 * d00$Phi + d_tm$Phi) / h2
  d2phi_tt <- (d_tp$phi - 2 * d00$phi + d_tm$phi) / h2

  if (nu_boundary) {
    d_np2 <- eval_at(t_half, nu + 2 * h, m)
    d2Phi_nn <- (d_np2$Phi - 2 * d_np$Phi + d00$Phi) / h2
    d2phi_nn <- (d_np2$phi - 2 * d_np$phi + d00$phi) / h2
  } else {
    d2Phi_nn <- (d_np$Phi - 2 * d00$Phi + d_nm$Phi) / h2
    d2phi_nn <- (d_np$phi - 2 * d00$phi + d_nm$phi) / h2
  }

  if (m_boundary) {
    d_mp2 <- eval_at(t_half, nu, m + 2 * h)
    d2Phi_mm_ <- (d_mp2$Phi - 2 * d_mp$Phi + d00$Phi) / h2
    d2phi_mm_ <- (d_mp2$phi - 2 * d_mp$phi + d00$phi) / h2
  } else {
    d2Phi_mm_ <- (d_mp$Phi - 2 * d00$Phi + d_mm$Phi) / h2
    d2phi_mm_ <- (d_mp$phi - 2 * d00$phi + d_mm$phi) / h2
  }

  # t_half / nu cross: corners (t_p, nu_p), (t_p, nu_m), (t_m, nu_p), (t_m, nu_m).
  # Under nu_boundary, use forward-in-nu / central-in-t_half instead (drops
  # the nu_m corners).
  if (nu_boundary) {
    # d_tp/d_tm (= eval_at(t_p/t_m, nu, m)) already computed above
    e_tn_pp <- eval_at(t_p, nu_p, m)
    e_tn_mp <- eval_at(t_m, nu_p, m)
    d2_tn_Phi <- ((e_tn_pp$Phi - d_tp$Phi) -
                    (e_tn_mp$Phi - d_tm$Phi)) / (2 * h2)
    d2_tn_phi <- ((e_tn_pp$phi - d_tp$phi) -
                    (e_tn_mp$phi - d_tm$phi)) / (2 * h2)
  } else {
    e_tn_pp <- eval_at(t_p, nu_p, m)
    e_tn_pm <- eval_at(t_p, nu_m, m)
    e_tn_mp <- eval_at(t_m, nu_p, m)
    e_tn_mm <- eval_at(t_m, nu_m, m)
    d2_tn_Phi <- (e_tn_pp$Phi - e_tn_pm$Phi - e_tn_mp$Phi + e_tn_mm$Phi) / (4 * h2)
    d2_tn_phi <- (e_tn_pp$phi - e_tn_pm$phi - e_tn_mp$phi + e_tn_mm$phi) / (4 * h2)
  }

  # t_half / m cross: corners (t_p, m_p), (t_p, m_m), (t_m, m_p), (t_m, m_m).
  # Under m_boundary, use forward-in-m / central-in-t_half instead (drops
  # the m_m corners).
  if (m_boundary) {
    # d_tp/d_tm (= eval_at(t_p/t_m, nu, m)) already computed above
    e_tm_pp <- eval_at(t_p, nu, m_p)
    e_tm_mp <- eval_at(t_m, nu, m_p)
    d2_tm_Phi <- ((e_tm_pp$Phi - d_tp$Phi) -
                    (e_tm_mp$Phi - d_tm$Phi)) / (2 * h2)
    d2_tm_phi <- ((e_tm_pp$phi - d_tp$phi) -
                    (e_tm_mp$phi - d_tm$phi)) / (2 * h2)
  } else {
    e_tm_pp <- eval_at(t_p, nu, m_p)
    e_tm_pm <- eval_at(t_p, nu, m_m)
    e_tm_mp <- eval_at(t_m, nu, m_p)
    e_tm_mm <- eval_at(t_m, nu, m_m)
    d2_tm_Phi <- (e_tm_pp$Phi - e_tm_pm$Phi - e_tm_mp$Phi + e_tm_mm$Phi) / (4 * h2)
    d2_tm_phi <- (e_tm_pp$phi - e_tm_pm$phi - e_tm_mp$phi + e_tm_mm$phi) / (4 * h2)
  }

  # nu / m cross: corners (nu_p, m_p), (nu_p, m_m), (nu_m, m_p), (nu_m, m_m).
  # Under m_boundary, use forward-in-m / central-in-nu (nu_m is valid here --
  # only m is boundary-restricted). Under nu_boundary, use forward-in-nu /
  # central-in-m (m_m is valid here -- only nu is boundary-restricted).
  if (m_boundary) {
    # d_np/d_nm (= eval_at(t_half, nu_p/nu_m, m)) already computed above
    e_nm_pp <- eval_at(t_half, nu_p, m_p)
    e_nm_mp <- eval_at(t_half, nu_m, m_p)
    d2_nm_Phi <- ((e_nm_pp$Phi - d_np$Phi) -
                    (e_nm_mp$Phi - d_nm$Phi)) / (2 * h2)
    d2_nm_phi <- ((e_nm_pp$phi - d_np$phi) -
                    (e_nm_mp$phi - d_nm$phi)) / (2 * h2)
  } else if (nu_boundary) {
    # d_mp/d_mm (= eval_at(t_half, nu, m_p/m_m)) already computed above
    e_nm_pp <- eval_at(t_half, nu_p, m_p)
    e_nm_pm <- eval_at(t_half, nu_p, m_m)
    d2_nm_Phi <- ((e_nm_pp$Phi - d_mp$Phi) -
                    (e_nm_pm$Phi - d_mm$Phi)) / (2 * h2)
    d2_nm_phi <- ((e_nm_pp$phi - d_mp$phi) -
                    (e_nm_pm$phi - d_mm$phi)) / (2 * h2)
  } else {
    e_nm_pp <- eval_at(t_half, nu_p, m_p)
    e_nm_pm <- eval_at(t_half, nu_p, m_m)
    e_nm_mp <- eval_at(t_half, nu_m, m_p)
    e_nm_mm <- eval_at(t_half, nu_m, m_m)
    d2_nm_Phi <- (e_nm_pp$Phi - e_nm_pm$Phi - e_nm_mp$Phi + e_nm_mm$Phi) / (4 * h2)
    d2_nm_phi <- (e_nm_pp$phi - e_nm_pm$phi - e_nm_mp$phi + e_nm_mm$phi) / (4 * h2)
  }

  list(
    d2Phi_dlog_thalf2    = d2Phi_tt,
    d2Phi_dnu2           = d2Phi_nn,
    d2Phi_dm2            = d2Phi_mm_,
    d2Phi_dlog_thalf_dnu = d2_tn_Phi,
    d2Phi_dlog_thalf_dm  = d2_tm_Phi,
    d2Phi_dnu_dm         = d2_nm_Phi,

    d2phi_dlog_thalf2    = d2phi_tt,
    d2phi_dnu2           = d2phi_nn,
    d2phi_dm2            = d2phi_mm_,
    d2phi_dlog_thalf_dnu = d2_tn_phi,
    d2phi_dlog_thalf_dm  = d2_tm_phi,
    d2phi_dnu_dm         = d2_nm_phi
  )
}

# ---------------------------------------------------------------------------
# G3 (late) phase second derivatives
# ---------------------------------------------------------------------------

#' Second derivatives of G3 phase shape functions
#'
#' Analogue of \code{.hzr_g3_phase_derivatives()} for second-order
#' shape derivatives over \code{\{log\_tau, gamma, alpha, eta\}}.
#'
#' @param time,tau,gamma,alpha,eta  G3 phase parameters.
#' @return List with \code{d2Phi_*} and \code{d2phi_*} entries for the ten
#'   unique pairs from \code{\{log\_tau, gamma, alpha, eta\}}.
#' @noRd
.hzr_g3_phase_second_derivatives <- function(time, tau, gamma, alpha, eta) {
  h <- .hzr_h2

  eval_at <- function(tau_p, gam_p, alp_p, eta_p) {
    # Guard against non-positive G3 parameters
    if (tau_p <= 0 || gam_p <= 0 || eta_p <= 0) return(NULL)
    tryCatch(
      {
        d <- hzr_decompos_g3(time, tau = tau_p, gamma = gam_p,
                              alpha = alp_p, eta = eta_p)
        list(Phi = d$G3, phi = d$g3)
      },
      error = function(e) NULL
    )
  }

  d00 <- eval_at(tau, gamma, alpha, eta)
  if (is.null(d00)) {
    n <- length(time)
    zeros <- rep(0, n)
    nms <- c("d2Phi_dlog_tau2", "d2Phi_dgamma2", "d2Phi_dalpha2",
             "d2Phi_deta2",
             "d2Phi_dlog_tau_dgamma", "d2Phi_dlog_tau_dalpha",
             "d2Phi_dlog_tau_deta",   "d2Phi_dgamma_dalpha",
             "d2Phi_dgamma_deta",     "d2Phi_dalpha_deta",
             "d2phi_dlog_tau2", "d2phi_dgamma2", "d2phi_dalpha2",
             "d2phi_deta2",
             "d2phi_dlog_tau_dgamma", "d2phi_dlog_tau_dalpha",
             "d2phi_dlog_tau_deta",   "d2phi_dgamma_dalpha",
             "d2phi_dgamma_deta",     "d2phi_dalpha_deta")
    return(setNames(replicate(length(nms), zeros, simplify = FALSE), nms))
  }

  # Perturbed values: log_tau, gamma, alpha, eta steps
  tau_p <- tau * exp(h)
  tau_m <- tau * exp(-h)
  gam_p <- gamma + h
  gam_m <- gamma - h
  alp_p <- alpha + h
  alp_m <- alpha - h
  eta_p <- eta   + h
  eta_m <- eta   - h

  # Fallback to d00 if perturbed eval returns NULL
  safe_eval <- function(tau_v, gam_v, alp_v, eta_v) {
    r <- eval_at(tau_v, gam_v, alp_v, eta_v)
    if (is.null(r)) d00 else r
  }

  d_Tp <- safe_eval(tau_p, gamma, alpha, eta)
  d_Tm <- safe_eval(tau_m, gamma, alpha, eta)
  d_Gp <- safe_eval(tau, gam_p, alpha, eta)
  d_Gm <- safe_eval(tau, gam_m, alpha, eta)
  d_Ap <- safe_eval(tau, gamma, alp_p, eta)
  d_Am <- safe_eval(tau, gamma, alp_m, eta)
  d_Ep <- safe_eval(tau, gamma, alpha, eta_p)
  d_Em <- safe_eval(tau, gamma, alpha, eta_m)

  h2 <- h * h

  # Diagonal second differences
  diag2 <- function(dp, dm, d0, which) {
    (dp[[which]] - 2 * d0[[which]] + dm[[which]]) / h2
  }

  # tau / gamma
  d_TpGp <- safe_eval(tau_p, gam_p, alpha, eta)
  d_TpGm <- safe_eval(tau_p, gam_m, alpha, eta)
  d_TmGp <- safe_eval(tau_m, gam_p, alpha, eta)
  d_TmGm <- safe_eval(tau_m, gam_m, alpha, eta)
  TG_Phi <- (d_TpGp$Phi - d_TpGm$Phi - d_TmGp$Phi + d_TmGm$Phi) / (4 * h2)
  TG_phi <- (d_TpGp$phi - d_TpGm$phi - d_TmGp$phi + d_TmGm$phi) / (4 * h2)

  # tau / alpha
  d_TpAp <- safe_eval(tau_p, gamma, alp_p, eta)
  d_TpAm <- safe_eval(tau_p, gamma, alp_m, eta)
  d_TmAp <- safe_eval(tau_m, gamma, alp_p, eta)
  d_TmAm <- safe_eval(tau_m, gamma, alp_m, eta)
  TA_Phi <- (d_TpAp$Phi - d_TpAm$Phi - d_TmAp$Phi + d_TmAm$Phi) / (4 * h2)
  TA_phi <- (d_TpAp$phi - d_TpAm$phi - d_TmAp$phi + d_TmAm$phi) / (4 * h2)

  # tau / eta
  d_TpEp <- safe_eval(tau_p, gamma, alpha, eta_p)
  d_TpEm <- safe_eval(tau_p, gamma, alpha, eta_m)
  d_TmEp <- safe_eval(tau_m, gamma, alpha, eta_p)
  d_TmEm <- safe_eval(tau_m, gamma, alpha, eta_m)
  TE_Phi <- (d_TpEp$Phi - d_TpEm$Phi - d_TmEp$Phi + d_TmEm$Phi) / (4 * h2)
  TE_phi <- (d_TpEp$phi - d_TpEm$phi - d_TmEp$phi + d_TmEm$phi) / (4 * h2)

  # gamma / alpha
  d_GpAp <- safe_eval(tau, gam_p, alp_p, eta)
  d_GpAm <- safe_eval(tau, gam_p, alp_m, eta)
  d_GmAp <- safe_eval(tau, gam_m, alp_p, eta)
  d_GmAm <- safe_eval(tau, gam_m, alp_m, eta)
  GA_Phi <- (d_GpAp$Phi - d_GpAm$Phi - d_GmAp$Phi + d_GmAm$Phi) / (4 * h2)
  GA_phi <- (d_GpAp$phi - d_GpAm$phi - d_GmAp$phi + d_GmAm$phi) / (4 * h2)

  # gamma / eta
  d_GpEp <- safe_eval(tau, gam_p, alpha, eta_p)
  d_GpEm <- safe_eval(tau, gam_p, alpha, eta_m)
  d_GmEp <- safe_eval(tau, gam_m, alpha, eta_p)
  d_GmEm <- safe_eval(tau, gam_m, alpha, eta_m)
  GE_Phi <- (d_GpEp$Phi - d_GpEm$Phi - d_GmEp$Phi + d_GmEm$Phi) / (4 * h2)
  GE_phi <- (d_GpEp$phi - d_GpEm$phi - d_GmEp$phi + d_GmEm$phi) / (4 * h2)

  # alpha / eta
  d_ApEp <- safe_eval(tau, gamma, alp_p, eta_p)
  d_ApEm <- safe_eval(tau, gamma, alp_p, eta_m)
  d_AmEp <- safe_eval(tau, gamma, alp_m, eta_p)
  d_AmEm <- safe_eval(tau, gamma, alp_m, eta_m)
  AE_Phi <- (d_ApEp$Phi - d_ApEm$Phi - d_AmEp$Phi + d_AmEm$Phi) / (4 * h2)
  AE_phi <- (d_ApEp$phi - d_ApEm$phi - d_AmEp$phi + d_AmEm$phi) / (4 * h2)

  list(
    d2Phi_dlog_tau2        = diag2(d_Tp, d_Tm, d00, "Phi"),
    d2Phi_dgamma2          = diag2(d_Gp, d_Gm, d00, "Phi"),
    d2Phi_dalpha2          = diag2(d_Ap, d_Am, d00, "Phi"),
    d2Phi_deta2            = diag2(d_Ep, d_Em, d00, "Phi"),
    d2Phi_dlog_tau_dgamma  = TG_Phi,
    d2Phi_dlog_tau_dalpha  = TA_Phi,
    d2Phi_dlog_tau_deta    = TE_Phi,
    d2Phi_dgamma_dalpha    = GA_Phi,
    d2Phi_dgamma_deta      = GE_Phi,
    d2Phi_dalpha_deta      = AE_Phi,

    d2phi_dlog_tau2        = diag2(d_Tp, d_Tm, d00, "phi"),
    d2phi_dgamma2          = diag2(d_Gp, d_Gm, d00, "phi"),
    d2phi_dalpha2          = diag2(d_Ap, d_Am, d00, "phi"),
    d2phi_deta2            = diag2(d_Ep, d_Em, d00, "phi"),
    d2phi_dlog_tau_dgamma  = TG_phi,
    d2phi_dlog_tau_dalpha  = TA_phi,
    d2phi_dlog_tau_deta    = TE_phi,
    d2phi_dgamma_dalpha    = GA_phi,
    d2phi_dgamma_deta      = GE_phi,
    d2phi_dalpha_deta      = AE_phi
  )
}

# ---------------------------------------------------------------------------
# Main analytic Hessian function
# ---------------------------------------------------------------------------
#' Analytic Hessian of the multiphase NLL
#'
#' Returns the \eqn{p \times p} Hessian of the **objective** (negative
#' log-likelihood) on the internal parameter scale, matching
#' \code{numDeriv::hessian(objective)}.  Coverage: event + right-censored
#' rows, with counting-process start-time corrections when
#' \code{time_lower > 0}. Returns \code{NULL} for data containing
#' left-censored (\code{status == -1}) or interval-censored
#' (\code{status == 2}) rows; the caller falls back to numDeriv in that case.
#'
#' **Structure** (see derivation in PLAN-hessian-analytic-multiphase.md):
#' \deqn{H = (\text{block-diagonal } A) + (\text{dense outer-product } B) + (\text{block-diagonal } C)}
#'
#' Term A: curvature of \eqn{\sum_i w_i H_i}; block-diagonal in phases.
#' Term B: \eqn{\sum_{i\in E} (w_i/h_i^2) \nabla h_i \nabla h_i^T}; dense.
#' Term C: curvature of \eqn{-\sum_{i\in E} w_i \log h_i}; block-diagonal.
#'
#' @noRd
.hzr_hessian_multiphase <- function(
    theta, time, status,
    time_lower = NULL, time_upper = NULL,
    x = NULL, weights = NULL,
    phases, covariate_counts, x_list) {

  n <- length(time)
  p <- length(theta)
  if (is.null(weights)) weights <- rep(1, n)

  # Coverage contract: decline for left/interval-censored rows
  if (any(status %in% c(-1L, 2L))) return(NULL)

  idx_event <- which(status == 1L)
  n_e       <- length(idx_event)

  # Counting-process start-time handling (mirrors gradient)
  need_start <- !is.null(time_lower) &&
                any(time_lower > 0 & time_lower < time & status %in% c(0L, 1L))
  start_vec <- if (need_start) {
    sv <- rep(0, n)
    epoch_idx <- status %in% c(0L, 1L) & time_lower < time
    sv[epoch_idx] <- time_lower[epoch_idx]
    sv
  } else {
    NULL
  }

  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)

  # --- Per-phase quantities --------------------------------------------------
  n_phases        <- length(phases)
  phase_mu        <- vector("list", n_phases)
  phase_Phi       <- vector("list", n_phases)
  phase_phi       <- vector("list", n_phases)
  phase_d1        <- vector("list", n_phases)   # .hzr_phase_derivatives output
  phase_d2        <- vector("list", n_phases)   # second-derivative output
  phase_Phi_start <- vector("list", n_phases)
  phase_d1_start  <- vector("list", n_phases)
  phase_d2_start  <- vector("list", n_phases)
  phase_x_tilde   <- vector("list", n_phases)   # [1 | x_j], n x (1+p_j)

  H_t <- rep(0, n)
  h_t <- rep(0, n)

  for (j in seq_along(phases)) {
    nm   <- names(phases)[j]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])

    # mu_j(i) and augmented covariate matrix
    if (length(pars$beta) > 0 && !is.null(x_list[[nm]])) {
      x_j  <- x_list[[nm]]
      eta_j <- pars$log_mu + as.numeric(x_j %*% pars$beta)
    } else {
      x_j   <- NULL
      eta_j <- rep(pars$log_mu, n)
    }
    mu_j <- exp(eta_j)
    phase_mu[[j]] <- mu_j

    x_tilde_j <- if (is.null(x_j)) matrix(1, n, 1L) else cbind(1, x_j)
    phase_x_tilde[[j]] <- x_tilde_j

    # Shape quantities
    if (phases[[nm]]$type == "constant") {
      phase_Phi[[j]] <- time
      phase_phi[[j]] <- rep(1, n)
      phase_d1[j]  <- list(NULL)
      phase_d2[j]  <- list(NULL)
      if (need_start) {
        phase_Phi_start[[j]] <- start_vec
        phase_d1_start[j]  <- list(NULL)
        phase_d2_start[j]  <- list(NULL)
      }
    } else if (phases[[nm]]$type == "g3") {
      tau_j <- exp(pars$log_tau)
      d3    <- hzr_decompos_g3(time, tau = tau_j, gamma = pars$gamma,
                                alpha = pars$alpha, eta = pars$eta)
      phase_Phi[[j]] <- d3$G3
      phase_phi[[j]] <- d3$g3
      phase_d1[[j]]  <- .hzr_g3_phase_derivatives(
        time, tau = tau_j, gamma = pars$gamma,
        alpha = pars$alpha, eta = pars$eta)
      phase_d2[[j]]  <- .hzr_g3_phase_second_derivatives(
        time, tau = tau_j, gamma = pars$gamma,
        alpha = pars$alpha, eta = pars$eta)
      if (need_start) {
        d3s <- hzr_decompos_g3(start_vec, tau = tau_j, gamma = pars$gamma,
                                alpha = pars$alpha, eta = pars$eta)
        phase_Phi_start[[j]] <- d3s$G3
        phase_d1_start[[j]]  <- .hzr_g3_phase_derivatives(
          start_vec, tau = tau_j, gamma = pars$gamma,
          alpha = pars$alpha, eta = pars$eta)
        phase_d2_start[[j]]  <- .hzr_g3_phase_second_derivatives(
          start_vec, tau = tau_j, gamma = pars$gamma,
          alpha = pars$alpha, eta = pars$eta)
      }
    } else {
      # CDF / hazard
      t_half_j <- exp(pars$log_t_half)
      d1 <- .hzr_phase_derivatives(time, t_half = t_half_j, nu = pars$nu,
                                    m = pars$m, type = phases[[nm]]$type)
      phase_Phi[[j]] <- d1$Phi
      phase_phi[[j]] <- d1$phi
      phase_d1[[j]]  <- d1
      phase_d2[[j]]  <- .hzr_phase_second_derivatives(
        time, t_half = t_half_j, nu = pars$nu, m = pars$m,
        type = phases[[nm]]$type)
      if (need_start) {
        d1s <- .hzr_phase_derivatives(start_vec, t_half = t_half_j,
                                       nu = pars$nu, m = pars$m,
                                       type = phases[[nm]]$type)
        phase_Phi_start[[j]] <- d1s$Phi
        phase_d1_start[[j]]  <- d1s
        phase_d2_start[[j]]  <- .hzr_phase_second_derivatives(
          start_vec, t_half = t_half_j, nu = pars$nu, m = pars$m,
          type = phases[[nm]]$type)
      }
    }

    H_t <- H_t + mu_j * phase_Phi[[j]]
    h_t <- h_t + mu_j * phase_phi[[j]]
    if (need_start) {
      H_t <- H_t - mu_j * phase_Phi_start[[j]]  # H(t) - H(start) net
    }
  }

  if (any(!is.finite(H_t)) || any(!is.finite(h_t))) {
    return(NULL)
  }

  # h_e = h at event times (guard against 0)
  h_e <- pmax(h_t[idx_event], .Machine$double.xmin)

  # w/h_e^2 and w/h_e vectors (length n_e)
  w_e     <- weights[idx_event]
  wh2_e   <- w_e / h_e^2     # for Term B
  wh1_e   <- w_e / h_e       # for Term C

  # --- Build the full p x p Hessian matrix ----------------------------------
  H_mat <- matrix(0, p, p)
  dimnames(H_mat) <- list(names(theta), names(theta))

  # We need a per-phase parameter offset so we know where in H_mat to write.
  # Build an offset table.
  phase_param_counts <- vapply(seq_along(phases), function(j) {
    nm <- names(phases)[j]
    if (phases[[nm]]$type == "constant") {
      1L + covariate_counts[[nm]]
    } else if (phases[[nm]]$type == "g3") {
      5L + covariate_counts[[nm]]  # log_mu, log_tau, gamma, alpha, eta
    } else {
      4L + covariate_counts[[nm]]  # log_mu, log_t_half, nu, m
    }
  }, integer(1))

  phase_offsets <- c(0L, cumsum(phase_param_counts))

  # --- Term B: dense outer product of h-gradient ----------------------------
  # Assemble X_h (n_e x p): row i is the gradient of h(t_i) w.r.t. theta.
  X_h <- matrix(0, n_e, p)

  for (j in seq_along(phases)) {
    nm   <- names(phases)[j]
    pos  <- phase_offsets[j]
    mu_j <- phase_mu[[j]]
    phi_j <- phase_phi[[j]]
    d1    <- phase_d1[[j]]

    mu_phi_e <- mu_j[idx_event] * phi_j[idx_event]  # length n_e

    # log_mu column: dh/d(log_mu_j) = mu_j phi_j
    X_h[, pos + 1L] <- mu_phi_e

    # Shape columns
    if (phases[[nm]]$type %in% c("cdf", "hazard")) {
      shape_dphi <- list(d1$dphi_dlog_thalf, d1$dphi_dnu, d1$dphi_dm)
      for (s in seq_along(shape_dphi)) {
        X_h[, pos + 1L + s] <- mu_j[idx_event] * shape_dphi[[s]][idx_event]
      }
      beta_start <- pos + 5L
    } else if (phases[[nm]]$type == "g3") {
      shape_dphi <- list(d1$dphi_dlog_tau, d1$dphi_dgamma,
                         d1$dphi_dalpha,   d1$dphi_deta)
      for (s in seq_along(shape_dphi)) {
        X_h[, pos + 1L + s] <- mu_j[idx_event] * shape_dphi[[s]][idx_event]
      }
      beta_start <- pos + 6L
    } else {
      # constant phase: no shape cols
      beta_start <- pos + 2L
    }

    # Beta columns: dh/d(beta_jk) = mu_j phi_j x_jk
    if (covariate_counts[[nm]] > 0 && !is.null(x_list[[nm]])) {
      x_j_e <- x_list[[nm]][idx_event, , drop = FALSE]
      for (k in seq_len(covariate_counts[[nm]])) {
        X_h[, beta_start + k - 1L] <- mu_phi_e * x_j_e[, k]
      }
    }
  }

  # Term B = X_h' diag(wh2_e) X_h (each row of X_h weighted by wh2_e)
  B_mat <- crossprod(X_h, wh2_e * X_h)
  H_mat <- H_mat + B_mat

  # --- Terms A and C: per-phase block contributions -------------------------
  for (j in seq_along(phases)) {
    nm   <- names(phases)[j]
    pos  <- phase_offsets[j]
    mu_j <- phase_mu[[j]]
    Phi_j <- phase_Phi[[j]]
    phi_j <- phase_phi[[j]]
    d1    <- phase_d1[[j]]
    d2    <- phase_d2[[j]]
    x_t   <- phase_x_tilde[[j]]

    wmu_Phi   <- weights * mu_j * Phi_j               # Term A weight (all rows)
    wmu_phi_e <- wh1_e * mu_j[idx_event] * phi_j[idx_event]  # Term C weight (events)

    # Counting-process start correction for Term A:
    if (need_start) {
      wmu_Phi_start <- weights * mu_j * phase_Phi_start[[j]]
      wmu_Phi_net <- wmu_Phi - wmu_Phi_start   # net for H(t) - H(start)
    } else {
      wmu_Phi_net <- wmu_Phi
    }

    # --- mu/beta sub-block (crossprod style) ---------------------------------
    # Term A mu/beta: X_t' diag(wmu_Phi_net) X_t
    # Term C mu/beta: - X_t_E' diag(wmu_phi_e) X_t_E
    x_t_e <- x_t[idx_event, , drop = FALSE]

    block_A_mub <- crossprod(x_t, wmu_Phi_net * x_t)
    block_C_mub <- crossprod(x_t_e, wmu_phi_e * x_t_e)

    # Number of mu/beta parameters: 1 + p_j
    n_mub <- 1L + covariate_counts[[nm]]

    if (phases[[nm]]$type == "constant") {
      # Only mu/beta block; no shape params
      idx_block <- pos + seq_len(n_mub)
      H_mat[idx_block, idx_block] <- H_mat[idx_block, idx_block] +
        block_A_mub - block_C_mub

    } else {
      # Shape params come AFTER log_mu, BEFORE beta (log_mu | shapes | beta)
      n_shape <- if (phases[[nm]]$type == "g3") 4L else 3L

      idx_mub   <- pos + 1L                             # log_mu index
      idx_shape <- pos + 1L + seq_len(n_shape)          # shape indices
      idx_beta  <- if (covariate_counts[[nm]] > 0)
                     pos + 1L + n_shape + seq_len(covariate_counts[[nm]])
                   else integer(0)

      # log_mu / log_mu (scalar)
      H_mat[idx_mub, idx_mub] <- H_mat[idx_mub, idx_mub] +
        sum(wmu_Phi_net) - sum(wmu_phi_e)

      # log_mu / beta cross terms (only if covariates present)
      if (length(idx_beta) > 0) {
        x_j <- x_list[[nm]]
        x_j_e <- x_j[idx_event, , drop = FALSE]
        cross_A <- colSums(wmu_Phi_net * x_j)
        cross_C <- colSums(wmu_phi_e * x_j_e)
        H_mat[idx_mub, idx_beta] <- H_mat[idx_mub, idx_beta] + cross_A - cross_C
        H_mat[idx_beta, idx_mub] <- H_mat[idx_beta, idx_mub] + cross_A - cross_C
      }

      # beta / beta (only if covariates)
      if (length(idx_beta) > 0) {
        x_j   <- x_list[[nm]]
        x_j_e <- x_j[idx_event, , drop = FALSE]
        bb_A  <- crossprod(x_j,   wmu_Phi_net * x_j)
        bb_C  <- crossprod(x_j_e, wmu_phi_e   * x_j_e)
        H_mat[idx_beta, idx_beta] <- H_mat[idx_beta, idx_beta] + bb_A - bb_C
      }

      # Shape parameters: extract first-derivative vectors
      if (phases[[nm]]$type %in% c("cdf", "hazard")) {
        shape_keys_d1_Phi <- c("dPhi_dlog_thalf", "dPhi_dnu", "dPhi_dm")
        shape_keys_d1_phi <- c("dphi_dlog_thalf", "dphi_dnu", "dphi_dm")
        shape_keys_d2_Phi <- list(
          c("d2Phi_dlog_thalf2",    "d2Phi_dlog_thalf2"),
          c("d2Phi_dnu2",           "d2Phi_dnu2"),
          c("d2Phi_dm2",            "d2Phi_dm2"),
          c("d2Phi_dlog_thalf_dnu", "d2Phi_dlog_thalf_dnu"),
          c("d2Phi_dlog_thalf_dm",  "d2Phi_dlog_thalf_dm"),
          c("d2Phi_dnu_dm",         "d2Phi_dnu_dm")
        )
        shape_keys_d2_phi <- list(
          c("d2phi_dlog_thalf2",    "d2phi_dlog_thalf2"),
          c("d2phi_dnu2",           "d2phi_dnu2"),
          c("d2phi_dm2",            "d2phi_dm2"),
          c("d2phi_dlog_thalf_dnu", "d2phi_dlog_thalf_dnu"),
          c("d2phi_dlog_thalf_dm",  "d2phi_dlog_thalf_dm"),
          c("d2phi_dnu_dm",         "d2phi_dnu_dm")
        )
        # Which pairs of shape indices correspond to the 6 unique pairs?
        shape_pairs <- list(c(1, 1), c(2, 2), c(3, 3), c(1, 2), c(1, 3), c(2, 3))
      } else {
        # G3: 4 shape params -> 10 unique pairs
        shape_keys_d1_Phi <- c("dPhi_dlog_tau", "dPhi_dgamma",
                                "dPhi_dalpha",   "dPhi_deta")
        shape_keys_d1_phi <- c("dphi_dlog_tau", "dphi_dgamma",
                                "dphi_dalpha",   "dphi_deta")
        shape_keys_d2_Phi <- list(
          c("d2Phi_dlog_tau2",       "d2Phi_dlog_tau2"),
          c("d2Phi_dgamma2",         "d2Phi_dgamma2"),
          c("d2Phi_dalpha2",         "d2Phi_dalpha2"),
          c("d2Phi_deta2",           "d2Phi_deta2"),
          c("d2Phi_dlog_tau_dgamma", "d2Phi_dlog_tau_dgamma"),
          c("d2Phi_dlog_tau_dalpha", "d2Phi_dlog_tau_dalpha"),
          c("d2Phi_dlog_tau_deta",   "d2Phi_dlog_tau_deta"),
          c("d2Phi_dgamma_dalpha",   "d2Phi_dgamma_dalpha"),
          c("d2Phi_dgamma_deta",     "d2Phi_dgamma_deta"),
          c("d2Phi_dalpha_deta",     "d2Phi_dalpha_deta")
        )
        shape_keys_d2_phi <- lapply(shape_keys_d2_Phi, function(k) {
          sub("Phi", "phi", k)
        })
        shape_pairs <- c(
          list(c(1, 1), c(2, 2), c(3, 3), c(4, 4)),
          list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4))
        )
      }

      n_shape <- length(shape_keys_d1_Phi)

      # log_mu / shape cross terms (Term A and C)
      for (s in seq_len(n_shape)) {
        dPhi_ds <- d1[[shape_keys_d1_Phi[s]]]
        dphi_ds <- d1[[shape_keys_d1_phi[s]]]
        if (need_start && !is.null(phase_d1_start[[j]])) {
          dPhi_ds_start <- phase_d1_start[[j]][[shape_keys_d1_Phi[s]]]
          a_val <- sum(weights * mu_j * dPhi_ds) -
                   sum(weights * mu_j * dPhi_ds_start)
        } else {
          a_val <- sum(weights * mu_j * dPhi_ds)
        }
        c_val <- sum(wh1_e * mu_j[idx_event] * dphi_ds[idx_event])
        col_s <- idx_shape[s]
        H_mat[idx_mub, col_s] <- H_mat[idx_mub, col_s] + a_val - c_val
        H_mat[col_s, idx_mub] <- H_mat[col_s, idx_mub] + a_val - c_val
      }

      # shape / shape block
      for (pair_idx in seq_along(shape_pairs)) {
        s_a <- shape_pairs[[pair_idx]][1]
        s_b <- shape_pairs[[pair_idx]][2]
        key_d2_Phi <- shape_keys_d2_Phi[[pair_idx]][1]
        key_d2_phi <- shape_keys_d2_phi[[pair_idx]][1]

        d2Phi_sasb <- d2[[key_d2_Phi]]
        d2phi_sasb <- d2[[key_d2_phi]]

        a_val <- sum(weights * mu_j * d2Phi_sasb)
        if (need_start && !is.null(phase_d2_start[[j]])) {
          a_val <- a_val - sum(weights * mu_j * phase_d2_start[[j]][[key_d2_Phi]])
        }
        c_val <- sum(wh1_e * mu_j[idx_event] * d2phi_sasb[idx_event])

        i_idx <- idx_shape[s_a]
        j_idx <- idx_shape[s_b]
        H_mat[i_idx, j_idx] <- H_mat[i_idx, j_idx] + a_val - c_val
        if (s_a != s_b) {
          H_mat[j_idx, i_idx] <- H_mat[j_idx, i_idx] + a_val - c_val
        }
      }

      # beta / shape cross terms
      if (length(idx_beta) > 0) {
        x_j   <- x_list[[nm]]
        x_j_e <- x_j[idx_event, , drop = FALSE]
        for (s in seq_len(n_shape)) {
          dPhi_ds <- d1[[shape_keys_d1_Phi[s]]]
          dphi_ds <- d1[[shape_keys_d1_phi[s]]]
          if (need_start && !is.null(phase_d1_start[[j]])) {
            dPhi_ds_start <- phase_d1_start[[j]][[shape_keys_d1_Phi[s]]]
            a_vec <- colSums((weights * mu_j * dPhi_ds -
                              weights * mu_j * dPhi_ds_start) * x_j)
          } else {
            a_vec <- colSums((weights * mu_j * dPhi_ds) * x_j)
          }
          c_vec <- colSums((wh1_e * mu_j[idx_event] *
                              dphi_ds[idx_event]) * x_j_e)
          col_s <- idx_shape[s]
          H_mat[idx_beta, col_s] <- H_mat[idx_beta, col_s] + a_vec - c_vec
          H_mat[col_s, idx_beta] <- H_mat[col_s, idx_beta] + a_vec - c_vec
        }
      }
    }
  }

  # Guard: any non-finite entry means the Hessian is unreliable; return NULL so
  # the caller falls back to numDeriv rather than silently passing a corrupted matrix.
  if (any(!is.finite(H_mat))) return(NULL)

  H_mat
}
