# test-multiphase-hessian.R
# Analytic Hessian tests for multiphase models.
# Validates shape second-derivative helpers and .hzr_hessian_multiphase()
# against numDeriv references.

# ---------------------------------------------------------------------------
# Task 1: Shape second-derivative helpers
# ---------------------------------------------------------------------------

test_that(".hzr_phase_second_derivatives diagonals match numDeriv", {
  skip_if_not_installed("numDeriv")
  t <- c(0.2, 0.5, 1.0, 2.0, 4.0)
  t_half <- 1.0
  nu <- 1.5
  m <- 1.0

  d2 <- .hzr_phase_second_derivatives(t, t_half = t_half, nu = nu, m = m,
                                       type = "cdf")

  # Reference: numDeriv Hessian of Phi_i w.r.t. (log_t_half, nu, m)
  par0 <- c(log(t_half), nu, m)
  for (i in seq_along(t)) {
    phi_fn_i <- function(par) {
      hzr_phase_cumhaz(t[i], t_half = exp(par[1]), nu = par[2], m = par[3],
                       type = "cdf")
    }
    H_nd <- numDeriv::hessian(phi_fn_i, par0)
    expect_equal(d2$d2Phi_dlog_thalf2[i],       H_nd[1, 1], tolerance = 1e-3,
                 label = paste0("d2Phi/d(log_thalf)^2 at t=", t[i]))
    expect_equal(d2$d2Phi_dnu2[i],              H_nd[2, 2], tolerance = 1e-3,
                 label = paste0("d2Phi/dnu^2 at t=", t[i]))
    expect_equal(d2$d2Phi_dm2[i],               H_nd[3, 3], tolerance = 1e-3,
                 label = paste0("d2Phi/dm^2 at t=", t[i]))
    expect_equal(d2$d2Phi_dlog_thalf_dnu[i],    H_nd[1, 2], tolerance = 1e-3,
                 label = paste0("d2Phi/d(log_thalf)dnu at t=", t[i]))
    expect_equal(d2$d2Phi_dlog_thalf_dm[i],     H_nd[1, 3], tolerance = 1e-3,
                 label = paste0("d2Phi/d(log_thalf)dm at t=", t[i]))
    expect_equal(d2$d2Phi_dnu_dm[i],            H_nd[2, 3], tolerance = 1e-3,
                 label = paste0("d2Phi/dnudm at t=", t[i]))
  }
})

test_that(".hzr_phase_second_derivatives phi second derivs match numDeriv", {
  skip_if_not_installed("numDeriv")
  t <- c(0.3, 0.8, 1.5, 3.0)
  t_half <- 0.5
  nu <- 2.0
  m <- 1.0
  d2 <- .hzr_phase_second_derivatives(t, t_half = t_half, nu = nu, m = m,
                                       type = "cdf")
  par0 <- c(log(t_half), nu, m)
  for (i in seq_along(t)) {
    haz_fn_i <- function(par) {
      hzr_phase_hazard(t[i], t_half = exp(par[1]), nu = par[2], m = par[3],
                       type = "cdf")
    }
    H_nd <- numDeriv::hessian(haz_fn_i, par0)
    expect_equal(d2$d2phi_dlog_thalf2[i], H_nd[1, 1], tolerance = 1e-3,
                 label = paste0("d2phi/d(log_thalf)^2 at t=", t[i]))
    expect_equal(d2$d2phi_dnu2[i],        H_nd[2, 2], tolerance = 1e-3,
                 label = paste0("d2phi/dnu^2 at t=", t[i]))
  }
})

# ---------------------------------------------------------------------------
# Task 1b: Case-boundary robustness (Case 3L: m=0,nu<0; Case 2L: nu=0,m<0)
# ---------------------------------------------------------------------------
# hzr_decompos() rejects m<0 && nu<0 as mathematically undefined. A shape
# sitting exactly at m=0 (with nu<0) or nu=0 (with m<0) is a legitimate,
# documented limiting case (Case 3L / Case 2L respectively), but the naive
# central-difference perturbation used below steps the OTHER side of the
# fixed value by -h, which crosses into that undefined region and used to
# raise "Decomposition undefined when both m < 0 and nu < 0 ...".

test_that(".hzr_phase_second_derivatives does not error at the Case 3L boundary (m=0, nu<0)", {
  t <- c(0.2, 0.5, 1.0, 2.0, 4.0)
  d2 <- .hzr_phase_second_derivatives(t, t_half = 1.666871, nu = -0.7462,
                                       m = 0, type = "cdf")
  expect_true(all(vapply(d2, function(v) all(is.finite(v)), logical(1))))
})

test_that(".hzr_phase_second_derivatives does not error at the Case 2L boundary (nu=0, m<0)", {
  t <- c(0.2, 0.5, 1.0, 2.0, 4.0)
  d2 <- .hzr_phase_second_derivatives(t, t_half = 1.0, nu = 0, m = -0.5,
                                       type = "cdf")
  expect_true(all(vapply(d2, function(v) all(is.finite(v)), logical(1))))
})

test_that(".hzr_phase_second_derivatives boundary values are continuous with the interior", {
  # m = 1e-3 is far enough from 0 (relative to h ~ 1.2e-4) that both m+h and
  # m-h stay positive, so this interior point still uses the ordinary
  # central-difference formula. Phi/phi are smooth in m here (m=0 is a
  # formula/case-dispatch artifact, not a true singularity), so the boundary
  # (forward-diff) value at m=0 should be close to the interior value.
  t <- c(0.5, 1.0, 2.0)
  d2_boundary <- .hzr_phase_second_derivatives(
    t, t_half = 1.666871, nu = -0.7462, m = 0, type = "cdf")
  d2_interior <- .hzr_phase_second_derivatives(
    t, t_half = 1.666871, nu = -0.7462, m = 1e-3, type = "cdf")
  expect_equal(d2_boundary$d2Phi_dm2, d2_interior$d2Phi_dm2, tolerance = 0.05)
  expect_equal(d2_boundary$d2Phi_dnu_dm, d2_interior$d2Phi_dnu_dm,
               tolerance = 0.05)
  expect_equal(d2_boundary$d2Phi_dlog_thalf_dm,
               d2_interior$d2Phi_dlog_thalf_dm, tolerance = 0.05)
})

test_that(".hzr_hessian_multiphase succeeds at the Case 3L fixed-shape boundary", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  phases <- list(
    early    = hzr_phase("cdf", t_half = 1.666871, nu = -0.7462, m = 0,
                         fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  theta <- c(early.log_mu = -1, early.log_t_half = log(1.666871),
             early.nu = -0.7462, early.m = 0, constant.log_mu = -1)
  cov_counts <- c(early = 0L, constant = 0L)
  H <- .hzr_hessian_multiphase(
    theta, time = avc$int_dead, status = avc$dead,
    phases = phases, covariate_counts = cov_counts,
    x_list = list(early = NULL, constant = NULL)
  )
  expect_true(is.matrix(H))
  expect_true(all(is.finite(H)))
})

test_that(".hzr_g3_phase_second_derivatives diagonals match numDeriv", {
  skip_if_not_installed("numDeriv")
  t <- c(0.5, 1.0, 3.0, 8.0)
  tau <- 2.0
  gamma <- 1.5
  alpha <- 0.8
  eta <- 0.6

  d2 <- .hzr_g3_phase_second_derivatives(t, tau = tau, gamma = gamma,
                                          alpha = alpha, eta = eta)
  par0 <- c(log(tau), gamma, alpha, eta)
  for (i in seq_along(t)) {
    g3_fn_i <- function(par) {
      hzr_decompos_g3(t[i], tau = exp(par[1]), gamma = par[2],
                       alpha = par[3], eta = par[4])$G3
    }
    H_nd <- numDeriv::hessian(g3_fn_i, par0)
    expect_equal(d2$d2Phi_dlog_tau2[i], H_nd[1, 1], tolerance = 1e-2,
                 label = paste0("d2G3/d(log_tau)^2 at t=", t[i]))
    expect_equal(d2$d2Phi_dgamma2[i],   H_nd[2, 2], tolerance = 1e-2,
                 label = paste0("d2G3/dgamma^2 at t=", t[i]))
    expect_equal(d2$d2Phi_dalpha2[i],   H_nd[3, 3], tolerance = 1e-2,
                 label = paste0("d2G3/dalpha^2 at t=", t[i]))
    expect_equal(d2$d2Phi_deta2[i],     H_nd[4, 4], tolerance = 1e-2,
                 label = paste0("d2G3/deta^2 at t=", t[i]))
  }
})

# ---------------------------------------------------------------------------
# Task 2: .hzr_hessian_multiphase() cross-check tests
# ---------------------------------------------------------------------------

test_that(".hzr_hessian_multiphase matches numDeriv (2-phase no covariates)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data    = avc,
    dist    = "multiphase",
    phases  = phases,
    fit     = TRUE,
    control = list(n_starts = 1, conserve = FALSE)
  )
  theta      <- fit$fit$par
  cov_counts <- fit$fit$covariate_counts
  x_list     <- fit$fit$x_list

  H_an <- .hzr_hessian_multiphase(
    theta, time = avc$int_dead, status = avc$dead,
    phases = phases, covariate_counts = cov_counts, x_list = x_list
  )
  expect_true(is.matrix(H_an))

  obj <- function(par) {
    -.hzr_logl_multiphase(par, avc$int_dead, avc$dead,
                           phases = phases,
                           covariate_counts = cov_counts,
                           x_list = x_list)
  }
  H_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-3,
               label = "2-phase no-cov NLL Hessian vs numDeriv")
})

test_that(".hzr_hessian_multiphase matches numDeriv (2-phase with covariates)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  phases <- list(
    early    = hzr_phase("cdf",
                          formula = ~ age + com_iv,
                          t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ early(age + com_iv),
    data    = avc,
    dist    = "multiphase",
    phases  = phases,
    fit     = TRUE,
    control = list(n_starts = 1, conserve = FALSE)
  )
  theta      <- fit$fit$par
  cov_counts <- fit$fit$covariate_counts
  x_list     <- fit$fit$x_list

  H_an <- .hzr_hessian_multiphase(
    theta, time = avc$int_dead, status = avc$dead,
    phases = phases, covariate_counts = cov_counts, x_list = x_list
  )
  obj <- function(par) {
    -.hzr_logl_multiphase(par, avc$int_dead, avc$dead,
                           phases = phases,
                           covariate_counts = cov_counts,
                           x_list = x_list)
  }
  H_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-3,
               label = "2-phase + covariates Hessian vs numDeriv")
})

test_that(".hzr_hessian_multiphase returns NULL for interval-censored rows", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  status_mixed <- avc$dead
  status_mixed[1:3] <- 2L

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  theta <- c(early.log_mu = 0, early.log_t_half = log(0.15), early.nu = 1.4,
             early.m = 1.0, constant.log_mu = 0)
  cov_counts <- c(early = 0L, constant = 0L)
  expect_null(.hzr_hessian_multiphase(
    theta, time = avc$int_dead, status = status_mixed,
    phases = phases, covariate_counts = cov_counts,
    x_list = list(early = NULL, constant = NULL)
  ))
})

test_that(".hzr_hessian_multiphase returns NULL for left-censored rows", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  status_left <- avc$dead
  status_left[1:3] <- -1L

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  theta <- c(early.log_mu = 0, early.log_t_half = log(0.15), early.nu = 1.4,
             early.m = 1.0, constant.log_mu = 0)
  cov_counts <- c(early = 0L, constant = 0L)
  expect_null(.hzr_hessian_multiphase(
    theta, time = avc$int_dead, status = status_left,
    phases = phases, covariate_counts = cov_counts,
    x_list = list(early = NULL, constant = NULL)
  ))
})

# ---------------------------------------------------------------------------
# Task 3: End-to-end wiring tests
# ---------------------------------------------------------------------------

test_that("multiphase fit vcov uses analytic Hessian (matches numDeriv to 1e-3)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = phases, fit = TRUE,
    control = list(n_starts = 1, conserve = FALSE)
  )
  theta      <- fit$fit$par
  cov_counts <- fit$fit$covariate_counts
  x_list     <- fit$fit$x_list
  free_idx   <- which(!fit$fit$fixed_mask)

  # numDeriv reference on the FREE-parameter objective (matches optimizer)
  obj_free <- function(par_free) {
    th <- theta
    th[free_idx] <- par_free
    -.hzr_logl_multiphase(th, avc$int_dead, avc$dead,
                           phases = phases,
                           covariate_counts = cov_counts,
                           x_list = x_list)
  }
  vcov_nd <- solve(numDeriv::hessian(obj_free, theta[free_idx]))

  # fit$fit$vcov rows/cols for free params only
  vcov_fit_free <- fit$fit$vcov[free_idx, free_idx]
  expect_equal(unname(vcov_fit_free), unname(vcov_nd),
               tolerance = 1e-3,
               label = "analytic vcov (free params) matches numDeriv reference")
})

test_that("CoE fit vcov (full-info) uses analytic Hessian", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  fit_coe <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = phases, fit = TRUE,
    control = list(n_starts = 1, conserve = TRUE)
  )
  # The conserved log_mu must now carry a finite SE (CoE full-info path)
  se <- sqrt(diag(fit_coe$fit$vcov))
  expect_true(all(is.finite(se[!is.na(se)])))
  expect_true(is.finite(fit_coe$fit$rcond))
})

# ---------------------------------------------------------------------------
# Task 4: Scale-invariance and 13-parameter anchor
# ---------------------------------------------------------------------------

test_that("multiphase SEs are invariant to covariate rescaling (log-mu z-stat)", {
  skip_on_cran()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  phases_base <- list(
    early    = hzr_phase("cdf", formula = ~ age,
                          t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  avc2 <- avc
  avc2$age <- avc$age / 10  # rescaled covariate

  f1 <- hazard(survival::Surv(int_dead, dead) ~ early(age),
               data = avc,  dist = "multiphase", phases = phases_base, fit = TRUE,
               control = list(n_starts = 1, conserve = FALSE))
  f2 <- hazard(survival::Surv(int_dead, dead) ~ early(age),
               data = avc2, dist = "multiphase", phases = phases_base, fit = TRUE,
               control = list(n_starts = 1, conserve = FALSE))

  # The beta z-statistic (estimate / SE) must be invariant under x -> x/10.
  se1 <- sqrt(diag(f1$fit$vcov))
  se2 <- sqrt(diag(f2$fit$vcov))
  beta_name <- grep("age", names(f1$fit$par), value = TRUE)[1]
  z1 <- f1$fit$par[beta_name] / se1[beta_name]
  # f2 beta for age is 10x larger (rescaled covariate), SE is 10x larger too
  beta_name2 <- grep("age", names(f2$fit$par), value = TRUE)[1]
  z2 <- f2$fit$par[beta_name2] / se2[beta_name2]
  expect_equal(unname(z1), unname(z2), tolerance = 1e-2,
               label = "beta z-stat invariant under covariate rescaling")
})

test_that("13-parameter multiphase anchor: stable SEs and rcond (supersedes placeholder)", {
  skip_on_cran()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  # 6 early-phase covariates + constant intercept = 10 free params
  # (log_mu_early + log_t_half + nu + 6 betas + log_mu_constant; m is fixed)
  # NOTE: AVC column names: age, com_iv, mal, opmos (NOT op_mos), op_age, status
  phases_hd <- list(
    early    = hzr_phase("cdf",
                          formula = ~ age + com_iv + mal + opmos + op_age + status,
                          t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  fit <- hazard(
    survival::Surv(int_dead, dead) ~
      early(age + com_iv + mal + opmos + op_age + status),
    data = avc, dist = "multiphase", phases = phases_hd, fit = TRUE,
    control = list(n_starts = 3, maxit = 800, conserve = TRUE)
  )

  expect_true(fit$fit$converged)
  free <- which(is.finite(diag(fit$fit$vcov)))
  expect_gte(length(free), 10L)
  se <- sqrt(diag(fit$fit$vcov)[free])
  expect_true(all(is.finite(se)))
  expect_true(all(se > 0))
  expect_true(is.finite(fit$fit$rcond))
  expect_true(isTRUE(fit$fit$pd),
              label = "10-free-param fit is positive-definite at optimum")
})
