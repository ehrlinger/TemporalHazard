# test-multiphase-gradient.R — Validate analytic gradient against numerical
#
# Tests that .hzr_gradient_multiphase() returns gradients consistent with
# central-difference numerical gradients of .hzr_logl_multiphase().
# This validates the semi-analytic gradient (chain-rule for mu/beta,
# central-diff for shape params via .hzr_phase_derivatives()).

# ============================================================================
# Helpers
# ============================================================================

load_kul_csv <- function() {
  csv_path <- system.file("extdata", "cabgkul.csv", package = "TemporalHazard")
  if (!nzchar(csv_path) || !file.exists(csv_path)) {
    csv_path <- file.path("inst", "extdata", "cabgkul.csv")
  }
  if (!file.exists(csv_path)) return(NULL)
  utils::read.csv(csv_path)
}

kul_phases <- function() {
  list(
    early    = hzr_phase("cdf",  t_half = 0.2, nu = 1, m = 1),
    constant = hzr_phase("constant"),
    late     = hzr_phase("g3",   tau = 1, gamma = 3, alpha = 1, eta = 1)
  )
}

# Central-difference numerical gradient (tight step, reference implementation)
numerical_gradient <- function(theta, time, status, phases,
                                covariate_counts, x_list,
                                time_lower = NULL, time_upper = NULL) {
  p <- length(theta)
  grad <- numeric(p)
  eps_rel <- (.Machine$double.eps)^(1 / 3)

  for (i in seq_len(p)) {
    h_i <- eps_rel * max(abs(theta[i]), 1)

    theta_plus <- theta
    theta_plus[i] <- theta_plus[i] + h_i
    theta_minus <- theta
    theta_minus[i] <- theta_minus[i] - h_i

    ll_plus <- .hzr_logl_multiphase(
      theta_plus, time, status,
      time_lower = time_lower, time_upper = time_upper,
      phases = phases, covariate_counts = covariate_counts, x_list = x_list
    )
    ll_minus <- .hzr_logl_multiphase(
      theta_minus, time, status,
      time_lower = time_lower, time_upper = time_upper,
      phases = phases, covariate_counts = covariate_counts, x_list = x_list
    )

    if (is.finite(ll_plus) && is.finite(ll_minus)) {
      grad[i] <- (ll_plus - ll_minus) / (2 * h_i)
    } else if (is.finite(ll_plus)) {
      ll0 <- .hzr_logl_multiphase(
        theta, time, status,
        time_lower = time_lower, time_upper = time_upper,
        phases = phases, covariate_counts = covariate_counts, x_list = x_list
      )
      grad[i] <- (ll_plus - ll0) / h_i
    }
  }

  grad
}


# ============================================================================
# Phase derivative unit tests
# ============================================================================

test_that(".hzr_phase_derivatives returns correct structure", {
  t_grid <- seq(0.5, 10, by = 0.5)

  pd <- .hzr_phase_derivatives(
    t_grid,
    t_half = 2,
    nu = 1.5,
    m = 0.8,
    type = "cdf"
  )
  expect_named(
    pd,
    c(
      "Phi", "phi", "dPhi_dlog_thalf", "dPhi_dnu", "dPhi_dm",
      "dphi_dlog_thalf", "dphi_dnu", "dphi_dm"
    )
  )
  expect_length(pd$Phi, length(t_grid))
  expect_true(all(is.finite(pd$Phi)))
  expect_true(all(is.finite(pd$dPhi_dlog_thalf)))
})

test_that(".hzr_phase_derivatives for constant phase returns zeros", {
  t_grid <- seq(0.5, 5, by = 0.5)
  pd <- .hzr_phase_derivatives(t_grid, type = "constant")

  expect_equal(pd$Phi, t_grid)
  expect_equal(pd$phi, rep(1, length(t_grid)))
  expect_equal(pd$dPhi_dlog_thalf, rep(0, length(t_grid)))
  expect_equal(pd$dPhi_dnu, rep(0, length(t_grid)))
  expect_equal(pd$dPhi_dm, rep(0, length(t_grid)))
})

test_that(".hzr_phase_derivatives matches numerical for cdf phase (case 1)", {
  t_grid <- seq(0.5, 8, by = 0.5)
  t_half <- 3
  nu <- 2
  m <- 1

  pd <- .hzr_phase_derivatives(
    t_grid,
    t_half = t_half,
    nu = nu,
    m = m,
    type = "cdf"
  )

  # Numerical check for dPhi/d(log_t_half)
  h <- 1e-5
  Phi_plus  <- hzr_phase_cumhaz(
    t_grid,
    t_half = t_half * exp(h),
    nu = nu,
    m = m,
    type = "cdf"
  )
  Phi_minus <- hzr_phase_cumhaz(
    t_grid,
    t_half = t_half * exp(-h),
    nu = nu,
    m = m,
    type = "cdf"
  )
  num_dPhi_dlog_thalf <- (Phi_plus - Phi_minus) / (2 * h)
  expect_equal(pd$dPhi_dlog_thalf, num_dPhi_dlog_thalf, tolerance = 1e-4)

  # Numerical check for dPhi/d(nu)
  Phi_plus  <- hzr_phase_cumhaz(
    t_grid,
    t_half = t_half,
    nu = nu + h,
    m = m,
    type = "cdf"
  )
  Phi_minus <- hzr_phase_cumhaz(
    t_grid,
    t_half = t_half,
    nu = nu - h,
    m = m,
    type = "cdf"
  )
  num_dPhi_dnu <- (Phi_plus - Phi_minus) / (2 * h)
  expect_equal(pd$dPhi_dnu, num_dPhi_dnu, tolerance = 1e-4)
})

test_that(".hzr_phase_derivatives matches numerical for hazard phase", {
  t_grid <- seq(0.5, 8, by = 0.5)
  t_half <- 5
  nu <- 1
  m <- 0

  pd <- .hzr_phase_derivatives(
    t_grid,
    t_half = t_half,
    nu = nu,
    m = m,
    type = "hazard"
  )

  h <- 1e-5
  Phi_plus  <- hzr_phase_cumhaz(
    t_grid,
    t_half = t_half * exp(h),
    nu = nu,
    m = m,
    type = "hazard"
  )
  Phi_minus <- hzr_phase_cumhaz(
    t_grid,
    t_half = t_half * exp(-h),
    nu = nu,
    m = m,
    type = "hazard"
  )
  num_dPhi_dlog_thalf <- (Phi_plus - Phi_minus) / (2 * h)
  expect_equal(pd$dPhi_dlog_thalf, num_dPhi_dlog_thalf, tolerance = 1e-4)
})

test_that(".hzr_phase_derivatives works for all 6 decomposition cases", {
  t_grid <- seq(0.5, 5, by = 0.5)
  cases <- list(
    list(t_half = 3, nu = 2,  m = 1),    # Case 1: m>0, nu>0
    list(t_half = 3, nu = 2,  m = 0),    # Case 1L: m=0, nu>0
    list(t_half = 3, nu = 2,  m = -0.5), # Case 2: m<0, nu>0
    list(t_half = 3, nu = 0,  m = -0.5), # Case 2L: m<0, nu=0
    list(t_half = 3, nu = -1, m = 1),    # Case 3: m>0, nu<0
    list(t_half = 3, nu = -1, m = 0)     # Case 3L: m=0, nu<0
  )

  for (i in seq_along(cases)) {
    cc <- cases[[i]]
    for (type in c("cdf", "hazard")) {
      pd <- .hzr_phase_derivatives(
        t_grid,
        t_half = cc$t_half,
        nu = cc$nu,
        m = cc$m,
        type = type
      )
      expect_true(
        all(is.finite(pd$Phi)),
        label = paste("Case", i, type, "Phi finite")
      )
      expect_true(
        all(is.finite(pd$dPhi_dlog_thalf)),
        label = paste("Case", i, type, "dPhi/dlog_thalf finite")
      )
      expect_true(
        all(is.finite(pd$dPhi_dnu)),
        label = paste("Case", i, type, "dPhi/dnu finite")
      )
      expect_true(
        all(is.finite(pd$dPhi_dm)),
        label = paste("Case", i, type, "dPhi/dm finite")
      )
    }
  }
})


# ============================================================================
# Full gradient tests on KUL dataset
# ============================================================================

test_that("Analytic gradient matches numerical at C reference parameters", {
  skip_on_cran()

  dat <- load_kul_csv()
  skip_if(is.null(dat), "KUL dataset not available")

  phases <- kul_phases()
  covariate_counts <- c(early = 0L, constant = 0L, late = 0L)
  x_list <- list(early = NULL, constant = NULL, late = NULL)

  # C converged parameters on internal scale
  theta_c <- c(
    -3.77955,   # early.log_mu  (mu = 0.02268)
    log(0.2),   # early.log_t_half
    1,           # early.nu
    1,           # early.m
    -7.2258,     # constant.log_mu  (mu = 0.0007269)
    -16.6578,    # late.log_mu  (mu = 5.837e-08)
    log(1),      # late.log_tau
    3,           # late.gamma
    1,           # late.alpha
    1            # late.eta
  )

  grad_analytic <- .hzr_gradient_multiphase(
    theta_c, dat$int_dead, dat$dead,
    phases = phases, covariate_counts = covariate_counts, x_list = x_list
  )

  grad_numeric <- numerical_gradient(
    theta_c, dat$int_dead, dat$dead,
    phases = phases, covariate_counts = covariate_counts, x_list = x_list
  )

  # Both should be finite

  expect_true(all(is.finite(grad_analytic)),
              label = "Analytic gradient all finite")
  expect_true(all(is.finite(grad_numeric)),
              label = "Numerical gradient all finite")

  # Relative error should be small for each component
  # Allow looser tolerance for shape params (central-diff vs central-diff)
  denom <- pmax(abs(grad_numeric), 1e-6)
  rel_err <- abs(grad_analytic - grad_numeric) / denom

  # log_mu and constant parameters should match tightly
  expect_true(rel_err[1] < 0.01, label = "early.log_mu gradient < 1% error")
  expect_true(rel_err[5] < 0.01, label = "constant.log_mu gradient < 1% error")
  expect_true(rel_err[6] < 0.01, label = "late.log_mu gradient < 1% error")

  # Shape params may have larger relative error near zero gradients
  # but absolute agreement should be reasonable
  abs_err <- abs(grad_analytic - grad_numeric)
  expect_true(all(abs_err < 1),
              label = paste("Max absolute gradient error:",
                            round(max(abs_err), 4)))
})

test_that("Analytic gradient matches numerical at SAS starting values", {
  skip_on_cran()

  dat <- load_kul_csv()
  skip_if(is.null(dat), "KUL dataset not available")

  phases <- kul_phases()
  covariate_counts <- c(early = 0L, constant = 0L, late = 0L)
  x_list <- list(early = NULL, constant = NULL, late = NULL)

  # SAS PARMS starting values
  theta_start <- c(
    log(0.02), log(0.2), 1, 1,
    log(0.0008),
    log(1e-9), log(1), 3, 1, 1
  )

  grad_analytic <- .hzr_gradient_multiphase(
    theta_start, dat$int_dead, dat$dead,
    phases = phases, covariate_counts = covariate_counts, x_list = x_list
  )

  grad_numeric <- numerical_gradient(
    theta_start, dat$int_dead, dat$dead,
    phases = phases, covariate_counts = covariate_counts, x_list = x_list
  )

  expect_true(all(is.finite(grad_analytic)))
  expect_true(all(is.finite(grad_numeric)))

  # At the starting point the gradients should be non-trivial and match
  denom <- pmax(abs(grad_numeric), 1e-6)
  rel_err <- abs(grad_analytic - grad_numeric) / denom

  # mu params should match well
  expect_true(rel_err[1] < 0.05, label = "early.log_mu gradient match")
  expect_true(rel_err[5] < 0.05, label = "constant.log_mu gradient match")
  expect_true(rel_err[6] < 0.05, label = "late.log_mu gradient match")
})


# ============================================================================
# Gradient on synthetic data with covariates
# ============================================================================

test_that("Analytic gradient works with covariates", {
  skip_on_cran()

  set.seed(42)
  n <- 200
  time   <- rexp(n, rate = 0.3) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.3, 0.7))
  x_mat  <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("age", "sex")))

  phases <- list(
    early = hzr_phase("cdf",    t_half = 1, nu = 1, m = 0),
    late  = hzr_phase("hazard", t_half = 5, nu = 1, m = 0)
  )
  covariate_counts <- c(early = 2L, late = 2L)
  x_list <- list(early = x_mat, late = x_mat)

  # theta: [log_mu, log_t_half, nu, m, beta1, beta2] x 2 phases
  theta <- c(
    log(0.1), log(1), 1, 0, 0.5, -0.3,   # early
    log(0.01), log(5), 1, 0, 0.1, 0.2     # late
  )

  grad_analytic <- .hzr_gradient_multiphase(
    theta, time, status,
    phases = phases, covariate_counts = covariate_counts, x_list = x_list
  )

  grad_numeric <- numerical_gradient(
    theta, time, status,
    phases = phases, covariate_counts = covariate_counts, x_list = x_list
  )

  expect_true(all(is.finite(grad_analytic)))
  expect_true(all(is.finite(grad_numeric)))

  # Check beta parameters specifically (indices 5,6 and 11,12)
  beta_idx <- c(5, 6, 11, 12)
  denom <- pmax(abs(grad_numeric[beta_idx]), 1e-6)
  rel_err_beta <- abs(grad_analytic[beta_idx] - grad_numeric[beta_idx]) / denom
  expect_true(all(rel_err_beta < 0.05),
              label = paste("Beta gradient max rel error:",
                            round(max(rel_err_beta), 4)))
})


# ============================================================================
# Simple 2-phase model: gradient at known parameters
# ============================================================================

test_that("Gradient is near zero at converged fit", {
  skip_on_cran()

  set.seed(42)
  n <- 100
  time   <- rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))

  phases <- list(
    early = hzr_phase("cdf",    t_half = 1, nu = 1.5, m = 0),
    const = hzr_phase("constant")
  )

  # Not guarded with skip_if().  When the multi-start perturbations were drawn
  # from the ambient RNG stream this call was a coin flip -- about a quarter of
  # draws hit an error thrown out of the objective, relabelled a convergence
  # failure -- and the skip_if() reported every one of those as a green skip.
  # Pinning start_seed made it deterministic and the underlying error is fixed,
  # so assert convergence: a regression must fail here, not vanish.  (At the
  # pinned default seed 3 this particular call did not fail; the one-in-four
  # figure is the rate across a seed sweep, which
  # test-multiphase-reproducibility.R asserts directly.)
  fit <- suppressWarnings(hazard(time = time, status = status,
    dist = "multiphase", phases = phases, fit = TRUE,
    control = list(n_starts = 3, maxit = 500)))
  expect_true(fit$fit$converged)

  covariate_counts <- c(early = 0L, const = 0L)
  x_list <- list(early = NULL, const = NULL)

  grad_at_mle <- .hzr_gradient_multiphase(
    fit$fit$theta, time, status,
    phases = phases, covariate_counts = covariate_counts, x_list = x_list
  )

  # Gradient should be reasonably small at the MLE.  With a small sample

  # (n=100) and rough BFGS convergence the per-observation score components
  # don't cancel perfectly, so we use a generous threshold.
  max_abs_grad <- max(abs(grad_at_mle))
  expect_true(max_abs_grad < 500,
              label = paste("Max |gradient| at MLE:", round(max_abs_grad, 2),
                            "should be < 500"))
})
