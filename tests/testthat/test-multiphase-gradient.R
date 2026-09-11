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
# The m derivative near m = 0
# ============================================================================
#
# hzr_decompos() changes formula at m = 0, and Case 2 (m < 0) meets the
# m >= 0 family in a |m|^nu cusp, so a difference stencil that straddles 0
# mixes two branches. numerical_gradient() above straddles too -- its step is
# eps^(1/3) * max(|theta|, 1), about 6e-6 for any |m| < 1 -- so it cannot be
# the reference here. This one keeps every stencil point on the side of `x`
# (the step is at most 1% of |x|) and Richardson-extrapolates; at x = 0 it
# differences forward, the m >= 0 side.
onside_derivative <- function(f, x, base = NULL) {
  if (x == 0) {
    h <- 1e-4 * 2^-(0:5)
    d <- sapply(h, function(hh) (f(hh) - f(0)) / hh)
    for (k in 1:4) d <- (2^k * d[-1] - d[-length(d)]) / (2^k - 1)
    return(d[length(d)])
  }
  if (is.null(base)) base <- 1e-2 * abs(x)
  h <- base * 2^-(0:3)
  d <- sapply(h, function(hh) (f(x + hh) - f(x - hh)) / (2 * hh))
  for (k in 1:3) d <- (4^k * d[-1] - d[-length(d)]) / (4^k - 1)
  d[length(d)]
}

test_that(".hzr_phase_derivatives dPhi/dm does not straddle m = 0", {
  t_grid <- c(0.05, 0.3, 2, 10)
  t_half <- 0.28
  # nu = -1 is Case 3L at m = 0 and Case 3 above it; m < 0 is undefined there,
  # so the old stencil already fell back to a forward difference and these
  # rows guard the new one-sided stencil rather than the original straddle.
  # Only the cdf type: the hazard type's Phi = -log(1 - G) loses digits to
  # cancellation at large t when nu < 0, which is a property of Phi, not of
  # the stencil.
  grid <- rbind(
    expand.grid(nu = c(0.5, 0.905, 2), m = c(-5e-6, 0, 5e-6),
                type = c("cdf", "hazard"), stringsAsFactors = FALSE),
    expand.grid(nu = -1, m = c(0, 5e-6), type = "cdf",
                stringsAsFactors = FALSE)
  )
  # Per entry, so a wrong small entry cannot hide behind a large one; the
  # floor stops an entry near zero demanding more than the reference can
  # deliver. A straddling stencil misses by 10% to 100000% on this grid.
  within <- function(v, r) {
    all(abs(v - r) <= 1e-3 * pmax(abs(r), 1e-2 * max(abs(r))))
  }
  for (k in seq_len(nrow(grid))) {
    nu <- grid$nu[k]
    m <- grid$m[k]
    type <- grid$type[k]
    pd <- .hzr_phase_derivatives(t_grid, t_half = t_half, nu = nu, m = m,
                                 type = type)
    cumhaz_at <- function(tt, z) {
      hzr_phase_cumhaz(tt, t_half = t_half, nu = nu, m = z, type = type)
    }
    hazard_at <- function(tt, z) {
      hzr_phase_hazard(tt, t_half = t_half, nu = nu, m = z, type = type)
    }
    ref_Phi <- sapply(t_grid, function(tt) {
      onside_derivative(function(z) cumhaz_at(tt, z), m)
    })
    ref_phi <- sapply(t_grid, function(tt) {
      onside_derivative(function(z) hazard_at(tt, z), m)
    })
    lab <- sprintf("nu = %g, m = %g, %s", nu, m, type)
    expect_true(within(pd$dPhi_dm, ref_Phi), label = paste("dPhi/dm,", lab))
    expect_true(within(pd$dphi_dm, ref_phi), label = paste("dphi/dm,", lab))
  }
})

test_that(".hzr_phase_derivatives dPhi/dm stays bounded as m -> 0 from below", {
  # For m < 0 the step is 1% of |m|, floored at 1e-10. Without the floor,
  # rounding grew without bound: dPhi/dm was 36 at m = -1e-15 where it should
  # be 0.151, and NA at -1e-300. At nu = 2 the derivative is flat in m this
  # close to 0, so it must hold its value. At nu = 0.9 it grows like
  # |m|^(nu - 1) toward the cusp, so it is only required to stay finite.
  d_at <- function(nu, m) {
    .hzr_phase_derivatives(c(0.05, 2), t_half = 0.28, nu = nu, m = m,
                           type = "cdf")$dPhi_dm
  }
  ref <- d_at(2, -1e-8)
  for (m in c(-1e-10, -1e-13, -1e-15)) {
    expect_true(all(abs(d_at(2, m) - ref) <= 1e-2 * abs(ref)),
                label = sprintf("nu = 2, m = %g", m))
  }
  for (m in c(-1e-13, -1e-15, -1e-300)) {
    expect_true(all(is.finite(d_at(0.9, m))),
                label = sprintf("nu = 0.9, m = %g", m))
  }
})

test_that("multiphase gradient is right for m within 1e-5 of 0, both signs", {
  set.seed(42)
  n <- 200
  time   <- rexp(n, rate = 0.3) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.3, 0.7))
  x_mat  <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("age", "sex")))
  phases <- list(
    early = hzr_phase("cdf",    t_half = 0.3, nu = 0.9, m = 0),
    late  = hzr_phase("hazard", t_half = 5,   nu = 1,   m = 1)
  )
  covariate_counts <- c(early = 2L, late = 2L)
  x_list <- list(early = x_mat, late = x_mat)
  ll <- function(th) {
    .hzr_logl_multiphase(
      th, time, status,
      phases = phases, covariate_counts = covariate_counts, x_list = x_list
    )
  }

  for (m in c(-5e-6, 5e-6)) {
    # [log_mu, log_t_half, nu, m, beta1, beta2] per phase; early m is index 4.
    theta <- c(log(0.1), log(0.3), 0.9, m, 0.5, -0.3,
               log(0.01), log(5), 1, 1, 0.1, 0.2)
    grad <- .hzr_gradient_multiphase(
      theta, time, status,
      phases = phases, covariate_counts = covariate_counts, x_list = x_list
    )
    ref <- vapply(seq_along(theta), function(i) {
      f <- function(z) {
        th <- theta
        th[i] <- z
        ll(th)
      }
      onside_derivative(f, theta[i],
                        base = if (i == 4L) NULL else 1e-3 * max(abs(theta[i]), 1))
    }, numeric(1))
    # Per component, so one wrong entry cannot hide behind the others.
    expect_true(all(abs(grad - ref) <= 1e-4 * pmax(abs(ref), 1)),
                label = sprintf("gradient at m = %g; worst component %d, error %.3g",
                                m, which.max(abs(grad - ref) / pmax(abs(ref), 1)),
                                max(abs(grad - ref))))
  }
})

test_that("multiphase Hessian m row does not straddle m = 0", {
  set.seed(42)
  n <- 200
  time   <- rexp(n, rate = 0.3) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.3, 0.7))
  x_mat  <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(NULL, c("age", "sex")))
  phases <- list(
    early = hzr_phase("cdf",    t_half = 0.3, nu = 0.9, m = 0),
    late  = hzr_phase("hazard", t_half = 5,   nu = 1,   m = 1)
  )
  covariate_counts <- c(early = 2L, late = 2L)
  x_list <- list(early = x_mat, late = x_mat)

  for (m in c(0, 5e-6)) {
    theta <- c(log(0.1), log(0.3), 0.9, m, 0.5, -0.3,
               log(0.01), log(5), 1, 1, 0.1, 0.2)
    H <- .hzr_hessian_multiphase(
      theta, time = time, status = status, time_lower = NULL,
      time_upper = NULL, x = NULL, weights = NULL, phases = phases,
      covariate_counts = covariate_counts, x_list = x_list
    )
    # H is the Hessian of the negative log-likelihood, so its m row is minus
    # the derivative of the score in m. m >= 0 is the smooth side, so the
    # reference differences the (corrected) score forward from m.
    score_at <- function(z) {
      th <- theta
      th[4] <- m + z
      .hzr_gradient_multiphase(
        th, time, status,
        phases = phases, covariate_counts = covariate_counts, x_list = x_list
      )
    }
    ref <- -vapply(seq_along(theta), function(j) {
      onside_derivative(function(z) score_at(z)[j], 0)
    }, numeric(1))
    expect_true(all(abs(H[4, ] - ref) <= 1e-2 * pmax(abs(ref), 1)),
                label = sprintf("Hessian m row at m = %g; worst error %.3g",
                                m, max(abs(H[4, ] - ref))))
  }
})

test_that("Hessian second derivatives step backward, not across, below m = 0", {
  # At m = -5e-5 the Hessian's 1.2e-4 step would cross 0, so it steps
  # backward (m_dir = -1). At nu = 2 the m < 0 side is smooth enough for that
  # to be accurate, which makes it the place to pin the sign of every m term:
  # dropping or flipping m_dir errs by about 200% here. (At nu < 2 the cusp
  # varies on the scale of |m| and a 1.2e-4 step understates the curvature;
  # that limit is documented, not tested.) The references difference Phi with
  # an m step of 0.1 * |m|, so they stay below 0.
  t_grid <- c(0.3, 2)
  t_half <- 0.28
  nu <- 2
  m <- -5e-5
  P <- function(th, n, mm) hzr_decompos(t_grid, th, n, mm)$G
  hm <- 0.1 * abs(m)
  ht <- 1e-3
  hn <- 1e-3
  ref_mm <- (P(t_half, nu, m + hm) - 2 * P(t_half, nu, m) +
               P(t_half, nu, m - hm)) / hm^2
  ref_tm <- (P(t_half * exp(ht), nu, m + hm) - P(t_half * exp(ht), nu, m - hm) -
               P(t_half * exp(-ht), nu, m + hm) +
               P(t_half * exp(-ht), nu, m - hm)) / (4 * ht * hm)
  ref_nm <- (P(t_half, nu + hn, m + hm) - P(t_half, nu + hn, m - hm) -
               P(t_half, nu - hn, m + hm) + P(t_half, nu - hn, m - hm)) /
    (4 * hn * hm)
  sd2 <- .hzr_phase_second_derivatives(t_grid, t_half = t_half, nu = nu,
                                       m = m, type = "cdf")
  within <- function(v, r) {
    all(abs(v - r) <= 0.05 * pmax(abs(r), 1e-2 * max(abs(r))))
  }
  expect_true(within(sd2$d2Phi_dm2, ref_mm), label = "d2Phi/dm2")
  expect_true(within(sd2$d2Phi_dlog_thalf_dm, ref_tm),
              label = "d2Phi/dlog_thalf dm")
  expect_true(within(sd2$d2Phi_dnu_dm, ref_nm), label = "d2Phi/dnu dm")

  # The same three terms for phi, which for the cdf type is the density g.
  g <- function(th, n, mm) hzr_decompos(t_grid, th, n, mm)$g
  gref_mm <- (g(t_half, nu, m + hm) - 2 * g(t_half, nu, m) +
                g(t_half, nu, m - hm)) / hm^2
  gref_tm <- (g(t_half * exp(ht), nu, m + hm) - g(t_half * exp(ht), nu, m - hm) -
                g(t_half * exp(-ht), nu, m + hm) +
                g(t_half * exp(-ht), nu, m - hm)) / (4 * ht * hm)
  gref_nm <- (g(t_half, nu + hn, m + hm) - g(t_half, nu + hn, m - hm) -
                g(t_half, nu - hn, m + hm) + g(t_half, nu - hn, m - hm)) /
    (4 * hn * hm)
  expect_true(within(sd2$d2phi_dm2, gref_mm), label = "d2phi/dm2")
  expect_true(within(sd2$d2phi_dlog_thalf_dm, gref_tm),
              label = "d2phi/dlog_thalf dm")
  expect_true(within(sd2$d2phi_dnu_dm, gref_nm), label = "d2phi/dnu dm")
})

test_that("the Hessian's m stencil never reaches 0 from below, m = -h included", {
  # At m = -h exactly a central stencil puts m + h on 0 itself, which belongs
  # to the m >= 0 family. Record every m the Hessian evaluates Phi at and
  # require all of them below 0. (A value comparison cannot see this: at
  # nu = 2 both stencils are accurate there.)
  h <- .hzr_h2
  real <- .hzr_phase_derivatives
  seen <- numeric(0)
  testthat::local_mocked_bindings(
    .hzr_phase_derivatives = function(time, t_half, nu, m, type) {
      seen <<- c(seen, m)
      real(time, t_half = t_half, nu = nu, m = m, type = type)
    }
  )
  for (m in c(-h, -h / 2, -2 * h)) {
    seen <- numeric(0)
    .hzr_phase_second_derivatives(c(0.3, 2), t_half = 0.28, nu = 2, m = m,
                                  type = "cdf")
    expect_true(length(seen) > 0 && all(seen < 0),
                label = sprintf("stencil points below 0 at m = %g; max %g",
                                m, max(seen)))
  }
})

test_that("Hessian second derivatives are finite with m and nu both at a boundary", {
  # m just below 0 with nu in [0, h): the nu / m cross term steps one-sided
  # in both. Only finiteness is checked -- the cusp below 0 is not resolved
  # at this step size (see above) -- but before this branch existed the same
  # point read a nu - h corner that does not exist.
  for (nu in c(0, .hzr_h2 / 2)) {
    sd2 <- .hzr_phase_second_derivatives(c(0.3, 2), t_half = 0.28, nu = nu,
                                         m = -5e-6, type = "cdf")
    expect_true(all(vapply(sd2, function(v) all(is.finite(v)), logical(1))),
                label = sprintf("all second derivatives finite at nu = %g", nu))
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
