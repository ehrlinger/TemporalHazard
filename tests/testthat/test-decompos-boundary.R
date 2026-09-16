# tests/testthat/test-decompos-boundary.R
#
# Phase 7d: boundary and limiting-case behavior of hzr_decompos().
#
# hzr_decompos() dispatches on the signs of (m, nu) into six branches plus two
# explicit error guards.  The branches include three "limiting" forms (m -> 0
# for nu > 0 and nu < 0, and nu -> 0 for m < 0) that must agree with the
# general branches as the parameter approaches the boundary.  These tests:
#   (1) confirm the previously-unhandled (nu == 0, m >= 0) combination now fails
#       loud instead of raising the cryptic "object 'G' not found";
#   (2) lock in continuity of the verified limiting branches;
#   (3) check internal consistency (g == dG/dt) and CDF sanity (monotone, [0,1]);
#   (4) record the known Case 3 (m > 0, nu < 0) <-> Case 3L discontinuity as a
#       skipped test pending validation against the C HAZARD G1 reference.

# ---------------------------------------------------------------------------
# (1) Unhandled-combination guard
# ---------------------------------------------------------------------------

test_that("nu == 0 with m >= 0 fails loud (not 'object G not found')", {
  # Previously these fell through every dispatch branch, leaving G unassigned
  # and raising a cryptic internal error.
  expect_error(
    hzr_decompos(c(1, 2), t_half = 3, nu = 0, m = 1),
    "Decomposition undefined for nu = 0 with m = 1"
  )
  expect_error(
    hzr_decompos(c(1, 2), t_half = 3, nu = 0, m = 0),
    "Decomposition undefined for nu = 0 with m = 0"
  )
})

test_that("m < 0 with nu < 0 still errors (existing guard)", {
  expect_error(
    hzr_decompos(c(1, 2), t_half = 3, nu = -1, m = -1),
    "Decomposition undefined when both m < 0 and nu < 0"
  )
})

# ---------------------------------------------------------------------------
# (2) Continuity of the verified limiting branches
# ---------------------------------------------------------------------------

test_that("Case 1 (m>0) limits to Case 1L (m=0) as m -> 0, nu > 0", {
  t <- c(0.3, 1, 2.5, 7)
  g_general <- hzr_decompos(t, t_half = 3, nu = 2, m = 1e-7)$G
  g_limit   <- hzr_decompos(t, t_half = 3, nu = 2, m = 0)$G
  expect_equal(g_general, g_limit, tolerance = 1e-5)
})

test_that("Case 2 (m<0) limits to Case 1L (m=0) as m -> 0, nu > 0", {
  t <- c(0.3, 1, 2.5, 7)
  g_general <- hzr_decompos(t, t_half = 3, nu = 2, m = -1e-7)$G
  g_limit   <- hzr_decompos(t, t_half = 3, nu = 2, m = 0)$G
  expect_equal(g_general, g_limit, tolerance = 1e-5)
})

test_that("Case 2 (nu>0) limits to Case 2L (nu=0) as nu -> 0, m < 0", {
  t <- c(0.3, 1, 2.5, 7)
  g_general <- hzr_decompos(t, t_half = 3, nu = 1e-7, m = -1)$G
  g_limit   <- hzr_decompos(t, t_half = 3, nu = 0,    m = -1)$G
  expect_equal(g_general, g_limit, tolerance = 1e-5)
})

# ---------------------------------------------------------------------------
# (3) Internal consistency (g == dG/dt) and CDF sanity
# ---------------------------------------------------------------------------

test_that("g equals the numerical derivative of G across the used branches", {
  check_deriv <- function(t_half, nu, m) {
    t0 <- 2.0
    h  <- 1e-6
    d  <- hzr_decompos(t0, t_half = t_half, nu = nu, m = m)
    dp <- hzr_decompos(t0 + h, t_half = t_half, nu = nu, m = m)
    dm <- hzr_decompos(t0 - h, t_half = t_half, nu = nu, m = m)
    num_g <- (dp$G - dm$G) / (2 * h)
    expect_equal(d$g, num_g, tolerance = 1e-4)
  }
  check_deriv(3,  2,  1)   # Case 1
  check_deriv(3,  2,  0)   # Case 1L
  check_deriv(3,  2, -1)   # Case 2
  check_deriv(3,  0, -1)   # Case 2L
})

test_that("G is a valid CDF (monotone, within [0, 1]) for the early-phase cases", {
  t <- seq(0.01, 30, length.out = 200)
  for (pars in list(c(nu = 2, m = 1), c(nu = 2, m = 0), c(nu = 2, m = -1))) {
    G <- hzr_decompos(t, t_half = 3, nu = pars[["nu"]], m = pars[["m"]])$G
    expect_true(all(G >= 0 & G <= 1),
                label = paste("G in [0,1] for nu=", pars[["nu"]], "m=", pars[["m"]]))
    expect_true(all(diff(G) >= -1e-9),
                label = paste("G monotone for nu=", pars[["nu"]], "m=", pars[["m"]]))
  }
})

# ---------------------------------------------------------------------------
# (4) Extreme t_half stability
# ---------------------------------------------------------------------------

test_that("extreme t_half values produce finite, valid G (nu > 0)", {
  t <- c(0.5, 1, 2, 5, 20)
  for (th in c(1e-6, 1e-3, 1e3, 1e6)) {
    G <- hzr_decompos(t, t_half = th, nu = 2, m = 1)$G
    expect_true(all(is.finite(G)), label = paste("finite G at t_half =", th))
    expect_true(all(G >= 0 & G <= 1), label = paste("G in [0,1] at t_half =", th))
    expect_true(all(diff(G) >= -1e-9), label = paste("G monotone at t_half =", th))
  }
})

# ---------------------------------------------------------------------------
# (5) Case 3 (m>0, nu<0) continuity and C parity
# ---------------------------------------------------------------------------
# Case 3's rho previously used a bare (2^m - 1)^nu (no /m divisor), leaving a
# spurious factor of m on the bt^(-1/nu) term. That diverged from the C G1
# evaluator (g1flag = 5) by up to ~0.2 and broke continuity with the m -> 0
# limit (Case 3L). The /m form matches C exactly and restores continuity.

test_that("Case 3 (m>0, nu<0) limits to Case 3L (m=0) as m -> 0", {
  t <- c(0.3, 1, 2.5, 7)
  g_general <- hzr_decompos(t, t_half = 3, nu = -2, m = 1e-7)$G
  g_limit   <- hzr_decompos(t, t_half = 3, nu = -2, m = 0)$G
  expect_equal(g_general, g_limit, tolerance = 1e-5)
})

test_that("Case 3 (m>0, nu<0) matches the C HAZARD g1flag=5 G1 evaluator", {
  # Direct port of hzd_set_rho() case 5 + hzd_ln_G1_and_SG1() case 5:
  #   rho_C = t_half * (2^m - 1)^nu;  G1 = 1 - ((t/rho_C)^(-1/nu) + 1)^(-1/m)
  c_g1flag5 <- function(time, t_half, nu, m) {
    rho_c <- t_half * (2^m - 1)^nu
    bt    <- time / rho_c
    1 - (bt^(-1 / nu) + 1)^(-1 / m)
  }
  t <- c(0.3, 1, 2.5, 7, 15)
  for (pars in list(c(nu = -2, m = 0.5), c(nu = -1, m = 2), c(nu = -3, m = 1))) {
    nu <- pars[["nu"]]
    m  <- pars[["m"]]
    expect_equal(
      hzr_decompos(t, t_half = 3, nu = nu, m = m)$G,
      c_g1flag5(t, t_half = 3, nu = nu, m = m),
      tolerance = 1e-10,
      label = paste("Case 3 vs C g1flag=5 for nu =", nu, "m =", m)
    )
  }
})

# G3 below the underflow of (t/tau)^gamma -----------------------------------
# With gamma large, (t/tau)^gamma underflows for t well below tau. The
# log-scale form used to clamp that to double.xmin, which froze G3 and put
# g3 wrong by more than 100 on the log scale, with no warning. Below the
# underflow the exact log forms reduce to linear functions of ln(t/tau), so
# compare against those.
test_that("G3 and g3 stay accurate where (t/tau)^gamma underflows", {
  tau <- exp(2.6238)
  gamma <- 220.08
  eta <- 0.0090875
  t <- c(0.1, 0.3, 0.5, 0.55)
  x <- gamma * log(t / tau)
  expect_true(all(x < -700))

  # alpha > 0: ln G3 = eta * (x - ln alpha)
  alpha <- 1.2174
  d <- hzr_decompos_g3(t, tau, gamma, alpha, eta)
  ln_scale <- log(gamma) + log(eta) - log(tau) - log(alpha)
  expect_equal(log(d$G3), eta * (x - log(alpha)), tolerance = 1e-12)
  expect_equal(
    log(d$g3),
    ln_scale + (eta - 1) * (x - log(alpha)) + (gamma - 1) * log(t / tau),
    tolerance = 1e-12
  )
  # G3 must still move with t below the underflow
  expect_true(all(diff(d$G3) > 0))

  # t / tau underflowing to exactly 0 (t = 0 is clamped to double.xmin).
  # With gamma = 2, alpha = 1, eta = 0.5, G3 = t / tau exactly, so g3 is
  # 1 / tau at every t. With alpha = 0, eta = 1, g3 = 2 t / tau^2 while
  # (t / tau)^2 is tiny.
  t_tiny <- c(0, .Machine$double.xmin, 1e-300, 1)
  d_lin <- hzr_decompos_g3(t_tiny, tau = 1e20, gamma = 2, alpha = 1, eta = 0.5)
  # Compare ratios to 1: expect_equal() goes absolute when the expected
  # value is below the tolerance, and 0 would pass against 1e-20.
  expect_equal(d_lin$g3 * 1e20, rep(1, 4), tolerance = 1e-12)
  # (t / tau)^2 underflows at t = 1e-150 while g3 = 2e-190 does not
  t_a0 <- c(1e-150, 1)
  d_a0 <- hzr_decompos_g3(t_a0, tau = 1e20, gamma = 2, alpha = 0, eta = 1)
  expect_equal(d_a0$g3 / (2 * t_a0 / 1e40), c(1, 1), tolerance = 1e-12)

  # alpha large enough that the division underflows while (t/tau)^gamma
  # does not
  t_near <- 10 * exp(-1 / 220)
  d_big <- hzr_decompos_g3(t_near, 10, 220, 1e308, 0.5)
  expect_equal(log(d_big$G3),
               0.5 * (log(log1p(exp(-1))) - log(1e308)),
               tolerance = 1e-12)

  # alpha = 0: ln G3 = eta * x
  d0 <- hzr_decompos_g3(t, tau, gamma, 0, eta)
  expect_equal(log(d0$G3), eta * x, tolerance = 1e-12)
  expect_equal(
    log(d0$g3),
    log(gamma) + log(eta) - log(tau) + (eta - 1) * x + (gamma - 1) * log(t / tau),
    tolerance = 1e-12
  )
})
