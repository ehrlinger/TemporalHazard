# tests/testthat/test-decompos-large-m.R
#
# Large-|m| numerical stability of hzr_decompos().
#
# The Case 1 (m > 0, nu > 0) and Case 3 (m > 0, nu < 0) branches formed
# `2^m` and `bt^(-1/nu)` directly.  Both overflow to Inf well inside the
# range an optimizer can reach, and `Inf^(-1/m)` then collapses silently to
# 0 (Case 1) or 1 (Case 3) with g and h going NaN.  No warning was raised --
# G = 0 is a legitimate-looking probability, which is exactly the
# "output that looks like a result and is not" failure mode AGENTS.md calls
# this package's signature defect.
#
# Measured breakdown of the pre-fix code (t_half = 1, t = 0.5):
#   Case 1, nu = 3/m : G correct to m = 700, G = 0 from m = 800
#   Case 1, nu = 2   : G already wrong by m = 700 (rho overflows sooner)
#   Case 3, nu = -3/m: G correct to m = 800, G = 1 from m = 1024
#
# Case 2 (m < 0, nu > 0) fails by a different mechanism: `1 - 2^m` rounds to
# exactly 1 once 2^m < eps, so (1 - 2^m)^(-nu) - 1 is 0 and rho is Inf.  That
# collapse is at m = -53, but the accuracy decays long before it -- at m = -20
# the pre-fix code was already wrong in the 10th digit.  Note this is a much
# more reachable region than the m ~ +750 of Case 1.
#
# Expected values below were computed independently of the R implementation,
# at 100 decimal digits with Python's `decimal` module, from
#   G(t) = [1 + (t_half/t)^(1/nu) (2^m - 1)]^(-1/m)          (Case 1)
#   G(t) = 1 - [1 + (t_half/t)^(1/nu) (2^m - 1)]^(-1/m)      (Case 3)
# with g obtained as a central difference of G at step t*1e-30 -- so the
# density is checked against the derivative of the CDF, not against itself.

tol <- 1e-12

# ---------------------------------------------------------------------------
# (1) Case 1 -- the reported defect
# ---------------------------------------------------------------------------

test_that("Case 1 G is correct for large m (was silently 0)", {
  ref <- list(
    list(t = 0.5,  nu = 3 / 1000, m = 1000, G = 3.96850262992049896e-01),
    list(t = 0.5,  nu = 3 / 2000, m = 2000, G = 3.96850262992049896e-01),
    list(t = 2.0,  nu = 3 / 1500, m = 1500, G = 6.29960524947436595e-01),
    list(t = 0.25, nu = 2,        m = 900,  G = 4.99615066482928194e-01),
    list(t = 1.5,  nu = 2,        m = 900,  G = 5.00112641882985209e-01),
    list(t = 0.5,  nu = 0.5,      m = 5000, G = 4.99861389780232535e-01)
  )
  for (r in ref) {
    got <- hzr_decompos(r$t, t_half = 1, nu = r$nu, m = r$m)$G
    expect_equal(got, r$G, tolerance = tol,
                 label = sprintf("G(t=%g, nu=%g, m=%g)", r$t, r$nu, r$m))
  }
})

test_that("Case 1 density g is correct for large m (was NaN)", {
  ref <- list(
    list(t = 0.5,  nu = 3 / 1000, m = 1000, g = 2.64566841994699931e-01),
    list(t = 2.0,  nu = 3 / 1500, m = 1500, g = 1.04993420824572761e-01),
    list(t = 0.25, nu = 2,        m = 900,  g = 1.11025570329539598e-03),
    list(t = 0.5,  nu = 0.5,      m = 5000, g = 3.99889111824186039e-04)
  )
  for (r in ref) {
    d <- hzr_decompos(r$t, t_half = 1, nu = r$nu, m = r$m)
    expect_equal(d$g, r$g, tolerance = 1e-9,
                 label = sprintf("g(t=%g, nu=%g, m=%g)", r$t, r$nu, r$m))
    expect_true(is.finite(d$h),
                label = sprintf("h(t=%g, nu=%g, m=%g) finite", r$t, r$nu, r$m))
  }
})

test_that("Case 1 converges to the analytic m -> Inf limit along m*nu = k", {
  # As m -> Inf with k = m*nu held fixed,
  #   G(t) -> 0.5 * (t/t_half)^(1/k)   for t <= t_half * 2^k.
  k <- 3
  t <- c(0.2, 0.5, 1, 3)
  analytic <- 0.5 * (t / 1)^(1 / k)
  for (m in c(1e3, 1e4, 1e5, 1e6)) {
    got <- hzr_decompos(t, t_half = 1, nu = k / m, m = m)$G
    expect_equal(got, analytic, tolerance = 1e-6,
                 label = sprintf("analytic limit at m=%g", m))
  }
})

# ---------------------------------------------------------------------------
# (2) Case 3 -- same overflow, opposite collapse (G -> 1)
# ---------------------------------------------------------------------------

test_that("Case 3 G and g are correct for large m (was G = 1, g = NaN)", {
  ref <- list(
    list(t = 0.5, nu = -3 / 1024, m = 1024,
         G = 3.70039475052563405e-01, g = 4.19973683298291045e-01),
    list(t = 0.5, nu = -3 / 2000, m = 2000,
         G = 3.70039475052563405e-01, g = 4.19973683298291045e-01),
    list(t = 2.0, nu = -1,        m = 1500,
         G = 5.00230995684740276e-01, g = 1.66589668105086554e-04)
  )
  for (r in ref) {
    d <- hzr_decompos(r$t, t_half = 1, nu = r$nu, m = r$m)
    expect_equal(d$G, r$G, tolerance = tol,
                 label = sprintf("Case 3 G(t=%g, nu=%g, m=%g)", r$t, r$nu, r$m))
    expect_equal(d$g, r$g, tolerance = 1e-9,
                 label = sprintf("Case 3 g(t=%g, nu=%g, m=%g)", r$t, r$nu, r$m))
    expect_true(is.finite(d$h))
  }
})

# ---------------------------------------------------------------------------
# (2b) Case 2 -- cancellation in 1 - 2^m, not overflow
# ---------------------------------------------------------------------------

test_that("Case 2 G and g are correct for large negative m (was 0 / NaN)", {
  ref <- list(
    list(t = 0.5, t_half = 1, nu = 2,   m = -60,
         G = 4.94257010176448075e-01, g = 1.64752336725482694e-02),
    list(t = 0.5, t_half = 1, nu = 2,   m = -100,
         G = 4.96546247718517963e-01, g = 9.93092495437035948e-03),
    list(t = 0.5, t_half = 1, nu = 2,   m = -1000,
         G = 4.99653546495226253e-01, g = 9.99307092990452421e-04),
    list(t = 2.0, t_half = 1, nu = 1.5, m = -500,
         G = 5.00693627855667289e-01, g = 5.00693627855667276e-04)
  )
  for (r in ref) {
    d <- hzr_decompos(r$t, t_half = r$t_half, nu = r$nu, m = r$m)
    lab <- sprintf("Case 2 (t=%g, nu=%g, m=%g)", r$t, r$nu, r$m)
    expect_equal(d$G, r$G, tolerance = tol, label = paste("G", lab))
    expect_equal(d$g, r$g, tolerance = 1e-9, label = paste("g", lab))
    expect_true(is.finite(d$h), label = paste("h finite", lab))
  }
})

test_that("Case 2 is accurate at moderate m where the old form had decayed", {
  # m = -20 is far from the m = -53 collapse, but the direct form was already
  # returning 0.42180716876655583 against a true 0.42180716891682779.
  ref <- list(
    list(t = 0.1, nu = 0.5, m = -20,
         G = 4.21807168916827790e-01, g = 2.10903579430077093e-01),
    list(t = 10,  nu = 0.5, m = -20,
         G = 5.31023701156430827e-01, g = 2.65511217549355008e-03)
  )
  for (r in ref) {
    d <- hzr_decompos(r$t, t_half = 3, nu = r$nu, m = r$m)
    expect_equal(d$G, r$G, tolerance = 1e-14,
                 label = sprintf("Case 2 G(t=%g, m=%g)", r$t, r$m))
    expect_equal(d$g, r$g, tolerance = 1e-9,
                 label = sprintf("Case 2 g(t=%g, m=%g)", r$t, r$m))
  }
})

test_that("Case 2 follows the m -> -Inf asymptote G = (1/2)(t/t_half)^(1/|m|)", {
  # rho -> Inf and bt -> 1, so G -> 1/2 for every t, approached as
  # 0.5 * (t/t_half)^(1/|m|).  The fixed branch reproduces that asymptote to
  # machine precision; the pre-fix code returned 0 for all of these.
  t <- c(0.2, 1, 5)
  for (m in c(-100, -500, -1000)) {
    expect_equal(hzr_decompos(t, t_half = 1, nu = 2, m = m)$G,
                 0.5 * t^(1 / abs(m)), tolerance = 1e-12,
                 label = sprintf("Case 2 asymptote at m=%g", m))
  }
})

test_that("Case 2 returns NA, not a plausible number, past the 2^m underflow", {
  # Below m ~ -1074 the term 2^m underflows to 0 outright and no rearrangement
  # recovers it in double precision.  The branch yields NA there.  That is a
  # deliberate line: NA propagates visibly, whereas the pre-fix code returned
  # 0 -- a perfectly plausible probability.  Pinned here so any future change
  # to this boundary is a decision rather than an accident.
  expect_true(all(is.na(hzr_decompos(c(0.5, 2), t_half = 1, nu = 2,
                                     m = -5000)$G)))
  # ...and the last decade before the cliff is still exact.
  expect_equal(hzr_decompos(0.5, t_half = 1, nu = 2, m = -1000)$G,
               4.99653546495226253e-01, tolerance = 1e-14)
})

# ---------------------------------------------------------------------------
# (3) No regression at moderate m, where the old code was already correct
# ---------------------------------------------------------------------------

test_that("moderate m still matches the high-precision reference", {
  ref <- list(
    list(t = 0.5, nu = 2,   m = 1,  G = 2.89897948556635643e-01,
         g = 2.05857127979289872e-01),
    list(t = 3.0, nu = 2,   m = 1,  G = 5.00000000000000000e-01,
         g = 4.16666666666666644e-02),
    list(t = 7.0, nu = 0.5, m = 5,  G = 6.83698216917979606e-01,
         g = 3.32320213449679991e-02),
    list(t = 1.0, nu = 1.5, m = 20, G = 4.82021015421515597e-01,
         g = 1.60673598141664598e-02)
  )
  for (r in ref) {
    d <- hzr_decompos(r$t, t_half = 3, nu = r$nu, m = r$m)
    expect_equal(d$G, r$G, tolerance = tol,
                 label = sprintf("G(t=%g, nu=%g, m=%g)", r$t, r$nu, r$m))
    expect_equal(d$g, r$g, tolerance = 1e-9,
                 label = sprintf("g(t=%g, nu=%g, m=%g)", r$t, r$nu, r$m))
  }
})

# ---------------------------------------------------------------------------
# (4) Structural sanity across the whole large-m range
# ---------------------------------------------------------------------------

test_that("G stays a monotone CDF across large m", {
  # Note on the upper end of the grid: along the ridge m*nu = k the family
  # saturates at t = t_half * 2^k, and G(t) = 1 beyond it -- this is the
  # analytic m -> Inf limit, not an overflow.  At k = 3, t_half = 1 the
  # 100-digit reference gives G(20) = 1 - 9.5e-110, which is exactly 1 in
  # double precision.  So G is asserted to be in (0, 1] and non-decreasing
  # over the full grid, and strictly inside (0, 1) and strictly increasing
  # over the sub-grid that lies below saturation.
  t     <- c(0.1, 0.5, 1, 2, 5, 20)
  inner <- t <= 5

  for (m in c(100, 800, 1024, 5000)) {
    for (nu in c(3 / m, 0.5, 2, -3 / m)) {
      d   <- hzr_decompos(t, t_half = 1, nu = nu, m = m)
      lab <- sprintf("m=%g nu=%g", m, nu)
      expect_true(all(is.finite(d$G)), label = paste("finite G,", lab))
      expect_true(all(is.finite(d$g)), label = paste("finite g,", lab))
      expect_true(all(d$G > 0 & d$G <= 1), label = paste("G in (0,1],", lab))
      expect_true(all(diff(d$G) >= 0), label = paste("G non-decreasing,", lab))
      expect_true(all(d$G[inner] < 1), label = paste("G < 1 below sat,", lab))
      expect_true(all(diff(d$G[inner]) > 0),
                  label = paste("G increasing below sat,", lab))
    }
  }
})

test_that("Case 1 saturates to G = 1 above t_half * 2^(m*nu), not below", {
  # Guards the tail assertion above: the saturation must be located where the
  # analytic limit puts it.  The pre-fix code returned G = 0 everywhere for
  # large m, which would fail the "below saturation" half of this test.
  k <- 3
  m <- 2000
  expect_true(all(hzr_decompos(c(0.1, 1, 5, 7.9), t_half = 1,
                               nu = k / m, m = m)$G < 1))
  expect_equal(hzr_decompos(c(50, 500), t_half = 1, nu = k / m, m = m)$G,
               c(1, 1))
})

test_that("g equals dG/dt at large m (central difference)", {
  for (m in c(900, 1024, 2000)) {
    for (nu in c(3 / m, 2, -3 / m)) {
      t0 <- 0.7
      h  <- 1e-5
      d  <- hzr_decompos(t0, t_half = 1, nu = nu, m = m)
      fd <- (hzr_decompos(t0 + h, t_half = 1, nu = nu, m = m)$G -
             hzr_decompos(t0 - h, t_half = 1, nu = nu, m = m)$G) / (2 * h)
      expect_equal(d$g, fd, tolerance = 1e-5,
                   label = sprintf("g vs dG/dt at m=%g nu=%g", m, nu))
    }
  }
})

# ---------------------------------------------------------------------------
# (5) Downstream: the multiphase likelihood becomes evaluable again
# ---------------------------------------------------------------------------

test_that("multiphase log-likelihood is finite for large m", {
  skip_on_cran()

  set.seed(42)
  n      <- 100
  time   <- rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))

  phases <- list(
    early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
    const = hzr_phase("constant")
  )
  cc <- c(early = 0L, const = 0L)
  xl <- list(early = NULL, const = NULL)

  k  <- 3
  ll <- vapply(c(100, 300, 500, 1000, 2000), function(m) {
    .hzr_logl_multiphase(
      c(log(0.5), log(1), k / m, m, log(0.5)),
      time, status, phases = phases, covariate_counts = cc, x_list = xl
    )
  }, numeric(1))

  expect_true(all(is.finite(ll)))
  # Along the ridge m*nu = k the likelihood approaches a finite supremum,
  # so successive values must be close rather than diverging.
  expect_lt(max(abs(diff(ll))), 0.1)
})
