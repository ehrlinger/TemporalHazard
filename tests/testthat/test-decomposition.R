# Too slow for the CRAN check time budget; runs on CI and locally.
# See "Where the check's time actually goes" in AGENTS.md.
skip_on_cran()

# test-decomposition.R — Tests for hzr_decompos() and phase helpers
#
# Validates all 6 parameter cases, mathematical identities, and golden
# reference values derived from the closed-form decomposition equations.

t_grid <- seq(0.1, 10, by = 0.1)

# ============================================================================
# Input validation
# ============================================================================

test_that("hzr_decompos rejects invalid inputs", {
  expect_error(
    hzr_decompos("a", t_half = 1, nu = 1, m = 0),
    "time must be numeric"
  )
  expect_error(
    hzr_decompos(1:5, t_half = -1, nu = 1, m = 0),
    "t_half must be a positive scalar"
  )
  expect_error(
    hzr_decompos(1:5, t_half = 0, nu = 1, m = 0),
    "t_half must be a positive scalar"
  )
  expect_error(
    hzr_decompos(1:5, t_half = 1, nu = NA, m = 0),
    "nu must be a numeric scalar"
  )
  expect_error(
    hzr_decompos(1:5, t_half = 1, nu = 1, m = Inf),
    "m must be a numeric scalar"
  )
})

test_that("hzr_decompos rejects m < 0 AND nu < 0", {
  expect_error(hzr_decompos(1:5, t_half = 1, nu = -1, m = -1),
               "both m < 0 and nu < 0")
})


# ============================================================================
# Case 1: m > 0, nu > 0 (standard sigmoidal)
# ============================================================================

test_that("Case 1 (m>0, nu>0): G is monotone increasing in [0,1]", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 1)
  expect_true(all(d$G >= 0 & d$G <= 1))
  expect_true(all(diff(d$G) >= -1e-12))  # monotone non-decreasing
})

test_that("Case 1: G(t_half) = 0.5", {
  d <- hzr_decompos(3, t_half = 3, nu = 2, m = 1)
  expect_equal(d$G, 0.5, tolerance = 1e-10)
})

test_that("Case 1: g >= 0 and h >= 0", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 1)
  expect_true(all(d$g >= -1e-15))
  expect_true(all(d$h >= -1e-15))
})


# ============================================================================
# Case 1L: m = 0, nu > 0 (Weibull-like)
# ============================================================================

test_that("Case 1L (m=0, nu>0): G(t_half) = 0.5", {
  d <- hzr_decompos(3, t_half = 3, nu = 2, m = 0)
  expect_equal(d$G, 0.5, tolerance = 1e-10)
})

test_that("Case 1L: G is monotone increasing in [0,1]", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
  expect_true(all(d$G >= 0 & d$G <= 1))
  expect_true(all(diff(d$G) >= -1e-12))
})


# ============================================================================
# Case 2: m < 0, nu > 0 (heavy-tailed)
# ============================================================================

test_that("Case 2 (m<0, nu>0): G(t_half) = 0.5", {
  d <- hzr_decompos(3, t_half = 3, nu = 2, m = -1)
  expect_equal(d$G, 0.5, tolerance = 1e-10)
})

test_that("Case 2: G is monotone increasing", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = -1)
  expect_true(all(diff(d$G) >= -1e-12))
  expect_true(all(d$g >= -1e-15))
})


# ============================================================================
# Case 2L: m < 0, nu = 0 (exponential decay)
# ============================================================================

test_that("Case 2L (m<0, nu=0): G(t_half) = 0.5", {
  d <- hzr_decompos(3, t_half = 3, nu = 0, m = -1)
  expect_equal(d$G, 0.5, tolerance = 1e-10)
})

test_that("Case 2L: G is monotone increasing", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 0, m = -1)
  expect_true(all(diff(d$G) >= -1e-12))
})


# ============================================================================
# Case 3: m > 0, nu < 0 (bounded cumulative)
# ============================================================================

test_that("Case 3 (m>0, nu<0): G(t_half) = 0.5", {
  d <- hzr_decompos(3, t_half = 3, nu = -2, m = 1)
  expect_equal(d$G, 0.5, tolerance = 1e-10)
})

test_that("Case 3: G is monotone increasing in [0,1]", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = -2, m = 1)
  expect_true(all(d$G >= 0 & d$G <= 1))
  expect_true(all(diff(d$G) >= -1e-12))
})


# ============================================================================
# Case 3L: m = 0, nu < 0 (bounded exponential)
# ============================================================================

test_that("Case 3L (m=0, nu<0): G(t_half) = 0.5", {
  d <- hzr_decompos(3, t_half = 3, nu = -2, m = 0)
  expect_equal(d$G, 0.5, tolerance = 1e-10)
})

test_that("Case 3L: G is monotone increasing in [0,1]", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = -2, m = 0)
  expect_true(all(d$G >= 0 & d$G <= 1))
  expect_true(all(diff(d$G) >= -1e-12))
})


# ============================================================================
# Mathematical identities
# ============================================================================

test_that("h(t) = g(t) / (1 - G(t)) identity holds across all cases", {
  cases <- list(
    list(t_half = 3, nu =  2, m =  1),   # Case 1
    list(t_half = 3, nu =  2, m =  0),   # Case 1L
    list(t_half = 3, nu =  2, m = -1),   # Case 2
    list(t_half = 3, nu =  0, m = -1),   # Case 2L
    list(t_half = 3, nu = -2, m =  1),   # Case 3
    list(t_half = 3, nu = -2, m =  0)    # Case 3L
  )

  for (pars in cases) {
    d <- hzr_decompos(
      t_grid,
      t_half = pars$t_half,
      nu = pars$nu,
      m = pars$m
    )
    # Only check where G(t) < 1 - eps (hazard blows up as G -> 1)
    valid <- d$G < 0.999
    if (sum(valid) > 0) {
      h_manual <- d$g[valid] / (1 - d$G[valid])
      expect_equal(
        d$h[valid],
        h_manual,
        tolerance = 1e-10,
        label = paste0("h = g/(1-G) for m=", pars$m,
                       ", nu=", pars$nu)
      )
    }
  }
})

test_that("G(t) integrates consistently with g(t) via numerical check", {
  # For a smooth case, check that diff(G) ~ g * dt
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
  dt <- diff(t_grid)
  dG <- diff(d$G)
  g_mid <- (d$g[-1] + d$g[-length(d$g)]) / 2  # midpoint rule
  ratio <- dG / (g_mid * dt)
  # Should be approximately 1 everywhere

  expect_true(
    all(abs(ratio - 1) < 0.05),
    label = "dG/dt approximately equals g(t)"
  )
})


# ============================================================================
# hzr_phase_cumhaz
# ============================================================================

test_that("hzr_phase_cumhaz('cdf') returns G(t)", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
  phi <- hzr_phase_cumhaz(t_grid, t_half = 3, nu = 2, m = 0, type = "cdf")
  expect_equal(phi, d$G)
})

test_that("hzr_phase_cumhaz('hazard') returns -log(1 - G(t))", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
  phi <- hzr_phase_cumhaz(t_grid, t_half = 3, nu = 2, m = 0, type = "hazard")
  expected <- -log(1 - d$G)
  expect_equal(phi, expected, tolerance = 1e-12)
})

test_that("hzr_phase_cumhaz('constant') returns t", {
  phi <- hzr_phase_cumhaz(t_grid, type = "constant")
  expect_equal(phi, t_grid)
})

test_that("hzr_phase_cumhaz('hazard') is monotone increasing", {
  phi <- hzr_phase_cumhaz(t_grid, t_half = 3, nu = 2, m = 0, type = "hazard")
  expect_true(all(diff(phi) >= -1e-12))
})

test_that("hzr_phase_cumhaz('cdf') is bounded [0, 1]", {
  phi <- hzr_phase_cumhaz(t_grid, t_half = 3, nu = 2, m = 0, type = "cdf")
  expect_true(all(phi >= 0 & phi <= 1))
})


# ============================================================================
# hzr_phase_hazard
# ============================================================================

test_that("hzr_phase_hazard('cdf') returns g(t)", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
  phi <- hzr_phase_hazard(t_grid, t_half = 3, nu = 2, m = 0, type = "cdf")
  expect_equal(phi, d$g)
})

test_that("hzr_phase_hazard('hazard') returns h(t)", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
  phi <- hzr_phase_hazard(t_grid, t_half = 3, nu = 2, m = 0, type = "hazard")
  expect_equal(phi, d$h)
})

test_that("hzr_phase_hazard('constant') returns vector of 1s", {
  phi <- hzr_phase_hazard(t_grid, type = "constant")
  expect_equal(phi, rep(1, length(t_grid)))
})


# ============================================================================
# Golden reference values (originally validated against mixhazard::decompos())
# ============================================================================
# These golden values were generated from the closed-form decomposition
# equations and cross-validated against the mixhazard package's decompos()
# function at 1e-12 tolerance.  They are now hardcoded so the test suite
# has no external dependency.

test_that("hzr_decompos matches golden reference values across all 6 cases", {
  snap_times <- c(0.1, 1.0, 3.0, 5.0, 10.0)

  golden <- list(
    # Case 1: m>0, nu>0 (standard sigmoidal)
    list(t_half = 3, nu =  2, m =  1,
         G = c(1.543870887948849e-01, 3.660254037844386e-01, 5.000000000000000e-01,
               5.635083268962916e-01, 6.461106321354770e-01),
         g = c(6.527585780416260e-01, 1.160254037844386e-01, 4.166666666666666e-02,
               2.459666924148338e-02, 1.143258415884856e-02),
         h = c(7.719354439744243e-01, 1.830127018922192e-01, 8.333333333333333e-02,
               5.635083268962916e-02, 3.230553160677384e-02)),

    # Case 1L: m=0, nu>0 (Weibull-like)
    list(t_half = 3, nu =  2, m =  0,
         G = c(2.244867998231077e-02, 3.010237439309285e-01, 5.000000000000000e-01,
               5.845520232324418e-01, 6.840991973811010e-01),
         g = c(4.261347015149480e-01, 1.806994562245466e-01, 5.776226504666210e-02,
               3.138515329720806e-02, 1.298599327498647e-02),
         h = c(4.359205422659926e-01, 2.585201638189699e-01, 1.155245300933242e-01,
               7.554532709824208e-02, 4.110781981979546e-02)),

    # Case 2: m<0, nu>0 (heavy-tailed)
    list(t_half = 3, nu =  2, m = -1,
         G = c(4.653741075440776e-02, 2.928932188134524e-01, 5.000000000000000e-01,
               5.917517095361370e-01, 6.984886554222364e-01),
         g = c(4.333920860207237e-01, 1.767766952966369e-01, 6.250000000000000e-02,
               3.402069087198858e-02, 1.370506111717107e-02),
         h = c(4.545454545454545e-01, 2.500000000000000e-01, 1.250000000000000e-01,
               8.333333333333333e-02, 4.545454545454546e-02)),

    # Case 2L: m<0, nu=0 (exponential decay)
    list(t_half = 3, nu =  0, m = -1,
         G = c(2.284003156575409e-02, 2.062994740159002e-01, 5.000000000000000e-01,
               6.850197375262816e-01, 9.007874342519875e-01),
         g = c(2.257718923587476e-01, 1.833837605982748e-01, 1.155245300933242e-01,
               7.277589362189646e-02, 2.292297007478435e-02),
         h = c(2.310490601866484e-01, 2.310490601866484e-01, 2.310490601866484e-01,
               2.310490601866484e-01, 2.310490601866483e-01)),

    # Case 3: m>0, nu<0 (bounded cumulative)
    list(t_half = 3, nu = -2, m =  1,
         G = c(1.543870887948848e-01, 3.660254037844386e-01, 5.000000000000000e-01,
               5.635083268962915e-01, 6.461106321354770e-01),
         g = c(6.527585780416260e-01, 1.160254037844387e-01, 4.166666666666666e-02,
               2.459666924148339e-02, 1.143258415884856e-02),
         h = c(7.719354439744243e-01, 1.830127018922193e-01, 8.333333333333333e-02,
               5.635083268962917e-02, 3.230553160677384e-02)),

    # Case 3L: m=0, nu<0 (bounded exponential)
    list(t_half = 3, nu = -2, m =  0,
         G = c(1.188705972421941e-01, 3.298064389861991e-01, 5.000000000000000e-01,
               5.913307632637709e-01, 7.179039946690142e-01),
         g = c(5.575380754920625e-01, 1.341019487465793e-01, 5.776226504666211e-02,
               3.656973241347520e-02, 1.784973505866168e-02),
         h = c(6.327539107729806e-01, 2.000943556421572e-01, 1.155245300933242e-01,
               8.948491622597644e-02, 6.327539107729806e-02))
  )

  for (ref in golden) {
    ours <- hzr_decompos(
      snap_times,
      t_half = ref$t_half,
      nu = ref$nu,
      m = ref$m
    )
    lbl <- paste0("m=", ref$m, ", nu=", ref$nu)
    expect_equal(
      ours$G,
      ref$G,
      tolerance = 1e-12,
      label = paste0("G for ", lbl)
    )
    expect_equal(
      ours$g,
      ref$g,
      tolerance = 1e-12,
      label = paste0("g for ", lbl)
    )
    expect_equal(
      ours$h,
      ref$h,
      tolerance = 1e-12,
      label = paste0("h for ", lbl)
    )
  }
})


# ============================================================================
# Numerical stability edge cases
# ============================================================================

test_that("hzr_decompos handles very small times without error", {
  d <- hzr_decompos(1e-10, t_half = 3, nu = 2, m = 0)
  expect_true(is.finite(d$G))
  expect_true(is.finite(d$g))
  expect_true(is.finite(d$h))
})

test_that("hzr_decompos handles very large times without error", {
  d <- hzr_decompos(1e6, t_half = 3, nu = 2, m = 0)
  expect_true(is.finite(d$G))
  # g and h may be 0 or Inf at extreme times, but should not be NaN
  expect_false(is.nan(d$G))
})

test_that("hzr_decompos handles time = 0 gracefully (clamped)", {
  d <- hzr_decompos(0, t_half = 3, nu = 2, m = 0)
  expect_true(is.finite(d$G))
  expect_true(d$G >= 0)
})
