test_that("a G3 constrained at alpha = eta = 1 is a Weibull cumulative hazard", {
  # PARMS ... WEIBULL is setopt(6) -> SETG3_weibull, g3flag += 2 in the
  # reference implementation. The R general form
  #   (((t/tau)^gamma + 1)^(1/alpha) - 1)^eta
  # collapses at alpha = eta = 1 to (t/tau)^gamma. This is the most common
  # production configuration: 76-88 blocks in every study profiled.
  t <- c(0.01, 0.5, 1, 5, 50, 180)
  tau <- 2.5
  gamma <- 1.4
  got <- hzr_decompos_g3(t, tau = tau, gamma = gamma, alpha = 1, eta = 1)
  expect_equal(got$G3, (t / tau)^gamma, tolerance = 1e-12)
})

test_that("the constrained hazard matches the Weibull hazard", {
  t <- c(0.5, 5, 50)
  tau <- 2.5
  gamma <- 1.4
  got <- hzr_decompos_g3(t, tau = tau, gamma = gamma, alpha = 1, eta = 1)
  expect_equal(got$g3, (gamma / tau) * (t / tau)^(gamma - 1),
               tolerance = 1e-12)
})

test_that("the constrained path survives extreme times", {
  # HAZARD branches on g3flag partly for numerical reasons: at the constraint
  # the general path computes expm1(log1p(exp(y))) where the answer is exp(y).
  # Measured worst case is 5e-13 relative, so the tolerance is set just wide
  # enough to admit it and no wider -- a looser one would stop detecting drift.
  t <- c(1e-6, 1e-3, 1e3, 1e6)
  got <- hzr_decompos_g3(t, tau = 1, gamma = 2, alpha = 1, eta = 1)
  expect_true(all(is.finite(got$G3)))
  expect_equal(got$G3, t^2, tolerance = 1e-11)
})

test_that("the comparison is not trivially satisfied", {
  # A max discrepancy of exactly zero across every case would more likely mean
  # nothing was compared than that the identity is exact. Assert the general
  # form differs from the Weibull form away from the constraint.
  t <- c(0.5, 5, 50)
  general <- hzr_decompos_g3(t, tau = 2.5, gamma = 1.4, alpha = 2, eta = 1.5)
  expect_false(isTRUE(all.equal(general$G3, (t / 2.5)^1.4)))
})
