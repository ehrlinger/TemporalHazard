# `.hzr_g3_phase_derivatives()` says, in its own comment, that it takes
# "central differences for log_tau". It stepped tau LINEARLY instead, with an
# absolute floor of 1e-10. Below tau = 1e-10 the minus step would reach 0, so
# it fell back to a forward difference across 100 times the parameter, and
# dPhi/dlog_tau came back 99.4% wrong (#352). It now steps log_tau directly.
#
# The reference here is ANALYTIC, from the closed form
# G3 = (((t/tau)^gamma + 1)^(1/alpha) - 1)^eta, so it carries no step size of
# its own and cannot be wrong about the thing under test. It was cross-checked
# against numDeriv Richardson at d = 1e-4, 1e-3 and 1e-2, agreeing to ~1e-10;
# numDeriv's DEFAULT d is 0.1 and is not truth (#332).

g3_ref_dlogtau <- function(time, tau, gamma, alpha, eta) {
  u <- gamma * (log(time) - log(tau))
  inner <- log1p(exp(u)) / alpha
  g3 <- exp(eta * log(expm1(inner)))
  g3 * eta * (-1 / expm1(-inner)) * (1 / (1 + exp(-u))) / alpha * (-gamma)
}

test_that("dPhi/dlog_tau is right where the step floor used to dominate", {
  # tau = 1e-12 is the shape that measured 0.9944 relative error: the floor
  # was 1e-10, a hundred times tau, and the fallback made it one-sided.
  for (tau in c(1e-12, 1e-11, 1e-10, 1e-9, 1e-8)) {
    got <- .hzr_g3_phase_derivatives(tau, tau = tau, gamma = 4, alpha = 1.5,
                                     eta = 0.5)$dPhi_dlog_tau
    want <- g3_ref_dlogtau(tau, tau, 4, 1.5, 0.5)
    expect_equal(got / want, 1, tolerance = 1e-6)
  }
})

test_that("ordinary tau is unchanged to the precision anyone relies on", {
  # A guard against the fix being a licence to move well-conditioned fits.
  # NOTE the margin: at tau = 1e-6 the OLD code was 1.3e-8 off, so this loop
  # fails on the previous implementation at that shape by a factor of only
  # 1.3. It is a discriminator there, not merely a guard, and widening the
  # tolerance past 2e-8 would silently stop it discriminating. The ~1e-10
  # regime the NEWS bullet describes begins around tau = 1e-5, not 1e-6.
  for (tau in c(1e-6, 1e-2, 1, 100)) {
    got <- .hzr_g3_phase_derivatives(tau, tau = tau, gamma = 4, alpha = 1.5,
                                     eta = 0.5)$dPhi_dlog_tau
    want <- g3_ref_dlogtau(tau, tau, 4, 1.5, 0.5)
    expect_equal(got / want, 1, tolerance = 1e-8)
  }
})

test_that("tau = 0 is handled by the base evaluation, not by a step fallback", {
  # The forward-difference fallback is DELETED. It was reachable only where
  # tau - max(|tau| * 1e-5, 1e-10) <= 0; under a log step tau * exp(-h) is
  # positive for every positive tau, down to and including denormals, so the
  # branch became unreachable. tau == 0 exactly is the one remaining boundary
  # and it is checked here rather than assumed: hzr_decompos_g3() already
  # returns G3 = Inf for tau <= 0, so the base evaluation is degenerate before
  # any step is taken and the deleted branch could not have rescued it.
  d <- .hzr_g3_phase_derivatives(1, tau = 0, gamma = 4, alpha = 1.5, eta = 0.5)
  expect_false(is.finite(d$Phi[[1L]]))
  expect_false(any(is.finite(d$dPhi_dlog_tau)))
})

test_that("a step that cannot move tau reports NaN, not a plausible zero", {
  # A MULTIPLICATIVE step needs tau to have bits left to move, and in the
  # denormal range it does not: 5e-324 IS the smallest positive double, and
  # both 5e-324 * exp(1e-5) and 5e-324 * exp(-1e-5) round back to it. The two
  # evaluation points coincide and the quotient is 0/(2h).
  #
  # An earlier version of this file asserted `expect_gt(5e-324 * exp(-1e-5),
  # 0)` and called that a boundary check. It is a statement about IEEE-754
  # that no change to this package can falsify, and it is true PRECISELY
  # BECAUSE the product equals tau -- so it certified the boundary while the
  # code returned 0 for a derivative whose true value is about -1.8e126.
  # This asserts the VALUE instead.
  expect_true(identical(5e-324 * exp(1e-5), 5e-324 * exp(-1e-5)))
  for (tau in c(5e-324, 1e-320)) {
    d <- .hzr_g3_phase_derivatives(tau, tau = tau, gamma = 4, alpha = 1.5,
                                   eta = 3)
    expect_true(is.nan(d$dPhi_dlog_tau[[1L]]))
    expect_true(is.nan(d$dphi_dlog_tau[[1L]]))
  }
  # tau = Inf coincides the same way and took the same silent-zero path.
  d <- .hzr_g3_phase_derivatives(1, tau = Inf, gamma = 4, alpha = 1.5, eta = 3)
  expect_true(is.nan(d$dPhi_dlog_tau[[1L]]))
  # A shape just above the collapse still differentiates normally, so the
  # guard cannot be passing by rejecting everything.
  d <- .hzr_g3_phase_derivatives(1e-300, tau = 1e-300, gamma = 4, alpha = 1.5,
                                 eta = 3)
  expect_true(is.finite(d$dPhi_dlog_tau[[1L]]))
})

test_that("the gamma and eta steps are left alone", {
  # DELIBERATELY UNCHANGED. The review of #332 reported dPhi/dgamma 3.1% off
  # at gamma = 5e-11 and dPhi/deta 0.5% off at eta = 5e-11. Against the
  # analytic reference both measure clean -- 1.1e-6 and 2.0e-7 at worst --
  # because G3 is very nearly LINEAR in gamma and eta there, so a 1e-10 step
  # is small relative to the scale on which the FUNCTION varies even when it
  # is 100% of the parameter. This test fails if the fix is widened to them.
  ref_g <- function(tm, tau, g, a, e) {
    u <- g * (log(tm) - log(tau))
    inner <- log1p(exp(u)) / a
    exp(e * log(expm1(inner))) * e * (-1 / expm1(-inner)) *
      (1 / (1 + exp(-u))) / a * (log(tm) - log(tau))
  }
  ref_e <- function(tm, tau, g, a, e) {
    u <- g * (log(tm) - log(tau))
    inner <- log1p(exp(u)) / a
    exp(e * log(expm1(inner))) * log(expm1(inner))
  }
  # NOT time == tau: there ln(t/tau) = 0 and dPhi/dgamma is identically zero,
  # so a probe placed there cannot see the quantity it is measuring.
  tm <- 3e-2
  for (g in c(5e-11, 1e-9, 1e-6)) {
    got <- .hzr_g3_phase_derivatives(tm, tau = 1e-2, gamma = g, alpha = 1.5,
                                     eta = 0.5)$dPhi_dgamma
    expect_equal(got / ref_g(tm, 1e-2, g, 1.5, 0.5), 1, tolerance = 1e-5)
  }
  for (e in c(5e-11, 1e-9, 1e-6)) {
    got <- .hzr_g3_phase_derivatives(tm, tau = 1e-2, gamma = 4, alpha = 1.5,
                                     eta = e)$dPhi_deta
    expect_equal(got / ref_e(tm, 1e-2, 4, 1.5, e), 1, tolerance = 1e-5)
  }
})
