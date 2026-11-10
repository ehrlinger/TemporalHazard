# test-coe-objective-consistency.R -- the reported objective must be the
# log-likelihood of the parameters the fit returns (#362).
#
# Under Conservation of Events the conserved phase's scale is adjusted after
# the optimizer returns, so `value` and `par` could describe different points.
# On the avc fit below they did: the fit reported -71.934410193 while the
# likelihood at its own returned theta is -78.1497959154, a gap of 6.215386.
# These tests fail on main, where the first assertion reads -71.934410193.
#
# The gap itself was a symptom. That fit stands on a likelihood discontinuity
# (#448): `early.nu` reaches -1.4e-16, the cdf phase's exponent 1/nu is
# -6.9e15, and `t_half` lands one ulp from five tied event times, so a one-ulp
# parameter change moves the log-likelihood by 6 to 30 units. This file pins
# the REPORT, not the soundness of the fit.

# #683: hazard() no longer returns that fit. Its nu ran to the step limit,
# and the limit layer refuses it. The fixture below is the same fit made
# reachable the way a user makes it: nu (and m) FIXED at the value the
# drifted fit reached, as SAS's FIXNU does, and every other parameter started
# from that fit's estimates. A fixed parameter never triggers the layer, so
# the fit returns ($fit$limit NULL) while still standing on the
# discontinuity, which is what opens the gap: measured, the reported
# objective is -75.5314090483 against the optimizer's -69.3233188349, a gap
# of 6.208, the same mechanism as the 6.215 above.
coe_448_theta <- c(-6.0546138418029578, -2.7225076583067618,
                   -1.4391540373781222e-16, 1, -1.2800810419808482,
                   6.3008368357595135, 1.1585115564266042,
                   9.7484245868043898e-02, 1.6829222863987869e-01)

coe_448_phases <- function() {
  # alpha = 1 is replaced by the value the constraint derives, and says so.
  # That warning is fixture construction, not the behaviour under test, and
  # it is silenced here for the same reason the fit call below is: so the
  # suite's warning multiset carries only warnings a test is actually about.
  suppressWarnings(list(
    early = hzr_phase("cdf", t_half = exp(coe_448_theta[[2L]]),
                      nu = coe_448_theta[[3L]], m = 1, fixed = c("nu", "m")),
    late  = hzr_phase("g3", tau = 5, gamma = 1, alpha = 1, eta = 1,
                      constraint = "alpha_gamma_eta")
  ))
}

coe_448_fit <- function(avc, phases) {
  suppressWarnings(hazard(
    time = avc$int_dead, status = avc$dead, dist = "multiphase",
    phases = phases, theta = coe_448_theta, fit = TRUE,
    control = list(n_starts = 1L)
  ))
}

coe_ll_at <- function(fit, time, status, phases) {
  .hzr_logl_multiphase(
    unname(fit$fit$theta), time, status,
    phases = .hzr_validate_phases(phases),
    covariate_counts = c(early = 0L, late = 0L),
    x_list = list(early = NULL, late = NULL)
  )
}

test_that("a CoE fit reports the log-likelihood of the theta it returns (#362)", {
  skip_on_cran()
  data(avc, package = "TemporalHazard")
  phases <- coe_448_phases()
  fit <- coe_448_fit(avc, phases)
  # Premises: without these the fixture is no longer the case under test.
  # The fixture is clean by the limit layer's reading, not by luck.
  expect_null(fit$fit$limit)
  expect_true(isTRUE(fit$fit$converged))
  expect_true(isTRUE(fit$spec$control$conserve_applied))

  independent <- coe_ll_at(fit, avc$int_dead, avc$dead, phases)
  expect_equal(fit$fit$objective, independent, tolerance = 1e-10)
  # The value this fixture reported before the fix, asserted as an inequality
  # so the test states what it is defending against.
  expect_false(isTRUE(all.equal(fit$fit$objective, -71.934410193,
                                tolerance = 1e-8)))
  # summary() and print() read the same field, so they move with it.
  expect_equal(summary(fit)$log_lik, independent, tolerance = 1e-10)
})

test_that("an ordinary CoE fit is unchanged: objective still equals its own likelihood", {
  # The normal path, where the adjustment was already at a fixed point. This
  # guards against the recomputation perturbing fits that were never wrong.
  set.seed(1)
  n <- 120
  tt <- stats::rexp(n) + 0.01
  st <- stats::rbinom(n, 1, 0.75)
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
                 late  = hzr_phase("constant"))
  fit <- suppressWarnings(hazard(
    time = tt, status = st, dist = "multiphase", phases = phases,
    fit = TRUE, control = list(n_starts = 1L)
  ))
  expect_true(isTRUE(fit$spec$control$conserve_applied))
  expect_equal(fit$fit$objective, coe_ll_at(fit, tt, st, phases),
               tolerance = 1e-10)
})

test_that("starts$objective keeps the optimizer's own value, and says so (#362)", {
  skip_on_cran()
  # The two fields answer different questions and can legitimately differ.
  # `starts$objective` records what each start's optimisation REACHED, on one
  # footing across rows, and `best` is the argmax of that column.
  # `fit$fit$objective` is the likelihood of the estimates actually returned,
  # after the conservation adjustment has moved the conserved scale. Rewriting
  # only the winning row to match would break the argmax invariant that
  # test-multiphase-reproducibility.R asserts, so the difference is kept and
  # pinned here rather than left to be discovered.
  data(avc, package = "TemporalHazard")
  fit <- coe_448_fit(avc, coe_448_phases())
  expect_null(fit$fit$limit)
  starts <- fit$fit$starts
  expect_identical(sum(starts$best), 1L)
  # The argmax invariant still holds over the column.
  expect_equal(starts$objective[starts$best], max(starts$objective))
  # And on this fit the winning row is the PRE-adjustment value, so it differs
  # from the reported objective by the size of the conservation adjustment.
  expect_false(isTRUE(all.equal(starts$objective[starts$best],
                                fit$fit$objective)))
  expect_gt(abs(starts$objective[starts$best] - fit$fit$objective), 1)
})
