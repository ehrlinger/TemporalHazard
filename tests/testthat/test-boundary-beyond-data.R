# #683: a G1 phase whose rise begins after the last observed time. Its t_half
# is not identified; the phase is a power law in t over the data (#729 drifted
# to log t_half 275.8 with converged = TRUE).

theta_for <- function(log_t_half, nu = 1, m = 1) {
  c(early.log_mu = 0, early.log_t_half = log_t_half, early.nu = nu, early.m = m)
}

test_that("a cdf phase whose rise lies past the data is recorded", {
  rec <- .hzr_phase_beyond_data_record("early", "cdf", theta_for(275.8),
                                       t_obs = c(0.1, 5, 17))
  expect_identical(rec$mechanism, "phase_beyond_data")
  expect_identical(rec$phase, "early")
  expect_identical(rec$parameter, "t_half")
  expect_lt(rec$G_max, .hzr_step_resolution)
  expect_match(rec$detail, "after the data ends")
})

test_that("a t_half past the data whose rise is visible is NOT recorded", {
  # t_half 2x beyond the last time: G(17) is far from 0, the shape is seen.
  expect_null(.hzr_phase_beyond_data_record("early", "cdf", theta_for(log(34)),
                                            t_obs = c(0.1, 5, 17)))
})

test_that("an interior t_half is not recorded, and nor is a constant phase", {
  expect_null(.hzr_phase_beyond_data_record("early", "cdf", theta_for(log(3)),
                                            t_obs = c(0.1, 5, 17)))
  expect_null(.hzr_phase_beyond_data_record(
    "c", "constant", c(c.log_mu = 0), t_obs = c(0.1, 5, 17)))
})

test_that("a t_half that overflows exp() is still recorded", {
  rec <- .hzr_phase_beyond_data_record("early", "cdf", theta_for(800),
                                       t_obs = c(0.1, 5, 17))
  expect_identical(rec$mechanism, "phase_beyond_data")
})

test_that("the boundary check carries the record into $boundary", {
  out <- .hzr_boundary_check_impl(
    theta = theta_for(275.8),
    phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 1)),
    time = c(0.1, 5, 17), fitted = TRUE)
  mech <- vapply(out$boundary, `[[`, "", "mechanism")
  expect_true("phase_beyond_data" %in% mech)
})
