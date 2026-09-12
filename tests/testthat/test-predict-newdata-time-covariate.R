# A covariate named `time` (#270).  In `newdata` the column `time` is the
# prediction time, so a model with a covariate of that name cannot be told
# its value there: the column was dropped from the covariates and the
# time-based types silently returned the baseline.  predict(newdata = ) now
# stops and says so.  predict() without newdata is unaffected.
#
# The check runs before any prediction, so unfitted objects with a supplied
# theta are enough for the error cases.

.tc_avc <- local({
  data(avc, package = "TemporalHazard")
  d <- na.omit(avc)
  d$time <- d$age          # a covariate that happens to be called `time`
  d
})

.tc_msg <- "covariate named 'time'"

test_that("a single-distribution covariate named time stops predict(newdata)", {
  w <- hazard(survival::Surv(int_dead, dead) ~ time, data = .tc_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, b = 0.002))
  for (type in c("cumulative_hazard", "survival", "hazard",
                 "linear_predictor")) {
    expect_error(predict(w, newdata = data.frame(time = 2), type = type),
                 .tc_msg, label = type)
  }
  # With another covariate beside it, too.
  w2 <- hazard(survival::Surv(int_dead, dead) ~ time + mal, data = .tc_avc,
               dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.002, 0.3))
  expect_error(predict(w2, newdata = data.frame(time = 2, mal = 1),
                       type = "cumulative_hazard"), .tc_msg)
})

test_that("a vector-interface x column named time stops predict(newdata)", {
  d <- .tc_avc
  v <- hazard(time = d$int_dead, status = d$dead, x = cbind(time = d$age),
              dist = "exponential", theta = c(log_lambda = -4, b = 0.002))
  expect_error(predict(v, newdata = data.frame(time = 2),
                       type = "cumulative_hazard"), .tc_msg)
})

test_that("a multiphase global or phase covariate named time stops it too", {
  # Multiphase needs a fit to have coefficients.  Its Hessian warnings
  # concern standard errors, which these predictions do not use.
  glob <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ time, data = .tc_avc,
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE))
  expect_error(predict(glob, newdata = data.frame(time = 2),
                       type = "cumulative_hazard"), .tc_msg)

  phase <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = .tc_avc,
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant", formula = ~ time)),
    fit = TRUE))
  expect_error(predict(phase, newdata = data.frame(time = 2),
                       type = "cumulative_hazard"), .tc_msg)
})

test_that("predict() without newdata still works for such a model", {
  w <- hazard(survival::Surv(int_dead, dead) ~ time, data = .tc_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, b = 0.002))
  d <- .tc_avc
  # Each subject at its own follow-up and its own covariate.
  want <- (0.01 * d$int_dead)^0.5 * exp(0.002 * d$time)
  expect_equal(unname(predict(w, type = "cumulative_hazard")), want,
               tolerance = 1e-12)
})

test_that("a model without a time covariate is unaffected", {
  w <- hazard(survival::Surv(int_dead, dead) ~ mal, data = .tc_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, b = 0.3))
  # A time-only newdata is still the baseline, and mal = 1 still moves it.
  expect_equal(predict(w, newdata = data.frame(time = c(1, 4)),
                       type = "cumulative_hazard"),
               (0.01 * c(1, 4))^0.5, tolerance = 1e-12)
  expect_equal(predict(w, newdata = data.frame(time = c(1, 4), mal = 1),
                       type = "cumulative_hazard"),
               (0.01 * c(1, 4))^0.5 * exp(0.3), tolerance = 1e-12)
})
