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

test_that("a single-distribution covariate named time stops the time-based types", {
  w <- hazard(survival::Surv(int_dead, dead) ~ time, data = .tc_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, b = 0.002))
  for (type in c("cumulative_hazard", "survival")) {
    expect_error(predict(w, newdata = data.frame(time = 2), type = type),
                 .tc_msg, label = type)
  }
  # With another covariate beside it, too.
  w2 <- hazard(survival::Surv(int_dead, dead) ~ time + mal, data = .tc_avc,
               dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.002, 0.3))
  expect_error(predict(w2, newdata = data.frame(time = 2, mal = 1),
                       type = "cumulative_hazard"), .tc_msg)
})

test_that("linear_predictor and hazard still read time as the covariate", {
  # These types have no prediction time, so `time` in newdata can only be
  # the covariate, and the answer was right before #270.  It must stay.
  w2 <- hazard(survival::Surv(int_dead, dead) ~ time + mal, data = .tc_avc,
               dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.002, 0.3))
  nd <- data.frame(time = c(2, 50), mal = c(1, 0))
  eta <- 0.002 * c(2, 50) + 0.3 * c(1, 0)
  expect_equal(predict(w2, newdata = nd, type = "linear_predictor"), eta,
               tolerance = 1e-12)
  expect_equal(predict(w2, newdata = nd, type = "hazard"), exp(eta),
               tolerance = 1e-12)
})

test_that("with time windows, time is the prediction time for every type", {
  w <- hazard(survival::Surv(int_dead, dead) ~ time + mal, data = .tc_avc,
              dist = "weibull", time_windows = 3,
              theta = c(mu = 0.01, nu = 0.5, 0.002, 0.3, 0.001, 0.2))
  expect_error(predict(w, newdata = data.frame(time = 2, mal = 1),
                       type = "linear_predictor"), .tc_msg)
})

test_that("a phase formula's constant named time is not a covariate", {
  # `time` here is a value in the formula's environment, not a data column.
  d <- .tc_avc
  d$time <- NULL
  time <- 50
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant", formula = ~ I(age > time))),
    fit = TRUE))
  got <- predict(fit, newdata = data.frame(time = c(1, 2), age = 60),
                 type = "cumulative_hazard")
  expect_true(all(is.finite(got) & got > 0))
})

test_that("hzr_deciles() and hzr_gof() stop on such a model, saying why", {
  # They evaluate the model at new rows through predict(newdata = ), where
  # follow-up time used to overwrite the covariate.  Stopping is the fix;
  # the message must name them, since the user passed no newdata.
  w <- hazard(survival::Surv(int_dead, dead) ~ time + mal, data = .tc_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0, 0),
              fit = TRUE)
  expect_error(hzr_deciles(w, time = 60), "hzr_deciles\\(\\)")
  expect_error(hzr_gof(w), "hzr_gof\\(\\)")
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
