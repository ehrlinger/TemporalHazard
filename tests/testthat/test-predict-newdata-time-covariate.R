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

.tc_msg <- "variable named 'time'"

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

test_that("a phase formula's constant named time is refused: newdata masks it", {
  # model.matrix() looks symbols up in newdata first, so newdata's `time`
  # would silently replace the formula constant: I(30 > 1) instead of
  # I(30 > 50).  Refusing is the only safe answer.
  d <- .tc_avc
  d$time <- NULL
  time <- 50
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant", formula = ~ I(age > time))),
    fit = TRUE))
  expect_error(predict(fit, newdata = data.frame(time = c(1, 2), age = 30),
                       type = "cumulative_hazard"), .tc_msg)
})

test_that("a global formula's constant named time is masked the same way", {
  d <- .tc_avc
  d$time <- NULL
  time <- 50
  th <- c(mu = 0.01, nu = 0.5, b_gt = 0.4, b_mal = 0.3)
  w <- hazard(survival::Surv(int_dead, dead) ~ I(age > time) + mal, data = d,
              dist = "weibull", theta = th)
  # Time-based: newdata's time is the prediction time, so refuse.
  expect_error(predict(w, newdata = data.frame(time = 2, age = 30, mal = 1),
                       type = "cumulative_hazard"), .tc_msg)
  # Eta-based with a `time` column: it would mask the constant, so refuse.
  expect_error(predict(w, newdata = data.frame(time = 2, age = 30, mal = 1),
                       type = "linear_predictor"), .tc_msg)
  # Eta-based without one: the constant is used, I(30 > 50) = FALSE.
  expect_equal(predict(w, newdata = data.frame(age = c(30, 60), mal = 1),
                       type = "linear_predictor"),
               c(0, 0.4) + 0.3, tolerance = 1e-12)
})

test_that("design-column newdata never re-evaluates a time constant", {
  # hzr_gof() and hzr_deciles() pass the fitted design columns, which are
  # taken by name, so a formula constant named `time` is never looked up
  # and nothing is masked.  They must match a fit written with the literal.
  d <- .tc_avc
  d$time <- NULL
  time <- 50
  th <- c(mu = 0.01, nu = 0.5, 0, 0)
  w <- hazard(survival::Surv(int_dead, dead) ~ I(age > time) + mal, data = d,
              dist = "weibull", theta = th, fit = TRUE)
  lit <- hazard(survival::Surv(int_dead, dead) ~ I(age > 50) + mal, data = d,
                dist = "weibull", theta = th, fit = TRUE)
  expect_equal(unname(w$fit$theta), unname(lit$fit$theta), tolerance = 1e-10)
  expect_equal(hzr_gof(w)$par_cumhaz, hzr_gof(lit)$par_cumhaz,
               tolerance = 1e-10)
  expect_equal(hzr_deciles(w, time = 60)$expected,
               hzr_deciles(lit, time = 60)$expected, tolerance = 1e-10)
})

test_that("hzr_gof() on log(time) uses the design column, not newdata's time", {
  # The design column is `log(time)`, taken by name and never evaluated,
  # so follow-up time cannot overwrite it: it must match a refit on the
  # precomputed column.
  d <- .tc_avc
  d$la <- log(d$time)
  th <- c(mu = 0.01, nu = 0.5, 0, 0)
  w <- hazard(survival::Surv(int_dead, dead) ~ log(time) + mal, data = d,
              dist = "weibull", theta = th, fit = TRUE)
  ref <- hazard(survival::Surv(int_dead, dead) ~ la + mal, data = d,
                dist = "weibull", theta = th, fit = TRUE)
  expect_equal(unname(w$fit$theta), unname(ref$fit$theta), tolerance = 1e-10)
  expect_equal(hzr_gof(w)$par_cumhaz, hzr_gof(ref)$par_cumhaz,
               tolerance = 1e-10)
})

test_that("a constant baked into predvars is not read from newdata", {
  # scale(age, center = time) stores the centre at fit time, so `time` is
  # never evaluated at newdata.
  d <- .tc_avc
  d$time <- NULL
  time <- 50
  th <- c(mu = 0.01, nu = 0.5, b = 0.01)
  w <- hazard(survival::Surv(int_dead, dead) ~
                scale(age, center = time, scale = FALSE),
              data = d, dist = "weibull", theta = th)
  expect_equal(unname(predict(w, newdata = data.frame(time = 2, age = 60),
                              type = "cumulative_hazard")),
               (0.01 * 2)^0.5 * exp(0.01 * (60 - 50)), tolerance = 1e-12)
})

test_that("a phase formula still counts when the global design is design-level", {
  # Global I(age > time) arrives as design columns and is not evaluated,
  # but the phase formula I(mal > time) is, and newdata's `time` would
  # mask its constant: this must stop (2e71588 returned masked values).
  d <- .tc_avc
  d$time <- NULL
  time <- 0.5
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ I(age > time), data = d,
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant", formula = ~ I(mal > time))),
    fit = TRUE))
  nd <- data.frame(time = c(1, 2), check.names = FALSE,
                   `I(age > time)TRUE` = 1, mal = 1)
  expect_error(predict(fit, newdata = nd, type = "cumulative_hazard"),
               .tc_msg)
})

test_that("a list element named time (cfg$time) is not the variable time", {
  d <- .tc_avc
  d$time <- NULL
  cfg <- list(time = 50)
  th <- c(mu = 0.01, nu = 0.5, b_gt = 0.4, b_mal = 0.3)
  w <- hazard(survival::Surv(int_dead, dead) ~ I(age > cfg$time) + mal,
              data = d, dist = "weibull", theta = th)
  nd <- data.frame(time = c(1, 2), age = 30, mal = 1)
  # I(30 > 50) is FALSE, so only mal contributes.
  expect_equal(predict(w, newdata = nd, type = "linear_predictor"),
               c(0.3, 0.3), tolerance = 1e-12)
  expect_equal(predict(w, newdata = nd, type = "cumulative_hazard"),
               (0.01 * c(1, 2))^0.5 * exp(0.3), tolerance = 1e-12)
})

test_that("time as the only covariate works for linear_predictor and hazard", {
  w <- hazard(survival::Surv(int_dead, dead) ~ time, data = .tc_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, b = 0.002))
  nd <- data.frame(time = c(2, 50))
  expect_equal(predict(w, newdata = nd, type = "linear_predictor"),
               0.002 * c(2, 50), tolerance = 1e-12)
  expect_equal(predict(w, newdata = nd, type = "hazard"),
               exp(0.002 * c(2, 50)), tolerance = 1e-12)
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
