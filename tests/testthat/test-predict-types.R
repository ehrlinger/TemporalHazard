# This file predicts from models built with fit = FALSE on purpose, so the
# warning that those numbers come from starting values is switched off for
# this file only (#398). A file that does not expect the warning sees it as
# an ordinary leaked warning.
withr::local_options(TemporalHazard.warn_unfitted_prediction = FALSE)

test_that("predict() supports survival type", {
  # Create a simple fitted model
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)

  fit <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(mu = 0.5, nu = 1.2),
    dist = "weibull",
    fit = FALSE
  )

  # Predict survival at original times
  surv <- predict(fit, type = "survival")

  expect_true(is.numeric(surv))
  expect_equal(length(surv), length(time))
  expect_true(all(surv > 0 & surv <= 1))  # Survival probability in (0, 1]
})

test_that("predict() supports cumulative_hazard type", {
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)

  fit <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(mu = 0.5, nu = 1.2),
    dist = "weibull",
    fit = FALSE
  )

  # Predict cumulative hazard at original times
  cumhaz <- predict(fit, type = "cumulative_hazard")

  expect_true(is.numeric(cumhaz))
  expect_equal(length(cumhaz), length(time))
  expect_true(all(cumhaz > 0))  # Cumulative hazard positive
})

test_that("survival and cumulative hazard are related correctly", {
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)

  fit <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(mu = 0.5, nu = 1.2),
    dist = "weibull",
    fit = FALSE
  )

  surv <- predict(fit, type = "survival")
  cumhaz <- predict(fit, type = "cumulative_hazard")

  # S(t) = exp(-H(t))
  # So cumhaz = -log(surv)
  expected_cumhaz <- -log(surv)

  expect_equal(cumhaz, expected_cumhaz, tolerance = 1e-12)
})

test_that("predict() works with covariates for survival", {
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)
  x <- matrix(c(0.5, 1.2, -0.3, 0.8, 1.5), ncol = 1)

  fit <- hazard(
    time = time,
    status = status,
    x = x,
    theta = c(mu = 0.5, nu = 1.2, beta = 0.3),
    dist = "weibull",
    fit = FALSE
  )

  # Predict survival at original times and data
  newdata <- data.frame(time = time, X1 = x[, 1])
  surv <- predict(fit, newdata = newdata, type = "survival")

  expect_true(is.numeric(surv))
  expect_equal(length(surv), length(time))
  expect_true(all(surv > 0 & surv <= 1))
})

test_that("predict() handles newdata time column correctly", {
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)

  fit <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(mu = 0.5, nu = 1.2),
    dist = "weibull",
    fit = FALSE
  )

  # Predict at different times
  new_time <- c(0.5, 1.0, 1.5, 2.0, 2.5)
  newdata <- data.frame(time = new_time)

  surv_new <- predict(fit, newdata = newdata, type = "survival")

  expect_equal(length(surv_new), length(new_time))
  expect_true(all(surv_new > 0 & surv_new <= 1))
})

test_that("predict() validates time is positive", {
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)

  fit <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(mu = 0.5, nu = 1.2),
    dist = "weibull",
    fit = FALSE
  )

  # Negative time should error
  newdata <- data.frame(time = c(1, -1, 2))
  expect_error(predict(fit, newdata = newdata, type = "survival"))
})

test_that("predict() errors if survival requested on non-Weibull/exponential", {
  time <- c(1, 2, 3)
  status <- c(1, 0, 1)

  fit <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(0.5, 0.3),
    dist = "gompertz",  # Not yet supported
    fit = FALSE
  )

  expect_error(predict(fit, type = "survival"), "only supported for")
})

test_that("predict() linear_predictor and hazard still work unchanged", {
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)
  x <- matrix(c(0.5, 1.2, -0.3, 0.8, 1.5), ncol = 1)

  # theta must have: shape params (2 for Weibull) + covariate coefs (1 for 1 covariate)
  fit <- hazard(
    time = time,
    status = status,
    x = x,
    theta = c(0.5, 0.3, 0.1),  # mu, nu, beta
    dist = "weibull",
    fit = FALSE
  )

  newdata <- data.frame(X1 = c(0.2, 0.5, 0.8))

  eta <- predict(fit, newdata = newdata, type = "linear_predictor")
  hz <- predict(fit, newdata = newdata, type = "hazard")

  # Check relationship: hz = exp(eta)
  expect_equal(hz, exp(eta), tolerance = 1e-12)
})

test_that("predict() works on fitted model with survival type", {
  time <- c(1, 2, 3, 4, 5, 6)
  status <- c(1, 0, 1, 1, 0, 1)

  fit <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(mu = 0.3, nu = 1.0),
    dist = "weibull",
    fit = TRUE,
    control = list(maxit = 100)
  )

  # Check if fit object has theta (may fail to converge)
  if (is.null(fit$fit$theta)) {
    skip("Model did not converge")
  }

  # Predict survival from fitted model
  surv <- predict(fit, type = "survival")

  expect_true(is.numeric(surv))
  expect_equal(length(surv), length(time))
  expect_true(all(surv > 0 & surv <= 1))
})

test_that("monotonicity: survival decreases with time", {
  fit <- hazard(
    time = c(1, 2),
    status = c(1, 0),
    x = NULL,
    theta = c(mu = 0.5, nu = 1.5),
    dist = "weibull",
    fit = FALSE
  )

  # Predict survival across increasing times
  new_times <- seq(0.1, 5, by = 0.5)
  newdata <- data.frame(time = new_times)
  surv <- predict(fit, newdata = newdata, type = "survival")

  # Survival should be monotonically decreasing
  expect_true(all(diff(surv) < 0))
})

test_that("monotonicity: cumulative hazard increases with time", {
  fit <- hazard(
    time = c(1, 2),
    status = c(1, 0),
    x = NULL,
    theta = c(mu = 0.5, nu = 1.5),
    dist = "weibull",
    fit = FALSE
  )

  # Predict cumulative hazard across increasing times
  new_times <- seq(0.1, 5, by = 0.5)
  newdata <- data.frame(time = new_times)
  cumhaz <- predict(fit, newdata = newdata, type = "cumulative_hazard")

  # Cumulative hazard should be monotonically increasing
  expect_true(all(diff(cumhaz) > 0))
})

test_that("covariate effects on survival are correct", {
  time <- c(1, 2, 3)
  status <- c(1, 0, 1)

  # The model carries the covariate its beta belongs to. It was built with
  # x = NULL, a 3-entry theta for a 2-parameter model, which hazard() now
  # refuses (#375); predict() had applied the stray beta to a newdata column
  # the model never had.
  fit <- hazard(
    time = time,
    status = status,
    x = matrix(c(0, 1, 0), ncol = 1, dimnames = list(NULL, "X1")),
    theta = c(mu = 0.5, nu = 1.2, beta = 0.5),  # 3 parameters
    dist = "weibull",
    fit = FALSE
  )

  # Compare survival with and without covariate effect
  t_eval <- 2.0

  # No covariate effect
  newdata_no_x <- data.frame(time = t_eval)
  # surv without covariate = exp(-(0.5 * 2)^1.2)

  # With covariate effect (x=1): exp(eta) where eta = 0.5*1 = 0.5
  newdata_x1 <- data.frame(time = t_eval, X1 = 1)
  surv_x1 <- predict(fit, newdata = newdata_x1, type = "survival")

  # Higher covariate effect should reduce survival (positive beta increases hazard)
  # So with positive covariate, survival should be lower
  expect_true(is.numeric(surv_x1))
  expect_true(surv_x1 > 0 && surv_x1 < 1)
  # What the comments above describe: a positive covariate lowers survival.
  surv_x0 <- predict(fit, newdata = data.frame(time = t_eval, X1 = 0),
                     type = "survival")
  expect_lt(surv_x1, surv_x0)
})
