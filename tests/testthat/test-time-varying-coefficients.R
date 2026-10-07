# This file predicts from models built with fit = FALSE on purpose, so the
# warning that those numbers come from starting values is switched off for
# this file only (#398). A file that does not expect the warning sees it as
# an ordinary leaked warning.
withr::local_options(TemporalHazard.warn_unfitted_prediction = FALSE)

library(testthat)

test_that("time-varying design expansion dimensions and names are correct", {
  x <- cbind(x1 = c(1, 2, 3, 4), x2 = c(2, 2, 2, 2))
  time <- c(0.5, 1.5, 2.5, 3.5)
  windows <- c(1, 2)

  x_tv <- .hzr_expand_time_varying_design(x, time, windows)

  expect_equal(dim(x_tv), c(4, 6))
  expect_equal(colnames(x_tv), c("x1_w1", "x2_w1", "x1_w2", "x2_w2", "x1_w3", "x2_w3"))

  expect_equal(unname(x_tv[1, ]), c(1, 2, 0, 0, 0, 0))
  expect_equal(unname(x_tv[2, ]), c(0, 0, 2, 2, 0, 0))
  expect_equal(unname(x_tv[3, ]), c(0, 0, 0, 0, 3, 2))
})

test_that("predict linear_predictor applies window-specific coefficients", {
  time <- c(0.5, 1.5, 2.5)
  x <- cbind(x1 = c(1, 1, 1))

  # exponential: theta = [log_lambda, beta_w1, beta_w2, beta_w3]
  fit <- hazard(
    time = time,
    status = c(1, 1, 1),
    x = x,
    time_windows = c(1, 2),
    theta = c(log(0.2), 0.0, 1.0, 2.0),
    dist = "exponential",
    fit = FALSE
  )

  pred <- predict(
    fit,
    newdata = data.frame(time = time, x1 = c(1, 1, 1)),
    type = "linear_predictor"
  )

  expect_equal(pred, c(0, 1, 2), tolerance = 1e-12)
})

test_that("predict hazard and cumulative_hazard remain coherent with time-varying coefficients", {
  time <- c(0.5, 1.5, 2.5)
  x <- cbind(x1 = c(1, 1, 1))

  fit <- hazard(
    time = time,
    status = c(1, 1, 1),
    x = x,
    time_windows = c(1, 2),
    theta = c(log(0.3), 0.0, 0.5, 1.0),
    dist = "exponential",
    fit = FALSE
  )

  lp <- predict(fit, newdata = data.frame(time = time, x1 = 1), type = "linear_predictor")
  hz_mult <- predict(fit, newdata = data.frame(time = time, x1 = 1), type = "hazard")
  ch <- predict(fit, newdata = data.frame(time = time, x1 = 1), type = "cumulative_hazard")
  s <- predict(fit, newdata = data.frame(time = time, x1 = 1), type = "survival")

  expect_equal(hz_mult, exp(lp), tolerance = 1e-12)

  lambda <- 0.3
  expect_equal(ch, lambda * time * exp(lp), tolerance = 1e-10)
  expect_equal(s, exp(-ch), tolerance = 1e-12)
})

test_that("time-varying predictor requires time at prediction for newdata", {
  fit <- hazard(
    time = c(1, 2, 3),
    status = c(1, 1, 1),
    x = cbind(x1 = c(1, 1, 1)),
    time_windows = c(1.5),
    theta = c(log(0.2), 0.1, 0.2),
    dist = "exponential",
    fit = FALSE
  )

  expect_error(
    predict(fit, newdata = data.frame(x1 = c(1, 1)), type = "linear_predictor"),
    "require a 'time' column"
  )
})

test_that("time-varying fit runs for Weibull and keeps expanded coefficient count", {
  set.seed(1001)
  n <- 80
  time <- rexp(n, rate = 0.4) + 0.05
  status <- rbinom(n, 1, 0.7)
  x <- cbind(x1 = rnorm(n))

  fit <- hazard(
    time = time,
    status = status,
    x = x,
    time_windows = c(1.5, 3.0),
    # Weibull: [mu, nu, beta_w1, beta_w2, beta_w3]
    theta = c(0.5, 1.1, 0.0, 0.0, 0.0),
    dist = "weibull",
    fit = TRUE,
    control = list(maxit = 200)
  )

  expect_true(is.finite(fit$fit$objective))
  expect_equal(length(fit$fit$theta), 5)
})

test_that("theta length validation reflects expanded time-varying design", {
  expect_error(
    hazard(
      time = c(1, 2, 3),
      status = c(1, 0, 1),
      x = cbind(x1 = c(1, 2, 3)),
      time_windows = c(1.5, 2.5),
      # needs at least 3 coefficients for expanded design (plus shape, depending dist)
      theta = c(0.1, 0.2),
      dist = "exponential",
      fit = FALSE
    ),
    "required coefficients"
  )
})
