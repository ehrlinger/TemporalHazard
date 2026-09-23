# This file predicts from models built with fit = FALSE on purpose, so the
# warning that those numbers come from starting values is switched off for
# this file only (#398). A file that does not expect the warning sees it as
# an ordinary leaked warning.
withr::local_options(TemporalHazard.warn_unfitted_prediction = FALSE)

test_that("exponential gradient computes correct shape for univariate case", {
  # Simple univariate exponential (no covariates)
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)

  # log(lambda) is the parameter
  theta <- c(log_lambda = log(0.5))

  # Compute gradient analytically
  grad <- .hzr_gradient_exponential(
    theta = theta,
    time = time,
    status = status,
    x = NULL
  )

  expect_length(grad, 1)
  expect_true(all(is.finite(grad)))
})

test_that("exponential gradient computes correct shape with covariates", {
  # Exponential with one covariate
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)
  x <- matrix(c(0.5, 1.2, -0.3, 0.8, 1.5), ncol = 1)

  theta <- c(log_lambda = log(0.5), beta = 0.3)

  # Compute gradient analytically
  grad <- .hzr_gradient_exponential(
    theta = theta,
    time = time,
    status = status,
    x = x
  )

  expect_length(grad, 2)
  expect_true(all(is.finite(grad)))
})

test_that("exponential gradient matches numerical gradient", {
  skip_if_not_installed("numDeriv")

  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)
  x <- matrix(c(0.5, 1.2, -0.3, 0.8, 1.5), ncol = 1)

  theta <- c(log_lambda = log(0.5), beta = 0.3)

  # Analytical gradient
  grad_analytical <- .hzr_gradient_exponential(
    theta = theta,
    time = time,
    status = status,
    x = x
  )

  # Numerical gradient via numDeriv
  objective <- function(th) {
    .hzr_logl_exponential(
      theta = th,
      time = time,
      status = status,
      x = x,
      return_gradient = FALSE
    )
  }

  grad_numerical <- numDeriv::grad(objective, theta)

  # Should match to machine precision
  expect_equal(grad_analytical, grad_numerical, tolerance = 1e-4)
})

test_that("exponential optimizer converges", {
  time <- c(1, 2, 3, 4, 5, 6)
  status <- c(1, 0, 1, 1, 0, 1)

  theta_start <- c(log_lambda = log(0.5))

  result <- .hzr_optim_exponential(
    time = time,
    status = status,
    x = NULL,
    theta_start = theta_start,
    control = list(maxit = 100)
  )

  # Check that optimization ran
  expect_type(result, "list")
  expect_true("par" %in% names(result))
  expect_true("value" %in% names(result))

  # Parameters should be finite
  expect_true(all(is.finite(result$par)))
  expect_true(is.finite(result$value))
})

test_that("hazard() accepts exponential distribution", {
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)

  fit <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(log_lambda = log(0.5)),
    dist = "exponential",
    fit = FALSE
  )

  expect_s3_class(fit, "hazard")
  expect_equal(fit$spec$dist, "exponential")
  expect_equal(fit$fit$theta, c(log_lambda = log(0.5)))
})

test_that("hazard() with fit=TRUE estimates exponential parameters", {
  time <- c(1, 2, 3, 4, 5, 6)
  status <- c(1, 0, 1, 1, 0, 1)

  fit <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(log_lambda = log(0.3)),
    dist = "exponential",
    fit = TRUE,
    control = list(maxit = 200)
  )

  # Should produce fitted model
  expect_true(!is.null(fit$fit$theta))
  expect_true(is.finite(fit$fit$objective))
})

test_that("predict() survival works for exponential", {
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)

  fit <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(log_lambda = log(0.5)),
    dist = "exponential",
    fit = FALSE
  )

  # Predict survival
  surv <- predict(fit, type = "survival")

  expect_true(is.numeric(surv))
  expect_equal(length(surv), length(time))
  expect_true(all(surv > 0 & surv <= 1))
})

test_that("predict() cumulative_hazard works for exponential", {
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)

  fit <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(log_lambda = log(0.5)),
    dist = "exponential",
    fit = FALSE
  )

  # Predict cumulative hazard
  cumhaz <- predict(fit, type = "cumulative_hazard")

  expect_true(is.numeric(cumhaz))
  expect_equal(length(cumhaz), length(time))
  expect_true(all(cumhaz > 0))
})

test_that("exponential: S(t) = exp(-H(t)) relationship holds", {
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)

  fit <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(log_lambda = log(0.5)),
    dist = "exponential",
    fit = FALSE
  )

  surv <- predict(fit, type = "survival")
  cumhaz <- predict(fit, type = "cumulative_hazard")

  # S(t) = exp(-H(t))
  expected_cumhaz <- -log(surv)

  expect_equal(cumhaz, expected_cumhaz, tolerance = 1e-12)
})

test_that("predict() works with covariates for exponential", {
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)
  x <- matrix(c(0.5, 1.2, -0.3, 0.8, 1.5), ncol = 1)

  fit <- hazard(
    time = time,
    status = status,
    x = x,
    theta = c(log_lambda = log(0.5), beta = 0.3),
    dist = "exponential",
    fit = FALSE
  )

  # Predict survival at original times and data
  newdata <- data.frame(time = time, X1 = x[, 1])
  surv <- predict(fit, newdata = newdata, type = "survival")

  expect_true(is.numeric(surv))
  expect_equal(length(surv), length(time))
  expect_true(all(surv > 0 & surv <= 1))
})

test_that("exponential: monotonicity of survival", {
  fit <- hazard(
    time = c(1, 2),
    status = c(1, 0),
    x = NULL,
    theta = c(log_lambda = log(0.5)),
    dist = "exponential",
    fit = FALSE
  )

  # Predict across increasing times
  new_times <- seq(0.1, 5, by = 0.5)
  newdata <- data.frame(time = new_times)
  surv <- predict(fit, newdata = newdata, type = "survival")

  # Survival should be monotonically decreasing
  expect_true(all(diff(surv) < 0))
})

test_that("exponential: monotonicity of cumulative hazard", {
  fit <- hazard(
    time = c(1, 2),
    status = c(1, 0),
    x = NULL,
    theta = c(log_lambda = log(0.5)),
    dist = "exponential",
    fit = FALSE
  )

  # Predict across increasing times
  new_times <- seq(0.1, 5, by = 0.5)
  newdata <- data.frame(time = new_times)
  cumhaz <- predict(fit, newdata = newdata, type = "cumulative_hazard")

  # Cumulative hazard should be monotonically increasing
  expect_true(all(diff(cumhaz) > 0))
})

test_that("exponential baseline hazard is constant", {
  # For exponential with no covariates, hazard should be constant
  # h(t) = lambda (does not depend on t)

  log_lambda <- log(0.5)
  lambda <- exp(log_lambda)

  time <- c(1, 2, 5, 10)
  fit <- hazard(
    time = time,
    status = c(1, 0, 1, 0),
    x = NULL,
    theta = c(log_lambda = log_lambda),
    dist = "exponential",
    fit = FALSE
  )

  # For exponential, cumulative hazard = lambda * t
  # So H(t) should be proportional to t
  cumhaz <- predict(fit, type = "cumulative_hazard")

  # Check proportionality: H(t) / t should be constant
  ratio <- cumhaz / time
  expect_true(all(abs(ratio - lambda) < 1e-12))
})
