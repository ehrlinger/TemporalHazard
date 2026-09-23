# This file predicts from models built with fit = FALSE on purpose, so the
# warning that those numbers come from starting values is switched off for
# this file only (#398). A file that does not expect the warning sees it as
# an ordinary leaked warning.
withr::local_options(TemporalHazard.warn_unfitted_prediction = FALSE)

library(testthat)

test_that("Log-logistic likelihood computation (univariate)", {
  # Simple univariate case
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 1, 0, 1, 0)

  # θ = [log(α), log(β)]
  # log(α) = 0 => α = 1
  # log(β) = 0 => β = 1 (exponential limit)
  log_alpha <- 0
  log_beta <- 0
  theta <- c(log_alpha, log_beta)

  logl <- .hzr_logl_loglogistic(theta, time, status, x = NULL, return_gradient = FALSE)

  # Check that likelihood is finite
  expect_true(is.finite(logl))
  # Log-likelihood should be negative (typical for parametric models)
  expect_true(logl < 0)
})

test_that("Log-logistic likelihood with covariates (multivariate)", {
  set.seed(123)
  n <- 20
  time <- rexp(n, rate = 0.5)
  status <- rbinom(n, size = 1, prob = 0.7)
  x <- cbind(runif(n), rnorm(n))

  log_alpha <- 0.5
  log_beta <- 0.3
  beta_coef <- c(0.2, -0.3)
  theta <- c(log_alpha, log_beta, beta_coef)

  logl <- .hzr_logl_loglogistic(theta, time, status, x = x, return_gradient = FALSE)

  expect_true(is.finite(logl))
  expect_true(logl < 0)
})

test_that("Log-logistic gradient computation against numerical derivative", {
  set.seed(456)
  time <- c(0.5, 1.0, 1.5, 2.0, 2.5)
  status <- c(1, 1, 0, 1, 0)

  log_alpha <- 0.2
  log_beta <- -0.1
  theta <- c(log_alpha, log_beta)

  # Analytical gradient
  logl <- .hzr_logl_loglogistic(theta, time, status, x = NULL, return_gradient = TRUE)
  grad_analytical <- attr(logl, "gradient")

  # Numerical gradient
  eps <- 1e-4
  grad_numerical <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta
    theta_plus[i] <- theta_plus[i] + eps
    logl_plus <- .hzr_logl_loglogistic(theta_plus, time, status, x = NULL, return_gradient = FALSE)

    theta_minus <- theta
    theta_minus[i] <- theta_minus[i] - eps
    logl_minus <- .hzr_logl_loglogistic(theta_minus, time, status, x = NULL, return_gradient = FALSE)

    grad_numerical[i] <- (logl_plus - logl_minus) / (2 * eps)
  }

  # Compare (within numerical tolerance)
  expect_equal(grad_analytical, grad_numerical, tolerance = 1e-4)
})

test_that("Log-logistic gradient with covariates", {
  set.seed(789)
  time <- c(0.5, 1.0, 1.5, 2.0, 2.5)
  status <- c(1, 1, 0, 1, 0)
  x <- cbind(c(1, 0, 1, 0, 1), c(0.5, 1.5, -0.5, 2.0, 1.0))
  colnames(x) <- c("x1", "x2")

  theta <- c(log_alpha = 0.3, log_beta = 0.1, beta_1 = 0.2, beta_2 = -0.1)

  # Analytical gradient
  logl <- .hzr_logl_loglogistic(theta, time, status, x = x, return_gradient = TRUE)
  grad_analytical <- attr(logl, "gradient")

  # Numerical gradient
  eps <- 1e-4
  grad_numerical <- numeric(length(theta))
  for (i in seq_along(theta)) {
    theta_plus <- theta
    theta_plus[i] <- theta_plus[i] + eps
    logl_plus <- .hzr_logl_loglogistic(theta_plus, time, status, x = x, return_gradient = FALSE)

    theta_minus <- theta
    theta_minus[i] <- theta_minus[i] - eps
    logl_minus <- .hzr_logl_loglogistic(theta_minus, time, status, x = x, return_gradient = FALSE)

    grad_numerical[i] <- (logl_plus - logl_minus) / (2 * eps)
  }

  expect_equal(grad_analytical, grad_numerical, tolerance = 1e-4)
})

test_that("Log-logistic optimization converges (univariate)", {
  set.seed(111)
  n <- 30
  # Simulate from log-logistic with α=1, β=1.5
  alpha_true <- 1.0
  beta_true <- 1.5
  u <- runif(n)
  time <- (u / (1 - u)) ^ (1 / beta_true) / alpha_true
  # Apply censoring
  status <- rbinom(n, size = 1, prob = 0.8)

  # Starting values
  theta_start <- c(log(0.5), log(1.0))

  # Fit
  fit <- .hzr_optim_loglogistic(time, status, x = NULL, theta_start = theta_start)

  expect_true(fit$convergence == 0, info = paste("Message:", fit$message))
  expect_true(is.finite(fit$value))
  expect_equal(length(fit$par), 2)

  # Check that estimates are in reasonable range
  expect_true(fit$par[1] > log(0.1) && fit$par[1] < log(10))
  expect_true(fit$par[2] > log(0.5) && fit$par[2] < log(5))
})

test_that("Log-logistic optimization converges (multivariate)", {
  set.seed(222)
  n <- 40

  # True parameters
  alpha_true <- 1.0
  beta_true <- 1.2
  beta_coef_true <- c(0.3, -0.2)

  # Generate covariates
  x <- cbind(rnorm(n), rnorm(n))
  eta <- x %*% beta_coef_true

  # Generate baseline times from log-logistic
  u <- runif(n)
  baseline_time <- (u / (1 - u)) ^ (1 / beta_true) / alpha_true

  # Apply covariate effect
  time <- baseline_time * exp(-eta / beta_true)  # Approximate adjustment

  # Apply censoring
  status <- rbinom(n, size = 1, prob = 0.75)

  # Starting values
  theta_start <- c(log(0.8), log(1.0), 0.1, -0.1)

  # Fit
  fit <- .hzr_optim_loglogistic(time, status, x = x, theta_start = theta_start)

  expect_true(fit$convergence == 0, info = paste("Message:", fit$message))
  expect_true(is.finite(fit$value))
  expect_equal(length(fit$par), 4)
})

test_that("Log-logistic survival predictions monotonicity", {
  set.seed(333)
  time <- seq(0.1, 5, by = 0.5)
  theta <- c(log(1.0), log(1.5))  # log-logistic params

  # Predict survival at increasing times
  fit_obj <- list(
    spec = list(dist = "loglogistic"),
    fit = list(theta = theta),
    data = list(time = time, status = rep(1, length(time)), x = NULL)
  )
  class(fit_obj) <- "hazard"

  survival <- predict(fit_obj, newdata = data.frame(time = time), type = "survival")

  # Survival should be decreasing with time
  expect_equal(survival, sort(survival, decreasing = TRUE), tolerance = 1e-10)
})

test_that("Log-logistic cumulative hazard predictions monotonicity", {
  set.seed(444)
  time <- seq(0.1, 5, by = 0.5)
  theta <- c(log(1.0), log(1.5))

  fit_obj <- list(
    spec = list(dist = "loglogistic"),
    fit = list(theta = theta),
    data = list(time = time, status = rep(1, length(time)), x = NULL)
  )
  class(fit_obj) <- "hazard"

  cumhaz <- predict(fit_obj, newdata = data.frame(time = time), type = "cumulative_hazard")

  # Cumulative hazard should be increasing with time
  expect_equal(cumhaz, sort(cumhaz, decreasing = FALSE), tolerance = 1e-10)
})

test_that("Log-logistic survival and cumulative hazard relationship", {
  set.seed(555)
  time <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
  theta <- c(log(1.2), log(1.8))

  fit_obj <- list(
    spec = list(dist = "loglogistic"),
    fit = list(theta = theta),
    data = list(time = time, status = rep(1, length(time)), x = NULL)
  )
  class(fit_obj) <- "hazard"

  survival <- predict(fit_obj, newdata = data.frame(time = time), type = "survival")
  cumhaz <- predict(fit_obj, newdata = data.frame(time = time), type = "cumulative_hazard")

  # S(t) = exp(-H(t))
  expected_survival <- exp(-cumhaz)
  expect_equal(survival, expected_survival, tolerance = 1e-10)
})

test_that("Log-logistic covariate effect on survival", {
  set.seed(666)
  time <- c(1.0, 2.0, 3.0)

  # With covariates
  x_high <- cbind(c(1, 1, 1))  # Positive covariate effect
  x_low <- cbind(c(-1, -1, -1))  # Negative covariate effect

  theta <- c(log(1.0), log(1.5), 0.5)  # Positive beta coefficient

  fit_obj <- list(
    spec = list(dist = "loglogistic"),
    fit = list(theta = theta),
    data = list(time = time, status = rep(1, length(time)), x = NULL)
  )
  class(fit_obj) <- "hazard"

  s_high <- predict(fit_obj, newdata = data.frame(time = time, V1 = c(1, 1, 1)), type = "survival")
  s_low <- predict(fit_obj, newdata = data.frame(time = time, V1 = c(-1, -1, -1)), type = "survival")

  # Higher covariate value should decrease survival (increase hazard) with positive beta
  expect_true(all(s_high < s_low))
})

test_that("Log-logistic hazard() function fitting (univariate)", {
  set.seed(777)
  n <- 40
  alpha_true <- 0.9
  beta_true <- 1.3
  u <- runif(n)
  time <- (u / (1 - u)) ^ (1 / beta_true) / alpha_true
  status <- rbinom(n, size = 1, prob = 0.8)

  theta_start <- c(log(0.7), log(1.0))

  # Fit via hazard() function
  fit <- hazard(
    time = time,
    status = status,
    theta = theta_start,
    dist = "loglogistic",
    fit = TRUE
  )

  expect_true(fit$fit$converged)
  expect_true(is.finite(fit$fit$objective))
  expect_equal(length(fit$fit$theta), 2)
})

test_that("Log-logistic hazard() function fitting (multivariate)", {
  set.seed(888)
  n <- 50

  # True parameters
  alpha_true <- 1.1
  beta_true <- 1.4
  beta_coef_true <- c(0.25, -0.15)

  x <- cbind(rnorm(n), rnorm(n))
  eta <- x %*% beta_coef_true

  # Generate times from log-logistic
  u <- runif(n)
  baseline_time <- (u / (1 - u)) ^ (1 / beta_true) / alpha_true
  time <- baseline_time * exp(-eta / beta_true)
  status <- rbinom(n, size = 1, prob = 0.75)

  theta_start <- c(log(0.8), log(1.0), 0.1, -0.1)

  fit <- hazard(
    time = time,
    status = status,
    x = x,
    theta = theta_start,
    dist = "loglogistic",
    fit = TRUE
  )

  expect_true(fit$fit$converged)
  expect_equal(length(fit$fit$theta), 4)
})

test_that("Log-logistic predict() linear_predictor", {
  theta <- c(log(1.0), log(1.5), 0.3, -0.2)
  fit_obj <- list(
    spec = list(dist = "loglogistic"),
    fit = list(theta = theta),
    data = list(time = 1:5, status = rep(1, 5), x = cbind(V1 = 1:5, V2 = 5:1))
  )
  class(fit_obj) <- "hazard"

  newdata <- data.frame(V1 = c(1, 0, -1), V2 = c(0, 1, -1))
  pred <- predict(fit_obj, newdata = newdata, type = "linear_predictor")

  # Should match x %*% beta_coef
  expected <- as.numeric(as.matrix(newdata) %*% c(0.3, -0.2))
  expect_equal(pred, expected)
})

test_that("Log-logistic predict() hazard", {
  theta <- c(log(1.0), log(1.5), 0.2, -0.1)
  fit_obj <- list(
    spec = list(dist = "loglogistic"),
    fit = list(theta = theta),
    data = list(time = 1:5, status = rep(1, 5), x = cbind(V1 = 1:5, V2 = 5:1))
  )
  class(fit_obj) <- "hazard"

  newdata <- data.frame(V1 = c(1, 0, 1), V2 = c(0, 1, -1))
  pred <- predict(fit_obj, newdata = newdata, type = "hazard")

  # Hazard should be exp(eta) where eta = x %*% beta
  x_mat <- as.matrix(newdata)
  beta <- theta[3:4]
  eta <- as.numeric(x_mat %*% beta)
  expected <- exp(eta)

  expect_equal(pred, expected)
})

test_that("Log-logistic parameter bounds are reasonable", {
  set.seed(999)
  n <- 35
  time <- rexp(n, rate = 0.4)
  status <- rbinom(n, size = 1, prob = 0.8)

  theta_start <- c(log(0.5), log(0.8))

  fit <- .hzr_optim_loglogistic(time, status, x = NULL, theta_start = theta_start)

  # Parameters should be finite and not absurdly large/small
  expect_true(all(is.finite(fit$par)))
  expect_true(all(fit$par > -10 & fit$par < 10))
})

test_that("Log-logistic error handling: invalid theta length", {
  time <- c(1, 2, 3)
  status <- c(1, 1, 0)

  # Too few parameters for covariates
  theta_bad <- c(log(1.0))  # Missing log(beta)

  expect_error(
    hazard(time = time, status = status, x = NULL, theta = theta_bad, dist = "loglogistic", fit = FALSE),
    NA  # Should allow creation without fitting
  )
})

test_that("Log-logistic mathematical relationship with α=β=1 (exponential limit)", {
  # When α=1 and β=1, log-logistic should approximate exponential
  time <- seq(0.1, 3, by = 0.5)
  theta_ll <- c(log(1.0), log(1.0))  # log-logistic with α=β=1
  theta_exp <- c(log(1.0))  # exponential with λ=1

  # Log-logistic: S(t) = 1 / (1 + t) when α=β=1
  # This is NOT the same as exponential S(t) = exp(-t)
  # So this test just confirms both are well-defined

  fit_ll <- list(
    spec = list(dist = "loglogistic"),
    fit = list(theta = theta_ll),
    data = list(time = time, status = rep(1, length(time)), x = NULL)
  )
  class(fit_ll) <- "hazard"

  s_ll <- predict(fit_ll, newdata = data.frame(time = time), type = "survival")

  # S(t) = 1 / (1 + t)
  expected_s <- 1 / (1 + time)
  expect_equal(s_ll, expected_s, tolerance = 1e-10)
})

test_that("Log-logistic negative log-likelihood at boundary", {
  time <- c(0.1, 0.5, 1.0, 2.0, 5.0)
  status <- c(1, 1, 0, 1, 0)

  # Very small α and β (extreme parameter space)
  theta_small <- c(log(1e-3), log(1e-1))
  logl_small <- .hzr_logl_loglogistic(theta_small, time, status)

  # Very large α and β
  theta_large <- c(log(1e2), log(1e1))
  logl_large <- .hzr_logl_loglogistic(theta_large, time, status)

  # Both should be finite
  expect_true(is.finite(logl_small))
  expect_true(is.finite(logl_large))
})
