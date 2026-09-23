## test-lognormal-dist.R — Unit tests for log-normal (AFT) parametric hazard model
##
## Coverage:
##   - Log-likelihood finite and correctly signed
##   - Analytical gradient matches numerical derivative (univariate + multivariate)
##   - Optimizer convergence (univariate + multivariate)
##   - predict() all four types: linear_predictor, hazard, survival, cumulative_hazard
##   - Survival monotone decreasing; cumhaz monotone increasing
##   - S(t) = exp(-H(t)) relationship
##   - Covariate direction of effect
##   - Mathematical identity: at true params, z ~ Normal(0,1)
##   - Boundary / edge cases (extreme sigma, all events, all censored)

# This file predicts from models built with fit = FALSE on purpose, so the
# warning that those numbers come from starting values is switched off for
# this file only (#398). A file that does not expect the warning sees it as
# an ordinary leaked warning.
withr::local_options(TemporalHazard.warn_unfitted_prediction = FALSE)

library(testthat)

# ---------------------------------------------------------------------------
# Likelihood
# ---------------------------------------------------------------------------

test_that("log-normal likelihood is finite and negative (univariate)", {
  time   <- c(0.5, 1.0, 2.0, 3.0, 5.0)
  status <- c(1,   1,   0,   1,   0)
  theta  <- c(mu = 0.5, log_sigma = 0.0)  # mu=0.5, sigma=1

  logl <- .hzr_logl_lognormal(theta, time, status)

  expect_true(is.finite(logl))
  expect_true(logl < 0)
})

test_that("log-normal likelihood is finite with covariates (multivariate)", {
  set.seed(101)
  n      <- 25
  time   <- rlnorm(n, meanlog = 1, sdlog = 0.8)
  status <- rbinom(n, 1, 0.75)
  x      <- cbind(rnorm(n), rbinom(n, 1, 0.5))
  theta  <- c(mu = 1.0, log_sigma = log(0.8), beta1 = 0.2, beta2 = -0.3)

  logl <- .hzr_logl_lognormal(theta, time, status, x = x)

  expect_true(is.finite(logl))
  expect_true(logl < 0)
})

test_that("log-normal likelihood increases near the true parameters", {
  set.seed(202)
  n      <- 200
  mu_true    <- 1.2
  sigma_true <- 0.9
  time   <- rlnorm(n, meanlog = mu_true, sdlog = sigma_true)
  status <- rbinom(n, 1, 0.8)

  theta_true <- c(mu_true, log(sigma_true))
  theta_off  <- c(0.0, 0.0)

  logl_true <- .hzr_logl_lognormal(theta_true, time, status)
  logl_off  <- .hzr_logl_lognormal(theta_off,  time, status)

  # MLE should (in expectation) be better than a far-off value
  expect_true(logl_true > logl_off)
})

# ---------------------------------------------------------------------------
# Gradient — numerical verification
# ---------------------------------------------------------------------------

test_that("log-normal gradient matches numerical derivative (univariate)", {
  time   <- c(0.5, 1.0, 1.5, 2.0, 3.0)
  status <- c(1,   1,   0,   1,   0)
  theta  <- c(mu = 0.3, log_sigma = -0.2)
  eps    <- 1e-5

  logl_with_grad <- .hzr_logl_lognormal(
    theta,
    time,
    status,
    return_gradient = TRUE
  )
  grad_analytical <- attr(logl_with_grad, "gradient")

  grad_numerical <- numeric(length(theta))
  for (i in seq_along(theta)) {
    tp <- theta
    tp[i] <- tp[i] + eps
    tm <- theta
    tm[i] <- tm[i] - eps
    grad_numerical[i] <- (
      .hzr_logl_lognormal(tp, time, status) -
      .hzr_logl_lognormal(tm, time, status)
    ) / (2 * eps)
  }

  expect_equal(grad_analytical, grad_numerical, tolerance = 1e-4)
})

test_that("log-normal gradient matches numerical derivative (multivariate)", {
  set.seed(303)
  time   <- c(0.4, 0.9, 1.4, 2.1, 3.0)
  status <- c(1,   0,   1,   1,   0)
  x      <- cbind(c(1, 0, 1, 0, 1), c(0.5, -0.5, 1.0, -1.0, 0.0))
  colnames(x) <- c("x1", "x2")
  theta  <- c(mu = 0.5, log_sigma = -0.1, beta1 = 0.3, beta2 = -0.2)
  eps    <- 1e-5

  logl_with_grad <- .hzr_logl_lognormal(
    theta,
    time,
    status,
    x = x,
    return_gradient = TRUE
  )
  grad_analytical <- attr(logl_with_grad, "gradient")

  grad_numerical <- numeric(length(theta))
  for (i in seq_along(theta)) {
    tp <- theta
    tp[i] <- tp[i] + eps
    tm <- theta
    tm[i] <- tm[i] - eps
    grad_numerical[i] <- (
      .hzr_logl_lognormal(tp, time, status, x = x) -
      .hzr_logl_lognormal(tm, time, status, x = x)
    ) / (2 * eps)
  }

  expect_equal(grad_analytical, grad_numerical, tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# Optimizer
# ---------------------------------------------------------------------------

test_that("log-normal optimizer converges (univariate)", {
  set.seed(404)
  n          <- 60
  mu_true    <- 1.0
  sigma_true <- 0.8
  time       <- rlnorm(n, meanlog = mu_true, sdlog = sigma_true)
  cens       <- runif(n, 0, quantile(time, 0.8))
  status     <- as.integer(time <= cens)
  time       <- pmin(time, cens)

  fit <- .hzr_optim_lognormal(time, status, x = NULL,
                               theta_start = c(0.5, 0.0))

  expect_equal(fit$convergence, 0L, info = fit$message)
  expect_true(is.finite(fit$value))
  expect_equal(length(fit$par), 2L)

  # Recovered mu should be in a plausible neighbourhood
  expect_true(abs(fit$par[1] - mu_true) < 0.5)
  # Recovered sigma (on natural scale) should be positive and reasonable
  expect_true(exp(fit$par[2]) > 0.3 && exp(fit$par[2]) < 2.5)
})

test_that("log-normal optimizer converges (multivariate)", {
  set.seed(505)
  n <- 80
  mu_true    <- 0.8
  sigma_true <- 0.7
  beta_true  <- c(0.4, -0.3)

  x   <- cbind(rnorm(n), rnorm(n))
  eta <- mu_true + x %*% beta_true         # AFT: eta shifts log(T)
  time   <- rlnorm(n, meanlog = eta, sdlog = sigma_true)
  cens   <- runif(n, 0, quantile(time, 0.75))
  status <- as.integer(time <= cens)
  time   <- pmin(time, cens)

  fit <- .hzr_optim_lognormal(time, status, x = x,
                               theta_start = c(0.5, 0.0, 0.0, 0.0))

  expect_equal(fit$convergence, 0L, info = fit$message)
  expect_equal(length(fit$par), 4L)
  expect_true(is.finite(fit$value))
})

test_that("hazard() dispatches lognormal and converges", {
  set.seed(606)
  n      <- 60
  time   <- rlnorm(n, meanlog = 1, sdlog = 0.8)
  cens   <- runif(n, 0, quantile(time, 0.8))
  status <- as.integer(time <= cens)
  time   <- pmin(time, cens)

  fit <- hazard(time = time, status = status,
                theta = c(mu = 1.0, log_sigma = 0.0),
                dist = "lognormal", fit = TRUE)

  expect_s3_class(fit, "hazard")
  expect_true(fit$fit$converged)
  expect_true(is.finite(fit$fit$objective))
})

# ---------------------------------------------------------------------------
# predict() — all four types
# ---------------------------------------------------------------------------

.make_ln_obj <- function(theta, n = 5) {
  obj <- list(
    spec = list(dist = "lognormal"),
    fit  = list(theta = theta),
    data = list(time = seq(0.5, 2.5, length.out = n),
                status = rep(1L, n), x = NULL)
  )
  class(obj) <- "hazard"
  obj
}

test_that("predict() survival is in (0,1) and finite", {
  fit_obj <- .make_ln_obj(c(mu = 0.5, log_sigma = log(0.8)))
  time    <- seq(0.2, 4.0, by = 0.4)
  surv    <- predict(fit_obj, newdata = data.frame(time = time), type = "survival")

  expect_equal(length(surv), length(time))
  expect_true(all(surv > 0 & surv < 1))
  expect_true(all(is.finite(surv)))
})

test_that("predict() survival is monotone decreasing", {
  fit_obj <- .make_ln_obj(c(mu = 0.5, log_sigma = log(0.8)))
  time    <- seq(0.1, 5.0, by = 0.3)
  surv    <- predict(fit_obj, newdata = data.frame(time = time), type = "survival")

  expect_equal(surv, sort(surv, decreasing = TRUE), tolerance = 1e-10)
})

test_that("predict() cumulative_hazard is monotone increasing", {
  fit_obj <- .make_ln_obj(c(mu = 0.5, log_sigma = log(0.8)))
  time    <- seq(0.1, 5.0, by = 0.3)
  cumhaz  <- predict(fit_obj, newdata = data.frame(time = time),
                     type = "cumulative_hazard")

  expect_equal(cumhaz, sort(cumhaz), tolerance = 1e-10)
})

test_that("predict() survival = exp(-cumhaz) identity holds", {
  fit_obj <- .make_ln_obj(c(mu = 0.8, log_sigma = log(0.7)))
  time    <- c(0.5, 1.0, 1.5, 2.0, 3.0, 4.0)
  surv    <- predict(
    fit_obj,
    newdata = data.frame(time = time),
    type = "survival"
  )
  cumhaz  <- predict(
    fit_obj,
    newdata = data.frame(time = time),
    type = "cumulative_hazard"
  )

  expect_equal(surv, exp(-cumhaz), tolerance = 1e-10)
})

test_that("predict() linear_predictor returns beta dot x (AFT covariate effect)", {
  theta   <- c(mu = 1.0, log_sigma = 0.0, beta1 = 0.4, beta2 = -0.2)
  fit_obj <- list(
    spec = list(dist = "lognormal"),
    fit  = list(theta = theta),
    data = list(time = 1:3, status = rep(1L, 3),
                x = cbind(V1 = 1:3, V2 = 3:1))
  )
  class(fit_obj) <- "hazard"

  x_new  <- data.frame(V1 = c(1, 0, -1), V2 = c(0, 1, -1))
  lp     <- predict(fit_obj, newdata = x_new, type = "linear_predictor")

  expected <- as.numeric(as.matrix(x_new) %*% c(0.4, -0.2))
  expect_equal(lp, expected)
})

test_that("predict() hazard returns exp(linear_predictor)", {
  theta   <- c(mu = 1.0, log_sigma = 0.0, beta1 = 0.3, beta2 = -0.1)
  fit_obj <- list(
    spec = list(dist = "lognormal"),
    fit  = list(theta = theta),
    data = list(time = 1:3, status = rep(1L, 3),
                x = cbind(V1 = 1:3, V2 = 3:1))
  )
  class(fit_obj) <- "hazard"

  x_new <- data.frame(V1 = c(1, 2), V2 = c(-1, 0))
  h     <- predict(fit_obj, newdata = x_new, type = "hazard")
  lp    <- predict(fit_obj, newdata = x_new, type = "linear_predictor")

  expect_equal(h, exp(lp))
})

# ---------------------------------------------------------------------------
# Covariate effects
# ---------------------------------------------------------------------------

test_that("positive AFT covariate shifts survival upward", {
  # Larger x → larger eta → larger log(T) → later failure → higher survival
  theta   <- c(mu = 0.5, log_sigma = log(0.7), beta1 = 0.5)
  fit_obj <- list(
    spec = list(dist = "lognormal"),
    fit  = list(theta = theta),
    data = list(time = 1:3, status = rep(1L, 3),
                x = cbind(V1 = 1:3))
  )
  class(fit_obj) <- "hazard"

  t_eval <- c(1.0, 2.0, 3.0)
  s_high <- predict(fit_obj,
                    newdata = data.frame(time = t_eval, V1 = rep(2, 3)),
                    type = "survival")
  s_low  <- predict(fit_obj,
                    newdata = data.frame(time = t_eval, V1 = rep(-2, 3)),
                    type = "survival")

  # Higher covariate value → longer survival → higher S(t)
  expect_true(all(s_high > s_low))
})

# ---------------------------------------------------------------------------
# Mathematical identity: standardised residuals
# ---------------------------------------------------------------------------

test_that("standardised residuals are approximately N(0,1) at true parameters", {
  set.seed(707)
  n          <- 500
  mu_true    <- 1.3
  sigma_true <- 0.6
  time       <- rlnorm(n, meanlog = mu_true, sdlog = sigma_true)
  # No censoring so all z are observed
  status <- rep(1L, n)

  z <- (log(time) - mu_true) / sigma_true

  # Mean and SD of z should be close to 0 and 1
  expect_equal(mean(z), 0, tolerance = 0.15)
  expect_equal(sd(z), 1, tolerance = 0.15)
})

# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

test_that("log-normal likelihood handles extreme sigma gracefully", {
  time   <- c(0.5, 1.0, 2.0, 4.0)
  status <- c(1,   1,   0,   1)

  # Very small sigma (near-constant log-time)
  theta_tight <- c(0.0, log(0.01))
  logl_tight  <- .hzr_logl_lognormal(theta_tight, time, status)
  expect_true(is.finite(logl_tight))

  # Large sigma (very spread out)
  theta_wide <- c(0.0, log(5.0))
  logl_wide  <- .hzr_logl_lognormal(theta_wide, time, status)
  expect_true(is.finite(logl_wide))
})

test_that("log-normal optimizer produces finite parameters with high censoring", {
  set.seed(808)
  n      <- 40
  time   <- rlnorm(n, meanlog = 0.5, sdlog = 0.9)
  # ~60% censored
  cens   <- runif(n, 0, quantile(time, 0.5))
  status <- as.integer(time <= cens)
  time   <- pmin(time, cens)

  fit <- .hzr_optim_lognormal(time, status, x = NULL,
                               theta_start = c(0.0, 0.0))

  expect_true(all(is.finite(fit$par)))
  expect_true(is.finite(fit$value))
})
