# Left truncation (counting-process entry) for dist = "loglogistic" -- issue #253.
#
# For status 0/1 rows `time_lower` is the entry time, so the row contributes
# H(time) - H(entry) rather than H(time).  The oracle below is written
# independently of the package from the documented parameterisation:
#   H(t) = log(1 + alpha t^beta exp(eta)),
#   h(t) = alpha beta t^(beta - 1) exp(eta) / (1 + alpha t^beta exp(eta)),
# with theta = (log alpha, log beta, b).

llogis_trunc_oracle <- function(theta, time, status, entry, x) {
  a <- exp(theta[1])
  b <- exp(theta[2])
  eta <- if (is.null(x)) rep(0, length(time)) else as.numeric(x %*% theta[-(1:2)])
  cumhaz <- function(t) log(1 + a * t^b * exp(eta))
  log_haz <- log(a) + log(b) + (b - 1) * log(time) + eta - cumhaz(time)
  start <- if (is.null(entry)) rep(0, length(time)) else
    ifelse(status %in% c(0, 1) & entry < time, entry, 0)
  sum(status * log_haz - (cumhaz(time) - cumhaz(start)))
}

# Simulate genuinely left-truncated log-logistic data: draw (entry, T) and
# keep only subjects with T > entry; about half enter at time 0.
sim_llogis_trunc <- function(n, seed, alpha = 0.4, beta = 1.6, b = 0.5) {
  set.seed(seed)
  time <- entry <- numeric(0)
  x <- numeric(0)
  while (length(time) < n) {
    xi <- rnorm(1)
    ei <- if (runif(1) < 0.5) 0 else runif(1, 0, 2)
    u <- runif(1)
    ti <- ((1 / u - 1) / (alpha * exp(b * xi)))^(1 / beta)
    if (ti > ei) {
      time <- c(time, ti)
      entry <- c(entry, ei)
      x <- c(x, xi)
    }
  }
  cens <- entry + rexp(n, 0.25)
  status <- as.integer(time <= cens)
  time <- pmin(time, cens)
  list(time = time, status = status, entry = entry, x = cbind(z = x))
}

theta_test <- c(log_alpha = log(0.5), log_beta = log(1.3), z = 0.3)

test_that("(a) loglogistic logl subtracts H(entry) on status 0/1 rows", {
  d <- sim_llogis_trunc(200, seed = 2531)
  # the fixture must actually carry truncation, or the test compares nothing
  expect_gt(sum(d$entry > 0), 50)
  expect_gt(sum(d$status == 0), 10)

  ll_pkg <- .hzr_logl_loglogistic(theta_test, d$time, d$status,
                                  time_lower = d$entry, x = d$x)
  ll_orc <- llogis_trunc_oracle(theta_test, d$time, d$status, d$entry, d$x)
  expect_equal(as.numeric(ll_pkg), ll_orc, tolerance = 1e-10)

  ll_null <- .hzr_logl_loglogistic(theta_test, d$time, d$status, x = d$x)
  expect_equal(as.numeric(ll_null),
               llogis_trunc_oracle(theta_test, d$time, d$status, NULL, d$x),
               tolerance = 1e-10)
  # the entry term is positive, so it must move the value materially
  expect_gt(as.numeric(ll_pkg) - as.numeric(ll_null), 1)
})

test_that("(b) loglogistic analytic gradient matches numDeriv with entry times", {
  skip_if_not_installed("numDeriv")
  d <- sim_llogis_trunc(200, seed = 2532)
  expect_true(any(d$entry == 0) && any(d$entry > 0))
  obj <- function(th) {
    llogis_trunc_oracle(th, d$time, d$status, d$entry, d$x)
  }
  g_nd <- numDeriv::grad(obj, theta_test)

  g_an <- .hzr_gradient_loglogistic(theta_test, d$time, d$status,
                                    time_lower = d$entry, x = d$x)
  expect_true(all(is.finite(g_an)))
  expect_equal(unname(g_an), g_nd, tolerance = 1e-6)

  g_attr <- attr(.hzr_logl_loglogistic(theta_test, d$time, d$status,
                                       time_lower = d$entry, x = d$x,
                                       return_gradient = TRUE), "gradient")
  expect_equal(unname(g_attr), g_nd, tolerance = 1e-6)

  # weighted: every per-row term scales by its weight
  w <- runif(length(d$time), 0.5, 2)
  obj_w <- function(th) {
    .hzr_logl_loglogistic(th, d$time, d$status, time_lower = d$entry,
                          x = d$x, weights = w)
  }
  g_an_w <- .hzr_gradient_loglogistic(theta_test, d$time, d$status,
                                      time_lower = d$entry, x = d$x,
                                      weights = w)
  expect_equal(unname(g_an_w), numDeriv::grad(obj_w, theta_test),
               tolerance = 1e-6)
})

test_that("(c) loglogistic analytic Hessian matches numDeriv with entry times", {
  skip_if_not_installed("numDeriv")
  d <- sim_llogis_trunc(200, seed = 2533)
  w <- runif(length(d$time), 0.5, 2)
  # .hzr_hessian_loglogistic() returns the Hessian of the NEGATIVE logl
  obj <- function(th) {
    -.hzr_logl_loglogistic(th, d$time, d$status, time_lower = d$entry,
                           x = d$x, weights = w)
  }
  h_an <- .hzr_hessian_loglogistic(theta_test, d$time, d$status,
                                   time_lower = d$entry, x = d$x, weights = w)
  h_nd <- numDeriv::hessian(obj, theta_test)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-5)
})

test_that("(d) hazard() loglogistic fit recovers the truncated-likelihood MLE", {
  d <- sim_llogis_trunc(200, seed = 2534)
  start <- c(log(0.3), log(1.2), 0)
  ref <- stats::optim(
    start,
    function(th) -llogis_trunc_oracle(th, d$time, d$status, d$entry, d$x),
    method = "BFGS", control = list(reltol = 1e-14, maxit = 1000)
  )
  expect_equal(ref$convergence, 0L)

  # The default reltol (1e-5) stops BFGS ~5e-10 log-lik short on this flat
  # ridge, about 1e-4 off in theta; tighten it so the comparison is about the
  # likelihood, not the stopping rule.
  ctl <- list(reltol = 1e-12)
  fit <- hazard(time = d$time, status = d$status, time_lower = d$entry,
                x = d$x, dist = "loglogistic", theta = start, fit = TRUE,
                control = ctl)
  expect_equal(unname(coef(fit)), ref$par, tolerance = 1e-4)

  fit0 <- hazard(time = d$time, status = d$status, x = d$x,
                 dist = "loglogistic", theta = start, fit = TRUE,
                 control = ctl)
  # ignoring the truncation must give a different answer
  expect_gt(max(abs(coef(fit) - coef(fit0))), 0.05)
})

test_that("(e) time_lower = 0 on every row is the same as time_lower = NULL", {
  d <- sim_llogis_trunc(150, seed = 2535)
  zero <- rep(0, length(d$time))
  ll_zero <- .hzr_logl_loglogistic(theta_test, d$time, d$status,
                                   time_lower = zero, x = d$x)
  ll_null <- .hzr_logl_loglogistic(theta_test, d$time, d$status, x = d$x)
  expect_identical(as.numeric(ll_zero), as.numeric(ll_null))

  g_zero <- .hzr_gradient_loglogistic(theta_test, d$time, d$status,
                                      time_lower = zero, x = d$x)
  g_null <- .hzr_gradient_loglogistic(theta_test, d$time, d$status, x = d$x)
  expect_true(all(is.finite(g_zero)))
  expect_equal(g_zero, g_null, tolerance = 1e-12)

  h_zero <- .hzr_hessian_loglogistic(theta_test, d$time, d$status,
                                     time_lower = zero, x = d$x)
  h_null <- .hzr_hessian_loglogistic(theta_test, d$time, d$status, x = d$x)
  expect_true(all(is.finite(h_zero)))
  expect_equal(h_zero, h_null, tolerance = 1e-12)
})

test_that("(f) status-2 rows still read time_lower as the interval lower bound", {
  d <- sim_llogis_trunc(150, seed = 2536)
  n <- length(d$time)
  status <- d$status
  lower <- d$entry
  upper <- d$time
  idx2 <- seq(1, n, by = 5)
  status[idx2] <- 2L
  lower[idx2] <- d$time[idx2] * 0.6
  upper[idx2] <- d$time[idx2] * 1.4

  a <- exp(theta_test[1])
  b <- exp(theta_test[2])
  eta <- as.numeric(d$x %*% theta_test[3])
  cdf <- function(t) {
    term <- a * t^b * exp(eta)
    term / (1 + term)
  }
  ll_int <- sum(log(cdf(upper) - cdf(lower))[idx2])
  ll_exact <- llogis_trunc_oracle(theta_test, d$time[-idx2], status[-idx2],
                                  lower[-idx2], d$x[-idx2, , drop = FALSE])

  ll_pkg <- .hzr_logl_loglogistic(theta_test, d$time, status,
                                  time_lower = lower, time_upper = upper,
                                  x = d$x)
  expect_equal(as.numeric(ll_pkg), ll_int + ll_exact, tolerance = 1e-10)
})
