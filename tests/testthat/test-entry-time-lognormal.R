## test-entry-time-lognormal.R -- counting-process entry (left truncation)
## for dist = "lognormal" (issue #253).
##
## For status 0/1 rows `time_lower` is the entry time. A row that enters at
## s > 0 contributes log f(t) - log S(s) (event) or log S(t) - log S(s)
## (right-censored). The oracle below is written from stats::dlnorm/plnorm,
## independently of the package's z / Mills-ratio algebra.

# Simulated left-truncated data: one covariate, a subset of rows entering
# late, the rest entering at 0, and independent right censoring.
.hzr_ln_entry_data <- function(n = 200, seed = 253) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n), ncol = 1L, dimnames = list(NULL, "age"))
  eta <- 0.4 + 0.5 * x[, 1]
  t_event <- exp(stats::rnorm(n, eta, 0.8))
  t_cens <- stats::rexp(n, 0.15)
  time <- pmin(t_event, t_cens)
  status <- as.integer(t_event <= t_cens)
  entry <- rep(0, n)
  late <- stats::runif(n) < 0.6
  entry[late] <- time[late] * stats::runif(sum(late), 0.1, 0.9)
  list(time = time, status = status, entry = entry, x = x, late = late)
}

# Hand-written truncated log-likelihood in the package's parameterisation:
# theta = (mu, log sigma, beta); log T ~ N(mu + x beta, sigma^2).
.hzr_ln_entry_oracle <- function(theta, time, status, entry, x) {
  eta <- theta[1] + as.numeric(x %*% theta[-(1:2)])
  sdlog <- exp(theta[2])
  ll <- ifelse(status == 1,
               stats::dlnorm(time, eta, sdlog, log = TRUE),
               stats::plnorm(time, eta, sdlog, lower.tail = FALSE,
                             log.p = TRUE))
  in_late <- entry > 0
  ll[in_late] <- ll[in_late] -
    stats::plnorm(entry[in_late], eta[in_late], sdlog, lower.tail = FALSE,
                  log.p = TRUE)
  sum(ll)
}

theta_test <- c(mu = 0.3, log_sigma = log(0.9), age = 0.4)  # not the MLE

test_that("(a) lognormal logl includes the entry term and matches the oracle", {
  d <- .hzr_ln_entry_data()
  expect_true(sum(d$late) > 50)          # the truncation is not vacuous
  expect_true(sum(!d$late) > 50)         # and some rows enter at 0
  ll <- .hzr_logl_lognormal(theta_test, d$time, d$status,
                            time_lower = d$entry, x = d$x)
  ll_oracle <- .hzr_ln_entry_oracle(theta_test, d$time, d$status, d$entry, d$x)
  expect_equal(as.numeric(ll), ll_oracle, tolerance = 1e-10)

  ll_null <- .hzr_logl_lognormal(theta_test, d$time, d$status, x = d$x)
  expect_gt(abs(as.numeric(ll) - as.numeric(ll_null)), 1)
})

test_that("(b) lognormal analytic gradient matches numDeriv with entry times", {
  skip_if_not_installed("numDeriv")
  d <- .hzr_ln_entry_data()
  grad_an <- .hzr_gradient_lognormal(theta_test, d$time, d$status,
                                     time_lower = d$entry, x = d$x)
  obj <- function(th) {
    .hzr_ln_entry_oracle(th, d$time, d$status, d$entry, d$x)
  }
  grad_nd <- numDeriv::grad(obj, theta_test)
  expect_equal(as.numeric(grad_an), grad_nd, tolerance = 1e-6)

  # The return_gradient path attaches the same vector.
  ll <- .hzr_logl_lognormal(theta_test, d$time, d$status,
                            time_lower = d$entry, x = d$x,
                            return_gradient = TRUE)
  expect_equal(as.numeric(attr(ll, "gradient")), grad_nd, tolerance = 1e-6)
})

test_that("(c) lognormal analytic Hessian matches numDeriv with entry times", {
  skip_if_not_installed("numDeriv")
  d <- .hzr_ln_entry_data()
  w <- seq(0.5, 2, length.out = length(d$time))
  obj <- function(th) {
    -.hzr_logl_lognormal(th, d$time, d$status, time_lower = d$entry,
                         x = d$x, weights = w)
  }
  h_an <- .hzr_hessian_lognormal(theta_test, d$time, d$status,
                                 time_lower = d$entry, x = d$x, weights = w)
  h_nd <- numDeriv::hessian(obj, theta_test)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-5)

  h_null <- .hzr_hessian_lognormal(theta_test, d$time, d$status,
                                   x = d$x, weights = w)
  expect_gt(max(abs(h_an - h_null)), 1e-2)
})

test_that("(d) hazard(dist = 'lognormal') fits the truncated likelihood", {
  d <- .hzr_ln_entry_data()
  start <- c(0, 0, 0)
  # Tight reltol: the default (1e-5) lets BFGS stop about 1e-4 short on this
  # likelihood, which would test the stopping rule rather than the MLE.
  fit <- hazard(
    time = d$time, status = d$status, time_lower = d$entry, x = d$x,
    dist = "lognormal", theta = start, fit = TRUE,
    control = list(reltol = 1e-12)
  )
  expect_true(fit$fit$converged)

  oracle_fit <- stats::optim(
    start,
    function(th) -.hzr_ln_entry_oracle(th, d$time, d$status, d$entry, d$x),
    method = "BFGS", control = list(reltol = 1e-14, maxit = 1000)
  )
  expect_equal(oracle_fit$convergence, 0L)
  expect_equal(unname(coef(fit)), oracle_fit$par, tolerance = 1e-4)

  fit_null <- hazard(
    time = d$time, status = d$status, x = d$x,
    dist = "lognormal", theta = start, fit = TRUE
  )
  expect_gt(max(abs(coef(fit) - coef(fit_null))), 0.05)
})

test_that("(e) time_lower = 0 rows equal time_lower = NULL, with no NaN", {
  d <- .hzr_ln_entry_data()
  zero <- rep(0, length(d$time))
  ll0 <- .hzr_logl_lognormal(theta_test, d$time, d$status,
                             time_lower = zero, x = d$x,
                             return_gradient = TRUE)
  lln <- .hzr_logl_lognormal(theta_test, d$time, d$status, x = d$x,
                             return_gradient = TRUE)
  expect_true(is.finite(ll0))
  expect_equal(as.numeric(ll0), as.numeric(lln), tolerance = 1e-12)

  g0 <- attr(ll0, "gradient")
  expect_false(anyNA(g0))
  expect_true(all(is.finite(g0)))
  expect_equal(g0, attr(lln, "gradient"), tolerance = 1e-12)

  h0 <- .hzr_hessian_lognormal(theta_test, d$time, d$status,
                               time_lower = zero, x = d$x)
  expect_false(anyNA(h0))
  expect_true(all(is.finite(h0)))
  expect_equal(h0, .hzr_hessian_lognormal(theta_test, d$time, d$status,
                                          x = d$x), tolerance = 1e-12)
})

test_that("(f) status-2 rows still read time_lower as the interval lower bound", {
  set.seed(2532)
  n <- 120
  time <- exp(stats::rnorm(n, 0.5, 0.7))
  status <- rep(c(1L, 0L, 2L), length.out = n)
  lower <- time               # status 0/1: time_lower = time, no entry
  upper <- time
  iv <- status == 2L
  lower[iv] <- time[iv] * 0.7
  upper[iv] <- time[iv] * 1.4
  theta <- c(mu = 0.4, log_sigma = log(0.8))

  ll <- .hzr_logl_lognormal(theta, time, status,
                            time_lower = lower, time_upper = upper)
  sdlog <- exp(theta[2])
  expected <-
    sum(stats::dlnorm(time[status == 1L], theta[1], sdlog, log = TRUE)) +
    sum(stats::plnorm(time[status == 0L], theta[1], sdlog,
                      lower.tail = FALSE, log.p = TRUE)) +
    sum(log(stats::plnorm(upper[iv], theta[1], sdlog) -
              stats::plnorm(lower[iv], theta[1], sdlog)))
  expect_equal(as.numeric(ll), unname(expected), tolerance = 1e-10)
})
