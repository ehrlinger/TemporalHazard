# `time_lower` as a counting-process entry time, dist = "exponential" --------
#
# For status 0/1 rows `time_lower` is the left-truncation (entry) time: the
# row contributes H(time) - H(time_lower), not H(time). "weibull" and
# "multiphase" honoured it; the exponential likelihood used `time_lower` only
# as the status-2 interval bound, so a left-truncated exponential fit silently
# ignored the truncation. Issue #253.
#
# Every assertion below is against an oracle that does not share the
# implementation's code path: a hand-written sum, the Weibull likelihood at
# nu = 1, numDeriv, or a separate stats::optim fit.

.exp_entry_data <- function() {
  set.seed(253)
  n <- 200
  z <- rnorm(n)
  entry <- runif(n, 0, 2)
  entry[1:40] <- 0                          # some rows enter at the origin
  rate <- 0.3 * exp(0.5 * z)
  t_event <- entry + stats::rexp(n, rate)   # memoryless: truncated draw
  t_cens <- entry + runif(n, 0, 6)
  list(
    time = pmin(t_event, t_cens),
    status = as.integer(t_event <= t_cens),
    entry = entry,
    z = z,
    x = matrix(z, ncol = 1, dimnames = list(NULL, "z")),
    theta = c(log(0.25), 0.4)               # deliberately not the MLE
  )
}

# Hand-written left-truncated exponential log-likelihood.
.exp_entry_oracle <- function(theta, d) {
  lambda <- exp(theta[1])
  eta <- theta[2] * d$z
  sum(d$status * (log(lambda) + eta) - lambda * exp(eta) * (d$time - d$entry))
}

test_that("exponential logl honours time_lower as entry time (hand oracle)", {
  d <- .exp_entry_data()
  expect_gt(sum(d$entry > 0), 100)
  expect_true(any(d$status == 0) && any(d$status == 1))

  ll <- .hzr_logl_exponential(d$theta, d$time, d$status,
                              time_lower = d$entry, x = d$x)
  ll_null <- .hzr_logl_exponential(d$theta, d$time, d$status, x = d$x)

  expect_equal(ll, .exp_entry_oracle(d$theta, d), tolerance = 1e-10)
  # The entry term is not negligible: it must move the likelihood.
  expect_gt(abs(ll - ll_null), 1)
})

test_that("exponential logl equals the Weibull logl at nu = 1 with entry", {
  d <- .exp_entry_data()
  ll_exp <- .hzr_logl_exponential(d$theta, d$time, d$status,
                                  time_lower = d$entry, x = d$x)
  # Weibull theta is on the natural scale: c(mu, nu, beta).
  ll_wei <- .hzr_logl_weibull(c(exp(d$theta[1]), 1, d$theta[2]),
                              d$time, d$status,
                              time_lower = d$entry, x = d$x)
  expect_equal(ll_exp, ll_wei, tolerance = 1e-10)
})

test_that("exponential gradient matches numDeriv with entry times", {
  skip_if_not_installed("numDeriv")
  d <- .exp_entry_data()
  obj <- function(p) {
    .hzr_logl_exponential(p, d$time, d$status, time_lower = d$entry, x = d$x)
  }
  g_nd <- numDeriv::grad(obj, d$theta)
  g_an <- .hzr_gradient_exponential(d$theta, d$time, d$status,
                                    time_lower = d$entry, x = d$x)
  expect_equal(g_an, g_nd, tolerance = 1e-6)

  # The return_gradient path passes precomputed pieces; it must agree too.
  g_attr <- attr(.hzr_logl_exponential(d$theta, d$time, d$status,
                                       time_lower = d$entry, x = d$x,
                                       return_gradient = TRUE), "gradient")
  expect_equal(g_attr, g_nd, tolerance = 1e-6)
})

test_that("exponential analytic Hessian matches numDeriv with entry times", {
  skip_if_not_installed("numDeriv")
  d <- .exp_entry_data()
  # .hzr_hessian_exponential() returns the Hessian of the objective (-logl).
  obj <- function(p) {
    -.hzr_logl_exponential(p, d$time, d$status, time_lower = d$entry, x = d$x)
  }
  h_nd <- numDeriv::hessian(obj, d$theta)
  h_an <- .hzr_hessian_exponential(d$theta, d$time, d$status,
                                   time_lower = d$entry, x = d$x)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-5)
})

test_that("hazard(dist = 'exponential') fit recovers the truncated MLE", {
  d <- .exp_entry_data()
  ref <- stats::optim(c(0, 0), function(p) -.exp_entry_oracle(p, d),
                      method = "BFGS",
                      control = list(reltol = 1e-14, maxit = 1000))
  expect_equal(ref$convergence, 0L)

  fit <- hazard(time = d$time, status = d$status, time_lower = d$entry,
                x = d$x, dist = "exponential", theta = c(log(0.3), 0),
                fit = TRUE)
  fit_null <- hazard(time = d$time, status = d$status,
                     x = d$x, dist = "exponential", theta = c(log(0.3), 0),
                     fit = TRUE)

  expect_equal(unname(coef(fit)), ref$par, tolerance = 1e-4)
  # Ignoring the entry times gives a materially different rate.
  expect_gt(abs(coef(fit)[[1]] - coef(fit_null)[[1]]), 0.05)
})

test_that("time_lower = 0 on status 0/1 rows equals time_lower = NULL", {
  d <- .exp_entry_data()
  ll_zero <- .hzr_logl_exponential(d$theta, d$time, d$status,
                                   time_lower = rep(0, length(d$time)),
                                   x = d$x)
  ll_null <- .hzr_logl_exponential(d$theta, d$time, d$status, x = d$x)
  expect_equal(ll_zero, ll_null, tolerance = 1e-12)
})

test_that("status-2 rows still use time_lower as the interval lower bound", {
  set.seed(2532)
  n <- 60
  z <- rnorm(n)
  x <- matrix(z, ncol = 1, dimnames = list(NULL, "z"))
  tl <- runif(n, 0.1, 2)
  tu <- tl + runif(n, 0.2, 3)
  tt <- (tl + tu) / 2
  th <- c(log(0.25), 0.4)

  # Status-2-only data: the entry logic cannot apply. Value recorded from
  # the implementation before the #253 change.
  ll_int <- .hzr_logl_exponential(th, tt, rep(2L, n),
                                  time_lower = tl, time_upper = tu, x = x)
  expect_equal(ll_int, -91.385091652688772, tolerance = 1e-12)

  # Mixed data: status 0/1 rows take time_lower as entry, status 2 rows as
  # the interval bound. Each block is checked against its own hand oracle.
  d <- .exp_entry_data()
  lambda <- exp(th[1])
  hz <- function(t, zz) lambda * t * exp(th[2] * zz)
  ll_mixed <- .hzr_logl_exponential(
    th,
    time = c(d$time, tt),
    status = c(d$status, rep(2L, n)),
    time_lower = c(d$entry, tl),
    time_upper = c(d$time, tu),
    x = rbind(d$x, x)
  )
  oracle_entry <- sum(d$status * (log(lambda) + th[2] * d$z) -
                        hz(d$time, d$z) + hz(d$entry, d$z))
  oracle_int <- sum(log(exp(-hz(tl, z)) - exp(-hz(tu, z))))
  expect_equal(ll_mixed, oracle_entry + oracle_int, tolerance = 1e-10)
})
