# test-weak-single-parameter.R
# The weak-direction detector also names ONE parameter the data do not
# determine (#415). The pairwise reading cannot: a correlation matrix
# normalises that parameter's variance away.

nm <- c("late.log_mu", "late.log_tau", "late.gamma", "late.alpha")
th <- c(-0.5, 1, 4, 0.8)
gated <- 1e-12  # below .hzr_rcond_tol, as on every fit this reading targets

# gamma's standard error, relative to gamma, is `rel`; everything else is
# determined and uncorrelated with it.
single_vcov <- function(rel, gamma = th[3]) {
  diag(c(0.01, 0.02, (rel * gamma)^2, 0.001))
}

test_that("a parameter the data do not determine is named on its own", {
  w <- .hzr_weak_direction(single_vcov(50), gated, nm, theta = th)
  expect_true(is.list(w))
  expect_true(isTRUE(w$single))
  expect_identical(w$params, "late.gamma")
  expect_equal(w$se_metric, 50, tolerance = 1e-8)
  expect_match(.hzr_weak_direction_message(w), "near-flat direction on its own",
               fixed = TRUE)
})

test_that("without the fitted values the reading does not run", {
  # Every pre-#415 caller passes no theta, so its answer is unchanged.
  expect_null(.hzr_weak_direction(single_vcov(50), gated, nm))
})

test_that("a parameter's own precision is not the criterion; the gate is", {
  # The same dominant direction, relative standard error 0.3: at gamma-hat =
  # 1e7 that is what a flat-to-infinity likelihood reports, so it is named
  # when the Hessian is ill-conditioned ...
  w <- .hzr_weak_direction(single_vcov(0.3, gamma = 1e7), gated, nm,
                           theta = replace(th, 3, 1e7))
  expect_identical(w$params, "late.gamma")
  expect_match(.hzr_weak_direction_message(w), "'late.gamma' = 1e+07",
               fixed = TRUE)
  # ... and not when it is not.
  expect_null(.hzr_weak_direction(single_vcov(0.3), 1e-4, nm, theta = th))
})

test_that("no single parameter dominating means nothing is named", {
  # The flattest direction is shared by log_mu and log_tau (loading 0.707
  # each, correlation 0.5, below the pairwise reading's bar too).
  v <- diag(c(1, 1, (0.1 * th[3])^2, 0.001))
  v[1, 2] <- v[2, 1] <- 0.5
  expect_null(.hzr_weak_direction(v, gated, nm, theta = th))
})

test_that("the rcond gate still decides whether anything is read", {
  expect_null(.hzr_weak_direction(single_vcov(50), 1e-4, nm, theta = th))
})

test_that("a log-scale estimate near 0 is not scaled by its own value", {
  # gamma is the flat one (relative standard error 50). log_tau = 1e-4 with a
  # standard error of 0.5 is determined to within a factor of e, but scaling
  # it by its own value, as #415's candidate did, reads a relative standard
  # error of 5000 and names log_tau instead.
  th0 <- th
  th0[2] <- 1e-4
  v <- diag(c(0.01, 0.25, (50 * th[3])^2, 0.001))
  w <- .hzr_weak_direction(v, gated, nm, theta = th0)
  expect_identical(w$params, "late.gamma")
  value_scaled <- v / outer(abs(th0), abs(th0))
  expect_identical(nm[which.max(abs(eigen(value_scaled)$vectors[, 1]))],
                   "late.log_tau")
})

test_that("a pair the correlation reading finds is still reported as the pair", {
  v <- matrix(c(1, 0.999, 0.999, 1), 2) * 100
  w <- .hzr_weak_direction(v, gated, c("early.m", "early.nu"),
                           theta = c(27, 0.027))
  expect_false(isTRUE(w$single))
  expect_setequal(w$params, c("early.m", "early.nu"))
})
