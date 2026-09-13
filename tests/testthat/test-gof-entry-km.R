# #254 follow-up: hzr_gof() subtracts H(entry) from each subject's expected
# events, so its Kaplan-Meier and risk-set columns must use the same risk set.
# A subject who enters late is not at risk before entry.

.gof_truncated_fit <- function(entry_share = 0.5) {
  set.seed(7)
  n <- 200
  t <- stats::rexp(n, 0.2) + 0.1
  e <- ifelse(stats::runif(n) < entry_share, t * stats::runif(n, 0.2, 0.9), 0)
  s <- stats::rbinom(n, 1, 0.7)
  f <- suppressWarnings(hazard(time = t, status = s, time_lower = e,
                               dist = "weibull", theta = c(0.2, 1),
                               fit = TRUE))
  list(fit = f, t = t, e = e, s = s)
}

test_that("hzr_gof's risk set and Kaplan-Meier honour entry times", {
  k <- .gof_truncated_fit()
  expect_gt(sum(k$e > 0), 50)
  g <- hzr_gof(k$fit)
  km <- survival::survfit(survival::Surv(k$e, k$t, k$s) ~ 1)
  # At the first exit time only subjects who have entered are at risk.
  expect_equal(g$n_risk[1], km$n.risk[1])
  expect_lt(g$n_risk[1], 200)
  idx <- match(km$time, g$time)
  expect_false(anyNA(idx))
  expect_equal(g$km_surv[idx], km$surv, tolerance = 1e-12)
  expect_equal(g$n_risk[idx], km$n.risk)
})

test_that("hzr_gof's Kaplan-Meier is unchanged when no row has an entry time", {
  k <- .gof_truncated_fit(entry_share = 0)
  expect_true(all(k$e == 0))
  g <- hzr_gof(k$fit)
  km <- survival::survfit(survival::Surv(k$t, k$s) ~ 1)
  idx <- match(km$time, g$time)
  expect_equal(g$km_surv[idx], km$surv, tolerance = 1e-12)
  expect_equal(g$n_risk[idx], km$n.risk)
})
