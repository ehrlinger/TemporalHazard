# A candidate holding NA is never tested by the score criterion (#686).
#
# The score is computed over the rows the base model was fitted on, so a
# candidate with a missing value on one of them has no score. That is a
# decline, as it should be. The defect was that it was declined as
# `non_numeric` and a completed run said nothing: a candidate that would have
# entered was passed over in favour of a weaker one, with no warning, while
# `criterion = "wald"` on the same data warned three times.

missing_candidate_data <- function(seed = 7, n = 300) {
  set.seed(seed)
  strong <- stats::rnorm(n)
  data.frame(
    tt = stats::rexp(n, rate = 0.2 * exp(0.8 * strong)) + 0.01,
    ev = stats::rbinom(n, 1, 0.8),
    strong = strong,
    noise = stats::rnorm(n)
  )
}

missing_candidate_run <- function(d) {
  base <- hazard(survival::Surv(tt, ev) ~ 1, data = d, dist = "weibull",
                 theta = c(mu = 0.2, nu = 1), fit = TRUE)
  warned <- character()
  res <- withCallingHandlers(
    hzr_stepwise(base, scope = ~ strong + noise, data = d,
                 direction = "forward", slentry = 0.3),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(res = res, warned = warned,
       entered = res$steps$variable[res$steps$action == "enter"])
}

test_that("a complete candidate enters, with no missing-value warning", {
  # The known positive: without it, the test below could pass because
  # `strong` was never strong enough to enter.
  out <- missing_candidate_run(missing_candidate_data())
  expect_true("strong" %in% out$entered)
  expect_false(any(grepl("missing", out$warned, ignore.case = TRUE)))
  expect_false("missing_values" %in% names(out$res$criteria$uncomputable_reasons))
})

test_that("a candidate holding NA is declined as missing_values, and the run warns", {
  d <- missing_candidate_data()
  d$strong[5L] <- NA
  out <- missing_candidate_run(d)

  expect_false("strong" %in% out$entered)
  reasons <- out$res$criteria$uncomputable_reasons
  expect_true(reasons[["missing_values"]] >= 1L)
  expect_false("non_numeric" %in% names(reasons))

  hit <- grepl("missing value", out$warned, fixed = TRUE)
  expect_equal(sum(hit), 1L)
  expect_match(out$warned[hit], "without testing", fixed = TRUE)
})
