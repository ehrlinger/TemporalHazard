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
                 direction = "forward", slentry = 0.3, trace = FALSE),
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

test_that("hzr_bootstrap() counts replicates that declined a candidate holding NA", {
  # Replicates run with their warnings muffled, so the stepwise warning above
  # cannot reach the user from inside one; the candidate was pooled as not
  # selected, and its frequency then tracked how often the NA rows were drawn.
  # Thirty NA rows: a resample missing all of them is ~exp(-30), so every
  # replicate declines `strong` and the count is exact.
  d <- missing_candidate_data()
  d$strong[1:30] <- NA
  base <- hazard(survival::Surv(tt, ev) ~ 1, data = d, dist = "weibull",
                 theta = c(mu = 0.2, nu = 1), fit = TRUE)
  warned <- character()
  boot <- withCallingHandlers(
    hzr_bootstrap(base, n_boot = 2, seed = 1, scope = ~ strong + noise,
                  direction = "forward", slentry = 0.3),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(boot$n_success, 2L)
  hit <- grepl("successful replicates completed after declining a candidate",
               warned, fixed = TRUE)
  expect_equal(sum(hit), 1L)
  expect_match(warned[hit], "^2 of 2 successful replicates")
  expect_match(warned[hit], "missing value", fixed = TRUE)
})
