# The Surv() response of a formula is evaluated in `data`, then in the
# formula's environment, as the right-hand side is and as stats::lm() does
# (#447). It used to skip the formula's environment: a name the response
# used resolved to a same-named global, or to the parser's own local
# variables, so a fit inside a function could silently use the wrong value.

set.seed(447)
d447 <- data.frame(tt = rexp(80) * 5, ss = rbinom(80, 1, 0.7),
                   age = rnorm(80, 60, 8))

fit447 <- function(formula) {
  hazard(formula, data = d447, dist = "weibull",
         theta = c(mu = 0.2, nu = 1, 0), fit = FALSE)
}

test_that("a caller-local wins over a same-named global in the response", {
  assign("cutoff_447", 100, envir = globalenv())
  withr::defer(rm("cutoff_447", envir = globalenv()))
  local_fit <- function() {
    cutoff_447 <- 3
    fit447(survival::Surv(pmin(tt, cutoff_447), ss) ~ age)
  }
  f <- local_fit()
  # Times above 3 exist, so a cap at 100 would leave them in place.
  expect_gt(max(d447$tt), 3)
  expect_equal(f$data$time, pmin(d447$tt, 3))
})

test_that("the response sees a caller function and a caller constant", {
  local_fit <- function() {
    helper <- function(x) x * 2
    shift <- 1
    list(
      fn = fit447(survival::Surv(helper(tt), ss) ~ age),
      const = fit447(survival::Surv(tt + shift, ss) ~ age)
    )
  }
  fits <- local_fit()
  expect_equal(fits$fn$data$time, d447$tt * 2)
  expect_equal(fits$const$data$time, d447$tt + 1)
})

test_that("the parser's own variables are not visible to the response", {
  # `data` is an argument of the parser. A caller's `data` must win.
  local_fit <- function() {
    data <- 2 # nolint: object_usage_linter.
    fit447(survival::Surv(tt * data, ss) ~ age)
  }
  expect_equal(local_fit()$data$time, d447$tt * 2)
})

test_that("a column of `data` still wins over the formula's environment", {
  local_fit <- function() {
    tt <- rep(1, nrow(d447)) # nolint: object_usage_linter.
    fit447(survival::Surv(tt, ss) ~ age)
  }
  expect_equal(local_fit()$data$time, d447$tt)
})

test_that("a formula without an environment still fits from `data`", {
  f <- survival::Surv(tt, ss) ~ age
  environment(f) <- NULL
  expect_null(environment(f))
  expect_equal(fit447(f)$data$time, d447$tt)
})

test_that("a stepwise refit resolves a response local as the fit did", {
  # Candidate refits rebuild the formula; each must keep the original
  # formula's environment, or the response would change between the fit
  # and its refits.
  assign("cutoff_447", 100, envir = globalenv())
  withr::defer(rm("cutoff_447", envir = globalenv()))
  d <- d447
  d$z <- stats::rnorm(nrow(d))
  local_sw <- function() {
    cutoff_447 <- 3
    base <- hazard(survival::Surv(pmin(tt, cutoff_447), ss) ~ age, data = d,
                   dist = "weibull", theta = c(mu = 0.2, nu = 1, 0),
                   fit = TRUE)
    suppressWarnings(hzr_stepwise(base, scope = ~ z, data = d,
                                  slentry = 0.99, criterion = "wald"))
  }
  # hzr_stepwise() reports its steps as it goes; keep that out of the log.
  utils::capture.output(sw <- local_sw())
  expect_gt(nrow(sw$steps), 0L)
  expect_equal(sw$data$time, pmin(d447$tt, 3))
})
