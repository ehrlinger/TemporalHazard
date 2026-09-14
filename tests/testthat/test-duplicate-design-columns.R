# tests/testthat/test-duplicate-design-columns.R
#
# A factor's dummy columns are named <factor><level>, so a factor `g` with
# level `b` and a numeric column `gb` both produce a design column `gb`.
# hazard() used to fit that design without a word: coef() carried two `gb`
# names, and predict() on a one-row newdata returned two values. Every
# name-based step downstream (coef names, predict by name, stepwise scope,
# bootstrap) assumes the names are unique, so hazard() now refuses the design
# and names the colliding columns.

.dup_data <- function(numeric_name = "gb") {
  set.seed(7)
  n <- 300
  d <- data.frame(time = rexp(n, 0.2), status = rbinom(n, 1, 0.7),
                  g = factor(sample(c("a", "b"), n, TRUE)), v = runif(n))
  names(d)[4L] <- numeric_name
  d
}

test_that("a single-distribution formula fit refuses duplicated design columns", {
  expect_error(
    hazard(survival::Surv(time, status) ~ g + gb, data = .dup_data(),
           dist = "weibull", theta = c(mu = 0.2, nu = 1, gb = 0, gb = 0),
           fit = TRUE),
    "duplicated column name.*'gb'"
  )
})

test_that("a multiphase phase formula refuses duplicated design columns", {
  expect_error(
    hazard(survival::Surv(time, status) ~ 1, data = .dup_data(),
           dist = "multiphase",
           phases = list(
             early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                                  fixed = "shapes", formula = ~ g + gb),
             constant = hzr_phase("constant")
           ),
           fit = TRUE, control = list(n_starts = 1L, maxit = 50L)),
    "phase 'early'.*duplicated column name.*'gb'"
  )
})

test_that("the vector interface refuses an x with duplicated colnames", {
  d <- .dup_data()
  x <- cbind(d$gb, runif(nrow(d)), runif(nrow(d)))
  colnames(x) <- c("gb", "u", "gb")
  expect_error(
    hazard(time = d$time, status = d$status, x = x, dist = "weibull",
           theta = c(mu = 0.2, nu = 1, 0, 0, 0), fit = TRUE),
    "duplicated column name.*'gb'"
  )
})

test_that("unnamed columns of a vector-interface x are not a collision", {
  # cbind(a = v1, v2, v3) names its columns c("a", "", ""); the empty names
  # are absent names, not a repeated one, and this fit ran on main.
  d <- .dup_data()
  x <- cbind(a = d$gb, runif(nrow(d)), runif(nrow(d)))
  expect_identical(colnames(x), c("a", "", ""))
  f <- hazard(time = d$time, status = d$status, x = x, dist = "weibull",
              theta = c(mu = 0.2, nu = 1, 0, 0, 0), fit = TRUE)
  expect_true(isTRUE(f$fit$converged))
  expect_length(stats::predict(f, type = "survival"), nrow(d))
})

test_that("NA column names of a vector-interface x are not a collision", {
  d <- .dup_data()
  x <- cbind(d$gb, runif(nrow(d)), runif(nrow(d)))
  colnames(x) <- c("a", NA_character_, NA_character_)
  f <- hazard(time = d$time, status = d$status, x = x, dist = "weibull",
              theta = c(mu = 0.2, nu = 1, 0, 0, 0), fit = TRUE)
  expect_true(isTRUE(f$fit$converged))
})

test_that("time_windows expansion of unnamed columns is refused", {
  # The expansion names each window's column <name>_w<k>, so absent names,
  # accepted above, become "_w1" or "NA_w1" twice: a duplicate the check on
  # `x` cannot see. The check runs on the expanded design as well.
  d <- .dup_data()
  x <- cbind(a = d$gb, runif(nrow(d)), runif(nrow(d)))
  expect_error(
    hazard(time = d$time, status = d$status, x = x, dist = "weibull",
           theta = c(mu = 0.2, nu = 1, rep(0, 6)), time_windows = 2,
           fit = TRUE),
    "duplicated column name.*'_w1'"
  )
  colnames(x) <- c("a", NA_character_, NA_character_)
  expect_error(
    hazard(time = d$time, status = d$status, x = x, dist = "weibull",
           theta = c(mu = 0.2, nu = 1, rep(0, 6)), time_windows = 2,
           fit = TRUE),
    "duplicated column name.*'NA_w1'"
  )
})

test_that("time_windows expansion of unique names is not refused", {
  d <- .dup_data()
  x <- cbind(a = d$gb, b = runif(nrow(d)), c = runif(nrow(d)))
  f <- hazard(time = d$time, status = d$status, x = x, dist = "weibull",
              theta = c(mu = 0.2, nu = 1, rep(0, 6)), time_windows = 2,
              fit = TRUE)
  expect_true(isTRUE(f$fit$converged))
})

test_that("a fit without a collision is unchanged", {
  # Coefficients pinned against origin/main 4b68020, before the refusal.
  f <- hazard(survival::Surv(time, status) ~ g + u, data = .dup_data("u"),
              dist = "weibull", theta = c(mu = 0.2, nu = 1, gb = 0, u = 0),
              fit = TRUE)
  expect_identical(colnames(f$data$x), c("gb", "u"))
  expect_equal(
    coef(f),
    c(mu = 0.138345781331580, nu = 1.123430810549730,
      gb = -0.290921889612059, u = 0.299436541170797),
    tolerance = 1e-8
  )
})
