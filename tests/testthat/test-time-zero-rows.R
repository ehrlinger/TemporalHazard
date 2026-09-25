# test-time-zero-rows.R
# A row at time 0 is dropped before fitting, whatever its status, as PROC
# HAZARD drops it at input (#374): hazard/src/hazard/readt.c:12-14.

tz_data <- function(n = 60, seed = 374) {
  withr::local_seed(seed)
  data.frame(time = stats::rexp(n, 0.4) + 0.01,
             status = stats::rbinom(n, 1, 0.7),
             x = stats::rnorm(n))
}

tz_theta <- list(exponential = c(0), weibull = c(0.3, 1.2),
                 lognormal = c(0.5, 1), loglogistic = c(0.5, 1))

test_that("an event at time 0 fits as if the row were absent, per family", {
  d <- tz_data()
  for (dist in names(tz_theta)) {
    ref <- suppressWarnings(hazard(time = d$time, status = d$status,
                                   dist = dist, theta = tz_theta[[dist]],
                                   fit = TRUE))
    w <- NULL
    got <- withCallingHandlers(
      hazard(time = c(0, d$time), status = c(1, d$status), dist = dist,
             theta = tz_theta[[dist]], fit = TRUE),
      hzr_time_zero_dropped = function(e) {
        w <<- e
        invokeRestart("muffleWarning")
      },
      warning = function(e) invokeRestart("muffleWarning")
    )
    expect_s3_class(w, "hzr_time_zero_dropped")
    expect_match(conditionMessage(w), "^1 row\\(s\\) with time = 0")
    # Not the clamp, not log(double.xmax): the likelihood of the other rows.
    expect_lt(abs(got$fit$objective), 1e6)
    expect_equal(got$fit$objective, ref$fit$objective, tolerance = 1e-8,
                 info = dist)
    expect_equal(coef(got), coef(ref), tolerance = 1e-6, info = dist)
    expect_identical(got$data$dropped_time_zero, 1L)
    expect_identical(length(got$data$time), nrow(d))
  }
})

test_that("the rule is status-blind and tests the upper bound only", {
  d <- tz_data()
  # Right- and left-censored rows at 0 are dropped too, as readt.c drops them.
  got <- suppressWarnings(hazard(
    time = c(0, 0, d$time), status = c(0, -1, d$status),
    dist = "weibull", theta = tz_theta$weibull, fit = TRUE))
  expect_identical(got$data$dropped_time_zero, 2L)
  # A lower bound of 0 is admissible (readct.c allows CT = 0): nothing drops.
  lo <- c(0, rep(0, nrow(d) - 1))
  st <- replace(d$status, 1, 2)
  kept <- suppressWarnings(hazard(time = d$time, status = st, time_lower = lo,
                                  dist = "weibull", theta = tz_theta$weibull,
                                  fit = TRUE))
  expect_identical(kept$data$dropped_time_zero, 0L)
  expect_identical(length(kept$data$time), nrow(d))
  # An interval row's `time` is its LOWER bound: (0, 2] is kept, and an
  # interval or left-censored row whose UPPER bound is 0 is dropped.
  iv <- suppressWarnings(hazard(
    time = c(0, 0, d$time), status = c(2, -1, d$status),
    time_upper = c(2, 0, d$time), dist = "weibull",
    theta = tz_theta$weibull, fit = TRUE))
  expect_identical(iv$data$dropped_time_zero_rows, 2L)
})

test_that("no warning and no drop when no row is at time 0", {
  d <- tz_data()
  expect_no_warning(f <- hazard(time = d$time, status = d$status,
                                dist = "exponential", theta = 0, fit = TRUE),
                    class = "hzr_time_zero_dropped")
  expect_identical(f$data$dropped_time_zero, 0L)
})

test_that("the formula interface drops the row from data too, so phases align", {
  skip_if_not_installed("numDeriv")
  d <- tz_data()
  d0 <- rbind(data.frame(time = 0, status = 1, x = 0.3), d)
  ph <- list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0,
                               formula = ~ x),
             const = hzr_phase("constant"))
  ref <- suppressWarnings(hazard(survival::Surv(time, status) ~ 1, data = d,
                                 dist = "multiphase", phases = ph, fit = TRUE,
                                 control = list(n_starts = 1L)))
  got <- suppressWarnings(hazard(survival::Surv(time, status) ~ 1, data = d0,
                                 dist = "multiphase", phases = ph, fit = TRUE,
                                 control = list(n_starts = 1L)))
  expect_identical(got$data$dropped_time_zero, 1L)
  expect_identical(nrow(got$data$frame), nrow(d))
  expect_equal(got$fit$objective, ref$fit$objective, tolerance = 1e-8)
  # A downstream consumer: residuals and predictions are per retained row.
  expect_length(stats::predict(got, type = "survival"), nrow(d))
})

test_that("hzr_stepwise() given the original frame is aligned with the fit", {
  # The package's own examples pass the frame given to hazard(); it still
  # holds the time-0 row the fit dropped. With the score criterion that was a
  # hard stop blaming NA removal.
  withr::local_seed(1)
  n <- 80
  d <- data.frame(time = c(0, stats::rexp(n, 0.4) + 0.01),
                  status = c(1, stats::rbinom(n, 1, 0.7)),
                  x1 = stats::rnorm(n + 1), x2 = stats::rnorm(n + 1))
  f <- suppressWarnings(hazard(survival::Surv(time, status) ~ 1, data = d,
                               dist = "weibull", theta = c(0.3, 1.2),
                               fit = TRUE))
  expect_identical(f$data$dropped_time_zero_rows, 1L)
  sw_orig <- suppressWarnings(hzr_stepwise(f, scope = c("x1", "x2"), data = d,
                                           trace = FALSE))
  sw_kept <- suppressWarnings(hzr_stepwise(f, scope = c("x1", "x2"),
                                           data = f$data$frame, trace = FALSE))
  expect_identical(sw_orig$steps, sw_kept$steps)
  expect_equal(sw_orig$fit$fit$objective, sw_kept$fit$fit$objective)
})

test_that("a data-dependent formula term is built without the dropped row", {
  # scale(age) centres on the data it is evaluated on; the dropped row
  # (age 150) must not move that centre.
  withr::local_seed(2)
  n <- 80
  d <- data.frame(time = stats::rexp(n, 0.4) + 0.01,
                  status = stats::rbinom(n, 1, 0.7),
                  age = stats::rnorm(n, 60, 10))
  d0 <- rbind(data.frame(time = 0, status = 1, age = 150), d)
  f <- function(dd) {
    suppressWarnings(hazard(survival::Surv(time, status) ~ scale(age),
                            data = dd, dist = "weibull",
                            theta = c(0.3, 1.2, 0), fit = TRUE))
  }
  with0 <- f(d0)
  without <- f(d)
  expect_equal(coef(with0), coef(without), tolerance = 1e-6)
  # Predictions on new data rebuild scale(age) from the stored design.
  nd <- data.frame(age = c(50, 70), time = 2)
  expect_equal(stats::predict(with0, newdata = nd, type = "survival"),
               stats::predict(without, newdata = nd, type = "survival"),
               tolerance = 1e-6)
})

test_that("hzr_stepwise() does not trim a frame that is not the fit's", {
  withr::local_seed(1)
  n <- 80
  d <- data.frame(time = c(0, stats::rexp(n, 0.4) + 0.01),
                  status = c(1, stats::rbinom(n, 1, 0.7)),
                  x1 = stats::rnorm(n + 1), x2 = stats::rnorm(n + 1))
  f <- suppressWarnings(hazard(survival::Surv(time, status) ~ 1, data = d,
                               dist = "weibull", theta = c(0.3, 1.2),
                               fit = TRUE))
  other <- d
  other$time <- stats::rexp(n + 1, 0.4) + 0.01  # same size, different cohort
  expect_error(suppressWarnings(hzr_stepwise(f, scope = c("x1", "x2"),
                                             data = other, trace = FALSE)),
               "rows but the fitted model used")
})
