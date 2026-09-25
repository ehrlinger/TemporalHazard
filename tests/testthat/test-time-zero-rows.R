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
  # A downstream consumer: predictions are per retained row, and equal the
  # reference fit's.
  expect_equal(stats::predict(got, type = "survival"),
               stats::predict(ref, type = "survival"), tolerance = 1e-6)
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

test_that("a list of columns given as data is filtered like a data frame", {
  d <- tz_data()
  lst <- list(tt = c(0, d$time), ss = c(1, d$status))
  f <- suppressWarnings(hazard(time = tt, status = ss, data = lst,
                               dist = "weibull", theta = tz_theta$weibull,
                               fit = TRUE))
  expect_identical(f$data$dropped_time_zero, 1L)
  expect_identical(lengths(f$data$frame), c(tt = nrow(d), ss = nrow(d)))
  expect_equal(f$data$frame$tt, f$data$time)
})

test_that("hzr_bootstrap() still refuses Surv() vectors from outside data", {
  # The #278 guard compared a vector's length with the stored frame, which
  # has lost the time-0 row; the caller's vector has not. It let the call
  # through and returned identical replicates as full success.
  # The Surv() vectors live where the formula interface looks for a name that
  # is not a column: the global environment (a local frame is not searched,
  # on main as here).
  withr::local_seed(3)
  tt <- stats::rexp(60, 0.4) + 0.01
  tt[4] <- 0
  assign("tz_tt", tt, envir = globalenv())
  assign("tz_ss", stats::rbinom(60, 1, 0.7), envir = globalenv())
  withr::defer(rm("tz_tt", "tz_ss", envir = globalenv()))
  dd <- data.frame(age = stats::rnorm(60))
  f <- suppressWarnings(hazard(survival::Surv(tz_tt, tz_ss) ~ 1, data = dd,
                               dist = "weibull", theta = c(0.5, 1),
                               fit = TRUE))
  expect_identical(f$data$dropped_time_zero, 1L)
  expect_error(suppressWarnings(hzr_bootstrap(f, n_boot = 3)),
               "not columns of it")
})

test_that("a rebuilt design is validated like any other", {
  # After the drop, scale() of a column that was constant on the retained
  # rows is NaN; the rebuilt design lost every row and the fit returned the
  # optimizer's clamp with converged = TRUE.
  withr::local_seed(4)
  dd <- data.frame(tm = c(0, stats::rexp(59, 0.4) + 0.01),
                   st = stats::rbinom(60, 1, 0.7),
                   age = stats::rnorm(60), redo = c(1, rep(0, 59)))
  expect_error(suppressWarnings(hazard(
    survival::Surv(tm, st) ~ age + scale(redo), data = dd, dist = "weibull",
    theta = c(0.5, 1, 0, 0), fit = TRUE)),
    "Predictor rows must match the length of 'time'", fixed = TRUE)
})

test_that("a term read from outside data falls back to subsetting, and says so", {
  withr::local_seed(5)
  n <- 60
  zz <- stats::rnorm(n + 1)
  dd <- data.frame(time = c(0, stats::rexp(n, 0.4) + 0.01),
                   status = c(1, stats::rbinom(n, 1, 0.7)))
  w <- NULL
  f <- withCallingHandlers(
    hazard(survival::Surv(time, status) ~ zz, data = dd, dist = "weibull",
           theta = c(0.5, 1, 0), fit = TRUE),
    hzr_time_zero_dropped = function(e) {
      w <<- e
      invokeRestart("muffleWarning")
    },
    warning = function(e) invokeRestart("muffleWarning"))
  expect_match(conditionMessage(w), "read values outside `data`", fixed = TRUE)
  expect_true(isTRUE(f$data$time_zero_design_subset))
  # The RIGHT rows were kept: the same fit as zz put in data without the row.
  ref <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ zz, dist = "weibull",
    data = data.frame(dd[-1, ], zz = zz[-1]), theta = c(0.5, 1, 0),
    fit = TRUE))
  expect_equal(unname(f$data$x[, "zz"]), zz[-1])
  expect_equal(unname(coef(f)), unname(coef(ref)), tolerance = 1e-6)
  # hzr_stepwise() cannot refit such a design against `data`, and says so.
  dd$x1 <- stats::rnorm(n + 1)
  expect_error(hzr_stepwise(f, scope = "x1", data = dd, trace = FALSE),
               "reads a per-row value from outside `data`", fixed = TRUE)
})

test_that("a Surv() response read from outside data falls back too", {
  # Worked on main; the rebuild must not turn it into an error.
  withr::local_seed(6)
  n <- 60
  gt <- c(0, stats::rexp(n, 0.4) + 0.01)
  assign("tz_gt", gt, envir = globalenv())
  withr::defer(rm("tz_gt", envir = globalenv()))
  dat <- data.frame(st = c(1, stats::rbinom(n, 1, 0.7)),
                    age = stats::rnorm(n + 1))
  f <- suppressWarnings(hazard(survival::Surv(tz_gt, st) ~ age, data = dat,
                               dist = "weibull", theta = c(0.5, 1, 0),
                               fit = TRUE))
  expect_identical(f$data$dropped_time_zero, 1L)
  expect_equal(f$data$time, gt[-1])
})

test_that("every row at time 0 reaches the no-observations stop", {
  d0 <- data.frame(tm = c(0, 0, 0), st = c(1, 0, 1), g = c("a", "b", "a"))
  expect_error(suppressWarnings(hazard(survival::Surv(tm, st) ~ g, data = d0,
                                       dist = "weibull", theta = c(0.5, 1, 0),
                                       fit = TRUE)),
               "no observations")
})

test_that("weights stay aligned with the rows that remain", {
  d <- tz_data()
  wts <- stats::runif(nrow(d), 0.5, 2)
  ref <- suppressWarnings(hazard(time = d$time, status = d$status,
                                 weights = wts, dist = "weibull",
                                 theta = tz_theta$weibull, fit = TRUE))
  got <- suppressWarnings(hazard(time = c(0, d$time), status = c(1, d$status),
                                 weights = c(5, wts), dist = "weibull",
                                 theta = tz_theta$weibull, fit = TRUE))
  expect_equal(got$data$weights, wts)
  expect_equal(got$fit$objective, ref$fit$objective, tolerance = 1e-8)
  # Formula path, weights as a column.
  d$w <- wts
  d0 <- rbind(data.frame(time = 0, status = 1, x = 0, w = 5), d)
  reff <- suppressWarnings(hazard(survival::Surv(time, status) ~ x, data = d,
                                  weights = w, dist = "weibull",
                                  theta = c(0.3, 1.2, 0), fit = TRUE))
  gotf <- suppressWarnings(hazard(survival::Surv(time, status) ~ x, data = d0,
                                  weights = w, dist = "weibull",
                                  theta = c(0.3, 1.2, 0), fit = TRUE))
  expect_equal(gotf$fit$objective, reff$fit$objective, tolerance = 1e-8)
})

test_that("a dropped row's own values are not checked", {
  d <- tz_data()
  # An NA weight, and an entry after exit, on the row that is dropped.
  expect_no_error(suppressWarnings(hazard(
    time = c(0, d$time), status = c(1, d$status),
    weights = c(NA, rep(1, nrow(d))), dist = "weibull",
    theta = tz_theta$weibull, fit = TRUE)))
  expect_no_error(suppressWarnings(hazard(
    time = c(0, d$time), status = c(1, d$status),
    time_lower = c(2, rep(0, nrow(d))), dist = "weibull",
    theta = tz_theta$weibull, fit = TRUE)))
})
