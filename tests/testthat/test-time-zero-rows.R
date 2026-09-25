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

test_that("hzr_stepwise() does not trim a frame whose dropped rows differ", {
  withr::local_seed(1)
  n <- 80
  d <- data.frame(time = c(0, stats::rexp(n, 0.4) + 0.01),
                  status = c(1, stats::rbinom(n, 1, 0.7)),
                  x1 = stats::rnorm(n + 1), x2 = stats::rnorm(n + 1))
  f <- suppressWarnings(hazard(survival::Surv(time, status) ~ 1, data = d,
                               dist = "weibull", theta = c(0.3, 1.2),
                               fit = TRUE))
  # Same retained rows, but the first row is a real observation elsewhere.
  other <- d
  other$time[1] <- 2
  expect_error(suppressWarnings(hzr_stepwise(f, scope = c("x1", "x2"),
                                             data = other, trace = FALSE)),
               "rows but the fitted model used")
})

test_that("a column named \"\" and a row at time 0 together (#470, #374)", {
  # Both drops happen in hazard()'s data preparation; each must happen, and
  # neither may undo the other.
  withr::local_seed(10)
  n <- 60
  d <- data.frame(t = stats::rexp(n) + 0.1, s = stats::rbinom(n, 1, 0.8),
                  age = stats::rnorm(n), mal = stats::rbinom(n, 1, 0.4),
                  x = stats::rnorm(n))
  d0 <- rbind(data.frame(t = 0, s = 1, age = 0.5, mal = 1, x = 0.1), d)
  names(d0)[names(d0) == "x"] <- ""
  expect_true(any(names(d0) == ""))
  clean <- d[, c("t", "s", "age", "mal")]
  fit <- function(dd) {
    suppressWarnings(hazard(survival::Surv(t, s) ~ age, data = dd,
                            dist = "weibull", theta = c(0.5, 1, 0),
                            fit = TRUE))
  }
  got <- fit(d0)
  ref <- fit(clean)
  expect_identical(got$data$dropped_time_zero, 1L)
  expect_equal(got$fit$objective, ref$fit$objective, tolerance = 1e-8)
  expect_equal(unname(coef(got)), unname(coef(ref)), tolerance = 1e-6)

  sw <- function(base, dd) {
    suppressWarnings(hzr_stepwise(base, scope = ~ age + mal, data = dd,
                                  direction = "forward", slentry = 0.9,
                                  trace = FALSE))
  }
  sw_got <- sw(got, d0)
  sw_ref <- sw(ref, clean)
  # The control must take a step, or equal screens prove nothing.
  expect_gt(length(sw_ref$fit$theta), length(ref$fit$theta))
  expect_equal(sw_got$fit$objective, sw_ref$fit$objective,
               tolerance = 1e-8)
  expect_equal(unname(sw_got$fit$theta), unname(sw_ref$fit$theta),
               tolerance = 1e-6)
})

test_that("a data-dependent term is built on the retained rows", {
  # The fit is built as if the dropped rows had not been given, so every
  # consumer that rebuilds from the stored frame agrees with it (John,
  # 2026-09-25).
  withr::local_seed(2)
  n <- 80
  d <- data.frame(time = stats::rexp(n, 0.4) + 0.01,
                  status = stats::rbinom(n, 1, 0.7),
                  age = stats::rnorm(n, 60, 10), bmi = stats::rnorm(n, 27, 4))
  d0 <- rbind(data.frame(time = 0, status = 1, age = 150, bmi = 40), d)
  fit <- function(dd) {
    suppressWarnings(hazard(survival::Surv(time, status) ~ scale(age),
                            data = dd, dist = "weibull",
                            theta = c(0.3, 1.2, 0), fit = TRUE))
  }
  with0 <- fit(d0)
  without <- fit(d)
  expect_equal(coef(with0), coef(without), tolerance = 1e-6)
  nd <- data.frame(age = c(50, 70), time = 2)
  expect_equal(stats::predict(with0, newdata = nd, type = "survival"),
               stats::predict(without, newdata = nd, type = "survival"),
               tolerance = 1e-6)
  # A consumer that rebuilds from the stored frame reproduces the fit.
  expect_equal(as.numeric(hzr_evaluate(with0, coef(with0))$logLik),
               with0$fit$objective, tolerance = 1e-8)
})

test_that("a phase formula is built on the retained rows too", {
  skip_if_not_installed("numDeriv")
  withr::local_seed(8)
  n <- 80
  d <- data.frame(time = stats::rexp(n, 0.4) + 0.01,
                  status = stats::rbinom(n, 1, 0.7),
                  age = stats::rnorm(n, 60, 10), bmi = stats::rnorm(n, 27, 4))
  d0 <- rbind(data.frame(time = 0, status = 1, age = 90, bmi = 40), d)
  ph <- list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0,
                               formula = ~ scale(age)),
             const = hzr_phase("constant"))
  fit <- function(dd) {
    suppressWarnings(hazard(survival::Surv(time, status) ~ 1, data = dd,
                            dist = "multiphase", fit = TRUE, phases = ph,
                            control = list(n_starts = 1L)))
  }
  got <- fit(d0)
  ref <- fit(d)
  expect_equal(got$fit$objective, ref$fit$objective, tolerance = 1e-8)
  expect_equal(unname(got$fit$x_list$early), unname(ref$fit$x_list$early))
  # The two consumers the review found disagreeing: hzr_evaluate() and the
  # multiphase score test, which rebuild the phase design from the frame.
  expect_equal(as.numeric(hzr_evaluate(got, got$fit$theta)$logLik),
               got$fit$objective, tolerance = 1e-8)
  sg <- suppressWarnings(hzr_stepwise(got, scope = list(early = ~ bmi),
                                      data = d0, direction = "forward",
                                      criterion = "score", slentry = 0.99,
                                      trace = FALSE))
  sr <- suppressWarnings(hzr_stepwise(ref, scope = list(early = ~ bmi),
                                      data = d, direction = "forward",
                                      criterion = "score", slentry = 0.99,
                                      trace = FALSE))
  expect_gt(nrow(sr$steps), 0L)
  expect_equal(sg$steps$p_value, sr$steps$p_value, tolerance = 1e-6)
})

test_that("a response that depends on the dropped rows is refused", {
  # Surv(time - min(time), status): without the row it makes 0, the
  # minimum changes, so there is no consistent fit.
  withr::local_seed(7)
  d <- data.frame(time = stats::rexp(60, 0.4) + 1,
                  status = stats::rbinom(60, 1, 0.7))
  expect_error(suppressWarnings(hazard(
    survival::Surv(time - min(time), status) ~ 1, data = d,
    dist = "weibull", theta = c(0.3, 1.2), fit = TRUE)),
    "depends on the 1 row(s) at time 0", fixed = TRUE)
})

test_that("a formula reading per-row values from outside data is refused", {
  withr::local_seed(5)
  n <- 60
  zz <- stats::rnorm(n + 1)
  dd <- data.frame(time = c(0, stats::rexp(n, 0.4) + 0.01),
                   status = c(1, stats::rbinom(n, 1, 0.7)))
  expect_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ zz, data = dd, dist = "weibull",
    theta = c(0.5, 1, 0), fit = TRUE)),
    "reads 'zz' from outside `data`", fixed = TRUE)
  # A phase formula too.
  dd$age <- stats::rnorm(n + 1)
  expect_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = dd, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", formula = ~ zz),
                  const = hzr_phase("constant")), fit = TRUE)),
    "reads 'zz' from outside `data`", fixed = TRUE)
  # Without a row at time 0 nothing changes: an outside vector still fits.
  zc <- zz[-1]
  expect_no_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ zc, data = dd[-1, ], dist = "weibull",
    theta = c(0.5, 1, 0), fit = TRUE)))
})
