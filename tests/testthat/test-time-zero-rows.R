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

test_that("a dropped row's weight and entry time are not checked", {
  # The drop runs before the weight and entry-time checks. (The bounds'
  # own type checks -- finite, non-negative -- still run on every row,
  # because a row's upper bound decides whether it is dropped.)
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

test_that("a helper reading an outside vector by position is refused (#590)", {
  # The check above reads the formula's own variables, and `g` is not one of
  # them: it is read inside `hf`. On the retained rows `hf` would pair them
  # with g[1:n], although `g` was aligned with every row given.
  withr::local_seed(5)
  n <- 50
  dd <- data.frame(time = c(0, 0, stats::rexp(n - 2, 0.3)),
                   status = stats::rbinom(n, 1, 0.7),
                   age = stats::rnorm(n, 50, 10))
  g <- stats::rnorm(n)
  hf <- function(a) a + g[seq_along(a)]
  expect_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ hf(age), data = dd, dist = "weibull",
    theta = c(0.3, 1, 0), fit = TRUE)),
    "the formula's 'hf(age)' depends on the order of the rows", fixed = TRUE)
  # A phase formula too.
  expect_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = dd, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                                    fixed = "shapes", formula = ~ hf(age)),
                  const = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1))),
    "the formula's 'hf(age)' depends on the order of the rows", fixed = TRUE)
  # Without a row at time 0 nothing is dropped, nothing is mis-paired, and
  # the helper still fits.
  dd1 <- dd
  dd1$time <- dd1$time + 1
  expect_no_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ hf(age), data = dd1, dist = "weibull",
    theta = c(0.3, 1, 0), fit = TRUE)))
})

test_that("a global formula with `.` beside the helper is checked (#590)", {
  # The global formula, with `.` beside the helper and a Surv() response
  # that reads a column through an expression.
  withr::local_seed(5)
  n <- 50
  dd <- data.frame(time = c(0, 0, stats::rexp(n - 2, 0.3)),
                   ev = "dead", z = stats::rnorm(n))
  g <- stats::rnorm(n)
  hf <- function(a) a + g[seq_along(a)]
  expect_error(suppressWarnings(hazard(
    survival::Surv(time, ev == "dead") ~ hf(z) + ., data = dd,
    dist = "weibull", theta = c(0.3, 1, 0, 0), fit = TRUE)),
    "the formula's 'hf(z)' depends on the order of the rows", fixed = TRUE)
})

test_that("the order check reads every row given, not the kept ones (#590)", {
  # Sorted by era, with the early rows at time 0: `late` is 0 over its first
  # 180 values, so a check of the 180 kept rows alone would see a constant,
  # while 20 of those rows are paired with the wrong era.
  withr::local_seed(5)
  n <- 200
  dd <- data.frame(time = c(rep(0, 20), stats::rexp(n - 20, 0.3)),
                   status = stats::rbinom(n, 1, 0.7),
                   age = stats::rnorm(n))
  late <- c(rep(0, 180), rep(1, 20))
  h <- function(a) a + late[seq_along(a)]
  expect_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ age + h(age), data = dd,
    dist = "weibull", theta = c(0.3, 1, 0, 0), fit = TRUE)),
    "the formula's 'h(age)' depends on the order of the rows", fixed = TRUE)
})

test_that("a helper that reads by position on some rows only is refused (#590)", {
  # Imputing an NA from an outside era indicator by position: the indicator
  # changes at row 101, which is not NA, so a check comparing neighbouring
  # rows sees no difference, but the 60 rows dropped before them move every
  # imputed row to another row's era.
  withr::local_seed(5)
  n <- 200
  gimp <- c(rep(0, 100), rep(1, 100))
  x <- stats::rnorm(n)
  x[setdiff(sample(61:n, 40), 101)] <- NA
  dd <- data.frame(time = c(rep(0, 60), stats::rexp(n - 60, 0.3)),
                   status = stats::rbinom(n, 1, 0.7), x = x)
  imp <- function(a) ifelse(is.na(a), gimp[seq_along(a)], a)
  expect_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ imp(x), data = dd, dist = "weibull",
    theta = c(0.1, 1, 0), fit = TRUE)),
    "the formula's 'imp(x)' depends on the order of the rows", fixed = TRUE)
  # On the fallback path too: poly() inside the helper refuses the NA only a
  # dropped row carries, and the era indicator is constant over its first
  # 160 values.
  dd <- data.frame(time = c(rep(0, 40), stats::rexp(160, 0.3)),
                   status = stats::rbinom(n, 1, 0.7), age = stats::rnorm(n))
  dd$age[1] <- NA
  late <- c(rep(0, 160), rep(1, 40))
  hp <- function(a) as.vector(poly(a, 1)) + late[seq_along(a)]
  expect_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = dd, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                                    fixed = "shapes", formula = ~ hp(age)),
                  const = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1))),
    "the formula's 'hp(age)' depends on the order of the rows", fixed = TRUE)
})

test_that("a column named \"\" does not switch the order check off (#590)", {
  # `.` would reach the column, the check's frame would fail on it, and the
  # check would pass while the fit, which drops the column, went ahead.
  withr::local_seed(5)
  n <- 50
  dd <- data.frame(time = c(0, 0, stats::rexp(n - 2, 0.3)),
                   status = stats::rbinom(n, 1, 0.7),
                   z = stats::rnorm(n), e = 1)
  names(dd)[4] <- ""
  g <- stats::rnorm(n)
  hf <- function(a) a + g[seq_along(a)]
  # The global formula.
  expect_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ hf(z) + ., data = dd,
    dist = "weibull", theta = c(0.3, 1, 0, 0), fit = TRUE)),
    "the formula's 'hf(z)' depends on the order of the rows", fixed = TRUE)
  # A phase formula.
  expect_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = dd, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                                    fixed = "shapes", formula = ~ hf(z) + .),
                  const = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1))),
    "the formula's 'hf(z)' depends on the order of the rows", fixed = TRUE)
})

test_that("one variable that fails on every row does not hide another (#590)", {
  # poly() refuses the NA that only a dropped row carries, so the check
  # evaluates it with that row patched; the helper beside it is still
  # checked.
  withr::local_seed(5)
  n <- 120
  dd <- data.frame(time = c(0, 0, 0, stats::rexp(n - 3, 0.3)),
                   status = stats::rbinom(n, 1, 0.7),
                   age = stats::rnorm(n, 50, 10), z = stats::rnorm(n))
  dd$age[1] <- NA
  g <- stats::rnorm(n)
  hf <- function(a) a + g[seq_along(a)]
  early <- function(f) {
    list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes", formula = f),
         const = hzr_phase("constant"))
  }
  expect_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = dd, dist = "multiphase",
    phases = early(~ poly(age, 2) + hf(z)),
    fit = TRUE, control = list(n_starts = 1))),
    "the formula's 'hf(z)' depends on the order of the rows", fixed = TRUE)
  # The helper inside poly() cannot be built on every row either; it is
  # evaluated with the dropped rows replaced by a kept one.
  expect_error(suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = dd, dist = "multiphase",
    phases = early(~ poly(hf(age), 2)),
    fit = TRUE, control = list(n_starts = 1))),
    "the formula's 'poly(hf(age), 2)' depends on the order of the rows",
    fixed = TRUE)
  # Without the helper, poly() over that NA still fits, as it did before.
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = dd, dist = "multiphase",
    phases = early(~ poly(age, 2)), fit = TRUE,
    control = list(n_starts = 1)))
  expect_identical(fit$data$dropped_time_zero, 3L)
})

test_that("terms computed from the rows as a set still fit after a drop (#590)", {
  # Known negatives for the order check: each of these follows its rows.
  withr::local_seed(5)
  n <- 50
  dd <- data.frame(time = c(0, 0, stats::rexp(n - 2, 0.3)),
                   status = stats::rbinom(n, 1, 0.7),
                   age = stats::rnorm(n, 50, 10),
                   grp = sample(c("a", "b", "c"), n, replace = TRUE))
  fits <- list(
    scale = list(~ scale(age), c(0.3, 1, 0)),
    poly_factor = list(~ poly(age, 2) + grp, c(0.3, 1, 0, 0, 0, 0)),
    ns = list(~ splines::ns(age, 3), c(0.3, 1, 0, 0, 0))
  )
  for (nm in names(fits)) {
    f <- stats::update(survival::Surv(time, status) ~ 1, fits[[nm]][[1]])
    fit <- suppressWarnings(hazard(f, data = dd, dist = "weibull",
                                   theta = fits[[nm]][[2]], fit = TRUE))
    expect_identical(fit$data$dropped_time_zero, 2L, info = nm)
    expect_true(isTRUE(fit$fit$converged), info = nm)
  }
  # scale(age) is still taken over the retained rows (#476).
  fit <- suppressWarnings(hazard(survival::Surv(time, status) ~ scale(age),
                                 data = dd, dist = "weibull",
                                 theta = c(0.3, 1, 0), fit = TRUE))
  expect_equal(unname(fit$data$x[, 1]),
               as.vector(scale(dd$age[dd$time > 0])))
})

test_that("the order check names only the columns that do not follow (#590)", {
  withr::local_seed(5)
  n <- 30
  d <- data.frame(age = stats::rnorm(n), z = stats::rnorm(n))
  g <- stats::rnorm(n + 5)
  hf <- function(a) a + g[seq_along(a)]
  expect_identical(
    .hzr_position_bound_columns(list(~ age + hf(z) + scale(age)), d, c(FALSE, rep(TRUE, n - 1))),
    "hf(z)")
  expect_identical(.hzr_position_bound_columns(list(~ age + z), d, c(FALSE, rep(TRUE, n - 1))),
                   character(0))
  # A list `data` with a scalar column, as hazard() accepts.
  expect_identical(
    .hzr_position_bound_columns(list(~ hf(age)), list(age = d$age, k = 1), c(FALSE, rep(TRUE, n - 1))),
    "hf(age)")
  # An outside vector that reads the same backwards over the kept rows: a
  # check that reversed the rows could not see it.
  gq <- rep(c(1, 2, 1), length.out = n)
  hq <- function(a) a + gq[seq_along(a)]
  expect_identical(.hzr_position_bound_columns(list(~ hq(age)), d, c(FALSE, rep(TRUE, n - 1))),
                   "hq(age)")
  # A factor with levels in order of appearance changes its reference level
  # when the rows are reordered, but its value on each row does not, so it
  # is accepted.
  d$ch <- sample(c("x", "y", "w"), n, replace = TRUE)
  expect_identical(.hzr_position_bound_columns(
    list(~ factor(ch, levels = unique(ch)) + age), d, c(FALSE, rep(TRUE, n - 1))),
    character(0))
  # A term that depends on the order of the rows is refused.
  expect_identical(.hzr_position_bound_columns(list(~ cumsum(z)), d, c(FALSE, rep(TRUE, n - 1))),
                   "cumsum(z)")
  # An NA in a column is compared as an NA, not as a moved value.
  d$age[3] <- NA
  expect_identical(
    .hzr_position_bound_columns(list(~ age + hf(z)), d, c(FALSE, rep(TRUE, n - 1))), "hf(z)")
  # Known negatives on every shape of `data`: own row names, a list with a
  # scalar column, a matrix column.
  rownames(d) <- paste0("r", seq_len(n))
  expect_identical(.hzr_position_bound_columns(list(~ age + ch), d, c(FALSE, rep(TRUE, n - 1))),
                   character(0))
  dl <- list(age = d$age, ch = d$ch, k = 1, m = cbind(d$z, d$z^2))
  expect_identical(.hzr_position_bound_columns(list(~ age + ch + m), dl, c(FALSE, rep(TRUE, n - 1))),
                   character(0))
  expect_identical(.hzr_position_bound_columns(list(~ age + hf(z)),
                                               c(dl, list(z = d$z)), c(FALSE, rep(TRUE, n - 1))),
                   "hf(z)")
  # An outside vector that alternates, or repeats with any period, still
  # is mis-paired, whatever the number of rows.
  for (m in c(n, n + 1)) {
    alt <- rep(c(0, 1), length.out = m)
    af <- function(a) alt[seq_along(a)]
    dm <- data.frame(age = stats::rnorm(m))
    expect_identical(.hzr_position_bound_columns(list(~ age + af(age)), dm, c(FALSE, rep(TRUE, m - 1))),
                     "af(age)", info = m)
    per <- rep(c("p", "q", "r"), length.out = m)
    pf <- function(a) factor(per[seq_along(a)])
    expect_identical(.hzr_position_bound_columns(list(~ pf(age)), dm, c(FALSE, rep(TRUE, m - 1))),
                     "pf(age)", info = m)
  }
  # A term that can be built in only one order of the rows is refused, and
  # named by its formula.
  bad <- d$z[2]  # first once the dropped row 1 is put last
  hord <- function(a) {
    if (isTRUE(all.equal(a[1], bad))) stop("order")
    a
  }
  expect_identical(.hzr_position_bound_columns(list(~ hord(z)), d, c(FALSE, rep(TRUE, n - 1))),
                   "hord(z)")
  # The tolerance follows the spread of the column, not the size of its
  # values: a large mean must not hide a mis-pairing.
  big <- 1.7e9 + stats::rnorm(n, 0, 100)
  gs <- stats::rnorm(n, 0, 1)
  hs <- function(a) a + gs[seq_along(a)]
  expect_identical(.hzr_position_bound_columns(list(~ hs(big)),
                                               data.frame(big = big),
                                               c(FALSE, rep(TRUE, n - 1))), "hs(big)")
  # With the dropped rows already last, a helper reading by position pairs
  # every kept row with its own value, and nothing is named.
  expect_identical(.hzr_position_bound_columns(
    list(~ hf(z) + cumsum(z)), d, c(rep(TRUE, n - 2), FALSE, FALSE)),
    character(0))
  # A variable that fails on every row and on the patched rows, but builds
  # on the kept rows alone, as a phase design is built, cannot be checked
  # and is refused.
  dk <- d
  dk$id <- c(NA, seq_len(n - 1))
  hid <- function(i) {
    if (anyNA(i) || anyDuplicated(i)) stop("ids")
    i + g[seq_along(i)]
  }
  expect_identical(.hzr_position_bound_columns(
    list(~ hid(id)), dk, c(FALSE, rep(TRUE, n - 1))), "hid(id)")
  # The documented limit: a helper that depends on how many rows it is
  # given sees n rows in both evaluations, so it is not seen.
  hm <- function(a) a * mean(g[seq_along(a)])
  expect_identical(.hzr_position_bound_columns(list(~ hm(z)), d, c(FALSE, rep(TRUE, n - 1))),
                   character(0))
})

test_that("an x, weights or data of the wrong length is refused, not misaligned", {
  # N2 (1.2.12 release review): a design k rows short was accepted when k
  # rows were at time 0, and paired with the wrong rows.
  withr::local_seed(3)
  n <- 60
  d <- data.frame(tm = stats::rexp(n, 0.3) + 0.05,
                  ev = stats::rbinom(n, 1, 0.7), age = stats::rnorm(n, 60, 10))
  d$age[40] <- NA
  d$tm[3] <- 0
  xx <- stats::model.matrix(~ age, d)[, -1, drop = FALSE]
  expect_equal(nrow(xx), n - 1)  # the premise: one row short
  expect_error(suppressWarnings(hazard(
    time = d$tm, status = d$ev, x = xx, dist = "weibull",
    theta = c(0.3, 1, 0), fit = TRUE)),
    "'x' has 59 row(s) but 'time' has 60", fixed = TRUE)
  expect_error(suppressWarnings(hazard(
    time = d$tm, status = d$ev, weights = rep(1, n - 1), dist = "weibull",
    theta = c(0.3, 1), fit = TRUE)),
    "'weights' has 59 row(s)", fixed = TRUE)
  expect_error(suppressWarnings(hazard(
    time = d$tm, status = d$ev, data = d[-60, ], dist = "multiphase",
    phases = list(early = hzr_phase("cdf", formula = ~ age),
                  const = hzr_phase("constant")), fit = TRUE)),
    "'data' has 59 row(s)", fixed = TRUE)
})

test_that("a list data with a scalar column is filtered; a ragged one is refused", {
  d <- tz_data()
  lst <- list(tt = c(0, d$time), ss = c(1, d$status), lab = "cohortA")
  f <- suppressWarnings(hazard(time = tt, status = ss, data = lst,
                               dist = "weibull", theta = tz_theta$weibull,
                               fit = TRUE))
  expect_identical(f$data$dropped_time_zero, 1L)
  expect_identical(f$data$frame$lab, "cohortA")
  lst$bad <- 1:3
  expect_error(suppressWarnings(hazard(time = tt, status = ss, data = lst,
                                       dist = "weibull",
                                       theta = tz_theta$weibull, fit = TRUE)),
               "'data$bad' has 3 row(s)", fixed = TRUE)
})

test_that("the length check refuses only inputs that would be misaligned", {
  withr::local_seed(11)
  n <- 40
  d1 <- data.frame(tm = stats::rexp(n, 0.3) + 0.05,
                   ev = stats::rbinom(n, 1, 0.7), age = stats::rnorm(n, 60, 10))
  d1$tm[3] <- 0
  d1$age[3] <- NA  # the dropped row is also the NA row
  # Formula path: model.matrix() drops the NA row, x is rebuilt anyway, and
  # the fit is the fit of the rows that remain.
  got <- suppressWarnings(hazard(survival::Surv(tm, ev) ~ age, data = d1,
                                 dist = "weibull", theta = c(0.3, 1, 0),
                                 fit = TRUE))
  ref <- suppressWarnings(hazard(survival::Surv(tm, ev) ~ age, data = d1[-3, ],
                                 dist = "weibull", theta = c(0.3, 1, 0),
                                 fit = TRUE))
  expect_equal(coef(got), coef(ref), tolerance = 1e-8)
  # Vector path: a longer `data` used only to look names up is left alone.
  lookup <- suppressWarnings(hazard(time = d1$tm[1:30], status = d1$ev[1:30],
                                    data = d1, dist = "weibull",
                                    theta = c(0.3, 1), fit = TRUE))
  lref <- suppressWarnings(hazard(time = d1$tm[c(1:2, 4:30)],
                                  status = d1$ev[c(1:2, 4:30)],
                                  dist = "weibull", theta = c(0.3, 1),
                                  fit = TRUE))
  expect_equal(coef(lookup), coef(lref), tolerance = 1e-8)
  # A scalar `weights` is not a row-aligned input, and is refused.
  expect_error(suppressWarnings(hazard(
    time = c(0, 1.5, 2), status = c(1, 1, 0), weights = 3,
    dist = "exponential", theta = 0, fit = TRUE)),
    "'weights' has 1 row(s)", fixed = TRUE)
})

test_that("an unnamed list column is named by position in the refusal", {
  d <- tz_data()
  lst <- list(tt = c(0, d$time), ss = c(1, d$status), 1:3)
  expect_error(suppressWarnings(hazard(time = tt, status = ss, data = lst,
                                       dist = "weibull",
                                       theta = tz_theta$weibull, fit = TRUE)),
               "'data[[3]]' has 3 row(s)", fixed = TRUE)
})

test_that("row numbers after a drop refer to the rows as given (#484)", {
  d <- tz_data()
  d0 <- rbind(data.frame(time = 0, status = 1, x = 0.1), d)
  # Original row 10 is the 9th row that remains.
  st <- d0$status
  st[10] <- NA
  expect_error(suppressWarnings(hazard(time = d0$time, status = st,
                                       dist = "weibull",
                                       theta = tz_theta$weibull, fit = TRUE)),
               "at index/indices 10\\.")
  st <- d0$status
  st[12] <- 7
  expect_error(suppressWarnings(hazard(time = d0$time, status = st,
                                       dist = "weibull",
                                       theta = tz_theta$weibull, fit = TRUE)),
               "at index/indices 12\\.")
  # objective = "sas": an interval row with upper <= lower, at original row 5.
  st <- d0$status
  st[5] <- 2
  lo <- d0$time
  up <- d0$time
  up[5] <- lo[5] / 2
  expect_error(suppressWarnings(hazard(
    time = d0$time, status = st, time_lower = lo, time_upper = up,
    dist = "multiphase", phases = list(c = hzr_phase("constant")),
    objective = "sas", fit = TRUE)),
    "at index/indices 5\\.")
  # With nothing dropped, positions are the rows given, as before.
  st <- d$status
  st[10] <- NA
  expect_error(suppressWarnings(hazard(time = d$time, status = st,
                                       dist = "weibull",
                                       theta = tz_theta$weibull, fit = TRUE)),
               "at index/indices 10\\.")
})

test_that("the re-parse after a drop does not repeat a term's warnings (#484)", {
  d <- tz_data()
  d0 <- rbind(data.frame(time = 0, status = 1, x = 0.1), d)
  noisy <- function(a) {
    warning("noisy term")
    a
  }
  count <- function(data) {
    k <- 0L
    withCallingHandlers(
      hazard(survival::Surv(time, status) ~ noisy(x), data = data,
             dist = "weibull", theta = c(tz_theta$weibull, 0), fit = TRUE),
      warning = function(w) {
        if (identical(conditionMessage(w), "noisy term")) k <<- k + 1L
        invokeRestart("muffleWarning")
      })
    k
  }
  k_kept <- count(d)
  expect_gt(k_kept, 0L)
  expect_identical(count(d0), k_kept)
  # A warning only the retained rows raise is not a repeat, and still shows.
  only_kept <- function(a) {
    if (length(a) == nrow(d)) warning("retained rows only")
    a
  }
  msgs <- character(0)
  withCallingHandlers(
    hazard(survival::Surv(time, status) ~ only_kept(x), data = d0,
           dist = "weibull", theta = c(tz_theta$weibull, 0), fit = TRUE),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_gt(sum(msgs == "retained rows only"), 0L)
})

test_that("hzr_stepwise() names the drop when given the data before it (#484)", {
  withr::local_seed(6)
  n <- 120
  dm <- data.frame(time = c(rep(0, 4), stats::rexp(n - 4, 0.3) + 0.01),
                   status = stats::rbinom(n, 1, 0.7),
                   x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  fit <- suppressWarnings(hazard(
    time = dm$time, status = dm$status, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1),
                  const = hzr_phase("constant")), fit = TRUE))
  expect_identical(fit$data$dropped_time_zero, 4L)
  msgs <- character(0)
  withCallingHandlers(
    hzr_stepwise(fit, data = dm, scope = list(const = ~ x1 + x2),
                 criterion = "wald", trace = FALSE),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  hit <- grep("it dropped at time 0", msgs, fixed = TRUE, value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, "has 120 rows: the fit's 116 plus the 4", fixed = TRUE)
  expect_match(hit, "(rows 1, 2, 3, 4 of", fixed = TRUE)
  # Given the rows that remain, there is nothing to say.
  msgs <- character(0)
  withCallingHandlers(
    hzr_stepwise(fit, data = dm[-(1:4), ], scope = list(const = ~ x1 + x2),
                 criterion = "wald", trace = FALSE),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("it dropped at time 0", msgs, fixed = TRUE)))
  # A formula fit given a frame whose dropped row differs is not trimmed
  # either, and stops on a row count; the drop is named first.
  df <- dm[-(2:4), ]
  ff <- suppressWarnings(hazard(survival::Surv(time, status) ~ 1, data = df,
                                dist = "weibull", theta = c(0.3, 1.2),
                                fit = TRUE))
  other <- df
  other$time[1] <- 2
  msgs <- character(0)
  expect_error(
    withCallingHandlers(
      hzr_stepwise(ff, scope = c("x1", "x2"), data = other, trace = FALSE),
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    "rows but the fitted model used")
  expect_true(any(grepl("plus the 1 it dropped at time 0 (rows 1 of", msgs,
                        fixed = TRUE)))
  # The row-order note compares `data` with the stored frame; it used to
  # compare the frame with the fit, and printed "116 rows, not the fit's 116".
  expect_true(any(grepl("`data` has 117 rows, not the 116 of the data frame",
                        msgs, fixed = TRUE)))
  expect_false(any(grepl("rows, not the fit's", msgs, fixed = TRUE)))
})

test_that("hzr_bootstrap() resamples the rows that remain (#484)", {
  d <- tz_data()
  d0 <- rbind(data.frame(time = 0, status = 1, x = 0.1), d)
  boot <- function(data) {
    f <- suppressWarnings(hazard(survival::Surv(time, status) ~ x,
                                 data = data, dist = "weibull",
                                 theta = c(tz_theta$weibull, 0), fit = TRUE))
    suppressWarnings(hzr_bootstrap(f, n_boot = 5, seed = 484))
  }
  with_zero <- boot(d0)
  without <- boot(d)
  # Premise: the replicates differ from one another, so equality below
  # compares resamples, not one repeated fit.
  est <- with_zero$replicates$estimate[with_zero$replicates$parameter == "x"]
  expect_gt(stats::sd(est), 0)
  expect_identical(with_zero$n_success, 5L)
  expect_equal(with_zero$replicates, without$replicates, tolerance = 1e-10)
})

test_that("a transformed response that does not depend on the dropped rows is fitted (#484)", {
  d <- tz_data()
  d0 <- rbind(data.frame(time = 0, status = 1, x = 0.1), d)
  # time / max(time) is the same on the rows that remain: max is not at 0.
  got <- suppressWarnings(hazard(
    survival::Surv(time / max(time), status) ~ x, data = d0,
    dist = "weibull", theta = c(tz_theta$weibull, 0), fit = TRUE))
  ref <- suppressWarnings(hazard(
    survival::Surv(time / max(time), status) ~ x, data = d,
    dist = "weibull", theta = c(tz_theta$weibull, 0), fit = TRUE))
  expect_identical(got$data$dropped_time_zero, 1L)
  expect_equal(coef(got), coef(ref), tolerance = 1e-8)
})
