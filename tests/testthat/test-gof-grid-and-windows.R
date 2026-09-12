# hzr_gof() items from Copilot's review of #285. Each test failed on 6538323.

.gof_grid_data <- function(n = 120) {
  set.seed(2851)
  stop_t <- round(stats::rexp(n, 0.3) + 0.1, 2)
  start_t <- ifelse(stats::runif(n) < 0.5,
                    round(stop_t * stats::runif(n, 0.1, 0.8), 2), 0)
  start_t[start_t >= stop_t] <- 0
  data.frame(start = start_t, stop = stop_t,
             event = stats::rbinom(n, 1, 0.7))
}

# At risk at t: followed from 0 and not yet out, or entered before t and not
# yet out. Counted subject by subject, as an oracle for hzr_gof()'s n_risk.
.brute_n_risk <- function(t, entry, time) {
  vapply(t, function(u) sum((entry < u | entry == 0) & time >= u), numeric(1))
}

test_that("n_risk at a custom time_grid is the risk set at that time", {
  d <- .gof_grid_data()
  kt <- sort(unique(d$stop))
  # One point before any exit, then points between exits.
  grid <- c(0.05, (kt[1:6] + kt[2:7]) / 2)
  for (entry in list(NULL, d$start)) {
    fit <- suppressWarnings(hazard(time = d$stop, status = d$event,
                                   time_lower = entry, dist = "weibull",
                                   theta = c(0.3, 1), fit = TRUE))
    e <- if (is.null(entry)) rep(0, nrow(d)) else entry
    expect_equal(hzr_gof(fit, time_grid = grid)$n_risk,
                 .brute_n_risk(grid, e, d$stop))
    # On the default grid the same count is survfit's own n.risk.
    g <- hzr_gof(fit)
    expect_equal(g$n_risk, .brute_n_risk(g$time, e, d$stop))
  }
})

test_that("n_risk at a seq() grid counts the exits tied at each grid time", {
  set.seed(2854)
  n <- 200
  stop_t <- pmax(round(stats::rexp(n, 0.8), 1), 0.1)
  fit <- suppressWarnings(hazard(time = stop_t,
                                 status = stats::rbinom(n, 1, 0.7),
                                 dist = "weibull", theta = c(0.3, 1),
                                 fit = TRUE))
  grid <- seq(0.1, 1, by = 0.1)
  # seq() lands some points a few ulps off the data (0.30000000000000004).
  expect_false(all(grid %in% stop_t))
  g <- hzr_gof(fit, time_grid = grid)
  expect_equal(g$n_risk, .brute_n_risk(round(grid, 10), rep(0, n), stop_t))
  expect_gt(sum(g$n_event), 0)
})

test_that("an unsorted or repeated time_grid is sorted and de-duplicated", {
  d <- .gof_grid_data()
  fit <- suppressWarnings(hazard(time = d$stop, status = d$event,
                                 dist = "weibull", theta = c(0.3, 1),
                                 fit = TRUE))
  kt <- sort(unique(d$stop))
  sorted <- hzr_gof(fit, time_grid = kt[c(5, 10)])
  shuffled <- hzr_gof(fit, time_grid = kt[c(10, 5, 10)])
  expect_equal(shuffled$time, kt[c(5, 10)])
  expect_equal(shuffled$cum_observed, sorted$cum_observed)
  expect_equal(shuffled$cum_expected, sorted$cum_expected)
  expect_error(hzr_gof(fit, time_grid = c(1, NA)), "time_grid")
  expect_error(hzr_gof(fit, time_grid = c(-1, 1)), "time_grid")
})

test_that("km_surv at time 0 is the Kaplan-Meier value when an event is at 0", {
  d <- .gof_grid_data()
  d$stop[1:2] <- 0
  d$event[1:2] <- c(1, 0)
  fit <- suppressWarnings(hazard(time = d$stop, status = d$event,
                                 dist = "exponential",
                                 theta = c(log_rate = log(0.3)), fit = TRUE))
  km <- survival::survfit(survival::Surv(d$stop, d$event) ~ 1)
  expect_equal(km$time[1], 0)
  expect_lt(km$surv[1], 1)
  expect_no_warning(g <- hzr_gof(fit))
  # The default grid is the Kaplan-Meier times, so the columns coincide.
  expect_equal(g$km_surv, km$surv)
})

test_that("with time_windows, H(entry) uses the exit-time design, as the likelihood does", {
  set.seed(2853)
  n <- 400
  x <- stats::rnorm(n)
  cut <- 2
  # Piecewise-exponential exits: x raises the hazard before the cut and
  # lowers it after, so the two window coefficients differ.
  r1 <- 0.3 * exp(1.5 * x)
  r2 <- 0.3 * exp(-0.5 * x)
  e <- stats::rexp(n)
  stop_t <- ifelse(e < r1 * cut, e / r1, cut + (e - r1 * cut) / r2)
  start_t <- ifelse(stats::runif(n) < 0.6,
                    stop_t * stats::runif(n, 0.1, 0.9), 0)
  expect_gt(sum(start_t > 0 & start_t <= cut & stop_t > cut), 20)
  fit <- suppressWarnings(hazard(
    time = stop_t, status = stats::rbinom(n, 1, 0.8), time_lower = start_t,
    x = matrix(x, ncol = 1, dimnames = list(NULL, "x")), time_windows = cut,
    dist = "weibull", theta = c(0.3, 1, 0, 0), fit = TRUE,
    control = list(reltol = 1e-12)))
  expect_true(fit$fit$converged)
  s <- attr(hzr_gof(fit), "summary")
  # A Weibull fit conserves events at its maximum, in the likelihood's terms.
  expect_equal(s$total_expected / s$total_observed, 1, tolerance = 1e-6)
})
