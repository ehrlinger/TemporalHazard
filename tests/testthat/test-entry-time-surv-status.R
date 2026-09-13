# #253 x #226: a counting-type Surv passed as `status =` (#283 reads it by
# type, so Surv(start, stop, event) sets time_lower = start) must reach the
# same entry-time likelihood as the other two routes to an entry time: the
# formula interface and an explicit `time_lower =`. A parity test only sees
# the routes diverging, so each family also asserts that the entry time
# changed the fit.

.surv_status_data <- function(n = 160) {
  set.seed(2532)
  x <- stats::rnorm(n)
  stop_t <- stats::rexp(n, 0.3 * exp(0.3 * x)) + 0.1
  start_t <- ifelse(stats::runif(n) < 0.5, stop_t * stats::runif(n, 0.1, 0.8), 0)
  data.frame(start = start_t, stop = stop_t, event = stats::rbinom(n, 1, 0.7),
             x = x)
}

.surv_status_starts <- list(
  weibull     = c(0.3, 1, 0),
  exponential = c(log_rate = 0, 0),
  loglogistic = c(log_alpha = 0, log_beta = 0, 0),
  lognormal   = c(mu = 0, log_sigma = 0, 0)
)

test_that("a counting Surv as `status` fits like time_lower = and the formula", {
  d <- .surv_status_data()
  expect_gt(sum(d$start > 0), 50)
  xm <- matrix(d$x, ncol = 1, dimnames = list(NULL, "x"))
  for (dist in names(.surv_status_starts)) {
    fit <- function(...) {
      suppressWarnings(hazard(..., dist = dist,
                              theta = .surv_status_starts[[dist]],
                              fit = TRUE, control = list(reltol = 1e-12)))
    }
    f_surv <- fit(time = d$stop, status = survival::Surv(d$start, d$stop, d$event),
                  x = xm)
    f_tl   <- fit(time = d$stop, status = d$event, time_lower = d$start, x = xm)
    f_form <- fit(survival::Surv(start, stop, event) ~ x, data = d)
    f_none <- fit(time = d$stop, status = d$event, x = xm)
    # The Surv's start column is what the fit stored as the entry time.
    expect_equal(f_surv$data$time_lower, d$start, info = dist)
    expect_equal(f_surv$fit$objective, f_tl$fit$objective, tolerance = 1e-8,
                 info = dist)
    expect_equal(f_surv$fit$objective, f_form$fit$objective, tolerance = 1e-8,
                 info = dist)
    expect_equal(unname(coef(f_surv)), unname(coef(f_tl)), tolerance = 1e-6,
                 info = dist)
    # The entry time matters, so a route that dropped it would fail above.
    expect_gt(abs(f_surv$fit$objective - f_none$fit$objective), 1,
              label = paste(dist, "entry-time effect"))
  }
})

test_that("a counting Surv as `status` honours the entry rule in multiphase", {
  d <- .surv_status_data()
  ph <- list(early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                                  fixed = "shapes"),
             constant = hzr_phase("constant"))
  fit <- function(...) {
    suppressWarnings(suppressMessages(hazard(..., dist = "multiphase",
                                             phases = ph, fit = TRUE)))
  }
  f_surv <- fit(time = d$stop, status = survival::Surv(d$start, d$stop, d$event))
  f_tl   <- fit(time = d$stop, status = d$event, time_lower = d$start)
  f_none <- fit(time = d$stop, status = d$event)
  expect_equal(f_surv$fit$objective, f_tl$fit$objective, tolerance = 1e-8)
  expect_lt(f_surv$fit$objective, 0)
  expect_gt(abs(f_surv$fit$objective - f_none$fit$objective), 1)
})
