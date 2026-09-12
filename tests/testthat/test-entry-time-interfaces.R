# #253: how `time_lower` reaches each family, across statuses and interfaces.
#
# On status 0/1 rows `time_lower` is the counting-process entry time when
# 0 < time_lower < time. On status -1 (left-censored) rows it is not used at
# all: the bound is `time_upper`. On status 2 rows it is the interval's lower
# bound (pinned per family in test-entry-time-<family>.R). These tests pin the
# status -1 half and the formula/vector parity of the entry time.

.entry_dists <- c("weibull", "exponential", "loglogistic", "lognormal")

.entry_starts <- list(
  weibull     = c(0.5, 1, 0),
  exponential = c(log_rate = 0, 0),
  loglogistic = c(log_alpha = 0, log_beta = 0, 0),
  lognormal   = c(mu = 0, log_sigma = 0, 0)
)

.mixed_status_data <- function(n = 160) {
  set.seed(2531)
  x <- stats::rnorm(n)
  time <- stats::rexp(n, 0.3 * exp(0.4 * x)) + 0.05
  status <- sample(c(1, 0, -1), n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  entry <- ifelse(status %in% c(0, 1) & stats::runif(n) < 0.5,
                  time * stats::runif(n, 0.1, 0.8), 0)
  list(time = time, status = status, x = x, entry = entry)
}

test_that("status -1 rows ignore time_lower in every family", {
  k <- .mixed_status_data()
  expect_gt(sum(k$status == -1), 20)
  expect_gt(sum(k$entry > 0), 20)
  left <- k$status == -1
  tl_a <- k$entry                    # 0 on the left-censored rows
  tl_b <- k$entry
  tl_b[left] <- k$time[left] * 0.5   # arbitrary values on those rows only
  for (dist in .entry_dists) {
    fit_with <- function(tl) {
      suppressWarnings(hazard(
        time = k$time, status = k$status, time_lower = tl,
        time_upper = k$time, x = matrix(k$x, ncol = 1), dist = dist,
        theta = .entry_starts[[dist]], fit = TRUE,
        control = list(reltol = 1e-12)))
    }
    f_a <- fit_with(tl_a)
    f_b <- fit_with(tl_b)
    expect_true(isTRUE(f_a$fit$converged), info = dist)
    expect_equal(f_b$fit$objective, f_a$fit$objective, tolerance = 1e-10,
                 info = dist)
    expect_equal(coef(f_b), coef(f_a), tolerance = 1e-8, info = dist)
  }
})

test_that("the formula and vector interfaces give the same entry-time fit", {
  # Surv(start, stop, event) is the formula route to an entry time; the
  # vector route is time_lower =. Both must reach the same likelihood.
  k <- .mixed_status_data()
  keep <- k$status %in% c(0, 1)
  d <- data.frame(start = k$entry[keep], stop = k$time[keep],
                  event = k$status[keep], x = k$x[keep])
  expect_gt(sum(d$start > 0), 20)
  for (dist in .entry_dists) {
    f_formula <- suppressWarnings(hazard(
      survival::Surv(start, stop, event) ~ x, data = d, dist = dist,
      theta = .entry_starts[[dist]], fit = TRUE,
      control = list(reltol = 1e-12)))
    f_vector <- suppressWarnings(hazard(
      time = d$stop, status = d$event, time_lower = d$start,
      x = matrix(d$x, ncol = 1, dimnames = list(NULL, "x")), dist = dist,
      theta = .entry_starts[[dist]], fit = TRUE,
      control = list(reltol = 1e-12)))
    expect_equal(f_formula$fit$objective, f_vector$fit$objective,
                 tolerance = 1e-8, info = dist)
    expect_equal(unname(coef(f_formula)), unname(coef(f_vector)),
                 tolerance = 1e-6, info = dist)
    # And the entry time matters: dropping it changes the fit.
    f_none <- suppressWarnings(hazard(
      time = d$stop, status = d$event,
      x = matrix(d$x, ncol = 1, dimnames = list(NULL, "x")), dist = dist,
      theta = .entry_starts[[dist]], fit = TRUE,
      control = list(reltol = 1e-12)))
    expect_gt(abs(f_none$fit$objective - f_vector$fit$objective), 1,
              label = paste(dist, "entry-time effect on the objective"))
  }
})
