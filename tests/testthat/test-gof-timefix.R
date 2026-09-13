# #286: survfit()'s timefix merges near-tied exit times into one
# Kaplan-Meier time. hzr_gof() matched each subject's raw exit time to the
# grid, so a subject merged onto a neighbour's time matched no point on the
# default grid and dropped out of both tallies, with no warning. The first
# five tests failed on 2daa3e1 (main after #285) apart from the custom-grid
# guard; main before #285 counted every event. The n_risk check, the gap
# sweep and the grid a few ulps off come from the #263/#264 session's
# parallel fix (09fb34a).

.near_tie_data <- function(t0, gap, n = 40) {
  set.seed(1)
  data.frame(
    time = c(t0, t0 * (1 + gap), round(stats::rexp(n - 2, 0.5) + 0.05, 2)),
    status = c(1, 1, stats::rbinom(n - 2, 1, 0.6))
  )
}

.gof_tf_avc <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  d$entry <- ifelse(seq_len(nrow(d)) %% 2 == 0, pmin(0.1, d$int_dead / 2), 0)
  d
}

test_that("near-tied exits merged by survfit stay in both tallies", {
  for (case in list(c(t0 = 1, gap = 1e-12), c(t0 = 200, gap = 1e-12))) {
    d <- .near_tie_data(case[["t0"]], case[["gap"]])
    km <- survival::survfit(survival::Surv(d$time, d$status) ~ 1)
    # The scenario: survfit merged the pair into one Kaplan-Meier time.
    expect_equal(sum(abs(km$time - case[["t0"]]) < 1e-6 * case[["t0"]]), 1)
    fit <- suppressWarnings(hazard(time = d$time, status = d$status,
                                   dist = "weibull", theta = c(0.5, 1),
                                   fit = TRUE,
                                   control = list(reltol = 1e-12)))
    g <- hzr_gof(fit)
    s <- attr(g, "summary")
    info <- paste("t0 =", case[["t0"]])
    expect_equal(s$total_observed, sum(d$status), info = info)
    expect_equal(sum(g$n_event), sum(d$status), info = info)
    # The risk set agrees with survfit()'s at its own (merged) times.
    expect_equal(g$n_risk, km$n.risk, info = info)
    # A Weibull fit conserves events, so E/O is 1 when every subject counts.
    expect_equal(s$total_expected / s$total_observed, 1, tolerance = 1e-5,
                 info = info)
  }
})

test_that("entry-time fits count every event (#286's case)", {
  d <- .gof_tf_avc()
  # The case that lost 13 of 68 events: survfit() merges some of these exit
  # times onto neighbours.
  km <- survival::survfit(survival::Surv(d$entry, d$int_dead, d$dead) ~ 1)
  expect_gt(sum(!(unique(d$int_dead) %in% km$time)), 0)
  fit <- suppressWarnings(hazard(
    time = d$int_dead, status = d$dead, time_lower = d$entry,
    x = cbind(age = d$age), dist = "weibull",
    theta = c(mu = 0.01, nu = 0.5, beta_age = 0), fit = TRUE))
  g <- hzr_gof(fit)
  s <- attr(g, "summary")
  expect_equal(s$total_observed, sum(d$dead))
  expect_equal(g$cum_observed[nrow(g)], sum(d$dead))
  expect_equal(s$total_expected / s$total_observed, 1, tolerance = 1e-5)
})

test_that("a multiphase fit with entry times and conservation counts every event", {
  skip_on_cran()
  d <- .gof_tf_avc()
  fit <- suppressWarnings(suppressMessages(hazard(
    time = d$int_dead, status = d$dead, time_lower = d$entry,
    x = cbind(age = d$age), dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = list(conserve = TRUE))))
  expect_true(isTRUE(fit$spec$control$conserve_applied))
  s <- attr(hzr_gof(fit), "summary")
  expect_equal(s$total_observed, sum(d$dead))
  expect_lt(abs(s$total_expected / s$total_observed - 1), 1e-6)
})

test_that("a custom grid of the raw exit times still counts every event", {
  d <- .near_tie_data(1, 1e-12)
  fit <- suppressWarnings(hazard(time = d$time, status = d$status,
                                 dist = "weibull", theta = c(0.5, 1),
                                 fit = TRUE))
  g <- hzr_gof(fit, time_grid = sort(unique(d$time)))
  expect_equal(attr(g, "summary")$total_observed, sum(d$status))
  # A grid holding the raw 1 + 1e-12 but not the merged time 1: that
  # subject is still found by its raw time. The subject exiting at exactly 1
  # has no grid point and is left out, as documented for custom grids.
  grid <- sort(unique(d$time[d$time != 1]))
  g <- hzr_gof(fit, time_grid = grid)
  expect_equal(attr(g, "summary")$total_observed, sum(d$status) - 1)
  # A grid holding only the merged time 1 counts both near-tied subjects,
  # once each, where survfit put them; n_event agrees. Matching the custom
  # grid on raw times alone would count 1 of these 2 events.
  g <- hzr_gof(fit, time_grid = 1)
  expect_equal(attr(g, "summary")$total_observed, 2)
  expect_equal(g$n_event, 2)
})

test_that("a subject the default grid cannot place is reported, not dropped silently", {
  # An exit at time 0 in entry-time data is Surv(0, 0), which survfit turns
  # into NA, so that subject has no Kaplan-Meier time.
  set.seed(2861)
  n <- 60
  stop_t <- stats::rexp(n, 0.4) + 0.1
  start_t <- ifelse(stats::runif(n) < 0.5, stop_t * stats::runif(n, 0.1, 0.8), 0)
  stop_t[1] <- 0
  start_t[1] <- 0
  status <- stats::rbinom(n, 1, 0.7)
  fit <- suppressWarnings(hazard(time = stop_t, status = status,
                                 time_lower = start_t, dist = "exponential",
                                 theta = c(log_rate = log(0.3)), fit = TRUE))
  g <- suppressWarnings(hzr_gof(fit))
  # Surv() warns about the NA it creates for that row; only hzr_gof()'s own
  # warning is under test.
  expect_warning(
    withCallingHandlers(hzr_gof(fit), warning = function(w) {
      if (grepl("NA created", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }),
    "1 of 60 subjects did not match")
  expect_equal(attr(g, "summary")$total_observed, sum(status[-1]))
})

test_that("no gap around survfit's merge threshold loses an event", {
  # Relative gaps from well inside to well outside survfit()'s merge
  # tolerance (about 1.5e-8 relative), at a small and a large time.
  for (t0 in c(1, 150)) {
    for (gap in 10^seq(-14, -6, by = 0.5)) {
      t <- c(t0, t0 * (1 + gap), t0 + seq_len(20))
      fit <- suppressWarnings(hazard(time = t, status = rep(1, 22),
                                     dist = "weibull",
                                     theta = c(mu = 1 / t0, nu = 1),
                                     fit = TRUE))
      expect_equal(attr(hzr_gof(fit), "summary")$total_observed, 22,
                   label = sprintf("t0 = %g, relative gap = %g", t0, gap))
    }
  }
})

test_that("a supplied grid a few ulps off the exit times counts them all", {
  d <- .gof_tf_avc()
  # The fit's Hessian is ill-conditioned this far from 0; only the counting
  # matters here.
  fit <- suppressWarnings(hazard(time = d$int_dead + 1e4, status = d$dead,
                                 dist = "weibull",
                                 theta = c(mu = 1e-4, nu = 1), fit = TRUE))
  exits <- sort(unique(d$int_dead + 1e4))
  grid <- exits * (1 + 2 * .Machine$double.eps)   # about two ulps above
  expect_true(all(grid != exits))
  expect_equal(attr(hzr_gof(fit, time_grid = grid), "summary")$total_observed,
               sum(d$dead))
})
