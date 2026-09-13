# `time_lower` as a counting-process entry time -------------------------------
#
# For status 0/1 rows the log-likelihood subtracts H(time_lower) from every
# row, unconditionally, whenever `time_lower` is supplied. The gradient and
# the Hessian used to define "start" differently -- they filtered on
# `time_lower < time` -- so a row entering the risk set at its own event or
# censoring time was differentiated as though it had no entry time at all,
# while its weight was still applied.
#
# The derivative was therefore taken of a different function from the one
# being evaluated. Measured on `avc` at fixed theta, the analytic gradient was
# out by 382 with time_lower = time, and by 126 on data where only SOME rows
# entered at their exit time -- which is ordinary left-truncated data, not a
# pathological input. Issue #136.
#
# The invariant is the guard: analytic derivatives must agree with numDeriv
# for EVERY arrangement of entry times, not merely the ones the old filter
# happened to admit.

.entry_case <- function() {
  data("avc", package = "TemporalHazard", envir = environment())
  ph <- list(early = hzr_phase("cdf", t_half = 0.19, nu = 1.44, m = 1,
                               fixed = "m"),
             constant = hzr_phase("constant"))
  list(
    time = avc$int_dead, status = avc$dead,
    theta = c(log(0.35), log(0.19), 1.44, 1, log(0.03)),
    phases = ph,
    weights = rep(1, nrow(avc)),
    counts = c(early = 0L, constant = 0L),
    x_list = list(early = NULL, constant = NULL)
  )
}

.entry_variants <- function(k) {
  mixed <- k$time
  half <- seq_len(floor(length(mixed) / 2))
  mixed[half] <- pmax(k$time[half] - 0.5, 0)     # some enter early, rest at exit
  list(
    "no entry times"      = NULL,
    "entry at zero"       = rep(0, length(k$time)),
    "entry at exit"       = k$time,
    "genuine truncation"  = pmax(k$time - 0.5, 0),
    "mixed entry times"   = mixed
  )
}

test_that("the multiphase gradient matches numDeriv for every entry-time layout", {
  skip_if_not_installed("numDeriv")
  k <- .entry_case()
  variants <- .entry_variants(k)
  for (nm in names(variants)) {
    tl <- variants[[nm]]
    f <- function(p) {
      .hzr_logl_multiphase(p, time = k$time, status = k$status, time_lower = tl,
                           time_upper = NULL, x = NULL, weights = k$weights,
                           phases = k$phases, covariate_counts = k$counts,
                           x_list = k$x_list)
    }
    analytic <- .hzr_gradient_multiphase(
      k$theta, time = k$time, status = k$status, time_lower = tl,
      time_upper = NULL, x = NULL, weights = k$weights, phases = k$phases,
      covariate_counts = k$counts, x_list = k$x_list)
    numeric_ <- numDeriv::grad(f, k$theta)
    expect_equal(analytic, numeric_, tolerance = 1e-5, info = nm)
  }
})

test_that("the multiphase Hessian matches numDeriv for every entry-time layout", {
  skip_if_not_installed("numDeriv")
  k <- .entry_case()
  variants <- .entry_variants(k)
  for (nm in names(variants)) {
    tl <- variants[[nm]]
    nll <- function(p) {
      -.hzr_logl_multiphase(p, time = k$time, status = k$status, time_lower = tl,
                            time_upper = NULL, x = NULL, weights = k$weights,
                            phases = k$phases, covariate_counts = k$counts,
                            x_list = k$x_list)
    }
    analytic <- .hzr_hessian_multiphase(
      k$theta, time = k$time, status = k$status, time_lower = tl,
      time_upper = NULL, x = NULL, weights = k$weights, phases = k$phases,
      covariate_counts = k$counts, x_list = k$x_list)
    expect_false(is.null(analytic), info = nm)
    expect_equal(unname(analytic), numDeriv::hessian(nll, k$theta),
                 tolerance = 1e-4, info = nm)
  }
})

# On status 0/1 rows every family reads `time_lower` as the entry time when
# 0 < time_lower < time. time_lower = 0 and time_lower = time both mean "no
# entry time"; the second is the mixed-interval layout, where only status-2
# rows carry a real lower bound. An entry AFTER exit (time_lower > time) is
# an error. SAS HAZARD refuses STIME >= TIME (setcoe_obs_loop.c, SETCOE960),
# but its STIME is only ever an entry time, while R's time_lower doubles as
# the interval bound, so here only the strict case is refused (#253).
.bad_entry_dists <- c("weibull", "exponential", "loglogistic", "lognormal",
                      "multiphase")

.bad_entry_call <- function(k, dist, tl, use_data = FALSE) {
  ph <- if (dist == "multiphase") k$phases else NULL
  th <- if (dist == "multiphase") k$theta else NULL
  if (use_data) {
    d <- data.frame(tt = k$time, st = k$status, lo = tl)
    hazard(time = tt, status = st, time_lower = lo, data = d,
           dist = dist, phases = ph, theta = th, fit = FALSE)
  } else {
    hazard(time = k$time, status = k$status, time_lower = tl,
           dist = dist, phases = ph, theta = th, fit = FALSE)
  }
}

test_that("entry after exit is an error for every dist", {
  k <- .entry_case()
  n <- length(k$time)
  one_bad <- pmax(k$time - 0.5, 0)
  j <- which(k$time > 0)[1L]
  one_bad[j] <- k$time[j] + 1          # a single row entering after it leaves
  for (dist in .bad_entry_dists) {
    for (use_data in c(FALSE, TRUE)) {
      info <- paste(dist, if (use_data) "data-masked" else "plain vectors")
      # Every row entering one unit after its own exit time.
      err <- tryCatch(suppressMessages(
        .bad_entry_call(k, dist, k$time + 1, use_data)),
        error = function(e) conditionMessage(e))
      expect_type(err, "character")
      expect_match(err, "counting-process entry time", fixed = TRUE,
                   info = info)
      # Derived, not hard-coded: the message counts every offending row.
      expect_match(err, paste0(n, " of ", n, " row(s)"), fixed = TRUE,
                   info = info)
      expect_match(err, "'time_lower' > 'time'", fixed = TRUE, info = info)
      expect_match(err, "cannot enter the risk set after it leaves",
                   fixed = TRUE, info = info)
      # The remedy has to be in the message, `time` included.
      expect_match(err, "set it to 0 or to 'time'", fixed = TRUE,
                   info = info)
      # One offending row among genuine entry times is still refused, and
      # the count says one.
      expect_error(suppressMessages(
        .bad_entry_call(k, dist, one_bad, use_data)),
        paste0("1 of ", n, " row(s)"), fixed = TRUE, info = info)
    }
  }
})

test_that("the formula interface refuses start >= stop through survival", {
  # Surv(start, stop, event) is the formula route to an entry time, and
  # survival refuses start >= stop itself, start == stop included: it sets
  # the entry to NA, which hazard() then rejects as a non-finite
  # 'time_lower'. So the formula path refuses a zero-length interval that
  # the vector path accepts as "no entry time"; the mixed-interval layout is
  # a vector-interface idiom. Pin both halves, so a change in survival that
  # lets such a row through is noticed here.
  expect_warning(s <- survival::Surv(c(0, 2), c(1, 2), c(1, 0)),
                 "Stop time must be > start time")
  expect_true(is.na(unclass(s)[2L, 1L]))
  d <- data.frame(a = c(0, 0.5, 2), b = c(1, 2, 2), e = c(1, 0, 1))
  expect_error(suppressWarnings(hazard(
    survival::Surv(a, b, e) ~ 1, data = d, dist = "weibull", fit = FALSE)),
    "'time_lower' must be a numeric vector of finite non-negative values",
    fixed = TRUE)
})

test_that("time_lower = 0 is 'no entry time' and matches time_lower = NULL", {
  set.seed(253)
  n <- 120
  tt <- stats::rexp(n, 0.4) + 0.05
  st <- stats::rbinom(n, 1, 0.7)
  starts <- list(weibull = c(0.5, 1), exponential = c(log_rate = 0),
                 loglogistic = c(log_alpha = 0, log_beta = 0),
                 lognormal = c(mu = 0, log_sigma = 0))
  for (dist in names(starts)) {
    f_null <- suppressWarnings(hazard(time = tt, status = st, dist = dist,
                                      theta = starts[[dist]], fit = TRUE))
    f_zero <- expect_no_error(suppressWarnings(hazard(
      time = tt, status = st, time_lower = rep(0, n), dist = dist,
      theta = starts[[dist]], fit = TRUE)))
    expect_equal(coef(f_zero), coef(f_null), tolerance = 1e-8, info = dist)
    expect_equal(f_zero$fit$objective, f_null$fit$objective,
                 tolerance = 1e-8, info = dist)
  }
  k <- .entry_case()
  mp <- function(tl) {
    suppressWarnings(suppressMessages(hazard(
      time = k$time, status = k$status, time_lower = tl, dist = "multiphase",
      phases = k$phases, theta = k$theta, fit = TRUE)))
  }
  f_null <- mp(NULL)
  f_zero <- expect_no_error(mp(rep(0, length(k$time))))
  expect_equal(coef(f_zero), coef(f_null), tolerance = 1e-8)
})

test_that("time_lower = time is 'no entry time' and matches NULL (#253)", {
  # The mixed-interval layout: status 0/1 rows carry their own time. Before
  # #253 multiphase read it as an entry at exit and returned +47915.76 with
  # converged = TRUE on avc; the three families below ignored time_lower.
  set.seed(253)
  n <- 120
  tt <- stats::rexp(n, 0.4) + 0.05
  st <- stats::rbinom(n, 1, 0.7)
  starts <- list(weibull = c(0.5, 1), exponential = c(log_rate = 0),
                 loglogistic = c(log_alpha = 0, log_beta = 0),
                 lognormal = c(mu = 0, log_sigma = 0))
  for (dist in names(starts)) {
    f_null <- suppressWarnings(hazard(time = tt, status = st, dist = dist,
                                      theta = starts[[dist]], fit = TRUE))
    f_same <- expect_no_error(suppressWarnings(hazard(
      time = tt, status = st, time_lower = tt, dist = dist,
      theta = starts[[dist]], fit = TRUE)))
    expect_equal(coef(f_same), coef(f_null), tolerance = 1e-8, info = dist)
    expect_equal(f_same$fit$objective, f_null$fit$objective,
                 tolerance = 1e-8, info = dist)
  }
  k <- .entry_case()
  mp <- function(tl) {
    suppressWarnings(suppressMessages(hazard(
      time = k$time, status = k$status, time_lower = tl, dist = "multiphase",
      phases = k$phases, theta = k$theta, fit = TRUE)))
  }
  f_null <- mp(NULL)
  f_same <- expect_no_error(mp(k$time))
  expect_equal(coef(f_same), coef(f_null), tolerance = 1e-8)
  expect_equal(f_same$fit$objective, f_null$fit$objective, tolerance = 1e-8)
  expect_lt(f_same$fit$objective, 0)
})

test_that("zero-length epochs among genuine entries are an error (#253)", {
  # Counting-process data where one row enters and exits at the same time
  # (a second event at the previous row's time; hzr_repeated_events() emits
  # these). Reading time_lower == time as "no entry" there would charge the
  # row its full H(0, t], so hazard() refuses the mix, as SAS refuses
  # STIME >= TIME. Found by r-reviewer; rows 3 of 4 is the offender.
  tt <- c(1, 2, 2, 3)
  st <- c(1, 1, 1, 0)
  tl <- c(0, 1, 2, 2)
  for (dist in .bad_entry_dists[.bad_entry_dists != "multiphase"]) {
    err <- tryCatch(hazard(time = tt, status = st, time_lower = tl,
                           dist = dist, fit = FALSE),
                    error = function(e) conditionMessage(e))
    expect_type(err, "character")
    expect_match(err, "1 of 4 row(s)", fixed = TRUE, info = dist)
    expect_match(err, "zero-length", fixed = TRUE, info = dist)
  }
  # Genuine entries only: fits.
  expect_no_error(hazard(time = tt, status = st,
                         time_lower = c(0, 1, 1.5, 2), dist = "weibull",
                         fit = FALSE))
  # The mixed-interval layout (no genuine entry on any status 0/1 row)
  # still reads time_lower == time as no entry.
  expect_no_error(hazard(time = c(1, 2, 3, 4), status = c(1, 0, 2, 1),
                         time_lower = c(1, 2, 2.5, 4),
                         time_upper = c(1, 2, 3.5, 4), dist = "weibull",
                         fit = FALSE))
})

test_that("a status-0 row at time 0 with time_lower 0 is not an error", {
  # time_lower == time == 0 is the SAS "no entry time" row, not an entry at
  # exit: STIME = 0 skips the check there, and so must hazard().
  tt <- c(0, 0.4, 1.1, 2.3, 3.0)
  st <- c(0, 1, 0, 1, 1)
  for (dist in .bad_entry_dists[.bad_entry_dists != "multiphase"]) {
    w <- expect_no_error(capture_warnings(hazard(
      time = tt, status = st, time_lower = rep(0, 5), dist = dist,
      fit = FALSE)))
    expect_false(any(grepl("entry time", w, fixed = TRUE)), info = dist)
  }
  # The same row through the formula interface: Surv(type = "interval")
  # gives status 0/1 rows an entry time of 0.
  d <- data.frame(t1 = tt, t2 = tt, e = c(0, 1, 0, 1, 1))
  w <- expect_no_error(capture_warnings(hazard(
    survival::Surv(t1, t2, e, type = "interval") ~ 1, data = d,
    dist = "weibull", fit = FALSE)))
  expect_false(any(grepl("entry time", w, fixed = TRUE)))
})

test_that("genuine left truncation does not warn", {
  k <- .entry_case()
  # Every row enters strictly before it leaves: an ordinary counting-process
  # dataset, which must stay silent or the warning is noise.
  w <- capture_warnings(suppressMessages(hazard(
    time = k$time, status = k$status, time_lower = pmax(k$time - 0.5, 0),
    dist = "multiphase", phases = k$phases, theta = k$theta, fit = FALSE)))
  expect_false(any(grepl("counting-process entry time", w)))
})
