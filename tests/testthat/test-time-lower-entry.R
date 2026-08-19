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
  for (nm in names(.entry_variants(k))) {
    tl <- .entry_variants(k)[[nm]]
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
  for (nm in names(.entry_variants(k))) {
    tl <- .entry_variants(k)[[nm]]
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

test_that("entry at or after exit warns, and names the NULL default", {
  k <- .entry_case()
  fit_with <- function(tl) {
    suppressMessages(hazard(
      time = k$time, status = k$status, time_lower = tl,
      dist = "multiphase", phases = k$phases, theta = k$theta, fit = FALSE))
  }
  w <- capture_warnings(fit_with(k$time))
  expect_match(w, "counting-process entry time", all = FALSE)
  expect_match(w, "310 of 310", all = FALSE)
  # The remedy has to be in the message: NULL is not the same as `time`.
  expect_match(w, "leave 'time_lower' as NULL", all = FALSE)
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
