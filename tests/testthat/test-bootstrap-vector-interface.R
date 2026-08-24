# hzr_bootstrap() must resample a fit built with the vector interface.
#
# hazard() accepts either a formula plus `data`, or bare time/status vectors.
# Resampling `data` alone is enough for the formula interface; for the vector
# one the stored call holds `time = d$col` as an EXPRESSION, so each replicate
# re-evaluated it against the original data and returned the original fit --
# n_success = n_boot, n_failed = 0, no warning, and n_boot identical
# replicates. A bootstrap summary that looked complete and contained nothing.

fixture <- function() {
  e <- new.env()
  utils::data("cabgkul", package = "TemporalHazard", envir = e)
  e$cabgkul
}

phases_fixed <- function() {
  list(
    early = hzr_phase("cdf", t_half = 0.19, nu = 1.4, m = 1,
                      fixed = c("t_half", "nu", "m")),
    late  = hzr_phase("g3", tau = 1, gamma = 1, alpha = 1, eta = 1.7,
                      fixed = c("tau", "gamma", "alpha", "eta"))
  )
}

sd_log_mu <- function(fit, n_boot = 20, seed = 1) {
  b <- suppressWarnings(hzr_bootstrap(fit, n_boot = n_boot, seed = seed))
  stats::sd(b$replicates$estimate[b$replicates$parameter == "early.log_mu"])
}

test_that("a vector-interface fit actually resamples", {
  d <- fixture()
  fit <- suppressWarnings(hazard(
    time = d$int_dead, status = as.integer(d$dead), data = d,
    dist = "multiphase", phases = phases_fixed(),
    fit = TRUE, control = list(n_starts = 1, maxit = 200)))

  # The regression: this was exactly 0.
  expect_gt(sd_log_mu(fit), 0)
})

test_that("vector and formula interfaces bootstrap identically", {
  # The strong form, and it has to compare the REPLICATES, not a summary of
  # them: two different sets of replicate estimates can share an SD, so
  # comparing sd() alone would pass on outputs that differ.
  d <- fixture()
  fv <- suppressWarnings(hazard(
    time = d$int_dead, status = as.integer(d$dead), data = d,
    dist = "multiphase", phases = phases_fixed(),
    fit = TRUE, control = list(n_starts = 1, maxit = 200)))
  ff <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d,
    dist = "multiphase", phases = phases_fixed(),
    fit = TRUE, control = list(n_starts = 1, maxit = 200)))

  bv <- suppressWarnings(hzr_bootstrap(fv, n_boot = 20, seed = 1))
  bf <- suppressWarnings(hzr_bootstrap(ff, n_boot = 20, seed = 1))

  # Same model, same data, same seed: every replicate estimate must match.
  expect_equal(bv$replicates, bf$replicates, tolerance = 1e-8)
  expect_equal(bv$n_success, bf$n_success)
  # And they must not be trivially equal by both being constant.
  expect_gt(stats::sd(bv$replicates$estimate[
    bv$replicates$parameter == "early.log_mu"]), 0)
})

test_that("time_lower and time_upper are resampled alongside time", {
  # Interval-censored rows carry their bounds in separate vectors; resampling
  # `time` but not the bounds would silently pair row i's time with row j's
  # interval.
  d <- fixture()
  status <- as.integer(d$dead)
  idx <- which(status == 0L)[1:20]
  status[idx] <- 2L
  lower <- pmax(d$int_dead - 0.1, 0)

  fit <- suppressWarnings(hazard(
    time = d$int_dead, status = status,
    time_lower = lower, time_upper = d$int_dead, data = d,
    dist = "multiphase", phases = phases_fixed(),
    fit = TRUE, control = list(n_starts = 1, maxit = 200)))

  expect_gt(sd_log_mu(fit), 0)
})

test_that("a partially-stored vector-interface fit is refused, not half-rewired", {
  # Rewiring some vector arguments and not others is worse than rewiring none:
  # the rewired ones follow the resample while the rest evaluate against the
  # original data, pairing row i's time with row j's status. Silent corruption
  # producing plausible numbers.
  d <- fixture()
  fit <- suppressWarnings(hazard(
    time = d$int_dead, status = as.integer(d$dead), data = d,
    dist = "multiphase", phases = phases_fixed(),
    fit = TRUE, control = list(n_starts = 1, maxit = 200)))

  # Simulate an object fitted by an older version that did not store `status`.
  fit$data$status <- NULL
  expect_error(suppressWarnings(hzr_bootstrap(fit, n_boot = 5, seed = 1)),
               "status")
})
