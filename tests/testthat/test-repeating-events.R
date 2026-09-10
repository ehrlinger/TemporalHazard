# tests/testthat/test-repeating-events.R
# Tests for Surv(start, stop, event) counting-process / epoch-decomposition
# support added in 0.9.4.
#
# The defining invariant of epoch decomposition is:
#
#   L([0, T], event=d)  ==  L([0, t1], 0) * L([t1, T], d)
#
# i.e. splitting a single observation into two contiguous epochs on a time
# s in (0, T) produces the same likelihood as the unsplit observation.
# This holds because each epoch contributes H(stop) - H(start), and the
# telescoping sum recovers H(T).
#
# These tests exercise that invariant on both Weibull and multiphase.
# They also verify the trivial case that Surv(0, t, event) matches
# Surv(t, event) exactly.

make_toy_rc <- function(seed = 91) {
  set.seed(seed)
  n <- 40
  x <- rnorm(n)
  u <- runif(n)
  t <- (-log(u) / exp(0.2 * x)) ^ (1 / 1.1)
  status <- as.integer(t < quantile(t, 0.6))
  data.frame(time = t, status = status, x = x)
}

# ---------------------------------------------------------------------------
# Trivial case: start = 0 matches right-censored Surv
# ---------------------------------------------------------------------------

test_that("Surv(0, t, d) fits identically to Surv(t, d) on Weibull", {
  skip_on_cran()
  df <- make_toy_rc()
  df$start <- 0

  # theta is required: without it hazard() used to skip the optimizer, and
  # this test compared NULL coefficients with NULL coefficients.
  fit_rc <- hazard(
    survival::Surv(time, status) ~ x,
    data = df, dist = "weibull", fit = TRUE, theta = c(0.5, 1, 0)
  )
  fit_cp <- hazard(
    survival::Surv(start, time, status) ~ x,
    data = df, dist = "weibull", fit = TRUE, theta = c(0.5, 1, 0)
  )
  expect_true(fit_rc$fit$converged)                 # guard: both fits ran
  expect_true(fit_cp$fit$converged)

  expect_equal(coef(fit_rc), coef(fit_cp), tolerance = 1e-3)
  expect_equal(fit_rc$fit$objective, fit_cp$fit$objective, tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# Epoch split invariance: splitting at an interior time preserves LL
# ---------------------------------------------------------------------------

split_epochs <- function(df, split_frac = 0.5) {
  # Each row [time_i, status_i] becomes two rows:
  #   [0,  s_i, 0]            -- left epoch, never an event
  #   [s_i, time_i, status_i] -- right epoch, preserves original status
  # where s_i = split_frac * time_i.
  s <- split_frac * df$time
  pre <- data.frame(
    start = 0, stop = s, status = 0L,
    x = df$x
  )
  post <- data.frame(
    start = s, stop = df$time, status = df$status,
    x = df$x
  )
  rbind(pre, post)
}

# The two split-invariance tests below are the canonical check that
# `H(stop) - H(start)` is actually applied to counting-process rows.
# Wired up in v0.9.6 (Phase 4f).

test_that("splitting each row into two epochs preserves Weibull log-lik", {
  df <- make_toy_rc()
  df_split <- split_epochs(df, split_frac = 0.3)

  theta <- c(mu = 0.8, nu = 1.1, beta = 0.2)

  ll_orig <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = df$time, status = df$status,
    x = matrix(df$x, ncol = 1)
  )

  ll_split <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta,
    time = df_split$stop, status = df_split$status,
    time_lower = df_split$start,
    x = matrix(df_split$x, ncol = 1)
  )

  expect_equal(ll_orig, ll_split, tolerance = 1e-8)
})

test_that("splitting each row into two epochs preserves multiphase log-lik", {
  df <- make_toy_rc()
  df_split <- split_epochs(df, split_frac = 0.4)

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  # Theta layout: [log_mu_e, log_thalf, nu, m, beta_e, log_mu_c, beta_c].
  # Match the phase spec's nu = m = 1 so the decomposition case dispatch
  # hits Case 1 (m > 0, nu > 0); otherwise nu = m = 0 lands in a gap.
  theta <- c(-4, log(0.3), 1, 1, -3, 0.1, 0.05)

  ll_orig <- TemporalHazard:::.hzr_logl_multiphase(
    theta = theta, time = df$time, status = df$status,
    phases = phases,
    covariate_counts = c(early = 1L, constant = 1L),
    x_list = list(early = matrix(df$x, ncol = 1),
                  constant = matrix(df$x, ncol = 1))
  )

  ll_split <- TemporalHazard:::.hzr_logl_multiphase(
    theta = theta,
    time = df_split$stop, status = df_split$status,
    time_lower = df_split$start,
    phases = phases,
    covariate_counts = c(early = 1L, constant = 1L),
    x_list = list(early = matrix(df_split$x, ncol = 1),
                  constant = matrix(df_split$x, ncol = 1))
  )

  expect_equal(ll_orig, ll_split, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Parser plumbing: Surv(start, stop, event) routes via the formula interface
# ---------------------------------------------------------------------------

test_that("Surv(start, stop, event) routes start -> time_lower in formula parser", {
  df <- data.frame(
    start = c(0, 0.5, 1.0),
    stop  = c(1.0, 1.5, 2.0),
    event = c(1, 0, 1)
  )
  parsed <- TemporalHazard:::.hzr_parse_formula(
    survival::Surv(start, stop, event) ~ 1,
    data = df
  )
  expect_equal(parsed$time_lower, df$start)
  expect_equal(parsed$time, df$stop)
  expect_equal(parsed$status, df$event)
})

# ---------------------------------------------------------------------------
# End-to-end fit: counting-process Surv converges and the split dataset
# yields the same MLE as the unsplit.
# ---------------------------------------------------------------------------

test_that("Weibull fit on split-epoch data matches unsplit fit", {
  df <- make_toy_rc()
  df_split <- split_epochs(df, split_frac = 0.35)

  # theta is required: without it hazard() used to skip the optimizer, and
  # this test compared NULL coefficients with NULL coefficients.
  fit_orig <- hazard(
    survival::Surv(time, status) ~ x,
    data = df, dist = "weibull", fit = TRUE, theta = c(0.5, 1, 0)
  )
  fit_split <- hazard(
    survival::Surv(start, stop, status) ~ x,
    data = df_split, dist = "weibull", fit = TRUE, theta = c(0.5, 1, 0)
  )
  expect_true(fit_orig$fit$converged)               # guard: both fits ran
  expect_true(fit_split$fit$converged)

  expect_equal(coef(fit_orig), coef(fit_split), tolerance = 1e-2)
  expect_equal(fit_orig$fit$objective, fit_split$fit$objective, tolerance = 1e-3)
})

# ---------------------------------------------------------------------------
# Counting-process data with nonzero starts is now accepted (v0.9.6)
# ---------------------------------------------------------------------------
# The v0.9.5 narrowing rejected `Surv(start, stop, event)` with any
# `start > 0`; v0.9.6 wires `H(stop) - H(start)` into the Weibull and
# multiphase likelihoods, so the guard is gone.

test_that("hazard() accepts counting-process Surv when every start is 0", {
  df <- data.frame(
    start = c(0, 0, 0, 0),
    stop  = c(1.5, 2.0, 3.0, 3.5),
    event = c(1, 0, 1, 0),
    x     = c(0.1, 0.2, 0.3, 0.4)
  )
  # Should not error -- start = 0 is equivalent to right-censored Surv.
  fit <- hazard(
    survival::Surv(start, stop, event) ~ x,
    data = df, dist = "weibull", fit = FALSE
  )
  expect_s3_class(fit, "hazard")
})
