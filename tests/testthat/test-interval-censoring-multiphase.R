# Interval-censoring under the multiphase additive model (roadmap 7c).
#
# The multiphase likelihood has a working interval-/left-censored code path
# (`.hzr_logl_multiphase()`, status %in% c(-1, 2)) but no real-data or SAS
# parity fixture, so it was previously untested in isolation.  These are
# R-only Tier-1 self-consistency invariants: every assertion follows from the
# mathematics of the additive cumulative hazard
# `H(t) = sum_j mu_j(x) Phi_j(t)` and needs no external reference.  They mirror
# the single-distribution checks in test-interval-censoring-weibull.R, anchored
# to the independently-exercised `.hzr_multiphase_cumhaz()` evaluator.

# Shared two-phase spec (early CDF + constant), shapes fixed so the only free
# parameters are the two intercepts -- keeps the manual cumulative hazard exact.
mp_ic_phases <- function() {
  list(
    early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
}
mp_ic_cc    <- c(early = 0L, constant = 0L)
mp_ic_xl    <- list(early = NULL, constant = NULL)
# theta layout: [log_mu_early, log_t_half, nu, m, log_mu_constant]
mp_ic_theta <- c(-1.0, log(0.3), 1, 1, -1.5)

# Convenience: the model's cumulative hazard at the shared theta.
mp_ic_H <- function(t) {
  TemporalHazard:::.hzr_multiphase_cumhaz(
    t, mp_ic_theta, mp_ic_phases(), mp_ic_cc, mp_ic_xl)
}

mp_ic_ll <- function(time, status, time_lower = NULL, time_upper = NULL,
                     weights = NULL) {
  TemporalHazard:::.hzr_logl_multiphase(
    theta = mp_ic_theta, time = time, status = status,
    time_lower = time_lower, time_upper = time_upper,
    phases = mp_ic_phases(), covariate_counts = mp_ic_cc,
    x_list = mp_ic_xl, weights = weights)
}

# ---------------------------------------------------------------------------

test_that("multiphase log-likelihood handles all four censoring statuses", {
  # One row of each: exact event, right-censored, left-censored, interval.
  time       <- c(1.0, 2.0, 1.5, 0.8)
  status     <- c(1L, 0L, -1L, 2L)
  time_lower <- c(1.0, 2.0, 0.0, 0.6)
  time_upper <- c(1.0, 2.0, 1.5, 1.4)

  ll <- mp_ic_ll(time, status, time_lower, time_upper)
  expect_true(is.finite(ll))
})

test_that("multiphase interval contribution equals log(S(lower) - S(upper))", {
  l <- 0.9
  u <- 1.7
  ll <- mp_ic_ll(time = 1.0, status = 2L, time_lower = l, time_upper = u)

  # -H(l) + log(1 - exp(-(H(u) - H(l))))  ==  log(S(l) - S(u)).
  manual_form  <- -mp_ic_H(l) + hzr_log1mexp(mp_ic_H(u) - mp_ic_H(l))
  manual_surv  <- log(exp(-mp_ic_H(l)) - exp(-mp_ic_H(u)))

  expect_equal(unname(ll), unname(manual_form), tolerance = 1e-10)
  expect_equal(unname(ll), unname(manual_surv), tolerance = 1e-10)
})

test_that("multiphase left-censored contribution equals log(1 - exp(-H(u)))", {
  u  <- 1.3
  ll <- mp_ic_ll(time = u, status = -1L)
  manual <- hzr_log1mexp(mp_ic_H(u))
  expect_equal(unname(ll), unname(manual), tolerance = 1e-10)
})

test_that("multiphase right-censored contribution is -(H(stop) - H(start))", {
  # Plain right-censored (no entry time): -H(stop).
  ll_plain <- mp_ic_ll(time = c(0.7, 1.2, 2.0), status = c(0L, 0L, 0L))
  expect_equal(unname(ll_plain), -sum(mp_ic_H(c(0.7, 1.2, 2.0))),
               tolerance = 1e-10)

  # Left-truncated (counting-process start): -(H(stop) - H(start)).
  stop_t  <- c(1.5, 2.2)
  start_t <- c(0.3, 0.8)
  ll_lt <- mp_ic_ll(time = stop_t, status = c(0L, 0L), time_lower = start_t)
  expect_equal(unname(ll_lt),
               -sum(mp_ic_H(stop_t) - mp_ic_H(start_t)),
               tolerance = 1e-10)
})

test_that("multiphase interval bounds are validated (lower > upper -> -Inf)", {
  # Specifically -Inf (an infeasible interval), not merely non-finite: a
  # likelihood of +Inf would be a bug, so assert the exact value.
  expect_equal(
    mp_ic_ll(time = 1.0, status = 2L, time_lower = 2.0, time_upper = 1.0),
    -Inf)
})

test_that("interval-censored multiphase: integer weights match row duplication", {
  # The interval branch of the weighted sum must scale by the row weight,
  # i.e. a weight-w interval row equals w identical interval rows.  Covers the
  # interval x weights x multiphase combination (the weights tests use only
  # event/right-censored multiphase rows).
  set.seed(91)
  n <- 24
  status <- rep(2L, n)
  lo <- runif(n, 0.1, 1.5)
  up <- lo + runif(n, 0.1, 1.0)
  w  <- sample(1:3, n, replace = TRUE)
  idx <- rep(seq_len(n), times = w)

  ll_w   <- mp_ic_ll(time = lo, status = status,
                     time_lower = lo, time_upper = up, weights = w)
  ll_dup <- mp_ic_ll(time = lo[idx], status = status[idx],
                     time_lower = lo[idx], time_upper = up[idx],
                     weights = rep(1, sum(w)))
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

test_that("multiphase fit converges on an interval-censored sample", {
  skip_on_cran()
  set.seed(92)
  n <- 200
  true_t <- rexp(n, rate = 0.6) + 0.02
  # Interval-censor everyone into 0.5-wide windows around the true time.
  brk        <- ceiling(true_t / 0.5)
  time_lower <- (brk - 1L) * 0.5
  time_upper <- brk * 0.5
  status     <- rep(2L, n)

  fit <- suppressWarnings(hazard(
    time = time_lower, status = status,
    time_lower = time_lower, time_upper = time_upper,
    dist = "multiphase", phases = mp_ic_phases(),
    fit = TRUE, control = list(n_starts = 1L, maxit = 500L, conserve = FALSE)
  ))

  # The optimizer must actually converge on interval-censored multiphase data
  # (not merely return a logical flag).
  expect_true(isTRUE(fit$fit$converged))
  expect_true(is.finite(fit$fit$objective))
  # Two intercepts estimated (early + constant), shapes held fixed.
  expect_length(fit$fit$theta, 5L)
})

# ---------------------------------------------------------------------------
# Formula interface parity
#
# The production job that surfaced this could not express its three
# interval-censored records through the formula interface -- the fit returned
# the optimizer's failure sentinel -- and fell back to the vector interface,
# which in turn could not be bootstrapped.  The two interfaces must describe
# the same data.
# ---------------------------------------------------------------------------

test_that("Surv(type = 'interval') fits identically to the vector interface", {
  set.seed(20260817)
  n  <- 400
  df <- data.frame(x = rnorm(n))
  t  <- rexp(n, rate = 0.25 * exp(0.6 * df$x))

  # Statuses: right-censored at 6, exact event, or observed in an interval.
  ev <- ifelse(t >= 6, 0L, 1L)
  t  <- pmin(t, 6)
  iv <- which(ev == 1L)[1:40]        # 40 events known only within a window
  lo <- t
  hi <- rep(NA_real_, n)
  lo[iv] <- pmax(t[iv] - 0.4, 1e-6)
  hi[iv] <- t[iv] + 0.4
  ev[iv] <- 3L                        # Surv "interval" code

  ph  <- mp_ic_phases()
  ctl <- list(n_starts = 1, maxit = 2000)

  df$lo <- lo
  df$hi <- hi
  df$ev <- ev
  fit_formula <- hazard(survival::Surv(lo, hi, ev, type = "interval") ~ x,
                        data = df, dist = "multiphase", phases = ph,
                        fit = TRUE, control = ctl)

  # Same data spelled out in TemporalHazard's own coding.
  fit_vector <- hazard(
    time        = lo,
    status      = ifelse(ev == 3L, 2, ev),
    time_lower  = ifelse(ev == 3L, lo, 0),
    time_upper  = ifelse(ev == 3L, hi, lo),
    x           = as.matrix(df["x"]),
    dist = "multiphase", phases = ph, fit = TRUE, control = ctl)

  # Assert both fits produced a real number before comparing them -- two
  # absent objectives compare equal and would report parity that was never
  # tested.
  expect_true(is.finite(fit_formula$fit$objective))
  expect_true(is.finite(fit_vector$fit$objective))
  expect_length(coef(fit_formula), length(coef(fit_vector)))

  expect_equal(fit_formula$fit$objective, fit_vector$fit$objective,
               tolerance = 1e-8)
  expect_equal(unname(coef(fit_formula)), unname(coef(fit_vector)),
               tolerance = 1e-6)
})
