# tests/testthat/test-weights.R
# Tests for the `weights` argument added in 0.9.4.
#
# The key invariant we check is that Fisher weighting is equivalent to
# row duplication for integer weights:
#
#   logL(theta; weights = w,  data = D)
#     ==  logL(theta; weights = 1,  data = expand(D, w))
#
# where expand(D, w) duplicates row i exactly w_i times.  This verifies
# that weights propagate through every likelihood term (event,
# right-censored, interval, left).  Running the same invariant through
# `hazard(fit = TRUE)` additionally exercises the optimizer, gradient,
# and Hessian paths because any mis-weighting there would show up as a
# different MLE or vcov between the two fits.

make_toy <- function(seed = 17) {
  set.seed(seed)
  n <- 30
  x <- rnorm(n)
  # Weibull-ish times
  u <- runif(n)
  t <- (-log(u) / exp(0.3 * x))^(1 / 1.2)
  # Mix of events, right-censored, and interval-censored
  status <- ifelse(t < 0.6, 1, 0)
  data.frame(time = t, status = status, x = x)
}

expand_rows <- function(df, weights) {
  df[rep(seq_len(nrow(df)), times = weights), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------------

test_that("weights must be numeric and the right length", {
  df <- make_toy()
  expect_error(
    hazard(survival::Surv(time, status) ~ x, data = df,
           dist = "weibull", weights = "a"),
    "numeric vector"
  )
  expect_error(
    hazard(survival::Surv(time, status) ~ x, data = df,
           dist = "weibull", weights = rep(1, nrow(df) - 1)),
    "numeric vector"
  )
})

test_that("weights must be non-negative and finite", {
  df <- make_toy()
  expect_error(
    hazard(survival::Surv(time, status) ~ x, data = df,
           dist = "weibull", weights = c(-1, rep(1, nrow(df) - 1))),
    "non-negative"
  )
  expect_error(
    hazard(survival::Surv(time, status) ~ x, data = df,
           dist = "weibull", weights = c(Inf, rep(1, nrow(df) - 1))),
    "non-negative"
  )
})

# ---------------------------------------------------------------------------
# Weighted log-likelihood equals duplicated-row log-likelihood (Weibull)
# ---------------------------------------------------------------------------

test_that("integer weights match row duplication at an arbitrary theta (Weibull)", {
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  df_exp <- expand_rows(df, w)
  theta <- c(mu = 0.8, nu = 1.2, beta = 0.25)

  ll_weighted <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta,
    time = df$time, status = df$status,
    x = matrix(df$x, ncol = 1),
    weights = w
  )
  ll_expanded <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta,
    time = df_exp$time, status = df_exp$status,
    x = matrix(df_exp$x, ncol = 1),
    weights = rep(1, nrow(df_exp))
  )

  expect_equal(ll_weighted, ll_expanded, tolerance = 1e-10)
})

test_that("integer weights match row duplication for all censoring flavours", {
  set.seed(3)
  n <- 20
  x <- rnorm(n)
  t <- runif(n, 0.1, 2)
  # Mix: 10 events, 5 right-censored, 5 left-censored
  status <- c(rep(1, 10), rep(0, 5), rep(-1, 5))
  # The parser's convention for left-censored data is time_upper = time for
  # every row (non-left rows ignore it). Mirror that here so the LL's
  # `any(upper <= 0)` guard doesn't trip on NAs.
  time_upper <- t
  w <- sample(1:3, n, replace = TRUE)

  theta <- c(mu = 0.7, nu = 1.3, beta = 0.1)
  ll_w <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = t, status = status,
    time_upper = time_upper, x = matrix(x, ncol = 1), weights = w
  )

  idx <- rep(seq_len(n), times = w)
  ll_dup <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = t[idx], status = status[idx],
    time_upper = time_upper[idx], x = matrix(x[idx], ncol = 1),
    weights = rep(1, sum(w))
  )

  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Weighted fit == duplicated-row fit (MLE + vcov)
# ---------------------------------------------------------------------------

test_that("fit with integer weights matches the duplicated-row fit (Weibull)", {
  skip_on_cran()
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  df_exp <- expand_rows(df, w)

  # Single-distribution fits only run the optimizer when start values are
  # supplied (hazard()'s `fit && !is.null(theta)` branch), so pass `theta`.
  fit_w <- hazard(
    survival::Surv(time, status) ~ x,
    data = df, dist = "weibull", theta = c(0.3, 1, 0),
    weights = w, fit = TRUE
  )
  fit_dup <- hazard(
    survival::Surv(time, status) ~ x,
    data = df_exp, dist = "weibull", theta = c(0.3, 1, 0),
    fit = TRUE
  )

  # Coefficients should match up to optimizer tolerance.
  expect_equal(coef(fit_w), coef(fit_dup), tolerance = 1e-3)
  # Log-likelihood should match tightly.
  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# Unit weights == unweighted (sanity)
# ---------------------------------------------------------------------------

test_that("unit weights are equivalent to no weights (Weibull)", {
  df <- make_toy()
  theta <- c(mu = 0.8, nu = 1.2, beta = 0.25)

  ll_noweight <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = df$time, status = df$status,
    x = matrix(df$x, ncol = 1)
  )
  ll_unit <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = df$time, status = df$status,
    x = matrix(df$x, ncol = 1), weights = rep(1, nrow(df))
  )
  expect_equal(ll_noweight, ll_unit, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Multiphase weighted likelihood equals duplicated-row likelihood
# ---------------------------------------------------------------------------

test_that("integer weights match duplication in multiphase likelihood", {
  skip_on_cran()
  set.seed(11)
  n <- 40
  t <- runif(n, 0.05, 3)
  status <- rbinom(n, 1, 0.4)
  w <- sample(1:2, n, replace = TRUE)
  idx <- rep(seq_len(n), times = w)

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  # Start values: [log_mu_early, t_half, nu, m, log_mu_constant]
  theta <- c(-4, log(0.3), 1, 1, -3)

  ll_w <- TemporalHazard:::.hzr_logl_multiphase(
    theta = theta, time = t, status = status,
    phases = phases,
    covariate_counts = c(early = 0L, constant = 0L),
    x_list = list(early = NULL, constant = NULL),
    weights = w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_multiphase(
    theta = theta, time = t[idx], status = status[idx],
    phases = phases,
    covariate_counts = c(early = 0L, constant = 0L),
    x_list = list(early = NULL, constant = NULL),
    weights = rep(1, sum(w))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Gradient is weighted (numerical check on Weibull)
# ---------------------------------------------------------------------------

test_that("weighted log-likelihood gradient equals numerical gradient", {
  skip_if_not_installed("numDeriv")
  df <- make_toy()
  w <- sample(1:3, nrow(df), replace = TRUE)
  theta <- c(mu = 0.8, nu = 1.2, beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  # Numerical gradient of the weighted LL
  num_g <- numDeriv::grad(
    function(th) {
      TemporalHazard:::.hzr_logl_weibull(
        theta = th, time = df$time, status = df$status,
        x = xm, weights = w
      )
    },
    theta
  )

  ana_g <- TemporalHazard:::.hzr_gradient_weibull(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )

  expect_equal(as.numeric(ana_g), num_g, tolerance = 1e-5)
})

# ---------------------------------------------------------------------------
# Zero weight excludes the row
# ---------------------------------------------------------------------------

test_that("a zero weight drops the row's contribution to the likelihood", {
  df <- make_toy()
  w <- rep(1, nrow(df))
  w[1] <- 0

  theta <- c(mu = 0.8, nu = 1.2, beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  ll_w <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  ll_drop <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = df$time[-1], status = df$status[-1],
    x = xm[-1, , drop = FALSE],
    weights = rep(1, nrow(df) - 1)
  )
  expect_equal(ll_w, ll_drop, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Per-distribution duplication-parity for exp / log-logistic / log-normal
# (wired up in v0.9.6 as part of Phase 4e).
# ---------------------------------------------------------------------------

test_that("integer weights match row duplication (exponential LL)", {
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  idx <- rep(seq_len(n), times = w)
  theta <- c(log_lambda = log(0.6), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  ll_w <- TemporalHazard:::.hzr_logl_exponential(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_exponential(
    theta = theta, time = df$time[idx], status = df$status[idx],
    x = matrix(df$x[idx], ncol = 1), weights = rep(1, sum(w))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

test_that("integer weights match row duplication (log-logistic LL)", {
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  idx <- rep(seq_len(n), times = w)
  theta <- c(log_alpha = log(0.6), log_beta = log(1.2), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  ll_w <- TemporalHazard:::.hzr_logl_loglogistic(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_loglogistic(
    theta = theta, time = df$time[idx], status = df$status[idx],
    x = matrix(df$x[idx], ncol = 1), weights = rep(1, sum(w))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

test_that("integer weights match row duplication (log-normal LL)", {
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  idx <- rep(seq_len(n), times = w)
  theta <- c(mu = 0, log_sigma = log(1), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  ll_w <- TemporalHazard:::.hzr_logl_lognormal(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_lognormal(
    theta = theta, time = df$time[idx], status = df$status[idx],
    x = matrix(df$x[idx], ncol = 1), weights = rep(1, sum(w))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

# Analytic gradients for exp / log-logistic / log-normal match numerical
# gradients of the weighted LL.  Catches the "accepts formal but never
# applies it" pattern in the analytic score path.

test_that("weighted exponential gradient equals numerical gradient", {
  skip_if_not_installed("numDeriv")
  df <- make_toy()
  w <- sample(1:3, nrow(df), replace = TRUE)
  theta <- c(log_lambda = log(0.6), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  num_g <- numDeriv::grad(
    function(th) {
      TemporalHazard:::.hzr_logl_exponential(
        theta = th, time = df$time, status = df$status,
        x = xm, weights = w
      )
    },
    theta
  )
  ana_g <- TemporalHazard:::.hzr_gradient_exponential(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  expect_equal(as.numeric(ana_g), num_g, tolerance = 1e-5)
})

test_that("weighted log-logistic gradient equals numerical gradient", {
  skip_if_not_installed("numDeriv")
  df <- make_toy()
  w <- sample(1:3, nrow(df), replace = TRUE)
  theta <- c(log_alpha = log(0.6), log_beta = log(1.2), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  num_g <- numDeriv::grad(
    function(th) {
      TemporalHazard:::.hzr_logl_loglogistic(
        theta = th, time = df$time, status = df$status,
        x = xm, weights = w
      )
    },
    theta
  )
  ana_g <- TemporalHazard:::.hzr_gradient_loglogistic(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  expect_equal(as.numeric(ana_g), num_g, tolerance = 1e-5)
})

test_that("weighted log-normal gradient equals numerical gradient", {
  skip_if_not_installed("numDeriv")
  df <- make_toy()
  w <- sample(1:3, nrow(df), replace = TRUE)
  theta <- c(mu = 0, log_sigma = log(1), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  num_g <- numDeriv::grad(
    function(th) {
      TemporalHazard:::.hzr_logl_lognormal(
        theta = th, time = df$time, status = df$status,
        x = xm, weights = w
      )
    },
    theta
  )
  ana_g <- TemporalHazard:::.hzr_gradient_lognormal(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  expect_equal(as.numeric(ana_g), num_g, tolerance = 1e-5)
})

# End-to-end: weighted fit == duplicated-row fit for each distribution.

test_that("fit with integer weights matches the duplicated-row fit (exponential)", {
  skip_on_cran()
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  df_exp <- expand_rows(df, w)

  fit_w <- hazard(
    survival::Surv(time, status) ~ x, data = df,
    dist = "exponential", theta = c(log(0.6), 0), weights = w, fit = TRUE
  )
  fit_dup <- hazard(
    survival::Surv(time, status) ~ x, data = df_exp,
    dist = "exponential", theta = c(log(0.6), 0), fit = TRUE
  )

  expect_equal(coef(fit_w), coef(fit_dup), tolerance = 1e-3)
  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-4)
})

test_that("fit with integer weights matches the duplicated-row fit (log-logistic)", {
  skip_on_cran()
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  df_exp <- expand_rows(df, w)

  fit_w <- hazard(
    survival::Surv(time, status) ~ x, data = df,
    dist = "loglogistic", theta = c(0, 0.2, 0), weights = w, fit = TRUE
  )
  fit_dup <- hazard(
    survival::Surv(time, status) ~ x, data = df_exp,
    dist = "loglogistic", theta = c(0, 0.2, 0), fit = TRUE
  )

  expect_equal(coef(fit_w), coef(fit_dup), tolerance = 1e-3)
  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-4)
})

test_that("fit with integer weights matches the duplicated-row fit (log-normal)", {
  skip_on_cran()
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  df_exp <- expand_rows(df, w)

  fit_w <- hazard(
    survival::Surv(time, status) ~ x, data = df,
    dist = "lognormal", theta = c(1, 0, 0), weights = w, fit = TRUE
  )
  fit_dup <- hazard(
    survival::Surv(time, status) ~ x, data = df_exp,
    dist = "lognormal", theta = c(1, 0, 0), fit = TRUE
  )

  expect_equal(coef(fit_w), coef(fit_dup), tolerance = 1e-3)
  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# Mixed-censoring (status == -1 and status == 2) duplication-parity for
# exp / log-logistic / log-normal.  The weighted LL now applies row
# weights to the left-censored term `log(1 - exp(-H(u)))` and the
# interval-censored term `log(F(u) - F(l))` as well; these tests
# exercise those newly weighted branches.
# ---------------------------------------------------------------------------

make_mixed <- function(seed = 23) {
  set.seed(seed)
  n <- 24
  x <- rnorm(n)
  t <- runif(n, 0.1, 2.0)
  # 10 events, 4 right-censored, 5 left-censored, 5 interval-censored
  status <- c(rep(1L, 10), rep(0L, 4), rep(-1L, 5), rep(2L, 5))
  # Interval bounds must be strictly positive with lower < upper.
  # Left-censored rows set time_upper = time; non-interval rows set
  # time_lower = time (ignored by the LL except on interval rows).
  time_lower <- t
  time_upper <- t
  idx_iv <- which(status == 2)
  time_lower[idx_iv] <- t[idx_iv]
  time_upper[idx_iv] <- t[idx_iv] + runif(length(idx_iv), 0.05, 0.3)
  list(time = t, status = status, x = x,
       time_lower = time_lower, time_upper = time_upper)
}

test_that("weighted LL matches duplicated-row LL under mixed censoring (exponential)", {
  m <- make_mixed()
  n <- length(m$time)
  w <- sample(1:3, n, replace = TRUE)
  idx <- rep(seq_len(n), times = w)
  theta <- c(log_lambda = log(0.5), beta = 0.2)
  xm <- matrix(m$x, ncol = 1)

  ll_w <- TemporalHazard:::.hzr_logl_exponential(
    theta = theta, time = m$time, status = m$status,
    time_lower = m$time_lower, time_upper = m$time_upper,
    x = xm, weights = w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_exponential(
    theta = theta,
    time = m$time[idx], status = m$status[idx],
    time_lower = m$time_lower[idx], time_upper = m$time_upper[idx],
    x = matrix(m$x[idx], ncol = 1), weights = rep(1, sum(w))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

test_that("weighted LL matches duplicated-row LL under mixed censoring (log-logistic)", {
  m <- make_mixed()
  n <- length(m$time)
  w <- sample(1:3, n, replace = TRUE)
  idx <- rep(seq_len(n), times = w)
  theta <- c(log_alpha = log(0.5), log_beta = log(1.1), beta = 0.2)
  xm <- matrix(m$x, ncol = 1)

  ll_w <- TemporalHazard:::.hzr_logl_loglogistic(
    theta = theta, time = m$time, status = m$status,
    time_lower = m$time_lower, time_upper = m$time_upper,
    x = xm, weights = w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_loglogistic(
    theta = theta,
    time = m$time[idx], status = m$status[idx],
    time_lower = m$time_lower[idx], time_upper = m$time_upper[idx],
    x = matrix(m$x[idx], ncol = 1), weights = rep(1, sum(w))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

test_that("weighted LL matches duplicated-row LL under mixed censoring (log-normal)", {
  m <- make_mixed()
  n <- length(m$time)
  w <- sample(1:3, n, replace = TRUE)
  idx <- rep(seq_len(n), times = w)
  theta <- c(mu = 0, log_sigma = log(1), beta = 0.2)
  xm <- matrix(m$x, ncol = 1)

  ll_w <- TemporalHazard:::.hzr_logl_lognormal(
    theta = theta, time = m$time, status = m$status,
    time_lower = m$time_lower, time_upper = m$time_upper,
    x = xm, weights = w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_lognormal(
    theta = theta,
    time = m$time[idx], status = m$status[idx],
    time_lower = m$time_lower[idx], time_upper = m$time_upper[idx],
    x = matrix(m$x[idx], ncol = 1), weights = rep(1, sum(w))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

test_that("weighted multiphase analytic gradient matches numerical gradient", {
  # `.hzr_gradient_multiphase()` must apply `weights` to the per-row score
  # (w_H, inv_h, and the interval-censored finite-difference correction).
  # Prior to 0.9.6 it ignored the argument, so weighted multiphase fits
  # optimized a weighted objective with an unweighted gradient -- the fit
  # still converged near the right MLE via line search but the final
  # gradient norm did not go to zero.  This test would have caught it.
  skip_if_not_installed("numDeriv")

  set.seed(31)
  n <- 40
  t <- runif(n, 0.1, 3)
  status <- rbinom(n, 1, 0.45)
  w <- sample(1:3, n, replace = TRUE)

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  cov_counts <- c(early = 0L, constant = 0L)
  x_list <- list(early = NULL, constant = NULL)
  # [log_mu_e, log_thalf, nu, m, log_mu_c]
  theta <- c(-4, log(0.3), 1, 1, -3)

  num_g <- numDeriv::grad(
    function(th) {
      TemporalHazard:::.hzr_logl_multiphase(
        theta = th, time = t, status = status,
        phases = phases, covariate_counts = cov_counts,
        x_list = x_list, weights = w
      )
    },
    theta
  )
  ana_g <- TemporalHazard:::.hzr_gradient_multiphase(
    theta = theta, time = t, status = status,
    weights = w,
    phases = phases, covariate_counts = cov_counts, x_list = x_list
  )
  expect_equal(as.numeric(ana_g), num_g, tolerance = 1e-4)
})

test_that("zero weight drops the row's contribution across mixed censoring (log-normal)", {
  # Targeted sanity check of the per-term weighting: zeroing one row of each
  # censoring flavour must equal dropping those rows entirely.
  m <- make_mixed()
  n <- length(m$time)
  w <- rep(1, n)
  # Zero out one event, one RC, one left, one interval row.
  w[c(1, 11, 15, 20)] <- 0
  keep <- w > 0
  theta <- c(mu = 0, log_sigma = log(1), beta = 0.2)
  xm <- matrix(m$x, ncol = 1)

  ll_w <- TemporalHazard:::.hzr_logl_lognormal(
    theta = theta, time = m$time, status = m$status,
    time_lower = m$time_lower, time_upper = m$time_upper,
    x = xm, weights = w
  )
  ll_drop <- TemporalHazard:::.hzr_logl_lognormal(
    theta = theta,
    time = m$time[keep], status = m$status[keep],
    time_lower = m$time_lower[keep], time_upper = m$time_upper[keep],
    x = xm[keep, , drop = FALSE], weights = rep(1, sum(keep))
  )
  expect_equal(ll_w, ll_drop, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Multiphase WITH covariates: weighted == duplicated (roadmap 7c)
# ---------------------------------------------------------------------------

# Shared fixture: multiphase data with a covariate in BOTH phases, integer
# weights, and the row-duplication index. covariate enters early + constant.
make_mp_cov <- function(seed = 41, n = 40) {
  set.seed(seed)
  t      <- runif(n, 0.05, 3)
  status <- rbinom(n, 1, 0.4)
  x      <- rnorm(n)
  w      <- sample(1:3, n, replace = TRUE)
  idx    <- rep(seq_len(n), times = w)
  list(t = t, status = status, x = x, w = w, idx = idx, n = n)
}

# Phases + theta shared by the direct-LL and gradient tests.
# Layout: [log_mu_e, log_thalf, nu, m, beta_e, log_mu_c, beta_c]
mp_cov_phases <- function() {
  list(
    early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
}
mp_cov_theta   <- c(-4, log(0.3), 1, 1, 0.5, -3, -0.3)
mp_cov_counts  <- c(early = 1L, constant = 1L)

test_that("integer weights match duplication in multiphase LL WITH covariates", {
  skip_on_cran()
  d <- make_mp_cov()
  x_list   <- list(early    = matrix(d$x, ncol = 1),
                   constant = matrix(d$x, ncol = 1))
  x_list_d <- list(early    = matrix(d$x[d$idx], ncol = 1),
                   constant = matrix(d$x[d$idx], ncol = 1))

  ll_w <- TemporalHazard:::.hzr_logl_multiphase(
    theta = mp_cov_theta, time = d$t, status = d$status,
    phases = mp_cov_phases(), covariate_counts = mp_cov_counts,
    x_list = x_list, weights = d$w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_multiphase(
    theta = mp_cov_theta, time = d$t[d$idx], status = d$status[d$idx],
    phases = mp_cov_phases(), covariate_counts = mp_cov_counts,
    x_list = x_list_d, weights = rep(1, length(d$idx))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

test_that("weighted multiphase analytic gradient matches numerical WITH covariates", {
  skip_if_not_installed("numDeriv")
  d <- make_mp_cov(seed = 42)
  x_list <- list(early    = matrix(d$x, ncol = 1),
                 constant = matrix(d$x, ncol = 1))

  num_g <- numDeriv::grad(
    function(th) {
      TemporalHazard:::.hzr_logl_multiphase(
        theta = th, time = d$t, status = d$status,
        phases = mp_cov_phases(), covariate_counts = mp_cov_counts,
        x_list = x_list, weights = d$w
      )
    },
    mp_cov_theta
  )
  ana_g <- TemporalHazard:::.hzr_gradient_multiphase(
    theta = mp_cov_theta, time = d$t, status = d$status, weights = d$w,
    phases = mp_cov_phases(), covariate_counts = mp_cov_counts, x_list = x_list
  )
  expect_equal(as.numeric(ana_g), num_g, tolerance = 1e-4)
})

test_that("weighted multiphase fit matches duplicated-row fit, covariates, CoE off", {
  skip_on_cran()
  d  <- make_mp_cov(seed = 43)
  df <- data.frame(time = d$t, status = d$status, x = d$x)
  df_exp <- expand_rows(df, d$w)

  mk_fit <- function(data, weights) {
    hazard(
      survival::Surv(time, status) ~ 1,
      data = data, dist = "multiphase",
      phases = list(
        early    = hzr_phase("cdf", t_half = 0.5, nu = 0, m = -0.4,
                             fixed = c("nu", "m"), formula = ~ x),
        constant = hzr_phase("constant", formula = ~ x)
      ),
      weights = weights, fit = TRUE,
      control = list(conserve = FALSE, n_starts = 1L, maxit = 500L)
    )
  }

  fit_w   <- mk_fit(df,     d$w)
  fit_dup <- mk_fit(df_exp, NULL)

  # Parity only means something over parameters the data determines. The
  # former spec pinned all three shapes (`fixed = "shapes"` is t_half, nu and
  # m), which left `early` contributing 2.7e-09 of the cumulative hazard: its
  # mu and beta drifted freely, and the two fits agreed on them only because
  # both took the identical deterministic path from the identical start
  # (#216). Freeing t_half revives the phase -- rcond 5.3e-12 -> 2.4e-04 --
  # so assert it is alive before trusting the comparison below.
  early_share <- function(f) {
    f$fit$phase_share$share[f$fit$phase_share$phase == "early"]
  }
  expect_gt(early_share(fit_w),   1e-3)
  expect_gt(early_share(fit_dup), 1e-3)

  # Observed agreement here is ~5e-17, far inside `tolerance`: the weighted and
  # duplicated objectives are algebraically the same function, so the optimizer
  # walks one path. That near-exactness is the expected result, not a sign that
  # nothing was compared -- the guard against that is `early_share` above. The
  # tolerance stays loose deliberately, because the one path is not guaranteed
  # to be bit-identical across the five CI platforms' BLAS implementations.
  expect_equal(coef(fit_w), coef(fit_dup), tolerance = 1e-3)
  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-4)
})

test_that(".hzr_logl_weibull gradient attribute respects weights", {
  # Regression: return_gradient = TRUE previously called .hzr_gradient_weibull
  # without forwarding `weights`, so the attached gradient was unit-weighted
  # even for weighted data.
  skip_if_not_installed("numDeriv")
  set.seed(314)
  n <- 120
  x <- cbind(z = rnorm(n))
  time <- rexp(n, 0.5) + 0.01
  status <- rbinom(n, 1, 0.75)
  w <- runif(n, 0.5, 2.5)
  theta <- c(mu = 0.8, nu = 1.4, z = 0.3)

  ll <- .hzr_logl_weibull(theta, time, status, x = x, weights = w,
                          return_gradient = TRUE)
  g_attr <- attr(ll, "gradient")

  obj <- function(th) {
    .hzr_logl_weibull(th, time, status, x = x, weights = w)
  }
  g_num <- numDeriv::grad(obj, theta)

  expect_equal(g_attr, g_num, tolerance = 1e-5)
})

# ---------------------------------------------------------------------------
# FRACTIONAL (non-integer) weights -- roadmap 7a.
#
# Every test above proves weighting via *integer* row duplication.  Fractional
# weights (e.g. inverse-probability or aggregated-cell weights) cannot be
# checked that way -- you cannot duplicate a row 2.5 times.  Instead we verify
# the two algebraic properties that DEFINE a correctly weighted, per-row
# additive log-likelihood  L(theta; w) = sum_i w_i * l_i(theta):
#
#   (A) ADDITIVE SPLIT.  Splitting one row of weight w_i into two identical
#       copies whose weights sum to w_i leaves L unchanged.  This is the
#       continuous generalisation of duplication -- it reduces to the integer
#       duplication tests when the split weights are 1, and exercises the
#       fractional arithmetic when they are not.  If any term failed to scale
#       linearly in w_i (e.g. a w_i applied inside a log, or a stray w_i^2),
#       the split would not balance.
#   (B) LINEAR SCALING.  L(theta; c*w) = c * L(theta; w) and the gradient
#       scales likewise; the MLE (argmax) is therefore invariant to a positive
#       scalar reweight while the objective scales by c.
#
# Together (A) and (B) pin the fractional-weight path without a SAS fixture --
# a SAS `WEIGHT`-statement parity run (roadmap 7a) remains the external check,
# but these establish internal correctness now.
# ---------------------------------------------------------------------------

# Replace each row i (weight w_i) with two identical copies carrying weights
# w_i * f_i and w_i * (1 - f_i).  Returns the doubled vectors + split weights.
.frac_split <- function(time, status, x, w, f) {
  list(
    time   = c(time, time),
    status = c(status, status),
    x      = rbind(x, x),
    w      = c(w * f, w * (1 - f))
  )
}

test_that("fractional weights: additive split leaves the LL unchanged (Weibull)", {
  df <- make_toy()
  n  <- nrow(df)
  set.seed(71)                      # seed the weight/split draws (make_toy() seeds its own data)
  w  <- runif(n, 0.3, 2.7)          # fractional
  f  <- runif(n, 0.1, 0.9)          # fractional split point per row
  theta <- c(mu = 0.8, nu = 1.2, beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  s <- .frac_split(df$time, df$status, xm, w, f)

  ll_w <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = df$time, status = df$status, x = xm, weights = w)
  ll_split <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = s$time, status = s$status, x = s$x, weights = s$w)

  expect_equal(ll_w, ll_split, tolerance = 1e-10)
})

test_that("fractional weights: additive split leaves the LL unchanged (exponential)", {
  df <- make_toy()
  n  <- nrow(df)
  set.seed(72)                      # seed the weight/split draws (make_toy() seeds its own data)
  w  <- runif(n, 0.3, 2.7)
  f  <- runif(n, 0.1, 0.9)
  theta <- c(log_lambda = log(0.6), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  s <- .frac_split(df$time, df$status, xm, w, f)

  ll_w <- TemporalHazard:::.hzr_logl_exponential(
    theta = theta, time = df$time, status = df$status, x = xm, weights = w)
  ll_split <- TemporalHazard:::.hzr_logl_exponential(
    theta = theta, time = s$time, status = s$status, x = s$x, weights = s$w)

  expect_equal(ll_w, ll_split, tolerance = 1e-10)
})

test_that("fractional weights: additive split leaves the LL unchanged (multiphase + covariates)", {
  d <- make_mp_cov(seed = 73)
  set.seed(73)                      # seed the weight/split draws (make_mp_cov() seeds its own data)
  w <- runif(d$n, 0.3, 2.7)
  f <- runif(d$n, 0.1, 0.9)
  xm <- matrix(d$x, ncol = 1)

  x_list   <- list(early = xm, constant = xm)
  x_list_s <- list(early = rbind(xm, xm), constant = rbind(xm, xm))

  ll_w <- TemporalHazard:::.hzr_logl_multiphase(
    theta = mp_cov_theta, time = d$t, status = d$status,
    phases = mp_cov_phases(), covariate_counts = mp_cov_counts,
    x_list = x_list, weights = w)
  ll_split <- TemporalHazard:::.hzr_logl_multiphase(
    theta = mp_cov_theta, time = c(d$t, d$t), status = c(d$status, d$status),
    phases = mp_cov_phases(), covariate_counts = mp_cov_counts,
    x_list = x_list_s, weights = c(w * f, w * (1 - f)))

  expect_equal(ll_w, ll_split, tolerance = 1e-10)
})

test_that("fractional weights: LL and gradient scale linearly in the weights (Weibull)", {
  # No numDeriv needed: this checks the analytic LL/gradient scale exactly in w.
  df <- make_toy()
  set.seed(74)                      # seed the weight draws (make_toy() seeds its own data)
  w  <- runif(nrow(df), 0.3, 2.7)
  cc <- 2.5                         # fractional scale factor
  theta <- c(mu = 0.8, nu = 1.2, beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  ll_1 <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = df$time, status = df$status, x = xm, weights = w)
  ll_c <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = df$time, status = df$status, x = xm, weights = cc * w)
  expect_equal(ll_c, cc * ll_1, tolerance = 1e-10)

  g_1 <- TemporalHazard:::.hzr_gradient_weibull(
    theta = theta, time = df$time, status = df$status, x = xm, weights = w)
  g_c <- TemporalHazard:::.hzr_gradient_weibull(
    theta = theta, time = df$time, status = df$status, x = xm, weights = cc * w)
  expect_equal(as.numeric(g_c), cc * as.numeric(g_1), tolerance = 1e-10)
})

test_that("fractional weights: the MLE is invariant to a positive scalar reweight (Weibull)", {
  skip_on_cran()
  df <- make_toy()
  set.seed(75)                      # seed the weight draws (make_toy() seeds its own data)
  w  <- runif(nrow(df), 0.3, 2.7)
  cc <- 2.0
  # Single-distribution fits run the optimizer only when start values are
  # supplied (the `fit && !is.null(theta)` branch in hazard()); pass `theta`
  # so this exercises a real fit rather than an unfitted object.
  theta0 <- c(0.3, 1, 0)

  fit_w  <- hazard(survival::Surv(time, status) ~ x, data = df,
                   dist = "weibull", theta = theta0, weights = w,      fit = TRUE)
  fit_cw <- hazard(survival::Surv(time, status) ~ x, data = df,
                   dist = "weibull", theta = theta0, weights = cc * w, fit = TRUE)

  # argmax unchanged; the maximised objective scales by the constant.
  expect_equal(coef(fit_cw), coef(fit_w), tolerance = 5e-3)
  expect_equal(fit_cw$fit$objective / fit_w$fit$objective, cc, tolerance = 1e-4)
})

test_that("fractional weights: additive-split fit matches the original fit (multiphase + covariates)", {
  skip_on_cran()
  d  <- make_mp_cov(seed = 76)
  set.seed(76)                      # seed the weight/split draws (make_mp_cov() seeds its own data)
  w  <- runif(d$n, 0.3, 2.7)
  f  <- runif(d$n, 0.1, 0.9)
  df <- data.frame(time = d$t, status = d$status, x = d$x)
  df_split <- rbind(df, df)
  w_split  <- c(w * f, w * (1 - f))

  mk_fit <- function(data, weights) {
    hazard(
      survival::Surv(time, status) ~ 1,
      data = data, dist = "multiphase",
      phases = list(
        early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                             fixed = "shapes", formula = ~ x),
        constant = hzr_phase("constant", formula = ~ x)
      ),
      weights = weights, fit = TRUE,
      control = list(conserve = FALSE, n_starts = 1L, maxit = 500L)
    )
  }

  fit_w     <- suppressWarnings(mk_fit(df,       w))
  fit_split <- suppressWarnings(mk_fit(df_split, w_split))

  expect_equal(coef(fit_w), coef(fit_split), tolerance = 1e-3)
  expect_equal(fit_w$fit$objective, fit_split$fit$objective, tolerance = 1e-4)
})
