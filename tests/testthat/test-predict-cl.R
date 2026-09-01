# tests/testthat/test-predict-cl.R
# Phase 4g -- predict.hazard(..., se.fit = TRUE, level = 0.95)
#
# Coverage:
# (1) API contract: output shape, level validation, decompose incompat
# (2) Jacobian parity: analytic vs numeric for Weibull + multiphase
# (3) Delta-method scale: log-CL for hazard/cumhaz, log-log for survival
# (4) Monotonicity: lower <= fit <= upper
# (5) Survival CLs stay in [0, 1]; cumhaz CLs stay in [0, Inf)
# (6) `vcov` NA (non-converged fit): returns NA CLs with a warning
# (7) Fixed shapes / CoE: free-parameter submatrix recovers meaningful CLs
# (8) Default of `se.fit = FALSE` is byte-identical to pre-0.9.8 predict()

make_toy <- function(seed = 17, n = 60) {
  set.seed(seed)
  data.frame(
    time = rexp(n, 0.3),
    status = rbinom(n, 1, 0.6),
    x = rnorm(n)
  )
}

# ---------------------------------------------------------------------------
# (1) API contract
# ---------------------------------------------------------------------------

test_that("se.fit returns a data frame with fit/se.fit/lower/upper", {
  skip_on_cran()
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)
  res <- predict(fit, newdata = data.frame(time = c(0.5, 1, 2), x = 0),
                 type = "survival", se.fit = TRUE)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("fit", "se.fit", "lower", "upper"))
  expect_equal(nrow(res), 3L)
  expect_true(all(is.finite(res$fit)))
  expect_true(all(is.finite(res$se.fit)))
  expect_true(all(res$se.fit >= 0))
})

test_that("invalid level triggers a clean error", {
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)
  expect_error(
    predict(fit, newdata = data.frame(time = 1, x = 0),
            type = "survival", se.fit = TRUE, level = 1.1),
    "level"
  )
  expect_error(
    predict(fit, newdata = data.frame(time = 1, x = 0),
            type = "survival", se.fit = TRUE, level = 0),
    "level"
  )
})

test_that("decompose + se.fit: cumulative_hazard works, survival errors", {
  skip_on_cran()
  df <- make_toy()
  phases <- list(early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                     fixed = "shapes"),
                 constant = hzr_phase("constant"))
  fit <- hazard(survival::Surv(time, status) ~ 1, data = df,
                dist = "multiphase", phases = phases, fit = TRUE,
                control = list(n_starts = 2))
  t_new <- c(0.5, 1, 2)

  # cumulative_hazard: long per-phase + total frame.
  res <- predict(fit, newdata = data.frame(time = t_new),
                 type = "cumulative_hazard", se.fit = TRUE, decompose = TRUE)
  expect_named(res, c("time", "component", "fit", "se.fit", "lower", "upper"))
  expect_equal(levels(res$component), c("total", "early", "constant"))
  expect_equal(nrow(res), length(t_new) * 3L)

  # Total rows reproduce the aggregate se.fit output exactly.
  agg <- predict(fit, newdata = data.frame(time = t_new),
                 type = "cumulative_hazard", se.fit = TRUE, decompose = FALSE)
  tot <- res[res$component == "total", ]
  expect_equal(tot$fit, agg$fit, tolerance = 1e-10)
  expect_equal(tot$se.fit, agg$se.fit, tolerance = 1e-10)
  expect_equal(tot$lower, agg$lower, tolerance = 1e-10)
  expect_equal(tot$upper, agg$upper, tolerance = 1e-10)

  # Per-phase point estimates reproduce the decompose=TRUE, se.fit=FALSE wide
  # cumhaz columns.
  wide <- predict(fit, newdata = data.frame(time = t_new),
                  type = "cumulative_hazard", se.fit = FALSE, decompose = TRUE)
  expect_equal(res[res$component == "early", "fit"], wide$early, tolerance = 1e-10)
  expect_equal(res[res$component == "constant", "fit"], wide$constant,
               tolerance = 1e-10)

  # component is an ordered factor (documented semantics).
  expect_true(is.ordered(res$component))

  # survival still rejected for multiphase.
  expect_error(
    predict(fit, newdata = data.frame(time = t_new),
            type = "survival", se.fit = TRUE, decompose = TRUE),
    "cumulative_hazard"
  )
})

test_that("decompose is ignored (not errored) for single-distribution se.fit", {
  # `decompose` applies only to multiphase models. A single-distribution fit
  # must not error on decompose = TRUE -- it returns the standard aggregate
  # se.fit frame regardless of type (decompose silently ignored).
  skip_on_cran()
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)
  res <- predict(fit, newdata = data.frame(time = c(0.5, 1, 2), x = 0),
                 type = "survival", se.fit = TRUE, decompose = TRUE)
  expect_named(res, c("fit", "se.fit", "lower", "upper"))
  expect_equal(nrow(res), 3L)
})

test_that("decompose + se.fit per-phase SE matches a numeric jacobian", {
  skip_if_not_installed("numDeriv")
  skip_on_cran()
  df <- make_toy()
  phases <- list(early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                     fixed = "shapes"),
                 constant = hzr_phase("constant"))
  fit <- hazard(survival::Surv(time, status) ~ 1, data = df,
                dist = "multiphase", phases = phases, fit = TRUE,
                control = list(n_starts = 2))
  t_new <- c(0.5, 1, 2)
  cov_counts <- c(early = 0L, constant = 0L)
  x_list <- list(early = NULL, constant = NULL)
  theta <- fit$fit$theta
  V <- fit$fit$vcov
  free_idx <- which(is.finite(diag(V)))

  res <- predict(fit, newdata = data.frame(time = t_new),
                 type = "cumulative_hazard", se.fit = TRUE, decompose = TRUE)

  # Numeric per-phase SE for the "early" component via numDeriv.
  H_early <- function(th) {
    TemporalHazard:::.hzr_multiphase_cumhaz(
      t_new, th, phases, cov_counts, x_list, per_phase = TRUE
    )$early
  }
  Jn <- numDeriv::jacobian(H_early, theta)[, free_idx, drop = FALSE]
  se_num <- sqrt(diag(Jn %*% V[free_idx, free_idx, drop = FALSE] %*% t(Jn)))
  se_ana <- res[res$component == "early", "se.fit"]
  expect_equal(se_ana, se_num, tolerance = 1e-5)
})

# ---------------------------------------------------------------------------
# (2) Jacobian parity: analytic vs numeric
# ---------------------------------------------------------------------------

test_that("Weibull analytic jacobian matches numDeriv (cumhaz)", {
  skip_if_not_installed("numDeriv")
  set.seed(3)
  t_new <- c(0.3, 0.8, 1.5, 2.2)
  x_new <- matrix(c(-0.5, 0, 0.4, 1.2), ncol = 1)
  theta <- c(mu = 0.6, nu = 1.2, beta = 0.25)

  cumhaz_fn <- function(th) {
    (th[1] * t_new) ^ th[2] * exp(as.numeric(x_new %*% th[3:length(th)]))
  }
  J_ana <- TemporalHazard:::.hzr_predict_jacobian_weibull(
    "cumulative_hazard", theta, t_new, x_new, length(theta)
  )
  J_num <- numDeriv::jacobian(cumhaz_fn, theta)
  expect_equal(J_ana, J_num, tolerance = 1e-6)
})

test_that("Weibull analytic jacobian matches numDeriv (hazard / relative)", {
  skip_if_not_installed("numDeriv")
  x_new <- matrix(c(-0.5, 0, 0.4), ncol = 1)
  theta <- c(mu = 0.6, nu = 1.2, beta = 0.25)
  haz_fn <- function(th) exp(as.numeric(x_new %*% th[3:length(th)]))
  J_ana <- TemporalHazard:::.hzr_predict_jacobian_weibull(
    "hazard", theta, time = NULL, x_new, length(theta)
  )
  J_num <- numDeriv::jacobian(haz_fn, theta)
  expect_equal(J_ana, J_num, tolerance = 1e-6)
})

test_that("Multiphase analytic jacobian matches numDeriv (cumhaz)", {
  skip_if_not_installed("numDeriv")
  skip_on_cran()
  t_new <- c(0.3, 0.8, 1.5, 2.2)
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  cov_counts <- c(early = 0L, constant = 0L)
  x_list <- list(early = NULL, constant = NULL)
  theta <- c(-3, log(0.3), 1, 1, -2)

  J_ana <- TemporalHazard:::.hzr_predict_jacobian_multiphase(
    theta, t_new, phases, cov_counts, x_list, length(theta)
  )
  J_num <- numDeriv::jacobian(
    function(th) {
      TemporalHazard:::.hzr_multiphase_cumhaz(
        t_new, th, phases, cov_counts, x_list
      )
    },
    theta
  )
  expect_equal(J_ana, J_num, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# (3, 4, 5) Scale, monotonicity, and range
# ---------------------------------------------------------------------------

test_that("survival CLs stay in [0, 1] and straddle the point estimate", {
  skip_on_cran()
  df <- make_toy()
  for (dist in c("weibull", "exponential", "loglogistic", "lognormal")) {
    theta_init <- switch(
      dist,
      weibull = c(0.5, 1, 0),
      exponential = c(log(0.3), 0),
      loglogistic = c(0, 0, 0),
      lognormal = c(0, 0, 0)
    )
    fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = dist,
                  theta = theta_init, fit = TRUE)
    res <- predict(fit, newdata = data.frame(time = c(0.5, 1, 2), x = 0),
                   type = "survival", se.fit = TRUE)
    expect_true(all(res$lower >= 0), info = dist)
    expect_true(all(res$upper <= 1), info = dist)
    expect_true(all(res$lower <= res$fit + 1e-8), info = dist)
    expect_true(all(res$fit <= res$upper + 1e-8), info = dist)
  }
})

test_that("cumulative_hazard CLs are strictly positive and log-symmetric", {
  skip_on_cran()
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)
  res <- predict(fit, newdata = data.frame(time = c(0.5, 1, 2), x = 0),
                 type = "cumulative_hazard", se.fit = TRUE)
  expect_true(all(res$lower > 0))
  expect_true(all(res$upper > res$lower))
  # log-symmetric: log(upper) - log(fit) ~= log(fit) - log(lower)
  log_halfwidth_upper <- log(res$upper) - log(res$fit)
  log_halfwidth_lower <- log(res$fit) - log(res$lower)
  expect_equal(log_halfwidth_upper, log_halfwidth_lower, tolerance = 1e-6)
})

test_that("linear_predictor CLs are symmetric on the natural scale", {
  skip_on_cran()
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)
  res <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)),
                 type = "linear_predictor", se.fit = TRUE)
  halfwidth_upper <- res$upper - res$fit
  halfwidth_lower <- res$fit - res$lower
  expect_equal(halfwidth_upper, halfwidth_lower, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# (6) Missing vcov -> NA CLs with a warning
# ---------------------------------------------------------------------------

test_that("missing vcov produces NA CLs with a warning", {
  skip_on_cran()
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)
  # Knock out vcov to simulate a non-converged / ill-conditioned fit
  fit$fit$vcov <- NULL
  expect_warning(
    res <- predict(fit, newdata = data.frame(time = 1, x = 0),
                   type = "survival", se.fit = TRUE),
    "Variance-covariance"
  )
  expect_true(all(is.na(res$se.fit)))
  expect_true(all(is.na(res$lower)))
  expect_true(all(is.na(res$upper)))
  expect_false(anyNA(res$fit))
})

# ---------------------------------------------------------------------------
# (7) Fixed-shapes multiphase recovers CLs from the free-parameter submatrix
# ---------------------------------------------------------------------------

test_that("multiphase CLs work with fixed shapes + CoE (free submatrix)", {
  skip_on_cran()
  df <- make_toy()
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  fit <- hazard(survival::Surv(time, status) ~ 1, data = df,
                dist = "multiphase", phases = phases, fit = TRUE,
                control = list(n_starts = 2))
  res <- predict(fit, newdata = data.frame(time = c(0.5, 1, 2)),
                 type = "survival", se.fit = TRUE)
  expect_true(all(is.finite(res$se.fit)))
  expect_true(all(res$se.fit > 0))
  expect_true(all(res$lower >= 0))
  expect_true(all(res$upper <= 1))
})

# ---------------------------------------------------------------------------
# per_phase Jacobian mode (decompose + se.fit support)
# ---------------------------------------------------------------------------

test_that("multiphase jacobian per_phase=TRUE returns blocks that sum to the total", {
  skip_on_cran()
  t_new <- c(0.3, 0.8, 1.5, 2.2)
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  cov_counts <- c(early = 0L, constant = 0L)
  x_list <- list(early = NULL, constant = NULL)
  theta <- c(-3, log(0.3), 1, 1, -2)
  p <- length(theta)

  J_total <- TemporalHazard:::.hzr_predict_jacobian_multiphase(
    theta, t_new, phases, cov_counts, x_list, p
  )
  J_list <- TemporalHazard:::.hzr_predict_jacobian_multiphase(
    theta, t_new, phases, cov_counts, x_list, p, per_phase = TRUE
  )

  expect_type(J_list, "list")
  expect_named(J_list, c("early", "constant"))
  expect_equal(dim(J_list$early), c(length(t_new), p))
  # Blocks sum to the aggregate Jacobian.
  expect_equal(Reduce(`+`, J_list), J_total, tolerance = 1e-12)
  # Each block touches only its own phase's columns:
  # early occupies cols 1:4 (log_mu, log_thalf, nu, m), constant col 5 (log_mu).
  expect_true(all(J_list$early[, 5] == 0))
  expect_true(all(J_list$constant[, 1:4] == 0))
})

test_that(".hzr_free_vcov restricts to finite-diagonal free parameters", {
  # All-free 2x2.
  V <- matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)
  res <- TemporalHazard:::.hzr_free_vcov(V, p = 2L)
  expect_equal(res$free_idx, 1:2)
  expect_equal(res$vcov_use, V)

  # One fixed parameter (NA row/col) -> restricted to the free submatrix.
  Vf <- matrix(NA_real_, 3, 3)
  Vf[c(1, 3), c(1, 3)] <- c(0.04, 0.00, 0.00, 0.09)
  res2 <- TemporalHazard:::.hzr_free_vcov(Vf, p = 3L)
  expect_equal(res2$free_idx, c(1L, 3L))
  expect_equal(dim(res2$vcov_use), c(2L, 2L))

  # Unusable vcov -> NULL with a warning.
  expect_warning(
    bad <- TemporalHazard:::.hzr_free_vcov(NULL, p = 3L),
    "unavailable"
  )
  expect_null(bad)
})

test_that(".hzr_predict_with_se_decomposed returns a tidy long frame", {
  skip_on_cran()
  df <- make_toy()
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  fit <- hazard(survival::Surv(time, status) ~ 1, data = df,
                dist = "multiphase", phases = phases, fit = TRUE,
                control = list(n_starts = 2))
  t_new <- c(0.5, 1, 2)
  cov_counts <- c(early = 0L, constant = 0L)
  x_list <- list(early = NULL, constant = NULL)

  res <- TemporalHazard:::.hzr_predict_with_se_decomposed(
    object = fit, time = t_new, x_list = x_list,
    cov_counts = cov_counts, phases = phases, level = 0.95
  )

  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "component", "fit", "se.fit", "lower", "upper"))
  expect_equal(levels(res$component), c("total", "early", "constant"))
  expect_equal(nrow(res), length(t_new) * 3L)  # total + 2 phases
  # Monotone, non-negative CLs everywhere SEs are finite.
  ok <- is.finite(res$se.fit)
  expect_true(all(res$lower[ok] <= res$fit[ok] + 1e-9))
  expect_true(all(res$fit[ok] <= res$upper[ok] + 1e-9))
  expect_true(all(res$lower[ok] >= -1e-12))
})

# ---------------------------------------------------------------------------
# (8) Backward compat: se.fit = FALSE reproduces the old scalar-vector return
# ---------------------------------------------------------------------------

test_that("se.fit = FALSE preserves pre-0.9.8 return shape", {
  skip_on_cran()
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)

  # Each type should return a bare numeric vector with se.fit = FALSE.
  for (type in c("survival", "cumulative_hazard", "hazard",
                 "linear_predictor")) {
    nd <- if (type %in% c("survival", "cumulative_hazard")) {
      data.frame(time = c(0.5, 1, 2), x = 0)
    } else {
      data.frame(x = c(-1, 0, 1))
    }
    out <- predict(fit, newdata = nd, type = type, se.fit = FALSE)
    expect_type(out, "double")
    expect_length(out, 3L)
    expect_false(is.data.frame(out))
  }
})

# ---------------------------------------------------------------------------
# conf.type for survival confidence limits (log-log default vs logit/HAZPRED)
# ---------------------------------------------------------------------------
test_that("predict(type='survival') conf.type selects the CL transform", {
  skip_if_not_installed("numDeriv")
  set.seed(7)
  n <- 400
  tev <- rexp(n, 0.3) + 0.01
  cens <- runif(n, 0.2, 6)
  d <- data.frame(time = pmin(tev, cens), status = as.integer(tev <= cens))
  set.seed(1)
  fit <- hazard(survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
                  constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1, maxit = 500, conserve = TRUE))
  grid <- data.frame(time = c(0.5, 1, 3))

  default <- predict(fit, newdata = grid, type = "survival", se.fit = TRUE)
  loglog  <- predict(fit, newdata = grid, type = "survival", se.fit = TRUE,
                     conf.type = "log-log")
  logit   <- predict(fit, newdata = grid, type = "survival", se.fit = TRUE,
                     conf.type = "logit")

  # Default is log-log.
  expect_equal(default$lower, loglog$lower)
  expect_equal(default$upper, loglog$upper)
  # logit differs from log-log but shares the point estimate and stays in (0,1).
  expect_equal(logit$fit, loglog$fit, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(logit$lower, loglog$lower)))
  expect_true(all(logit$lower > 0 & logit$upper < 1))
  expect_true(all(logit$lower <= logit$fit & logit$fit <= logit$upper))

  # Invalid conf.type is rejected.
  expect_error(predict(fit, newdata = grid, type = "survival", se.fit = TRUE,
                       conf.type = "bogus"))
})

# The SAS band level documented on predict() and hzr_nelson() (#204).
#
# SAS %KAPLAN / %NELSONT / HAZPRED default to CLEVEL = 0.68268948, documented
# in the macro source as "(1 sd)", which makes T_ALPHA 1 to seven decimals
# (the literal is truncated): a 68.3% band, not 95%. The help pages now tell
# a caller to pass `level = 2 * pnorm(1) - 1` to reproduce a SAS figure.
# These pin that recipe, so the documented incantation cannot quietly stop
# being the right one.

test_that("the documented SAS level is the 1-SD level", {
  lvl <- 2 * stats::pnorm(1) - 1
  # This is what SAS computes as T_ALPHA=PROBIT(0.5+(0.5*&CLEVEL)).
  expect_equal(stats::qnorm(0.5 + 0.5 * lvl), 1, tolerance = 1e-12)
  # The macro spells CLEVEL as a truncated literal, so its T_ALPHA is 1 only
  # to about seven decimals -- close enough to reproduce a plot, and worth
  # recording so "exactly 1" is not carried around as a fact.
  expect_equal(stats::qnorm(0.5 + 0.5 * 0.68268948), 1, tolerance = 1e-7)
  expect_lt(abs(lvl - 0.68268948), 1e-7)
})

test_that("the SAS level narrows a band by the z ratio, on the CL's own scale", {
  skip_on_cran()
  # NOT (upper - fit) / se.fit: predict() builds cumulative-hazard limits on
  # the LOG scale (and survival limits on log-log), so the half-width is not
  # z * se on the natural scale and asserting that fails against correct code.
  # The invariant that does hold is that the two levels differ by exactly the
  # ratio of their z values, measured on the scale the limits are built on.
  set.seed(4)
  n <- 200
  fit <- hazard(time = stats::rweibull(n, 1.3, 2), status = rep(1L, n),
                theta = c(0.4, 1.0), dist = "weibull", fit = TRUE)
  nd <- data.frame(time = c(0.5, 1, 2))

  sas <- predict(fit, newdata = nd, type = "cumulative_hazard",
                 se.fit = TRUE, level = 2 * stats::pnorm(1) - 1)
  d95 <- predict(fit, newdata = nd, type = "cumulative_hazard", se.fit = TRUE)

  # Guard the guard: a degenerate band makes every ratio below 0/0 or 1.
  expect_true(all(is.finite(sas$se.fit) & sas$se.fit > 0))
  expect_true(all(sas$upper > sas$fit))

  ratio <- log(d95$upper / d95$fit) / log(sas$upper / sas$fit)
  expect_equal(ratio, rep(stats::qnorm(0.975), nrow(nd)), tolerance = 1e-6)

  # And the plain statement the documentation makes: the default band is
  # materially wider than the SAS one, so reproducing a figure at level = 0.95
  # does not match it.
  expect_true(all(d95$upper > sas$upper))
})
