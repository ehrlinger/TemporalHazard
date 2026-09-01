# `.fit_weibull()` is defined in `helper-stepwise.R` and auto-sourced
# by testthat before any test-*.R files run.

# .hzr_aic -----------------------------------------------------------------

test_that(".hzr_aic matches -2*logLik + 2k for a single-distribution fit", {
  fit <- .fit_weibull()
  k <- length(fit$fit$theta)
  expected <- -2 * fit$fit$objective + 2 * k
  expect_equal(.hzr_aic(fit), expected)
})

test_that(".hzr_aic returns NA when logLik is not available", {
  fit <- .fit_weibull()
  fit$fit$objective <- NULL
  expect_true(is.na(.hzr_aic(fit)))
})

test_that(".hzr_aic counts only FREE parameters when fixed_mask is set", {
  # Multiphase fits expose fixed_mask where TRUE == FIXED (see
  # .hzr_optim_multiphase).  Simulate a fit with 8 theta entries of
  # which 3 are held fixed; k should be 5, not 8 (and not 3).
  fake <- structure(
    list(
      fit = list(
        theta      = setNames(seq_len(8), paste0("p", 1:8)),
        objective  = -100,
        fixed_mask = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
      )
    ),
    class = "hazard"
  )
  # 5 free params → AIC = -2*(-100) + 2*5 = 210
  expect_equal(.hzr_aic(fake), 210)
})

# Wald / entry -------------------------------------------------------------

test_that("Wald entry score equals candidate's p-value", {
  current <- hazard(
    time = rexp(200L, 0.5), status = rep(1L, 200L),
    theta = c(0.5, 1.0), dist = "weibull", fit = TRUE
  )
  candidate <- .fit_weibull(n = 200L, betas = 0.3, seed = 42L)

  covariate_row <- setdiff(
    rownames(summary(candidate)$coefficients), c("mu", "nu"))

  s <- .hzr_candidate_score(
    criterion = "wald", mode = "entry",
    current = current, candidate = candidate,
    names = covariate_row
  )

  wald <- .hzr_wald_p(candidate, covariate_row)
  expect_equal(s$score, wald$p_value)
  expect_equal(s$p_value, wald$p_value)
  expect_identical(s$criterion, "wald")
  expect_identical(s$mode, "entry")
  expect_true(is.na(s$delta_aic))   # no AIC produced in Wald-entry
})

# Wald / drop --------------------------------------------------------------

test_that("Wald drop score equals 1 - p_value (larger p => smaller score)", {
  fit <- .fit_weibull()
  covariate_row <- setdiff(
    rownames(summary(fit)$coefficients), c("mu", "nu"))

  s <- .hzr_candidate_score(
    criterion = "wald", mode = "drop",
    current = fit, names = covariate_row
  )

  wald <- .hzr_wald_p(fit, covariate_row)
  expect_equal(s$score, 1 - wald$p_value)
  expect_equal(s$p_value, wald$p_value)
  # ΔAIC is computed as a bonus for Wald-drop; check the formula
  expect_equal(s$delta_aic, wald$stat^2 - 2 * wald$df)
})

# AIC / entry --------------------------------------------------------------

test_that("AIC entry score equals AIC(candidate) - AIC(current)", {
  current <- hazard(
    time = rexp(200L, 0.5), status = rep(1L, 200L),
    theta = c(0.5, 1.0), dist = "weibull", fit = TRUE
  )
  candidate <- .fit_weibull(n = 200L, betas = 0.3, seed = 42L)
  covariate_row <- setdiff(
    rownames(summary(candidate)$coefficients), c("mu", "nu"))

  s <- .hzr_candidate_score(
    criterion = "aic", mode = "entry",
    current = current, candidate = candidate,
    names = covariate_row
  )

  expected_delta <- .hzr_aic(candidate) - .hzr_aic(current)
  expect_equal(s$score, expected_delta)
  expect_equal(s$delta_aic, expected_delta)
  expect_identical(s$criterion, "aic")
  # p_value is still populated for the candidate's new coef
  expect_true(is.finite(s$p_value))
})

# AIC / drop ---------------------------------------------------------------

test_that("AIC drop score uses Wald->LR approximation W - 2*df", {
  fit <- .fit_weibull()
  covariate_row <- setdiff(
    rownames(summary(fit)$coefficients), c("mu", "nu"))

  s <- .hzr_candidate_score(
    criterion = "aic", mode = "drop",
    current = fit, names = covariate_row
  )

  wald <- .hzr_wald_p(fit, covariate_row)
  expect_equal(s$score, wald$stat^2 - 2 * wald$df)
  expect_equal(s$delta_aic, s$score)
})

test_that("AIC drop: scalar case matches z^2 - 2 directly", {
  fit <- .fit_weibull(n = 500L, betas = 0.5)   # strong effect → large z
  covariate_row <- setdiff(
    rownames(summary(fit)$coefficients), c("mu", "nu"))
  wald <- .hzr_wald_p(fit, covariate_row)
  expect_equal(wald$df, 1L)

  s <- .hzr_candidate_score(
    criterion = "aic", mode = "drop",
    current = fit, names = covariate_row
  )
  expect_equal(s$score, wald$stat^2 - 2)
  # Strong effect should make dropping costly (score > 0)
  expect_gt(s$score, 0)
})

test_that("AIC drop: weak effect yields negative score (drop improves AIC)", {
  set.seed(123)
  n <- 300L
  x <- rnorm(n)
  # No effect: beta = 0
  time <- rexp(n, rate = 1)
  status <- rep(1L, n)
  fit <- hazard(
    time   = time,
    status = status,
    x      = matrix(x, ncol = 1L, dimnames = list(NULL, "x1")),
    theta  = c(0.5, 1.0, 0),
    dist   = "weibull",
    fit    = TRUE
  )

  covariate_row <- setdiff(
    rownames(summary(fit)$coefficients), c("mu", "nu"))

  s <- .hzr_candidate_score(
    criterion = "aic", mode = "drop",
    current = fit, names = covariate_row
  )
  # With no true effect and n = 300, z^2 << 2 on average; ΔAIC < 0.
  expect_lt(s$score, 0)
})

# Multi-df joint test ------------------------------------------------------

test_that("joint (df = 2) Wald drop uses chi^2 directly, not z^2", {
  fit <- .fit_weibull(n = 300L, betas = c(0.4, -0.3), seed = 7L)
  covariate_rows <- setdiff(
    rownames(summary(fit)$coefficients), c("mu", "nu"))
  expect_length(covariate_rows, 2L)

  s <- .hzr_candidate_score(
    criterion = "aic", mode = "drop",
    current = fit, names = covariate_rows
  )

  wald <- .hzr_wald_p(fit, covariate_rows)
  expect_equal(wald$df, 2L)
  # df > 1 path: W = stat (already chi^2), not stat^2
  expect_equal(s$score, wald$stat - 2 * 2)
})

# Error surface ------------------------------------------------------------

test_that("entry mode requires a candidate", {
  fit <- .fit_weibull()
  expect_error(
    .hzr_candidate_score(
      criterion = "wald", mode = "entry",
      current = fit, names = "beta1"
    ),
    "candidate"
  )
})

test_that("non-hazard current is rejected", {
  expect_error(
    .hzr_candidate_score(
      criterion = "wald", mode = "drop",
      current = list(), names = "x1"
    ),
    "`current` must be a fitted `hazard` object"
  )
})

test_that("unknown names propagate the Wald helper's error", {
  fit <- .fit_weibull()
  expect_error(
    .hzr_candidate_score(
      criterion = "wald", mode = "drop",
      current = fit, names = "nonexistent"
    ),
    "Unknown coefficient name"
  )
})

# NA propagation -----------------------------------------------------------

test_that("non-invertible vcov yields NA score without erroring", {
  fit <- .fit_weibull()
  fit$fit$vcov <- NULL    # simulate non-invertible Hessian
  covariate_row <- setdiff(
    rownames(summary(fit)$coefficients), c("mu", "nu"))

  s <- .hzr_candidate_score(
    criterion = "wald", mode = "drop",
    current = fit, names = covariate_row
  )
  expect_true(is.na(s$score))
  expect_true(is.na(s$p_value))
  expect_true(is.na(s$delta_aic))
})

# stat_type ----------------------------------------------------------------

test_that("stat_type names the distribution `stat` belongs to", {
  # `stat` is a Wald z for a scalar candidate and a chi-square for a joint
  # one, and a score Q is a chi-square too -- all three at df = 1 for the
  # scalar cases. `df` therefore cannot tell a caller which reference
  # distribution recomputes the p-value, which is what stat_type exists for.
  # Pinned here, against the scorer, rather than through a stepwise run: the
  # Wald path refits each candidate, and a multiphase refit is not stable
  # enough across platforms to hang a reporting assertion on.
  fit <- .fit_weibull()
  covariate_row <- setdiff(
    rownames(summary(fit)$coefficients), c("mu", "nu"))

  w <- .hzr_candidate_score(
    criterion = "wald", mode = "drop", current = fit, names = covariate_row)
  expect_identical(w$stat_type, "wald_z")
  expect_identical(w$df, 1L)
  # A scalar Wald is reported as a z, NOT as its square -- which is exactly
  # why it is indistinguishable from a 1-df Q on `df` alone.
  expect_equal(2 * stats::pnorm(-abs(w$stat)), w$p_value, tolerance = 1e-8)

  # The score branch, from a .hzr_score_q()-shaped result.
  q <- list(stat = 3.84, df = 1L, p_value = 0.05)
  sc <- .hzr_candidate_score(
    criterion = "score", mode = "entry", current = fit,
    score = q, names = covariate_row)
  expect_identical(sc$stat_type, "score_q")
  expect_identical(sc$df, 1L)

  # The two branches must not agree, or the field distinguishes nothing.
  expect_false(identical(w$stat_type, sc$stat_type))
})

test_that("a multi-df Wald is labelled chi-square, not z", {
  # The joint case: .hzr_candidate_score() reports the chi-square directly
  # rather than a z, so stat_type has to switch with it.
  fit <- .fit_weibull(betas = c(0.3, -0.2), seed = 7L)
  covariate_rows <- setdiff(
    rownames(summary(fit)$coefficients), c("mu", "nu"))
  skip_if(length(covariate_rows) < 2L, "fixture has fewer than 2 covariates")

  w <- .hzr_candidate_score(
    criterion = "wald", mode = "drop", current = fit, names = covariate_rows)
  expect_gt(w$df, 1L)
  expect_identical(w$stat_type, "wald_chisq")
})
