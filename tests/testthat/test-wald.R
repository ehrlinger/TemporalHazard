# `.fit_weibull()` is defined in `helper-stepwise.R` and auto-sourced
# by testthat before any test-*.R files run.

test_that("scalar Wald p-value matches summary.hazard", {
  fit <- .fit_weibull()

  s <- summary(fit)
  # Find the single covariate row in the summary table
  covariate_row <- setdiff(rownames(s$coefficients), c("mu", "nu"))
  expect_length(covariate_row, 1L)

  wald <- .hzr_wald_p(fit, covariate_row)

  expect_equal(wald$df, 1L)
  expect_equal(wald$stat,     s$coefficients[covariate_row, "z_stat"],
               tolerance = 1e-10)
  expect_equal(wald$p_value,  s$coefficients[covariate_row, "p_value"],
               tolerance = 1e-10)
  expect_equal(wald$estimate, unname(s$coefficients[covariate_row, "estimate"]),
               tolerance = 1e-10)
})

test_that("joint Wald chi-square matches hand calculation on 2 coefficients", {
  fit <- .fit_weibull(n = 300L, betas = c(0.4, -0.3), seed = 7L)

  s <- summary(fit)
  covariate_rows <- setdiff(rownames(s$coefficients), c("mu", "nu"))
  expect_length(covariate_rows, 2L)

  w <- .hzr_wald_p(fit, covariate_rows)

  # Hand calc: W = b' V^{-1} b.  fit$fit$theta / vcov are unnamed for
  # single-distribution fits, so index by position via match() against
  # the canonical names summary() prints.
  all_names <- rownames(s$coefficients)
  idx <- match(covariate_rows, all_names)
  coef_vec <- fit$fit$theta[idx]
  V <- fit$fit$vcov[idx, idx]
  expected_stat <- as.numeric(crossprod(coef_vec, solve(V) %*% coef_vec))
  expected_p <- pchisq(expected_stat, df = 2, lower.tail = FALSE)

  expect_equal(w$df, 2L)
  expect_equal(w$stat, expected_stat, tolerance = 1e-10)
  expect_equal(w$p_value, expected_p, tolerance = 1e-10)
})

test_that("scalar z^2 equals the df=1 chi-square via pchisq", {
  fit <- .fit_weibull(n = 150L, betas = 0.5, seed = 99L)

  covariate_row <- setdiff(rownames(summary(fit)$coefficients),
                           c("mu", "nu"))

  w <- .hzr_wald_p(fit, covariate_row)
  expect_equal(
    pchisq(w$stat^2, df = 1L, lower.tail = FALSE),
    w$p_value,
    tolerance = 1e-10
  )
})

test_that("unknown coefficient name produces an informative error", {
  fit <- .fit_weibull(n = 100L, seed = 11L)
  expect_error(
    .hzr_wald_p(fit, "not_in_model"),
    "Unknown coefficient name"
  )
})

test_that("no vcov on the fit returns p_value = NA", {
  fit <- .fit_weibull(n = 80L, seed = 5L)
  covariate_row <- setdiff(rownames(summary(fit)$coefficients),
                           c("mu", "nu"))

  fit$fit$vcov <- NULL

  w <- .hzr_wald_p(fit, covariate_row)
  expect_true(is.na(w$p_value))
  expect_equal(w$df, 1L)
  expect_identical(w$names, covariate_row)
})

test_that("non-finite variance (boundary parameter) returns p_value = NA", {
  fit <- .fit_weibull(n = 80L, seed = 5L)
  s <- summary(fit)
  covariate_row <- setdiff(rownames(s$coefficients), c("mu", "nu"))

  # vcov for single-dist fits has no dimnames; poison by positional
  # index
  idx <- match(covariate_row, rownames(s$coefficients))
  V <- fit$fit$vcov
  V[idx, idx] <- NA_real_
  fit$fit$vcov <- V

  w <- .hzr_wald_p(fit, covariate_row)
  expect_true(is.na(w$p_value))
})

test_that("works on a multiphase fit with phase-prefixed coefficient names", {
  data(avc)
  set.seed(1)
  fit <- hazard(
    Surv(int_dead, dead) ~ age,
    data = avc,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE,
    control = list(n_starts = 3L, maxit = 1000L)
  )

  coef_names <- names(coef(fit))
  expect_true(!is.null(coef_names))
  age_name <- grep("age", coef_names, value = TRUE)
  expect_true(length(age_name) >= 1L)

  w <- .hzr_wald_p(fit, age_name[1])
  expect_equal(w$df, 1L)
  # p_value is either a finite probability in [0, 1] or NA
  expect_true(is.na(w$p_value) ||
                (is.finite(w$p_value) && w$p_value >= 0 && w$p_value <= 1))
})

test_that("non-hazard object is rejected", {
  expect_error(.hzr_wald_p(list(), "x1"), "must be a `hazard` object")
})

test_that("empty or invalid names argument is rejected", {
  fit <- .fit_weibull(n = 50L, seed = 1L)
  expect_error(.hzr_wald_p(fit, character()), "non-empty character vector")
  expect_error(.hzr_wald_p(fit, c("mu", "")), "non-empty character vector")
  expect_error(.hzr_wald_p(fit, 1), "non-empty character vector")
})

test_that("a theta/design mismatch agrees its article with the dist name", {
  # The message read "but a exponential fit" (ledger item 2). It can become a
  # replicate's failure-reason key under hzr_bootstrap(scope = ), so the fix is
  # vowel-aware: a consonant-initial dist keeps its exact text.
  set.seed(1)
  d <- data.frame(time = rexp(60, 0.3) + 0.01, status = rbinom(60, 1, 0.7),
                  age = rnorm(60))
  msg <- function(dist, theta) {
    fit <- suppressWarnings(hazard(survival::Surv(time, status) ~ age, data = d,
                                   dist = dist, theta = theta, fit = TRUE))
    fit$fit$theta <- c(fit$fit$theta, 0)
    tryCatch(.hzr_wald_p(fit, "age"), error = conditionMessage)
  }
  expect_match(msg("exponential", c(0.3, 0)),
               "but an exponential fit with 1 covariate column(s) has 2 coefficients",
               fixed = TRUE)
  expect_match(msg("weibull", c(0.3, 1, 0)),
               "but a weibull fit with 1 covariate column(s) has 3 coefficients",
               fixed = TRUE)
})
