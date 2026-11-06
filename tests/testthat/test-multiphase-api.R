# Too slow for the CRAN check time budget; runs on CI and locally.
# See "Where the check's time actually goes" in AGENTS.md.
skip_on_cran()

# test-multiphase-api.R — End-to-end tests for multiphase hazard via hazard()
#
# Tests construction, fitting, predict, summary, and print for dist="multiphase".

# ============================================================================
# Helpers
# ============================================================================

make_multiphase_data <- function(n = 200, seed = 42) {
  set.seed(seed)
  data.frame(
    time   = rexp(n, rate = 0.3) + 0.01,
    status = rbinom(n, 1, prob = 0.7),
    age    = rnorm(n, mean = 60, sd = 10),
    nyha   = sample(1:4, n, replace = TRUE)
  )
}


# ============================================================================
# hazard() construction with phases
# ============================================================================

test_that("hazard() accepts dist='multiphase' with phases", {
  dat <- make_multiphase_data(n = 50)

  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ age,
    data   = dat,
    dist   = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
      bg    = hzr_phase("constant")
    ),
    fit = TRUE
  ))

  expect_s3_class(fit, "hazard")
  expect_equal(fit$spec$dist, "multiphase")
  expect_true(!is.null(fit$spec$phases))
  expect_equal(length(fit$spec$phases), 2)
  expect_equal(names(fit$spec$phases), c("early", "bg"))
})

test_that("hazard() errors when dist='multiphase' without phases", {
  dat <- make_multiphase_data(n = 50)
  expect_error(
    hazard(survival::Surv(time, status) ~ age, data = dat,
           dist = "multiphase", fit = TRUE),
    "phases.*required"
  )
})

test_that("hazard() warns when phases given with non-multiphase dist", {
  dat <- make_multiphase_data(n = 50)
  expect_warning(
    hazard(survival::Surv(time, status) ~ age, data = dat,
           dist = "weibull", theta = c(0.3, 1, 0),
           phases = list(hzr_phase("cdf")),
           fit = FALSE),
    "phases.*ignored"
  )
})


# ============================================================================
# Single-phase (constant) multiphase — should recover exponential-like fit
# ============================================================================

test_that("single constant phase fits and converges", {
  skip_on_cran()

  set.seed(99)
  dat <- data.frame(
    time   = rexp(200, rate = 0.3) + 0.01,
    status = rep(1L, 200)
  )

  fit <- hazard(
    time = dat$time, status = dat$status,
    dist = "multiphase",
    phases = list(bg = hzr_phase("constant")),
    fit = TRUE
  )

  expect_true(fit$fit$converged)
  expect_true(is.finite(fit$fit$objective))

  # theta should have bg.log_mu
  expect_true("bg.log_mu" %in% names(fit$fit$theta))
})


# ============================================================================
# Two-phase fit
# ============================================================================

test_that("two-phase (cdf + constant) fits and stores phases", {
  skip_on_cran()

  dat <- make_multiphase_data(n = 150)

  fit <- hazard(
    survival::Surv(time, status) ~ 1,
    data   = dat,
    dist   = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
      bg    = hzr_phase("constant")
    ),
    fit = TRUE,
    control = list(n_starts = 2, maxit = 500)
  )

  expect_s3_class(fit, "hazard")
  expect_true(is.finite(fit$fit$objective))

  # Phases metadata should be stored
  expect_true(!is.null(fit$fit$phases))
  expect_equal(names(fit$fit$phases), c("early", "bg"))
  expect_true(!is.null(fit$fit$covariate_counts))
})


# ============================================================================
# predict.hazard() for multiphase
# ============================================================================

test_that("predict survival and cumhaz for multiphase model", {
  skip_on_cran()

  dat <- make_multiphase_data(n = 100)

  fit <- hazard(
    survival::Surv(time, status) ~ 1,
    data   = dat,
    dist   = "multiphase",
    phases = list(bg = hzr_phase("constant")),
    fit = TRUE,
    control = list(n_starts = 2, maxit = 500)
  )

  # Predict on training data
  surv <- predict(fit, type = "survival")
  cumhaz <- predict(fit, type = "cumulative_hazard")

  expect_true(all(is.finite(surv)))
  expect_true(all(surv >= 0 & surv <= 1))
  expect_true(all(is.finite(cumhaz)))
  expect_true(all(cumhaz >= 0))
  expect_equal(surv, exp(-cumhaz), tolerance = 1e-10)
})

test_that("predict with newdata for multiphase model", {
  skip_on_cran()

  dat <- make_multiphase_data(n = 100)

  fit <- hazard(
    survival::Surv(time, status) ~ 1,
    data   = dat,
    dist   = "multiphase",
    phases = list(bg = hzr_phase("constant")),
    fit = TRUE,
    control = list(n_starts = 2, maxit = 500)
  )

  nd <- data.frame(time = c(0.5, 1.0, 2.0, 5.0))
  surv <- predict(fit, newdata = nd, type = "survival")

  expect_equal(length(surv), 4)
  expect_true(all(surv >= 0 & surv <= 1))
  # Survival should be monotone decreasing
  expect_true(all(diff(surv) <= 0))
})

test_that("predict decompose=TRUE returns data frame", {
  skip_on_cran()

  dat <- make_multiphase_data(n = 100)

  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1,
    data   = dat,
    dist   = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
      bg    = hzr_phase("constant")
    ),
    fit = TRUE,
    control = list(n_starts = 2, maxit = 500)
  ))

  nd <- data.frame(time = c(0.5, 1.0, 2.0))
  decomp <- predict(fit, newdata = nd, type = "cumulative_hazard",
                    decompose = TRUE)

  expect_true(is.data.frame(decomp))
  expect_true("time" %in% names(decomp))
  expect_true("total" %in% names(decomp))
  expect_true("early" %in% names(decomp))
  expect_true("bg" %in% names(decomp))

  # Total should equal sum of components
  expect_equal(decomp$total, decomp$early + decomp$bg, tolerance = 1e-10)
})

test_that("predict linear_predictor errors for multiphase", {
  skip_on_cran()

  dat <- make_multiphase_data(n = 50)
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1,
    data   = dat,
    dist   = "multiphase",
    phases = list(bg = hzr_phase("constant")),
    fit = TRUE,
    control = list(n_starts = 1, maxit = 300)
  ))

  expect_error(predict(fit, type = "linear_predictor"), "not supported.*multiphase")
})


# ============================================================================
# print and summary
# ============================================================================

test_that("print.hazard shows phase info for multiphase", {
  dat <- make_multiphase_data(n = 50)
  fit <- suppressWarnings(hazard(
    time = dat$time, status = dat$status,
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf"),
      bg    = hzr_phase("constant")
    ),
    fit = TRUE,
    control = list(n_starts = 1, maxit = 300)
  ))

  out <- capture.output(print(fit))
  expect_true(any(grepl("multiphase", out)))
  expect_true(any(grepl("phases.*2", out)))
})

test_that("summary.hazard shows per-phase coefficients", {
  skip_on_cran()

  dat <- make_multiphase_data(n = 100)
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1,
    data   = dat,
    dist   = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
      bg    = hzr_phase("constant")
    ),
    fit = TRUE,
    control = list(n_starts = 2, maxit = 500)
  ))

  s <- summary(fit)
  expect_s3_class(s, "summary.hazard")
  expect_equal(s$dist, "multiphase")
  expect_true(!is.null(s$phases))
  expect_true(!is.null(s$coefficients))

  # Should have rows for both phases
  rn <- rownames(s$coefficients)
  expect_true(any(grepl("^early\\.", rn)))
  expect_true(any(grepl("^bg\\.", rn)))

  # Print should work without error
  out <- capture.output(print(s))
  expect_true(any(grepl("Multiphase", out)))
  expect_true(any(grepl("Phase.*early", out)))
  expect_true(any(grepl("Phase.*bg", out)))
})


# ============================================================================
# coef and vcov
# ============================================================================

test_that("coef returns named theta for multiphase", {
  dat <- make_multiphase_data(n = 50)
  fit <- suppressWarnings(hazard(
    time = dat$time, status = dat$status,
    dist = "multiphase",
    phases = list(bg = hzr_phase("constant")),
    fit = TRUE,
    control = list(n_starts = 1, maxit = 300)
  ))

  cf <- coef(fit)
  expect_true(!is.null(cf))
  expect_true(!is.null(names(cf)))
  expect_true("bg.log_mu" %in% names(cf))
})


# ============================================================================
# Three-phase classic pattern
# ============================================================================

test_that("three-phase (cdf + constant + cdf) constructs and fits", {
  skip_on_cran()

  dat <- make_multiphase_data(n = 200)

  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1,
    data   = dat,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf",      t_half = 0.5, nu = 2, m = 0),
      constant = hzr_phase("constant"),
      late     = hzr_phase("cdf",      t_half = 5,   nu = 1, m = 0)
    ),
    fit = TRUE,
    control = list(n_starts = 2, maxit = 500)
  ))

  expect_s3_class(fit, "hazard")
  expect_true(is.finite(fit$fit$objective))
  expect_equal(length(fit$spec$phases), 3)

  # Predictions should work
  nd <- data.frame(time = c(0.5, 2, 8))
  surv <- predict(fit, newdata = nd, type = "survival")
  expect_equal(length(surv), 3)
  expect_true(all(surv >= 0 & surv <= 1))
})
