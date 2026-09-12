# tests/testthat/test-bootstrap-data-columns.R
# hzr_bootstrap() resamples the rows of the fit's `data`. A formula variable
# that is not a column of `data` was never resampled, so each replicate paired
# the original vector with resampled rows: a wrong interval, with every
# replicate reported as a success and no warning (#278). It is now refused.

avc_278 <- na.omit(avc[, c("int_dead", "dead", "age", "mal")])

weibull_278 <- function(formula, theta = c(0.3, 1, 0)) {
  hazard(formula, data = avc_278, dist = "weibull", theta = theta,
         fit = TRUE)
}

test_that("a function-local covariate is refused, not bootstrapped", {
  # The 40-replicate case from the #273 review: 40 successes, an interval
  # that excluded its own estimate, and no warning.
  local_fit <- function() {
    zz <- avc_278$age
    weibull_278(survival::Surv(int_dead, dead) ~ zz)
  }
  fit <- local_fit()
  expect_true(fit$fit$converged)
  expect_error(hzr_bootstrap(fit, n_boot = 40, seed = 1),
               "uses 'zz', which is not a column.*Add it to `data`")
})

test_that("a workspace covariate is refused, and only it is named", {
  assign("zz_ws_278", avc_278$age, envir = globalenv())
  withr::defer(rm("zz_ws_278", envir = globalenv()))
  fit <- weibull_278(survival::Surv(int_dead, dead) ~ zz_ws_278 + mal,
                     theta = c(0.3, 1, 0, 0))
  expect_error(hzr_bootstrap(fit, n_boot = 5, seed = 1),
               "uses 'zz_ws_278', which is not a column")
})

test_that("a phase-formula covariate outside `data` is refused", {
  local_mp <- function() {
    zz <- avc_278$age
    hazard(survival::Surv(int_dead, dead) ~ 1, data = avc_278,
           dist = "multiphase",
           phases = list(
             early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                  fixed = "m", formula = ~ zz + mal),
             constant = hzr_phase("constant")
           ),
           fit = TRUE, control = list(n_starts = 1L, conserve = FALSE))
  }
  expect_error(hzr_bootstrap(suppressWarnings(local_mp()), n_boot = 2,
                             seed = 1),
               "uses 'zz', which is not a column")
})

test_that("a response variable outside `data` is refused", {
  # The times held fixed while `dead` and `age` were resampled. A workspace
  # variable: the Surv() term is evaluated with `data` in front of the
  # search path, so a function-local one does not reach the fit at all.
  assign("tt_ws_278", avc_278$int_dead, envir = globalenv())
  withr::defer(rm("tt_ws_278", envir = globalenv()))
  fit <- weibull_278(survival::Surv(tt_ws_278, dead) ~ age)
  expect_error(hzr_bootstrap(fit, n_boot = 5, seed = 1),
               "uses 'tt_ws_278', which is not a column")
})

test_that("the response is checked where the fit evaluates it", {
  # The Surv() term is evaluated in `data`, then the search path, never in
  # the formula's environment, so a local of the same name must not hide
  # the per-row vector the fit actually used.
  assign("tt_g_278", avc_278$int_dead, envir = globalenv())
  withr::defer(rm("tt_g_278", envir = globalenv()))
  shadowed <- function() {
    tt_g_278 <- 5 # nolint: object_usage_linter.
    weibull_278(survival::Surv(tt_g_278, dead) ~ age)
  }
  expect_error(hzr_bootstrap(shadowed(), n_boot = 3, seed = 1),
               "uses 'tt_g_278', which is not a column")
})

test_that("a select-mode scope variable outside `data` is refused", {
  # Candidate refits resolve a scope's variables where the base model's
  # formula was written, so `zz` here does reach the refits.
  local_base <- function() {
    zz <- avc_278$age # nolint: object_usage_linter.
    weibull_278(survival::Surv(int_dead, dead) ~ 1, theta = c(0.3, 1))
  }
  base <- local_base()
  for (sc in list(~ zz + mal, "zz", "log(zz)")) {
    expect_error(hzr_bootstrap(base, n_boot = 5, seed = 1, scope = sc,
                               criterion = "wald"),
                 "uses 'zz', which is not a column")
  }
})

test_that("a multiphase scope is checked where its phase refit reads it", {
  # A phase's scope is refit into that phase's own formula, so a per-row
  # vector written beside the phase formula reaches the refit even though
  # the base formula cannot see it. A phase with no formula gets a fresh
  # one whose lookups reach the search path.
  mk_phases <- function() {
    zz <- avc_278$age + sin(seq_len(nrow(avc_278))) # nolint: object_usage_linter.
    list(
      early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                           fixed = "m", formula = ~ mal),
      constant = hzr_phase("constant")
    )
  }
  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc_278, dist = "multiphase",
    phases = mk_phases(), fit = TRUE,
    control = list(n_starts = 1L, conserve = FALSE)
  ))
  expect_error(
    hzr_bootstrap(base, n_boot = 3, seed = 1, scope = list(early = ~ zz),
                  criterion = "wald", slentry = 0.99),
    "uses 'zz', which is not a column"
  )
  assign("zz_mp_278", avc_278$age, envir = globalenv())
  withr::defer(rm("zz_mp_278", envir = globalenv()))
  expect_error(
    hzr_bootstrap(base, n_boot = 3, seed = 1,
                  scope = list(constant = ~ zz_mp_278),
                  criterion = "wald", slentry = 0.99),
    "uses 'zz_mp_278', which is not a column"
  )
})

test_that("a scope variable only the scope's own frame can see is refused", {
  # The refit cannot resolve `zl`, so the screen could never test it: after
  # one warning up front, it silently dropped out of every replicate's
  # candidates.
  mk_scope <- function() {
    zl <- avc_278$age # nolint: object_usage_linter.
    ~ zl + mal
  }
  base <- weibull_278(survival::Surv(int_dead, dead) ~ 1, theta = c(0.3, 1))
  expect_error(hzr_bootstrap(base, n_boot = 3, seed = 1, scope = mk_scope(),
                             criterion = "wald"),
               "uses 'zl', which is not a column")

  # The same for a phase with no formula of its own.
  mp <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc_278, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                           fixed = "m"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = list(n_starts = 1L, conserve = FALSE)
  ))
  mk_mscope <- function() {
    zl <- avc_278$age # nolint: object_usage_linter.
    list(constant = ~ zl)
  }
  expect_error(hzr_bootstrap(mp, n_boot = 2, seed = 1, scope = mk_mscope(),
                             criterion = "wald", slentry = 0.99),
               "uses 'zl', which is not a column")
})

test_that("a vector-interface fit with a direct `x` is refused", {
  # `x` was re-evaluated unresampled in every replicate: 20 of 20 succeeded,
  # and the interval for `age` excluded its own estimate.
  fit <- hazard(data = avc_278, time = int_dead, status = dead,
                x = as.matrix(avc_278["age"]), dist = "weibull",
                theta = c(0.3, 1, 0), fit = TRUE)
  # Refused before seeding, so the caller's random number stream is left
  # alone.
  set.seed(42)
  before <- .Random.seed
  expect_error(hzr_bootstrap(fit, n_boot = 3, seed = 1),
               "passed directly as `x`")
  expect_identical(.Random.seed, before)

  # Select mode too: the base refits reused the original `x` and the
  # candidate refits dropped it, so the `age` terms vanished from every
  # replicate while all of them reported success.
  mp <- suppressWarnings(hazard(
    data = avc_278, time = int_dead, status = dead,
    x = as.matrix(avc_278["age"]), dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                           fixed = "m"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = list(n_starts = 1L, conserve = FALSE)
  ))
  expect_error(hzr_bootstrap(mp, n_boot = 3, seed = 1,
                             scope = list(early = ~ mal),
                             criterion = "wald", slentry = 0.99),
               "passed directly as `x`")
})

test_that("constants outside `data` are not refused", {
  # One value, or a few, is the same in every replicate by design: only a
  # per-row vector is data the resample would miss.
  cutoff <- 60
  kn <- c(40, 60)
  fits <- list(
    weibull_278(survival::Surv(int_dead, dead) ~ I(age * pi)),
    weibull_278(survival::Surv(int_dead, dead) ~ I(age > cutoff)),
    weibull_278(survival::Surv(int_dead, dead) ~ splines::ns(age, knots = kn),
                theta = c(0.3, 1, 0, 0, 0))
  )
  for (fit in fits) {
    b <- suppressWarnings(hzr_bootstrap(fit, n_boot = 4, seed = 1))
    expect_equal(b$n_success, 4)
  }
})

test_that("a bootstrap over columns of `data` still runs and varies", {
  # `log(age)` reads the column `age`: the check works on the variables of
  # the terms, not on the term labels.
  fit <- weibull_278(survival::Surv(int_dead, dead) ~ log(age) + mal,
                     theta = c(0.3, 1, 0, 0))
  b <- suppressWarnings(hzr_bootstrap(fit, n_boot = 10, seed = 1))
  expect_equal(b$n_success, 10)
  expect_equal(b$n_failed, 0)
  row <- b$summary[b$summary$parameter == "log(age)", ]
  expect_equal(nrow(row), 1L)
  expect_gt(row$sd, 0)
})
