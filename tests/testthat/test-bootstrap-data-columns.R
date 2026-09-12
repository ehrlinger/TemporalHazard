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
