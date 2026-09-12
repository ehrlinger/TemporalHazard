# tests/testthat/test-gof-phase-covariates.R
# hzr_gof()'s par_cumhaz / par_surv curve is documented to sit at the
# covariate means.  A multiphase fit whose covariates enter only through the
# phase formulas stores no global design matrix (object$data$x is NULL), so a
# time-only newdata evaluated that curve at covariates = 0 instead.

.gof_pc_avc <- local({
  data(avc, package = "TemporalHazard")
  na.omit(avc)
})

test_that("phase-formula fit: par_cumhaz is at the covariate means, not at 0", {
  d <- .gof_pc_avc
  # The fit warns that its Hessian is singular.  That concerns standard
  # errors only, which hzr_gof() does not use.
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = d,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ age + mal),
      constant = hzr_phase("constant", formula = ~ age)
    ),
    fit    = TRUE
  ))
  expect_true(isTRUE(fit$fit$converged))
  expect_null(fit$data$x)

  gof <- hzr_gof(fit)
  nd_mean <- data.frame(time = gof$time, age = mean(d$age), mal = mean(d$mal))
  nd_zero <- data.frame(time = gof$time)
  at_mean <- predict(fit, newdata = nd_mean, type = "cumulative_hazard")
  at_zero <- predict(fit, newdata = nd_zero, type = "cumulative_hazard")

  # The two candidate curves are far apart, so the comparison can fail.
  expect_gt(max(abs(at_zero / at_mean - 1)), 0.5)

  expect_equal(unname(gof$par_cumhaz), unname(at_mean), tolerance = 1e-10)
  expect_equal(unname(gof$par_surv), unname(exp(-at_mean)), tolerance = 1e-10)

  decomp <- predict(fit, newdata = nd_mean, type = "cumulative_hazard",
                    decompose = TRUE)
  expect_equal(gof$par_cumhaz_early, decomp$early, tolerance = 1e-10)
  expect_equal(gof$par_cumhaz_constant, decomp$constant, tolerance = 1e-10)
})

test_that("global-covariate fit: par_cumhaz is at the design-matrix means", {
  d <- .gof_pc_avc
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ age + mal,
    data  = d,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5, beta_age = 0, beta_mal = 0),
    fit   = TRUE
  )
  expect_true(isTRUE(fit$fit$converged))

  gof <- hzr_gof(fit)
  nd_mean <- data.frame(time = gof$time, age = mean(d$age), mal = mean(d$mal))
  at_mean <- predict(fit, newdata = nd_mean, type = "cumulative_hazard")
  at_zero <- predict(fit, newdata = data.frame(time = gof$time, age = 0,
                                               mal = 0),
                     type = "cumulative_hazard")
  expect_gt(max(abs(at_zero / at_mean - 1)), 0.5)
  expect_equal(unname(gof$par_cumhaz), unname(at_mean), tolerance = 1e-10)
})

test_that("intercept-only multiphase fit: par_cumhaz is the one shared curve", {
  d <- .gof_pc_avc
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = d,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit    = TRUE
  )
  expect_true(isTRUE(fit$fit$converged))

  gof <- hzr_gof(fit)
  at_t <- predict(fit, newdata = data.frame(time = gof$time),
                  type = "cumulative_hazard")
  expect_equal(unname(gof$par_cumhaz), unname(at_t), tolerance = 1e-10)
  expect_gt(max(gof$par_cumhaz), 0.1)
})
