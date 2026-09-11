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

test_that("multiphase fit: one par_cumhaz_<phase> column per phase, no others", {
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
  gof <- hzr_gof(fit)
  # decompose = TRUE also returns `time` and `total`; neither is a phase.
  expect_identical(grep("^par_cumhaz_", names(gof), value = TRUE),
                   c("par_cumhaz_early", "par_cumhaz_constant"))
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

test_that("multiphase fit with only global covariates: curve at their means", {
  d <- .gof_pc_avc
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age,
    data   = d,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit    = TRUE
  ))
  expect_true(isTRUE(fit$fit$converged))

  gof <- hzr_gof(fit)
  nd_mean <- data.frame(time = gof$time, age = mean(d$age))
  at_mean <- predict(fit, newdata = nd_mean, type = "cumulative_hazard")
  at_zero <- predict(fit, newdata = data.frame(time = gof$time, age = 0),
                     type = "cumulative_hazard")
  expect_gt(max(abs(at_zero / at_mean - 1)), 0.5)
  expect_equal(unname(gof$par_cumhaz), unname(at_mean), tolerance = 1e-10)
})

test_that("multiphase fit with time_windows: curve switches windows with time", {
  # time_windows expands the global x into one column per window, each on
  # only while the row's time is in that window. Column means of the
  # expanded matrix would switch every window on at once. The mean patient
  # has mean(age) in the window active at each grid time.
  d <- .gof_pc_avc
  fit <- suppressWarnings(hazard(
    time = d$int_dead, status = d$dead, x = cbind(age = d$age),
    time_windows = 1,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit    = TRUE
  ))
  expect_true(isTRUE(fit$fit$converged))
  expect_identical(colnames(fit$fit$x_list$constant), c("age_w1", "age_w2"))

  gof <- hzr_gof(fit)
  n <- length(gof$time)
  x_mean <- .hzr_expand_time_varying_design(
    x = matrix(mean(d$age), nrow = n, ncol = 1,
               dimnames = list(NULL, "age")),
    time = gof$time, time_windows = 1
  )
  # Both windows are visited, so a wrong window would show.
  expect_true(any(gof$time <= 1) && any(gof$time > 1))
  ref <- .hzr_multiphase_cumhaz(
    gof$time, fit$fit$theta, fit$fit$phases, fit$fit$covariate_counts,
    list(early = x_mean, constant = x_mean)
  )
  expect_equal(unname(gof$par_cumhaz), ref, tolerance = 1e-10)
})

test_that("global and phase covariates together: curve at both sets of means", {
  # age enters globally (so data$x is set and reaches the constant phase);
  # mal enters only through the early phase's formula.  A newdata built from
  # data$x alone has no mal column, and predict() needs it for that formula.
  d <- .gof_pc_avc
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age,
    data   = d,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ mal),
      constant = hzr_phase("constant")
    ),
    fit     = TRUE,
    control = list(conserve = TRUE)
  ))
  expect_true(isTRUE(fit$fit$converged))
  expect_true(isTRUE(fit$spec$control$conserve_applied))
  expect_identical(colnames(fit$data$x), "age")
  expect_identical(colnames(fit$fit$x_list$early), "mal")

  gof <- hzr_gof(fit)
  # The per-subject tally is reached, and conserves events.
  s <- attr(gof, "summary")
  expect_equal(s$total_observed, sum(d$dead))
  expect_equal(s$total_expected / s$total_observed, 1, tolerance = 1e-6)

  # predict() with a newdata cannot evaluate this fit, so build the reference
  # from the model's structure instead: within a phase the covariates scale
  # the baseline, H_j(t | x) = exp(x beta_j) H0_j(t). A time-only newdata
  # gives each phase's baseline H0_j.
  base <- predict(fit, newdata = data.frame(time = gof$time),
                  type = "cumulative_hazard", decompose = TRUE)
  th <- fit$fit$theta
  early_ref <- base$early * exp(th[["early.mal"]] * mean(d$mal))
  constant_ref <- base$constant * exp(th[["constant.age"]] * mean(d$age))
  expect_gt(max(abs(base$total / (early_ref + constant_ref) - 1)), 0.1)

  expect_equal(gof$par_cumhaz_early, early_ref, tolerance = 1e-10)
  expect_equal(gof$par_cumhaz_constant, constant_ref, tolerance = 1e-10)
  expect_equal(unname(gof$par_cumhaz), early_ref + constant_ref,
               tolerance = 1e-10)
  expect_identical(grep("^par_cumhaz_", names(gof), value = TRUE),
                   c("par_cumhaz_early", "par_cumhaz_constant"))
})
