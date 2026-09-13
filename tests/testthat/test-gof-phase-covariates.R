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

# The mean patient's curve under time_windows = 1, built by hand: each phase's
# baseline (a time-only newdata) times exp(beta * mean(age)), using the w1
# coefficient up to time 1 and the w2 coefficient after it.
.gof_pc_window_ref <- function(fit, gof, d) {
  base <- predict(fit, newdata = data.frame(time = gof$time),
                  type = "cumulative_hazard", decompose = TRUE)
  th <- fit$fit$theta
  beta <- function(ph) {
    ifelse(gof$time <= 1, th[[paste0(ph, ".age_w1")]],
           th[[paste0(ph, ".age_w2")]])
  }
  base$early * exp(beta("early") * mean(d$age)) +
    base$constant * exp(beta("constant") * mean(d$age))
}

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
  # Both windows are visited, so a wrong window would show.
  expect_true(any(gof$time <= 1) && any(gof$time > 1))
  expect_equal(unname(gof$par_cumhaz), .gof_pc_window_ref(fit, gof, d),
               tolerance = 1e-10)
})

test_that("time_windows: a phase formula ignored at fit time gets the windows", {
  # Without `data`, the fit cannot evaluate a phase formula, so the phase
  # takes the window-expanded global x like any other. The curve has to
  # follow what the fit built, not whether a formula was written.
  d <- .gof_pc_avc
  fit <- suppressWarnings(hazard(
    time = d$int_dead, status = d$dead, x = cbind(age = d$age),
    time_windows = 1,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ age),
      constant = hzr_phase("constant")
    ),
    fit    = TRUE
  ))
  expect_identical(colnames(fit$fit$x_list$early), c("age_w1", "age_w2"))
  expect_false(is.null(fit$fit$phases$early$formula))

  gof <- hzr_gof(fit)
  expect_equal(unname(gof$par_cumhaz), .gof_pc_window_ref(fit, gof, d),
               tolerance = 1e-10)
})

test_that("a phase covariate with missing values is refused, not recycled", {
  # The fit drops rows with a missing phase covariate from that phase's
  # design matrix but keeps every row's time and status, so a per-subject
  # prediction would recycle the shorter design across the longer data.
  d <- .gof_pc_avc
  d$mal[1:5] <- NA
  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                         fixed = "shapes", formula = ~ mal),
    constant = hzr_phase("constant")
  )
  for (f in list(survival::Surv(int_dead, dead) ~ age,
                 survival::Surv(int_dead, dead) ~ 1)) {
    fit <- suppressWarnings(hazard(f, data = d, dist = "multiphase",
                                   phases = phases, fit = TRUE))
    expect_identical(nrow(fit$fit$x_list$early), nrow(d) - 5L)
    expect_error(hzr_gof(fit), "one design row per subject")
  }
})

test_that("time_windows: a phase formula whose columns share the window names", {
  # The global age expands to age_w1, age_w2; this phase formula uses data
  # columns of the same names. The fit evaluates the formula, so the phase
  # sits at its own columns' means, not at the windowed global age.
  d <- .gof_pc_avc
  set.seed(1)
  d$age_w1 <- stats::rnorm(nrow(d), 5, 1)
  d$age_w2 <- stats::rnorm(nrow(d), 50, 5)
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age,
    data = d, time_windows = 1, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ age_w1 + age_w2),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  # The collision: both phases' designs carry the same column names.
  expect_identical(colnames(fit$fit$x_list$early),
                   colnames(fit$fit$x_list$constant))

  gof <- hzr_gof(fit)
  base <- predict(fit, newdata = data.frame(time = gof$time),
                  type = "cumulative_hazard", decompose = TRUE)
  th <- fit$fit$theta
  early_ref <- base$early * exp(th[["early.age_w1"]] * mean(d$age_w1) +
                                  th[["early.age_w2"]] * mean(d$age_w2))
  expect_equal(gof$par_cumhaz_early, early_ref, tolerance = 1e-10)
  # The inherited phase still switches windows with time.
  beta_c <- ifelse(gof$time <= 1, th[["constant.age_w1"]],
                   th[["constant.age_w2"]])
  expect_equal(gof$par_cumhaz_constant,
               base$constant * exp(beta_c * mean(d$age)), tolerance = 1e-10)
})

test_that("#263: a phase literally named `time` keeps its column", {
  # predict(decompose = TRUE) names its grid column `time`, so a phase named
  # `time` collides with it there (the phase's values currently overwrite
  # the grid column). Selecting phases by excluding "time" would drop this
  # real phase; hzr_gof() selects by phase name and keeps its own grid.
  d <- .gof_pc_avc
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes"),
      time  = hzr_phase("constant")
    ),
    fit = TRUE
  )
  gof <- hzr_gof(fit)
  expect_identical(grep("^par_cumhaz_", names(gof), value = TRUE),
                   c("par_cumhaz_early", "par_cumhaz_time"))
  # The phases add up to the total, so the `time` phase's own contribution
  # is the total less the early phase, independent of decompose's columns.
  expect_equal(gof$par_cumhaz_time, gof$par_cumhaz - gof$par_cumhaz_early,
               tolerance = 1e-12)
  expect_gt(max(gof$par_cumhaz_time), 0.01)
  # hzr_gof()'s time column is its Kaplan-Meier grid, not the collided
  # column of decompose's output.
  km <- survival::survfit(survival::Surv(d$int_dead, d$dead) ~ 1)
  expect_equal(gof$time, km$time)
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
