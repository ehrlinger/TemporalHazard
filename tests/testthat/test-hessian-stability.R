test_that("fit object carries rcond and pd diagnostics (weibull)", {
  set.seed(2)
  n <- 250L
  dat <- data.frame(t = rweibull(n, shape = 1.4, scale = 3),
                    d = rep(1L, n))
  fit <- hazard(survival::Surv(t, d) ~ 1, data = dat, dist = "weibull",
                theta = c(mu = 1, nu = 1), fit = TRUE)
  expect_true(is.finite(fit$fit$rcond))
  expect_true(isTRUE(fit$fit$pd))
})

# NOTE: A more comprehensive version of this test (13-param, CoE, scale-invariance)
# lives in test-multiphase-hessian.R (Task 4). This simpler 2-param fit
# remains useful for the Layer-1 diagnostics contract.
test_that("fit object carries rcond and pd diagnostics (multiphase)", {
  skip_on_cran()
  set.seed(101)
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                    fixed = "m"),
                  constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 3, maxit = 500, conserve = TRUE)
  )
  expect_true(is.finite(fit$fit$rcond))
  expect_false(is.null(fit$fit$pd))
})

test_that("summary prints no Hessian note for a clean fit", {
  set.seed(3)
  dat <- data.frame(t = rweibull(250, 1.4, 3), d = rep(1L, 250))
  fit <- hazard(survival::Surv(t, d) ~ 1, data = dat, dist = "weibull",
                fit = TRUE, theta = c(mu = 1, nu = 1))
  out <- capture.output(print(summary(fit)))
  expect_false(any(grepl("ill-conditioned|not positive-definite", out)))
})

test_that("summary prints a note when the fit is flagged ill-conditioned", {
  s <- summary(structure(
    list(spec = list(dist = "weibull", phases = NULL),
         engine = "test",
         fit = list(theta = c(a = 1), vcov = matrix(1, 1, 1),
                    converged = TRUE, objective = -1,
                    rcond = 1e-12, pd = TRUE),
         data = list(time = 1:5, x = NULL),
         call = quote(hazard())),
    class = "hazard"))
  out <- capture.output(print(s))
  expect_true(any(grepl("ill-conditioned", out)))
})

test_that("summary prints a note when the fit is not positive-definite", {
  s <- summary(structure(
    list(spec = list(dist = "weibull", phases = NULL),
         engine = "test",
         fit = list(theta = c(a = 1), vcov = matrix(1, 1, 1),
                    converged = TRUE, objective = -1,
                    rcond = 0.5, pd = FALSE),
         data = list(time = 1:5, x = NULL),
         call = quote(hazard())),
    class = "hazard"))
  out <- capture.output(print(s))
  expect_true(any(grepl("not positive-definite", out)))
  expect_false(any(grepl("ill-conditioned", out)))  # rcond=0.5 is healthy
})

test_that("13-parameter multiphase deciles fit is stable (anchor)", {
  skip_on_cran()
  set.seed(102)
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  # Continuous covariates in `avc` span wildly different scales (age up to
  # ~286, op_age into the thousands). Unscaled, the multiphase Hessian picks
  # up non-finite entries and collapses to < 12 identifiable directions, so
  # scale them as any analyst fitting this model would. This is what makes
  # the high-dimensional inversion path a *meaningful* (not degenerate) test.
  for (cl in c("age", "opmos", "op_age", "inc_surg")) {
    avc[[cl]] <- as.numeric(scale(avc[[cl]]))
  }

  # A high-dimensional multiphase fit with covariates in both phases,
  # exercising the 12+-parameter inversion path named in DEVELOPMENT-PLAN
  # §7c ("hm.death.AVC.deciles"). Six covariates per phase plus the phase
  # shape parameters drive the free-parameter count past 12.
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf",
                           formula = ~ age + mal + inc_surg + opmos +
                             op_age + status,
                           t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
      constant = hzr_phase("constant",
                           formula = ~ age + mal + inc_surg + opmos +
                             op_age + status)
    ),
    fit = TRUE, control = list(n_starts = 3, maxit = 800, conserve = TRUE)
  )

  expect_true(fit$fit$converged)

  free <- which(is.finite(diag(fit$fit$vcov)))
  expect_gte(length(free), 12L)               # genuinely high-dimensional
  se <- sqrt(diag(fit$fit$vcov)[free])
  expect_true(all(is.finite(se)))             # no NaN/NA SEs on free params
  expect_true(all(se > 0))                    # positive variances
  expect_true(is.finite(fit$fit$rcond))       # conditioning was measured
})

test_that("summary reports when standard errors are unavailable", {
  # A fitted object whose Hessian could not be inverted: vcov is NA (scalar),
  # so has_vcov is FALSE while converged is TRUE.
  s <- summary(structure(
    list(spec = list(dist = "weibull", phases = NULL),
         engine = "test",
         fit = list(theta = c(a = 1), vcov = NA,
                    converged = TRUE, objective = -1,
                    rcond = NA_real_, pd = NA),
         data = list(time = 1:5, x = NULL),
         call = quote(hazard()),
         degraded = c("standard_errors", "weak_direction_check"),
         degraded_causes = c(standard_errors = "Hessian not invertible",
                             weak_direction_check = "standard errors unavailable")),
    class = "hazard"))
  out <- capture.output(print(s))
  expect_true("    standard_errors: Hessian not invertible" %in% out)
})

test_that("clean fit reports no standard-errors-unavailable note", {
  set.seed(7)
  dat <- data.frame(t = rweibull(250, 1.4, 3), d = rep(1L, 250))
  fit <- hazard(survival::Surv(t, d) ~ 1, data = dat, dist = "weibull",
                fit = TRUE, theta = c(mu = 1, nu = 1))
  out <- capture.output(print(summary(fit)))
  expect_false(any(grepl("standard errors unavailable", out)))
})
