# .hzr_formula_rhs_terms ----------------------------------------------------

test_that("rhs_terms extracts labels from common formula shapes", {
  expect_identical(.hzr_formula_rhs_terms(~ 1), character())
  expect_identical(.hzr_formula_rhs_terms(~ x), "x")
  expect_identical(.hzr_formula_rhs_terms(~ x + y), c("x", "y"))
  expect_identical(.hzr_formula_rhs_terms(Surv(time, status) ~ 1),
                   character())
  expect_identical(.hzr_formula_rhs_terms(Surv(time, status) ~ x + y),
                   c("x", "y"))
  expect_identical(.hzr_formula_rhs_terms(NULL), character())
})

# .hzr_formula_update -------------------------------------------------------

test_that("formula_update adds a variable to a one-sided formula", {
  f <- ~ x
  out <- .hzr_formula_update(f, "add", "y")
  expect_identical(.hzr_formula_rhs_terms(out), c("x", "y"))
})

test_that("formula_update adds to intercept-only formula", {
  out <- .hzr_formula_update(~ 1, "add", "x")
  expect_identical(.hzr_formula_rhs_terms(out), "x")
})

test_that("formula_update drops an existing variable", {
  out <- .hzr_formula_update(~ x + y + z, "drop", "y")
  expect_identical(.hzr_formula_rhs_terms(out), c("x", "z"))
})

test_that("dropping the only term leaves ~ 1", {
  out <- .hzr_formula_update(~ x, "drop", "x")
  expect_identical(.hzr_formula_rhs_terms(out), character())
  expect_identical(deparse(out), "~1")
})

test_that("adding a duplicate is a no-op", {
  f <- ~ x + y
  out <- .hzr_formula_update(f, "add", "x")
  expect_identical(.hzr_formula_rhs_terms(out), c("x", "y"))
})

test_that("dropping a missing variable is a no-op", {
  f <- ~ x + y
  out <- .hzr_formula_update(f, "drop", "z")
  expect_identical(.hzr_formula_rhs_terms(out), c("x", "y"))
})

test_that("two-sided formula preserves LHS across updates", {
  f <- Surv(time, status) ~ x + y
  added <- .hzr_formula_update(f, "add", "z")
  expect_identical(deparse(added[[2L]]), "Surv(time, status)")
  expect_identical(.hzr_formula_rhs_terms(added), c("x", "y", "z"))

  dropped <- .hzr_formula_update(f, "drop", "x")
  expect_identical(deparse(dropped[[2L]]), "Surv(time, status)")
  expect_identical(.hzr_formula_rhs_terms(dropped), "y")
})

test_that("formula environment is preserved across updates", {
  env <- new.env()
  f <- stats::as.formula("~ x", env = env)
  out <- .hzr_formula_update(f, "add", "y")
  expect_identical(environment(out), env)
})

test_that("invalid inputs are rejected", {
  expect_error(.hzr_formula_update("not a formula", "add", "x"),
               "must be a `formula` object")
  expect_error(.hzr_formula_update(~ 1, "add", character()),
               "non-empty character scalar")
  expect_error(.hzr_formula_update(~ 1, "add", ""),
               "non-empty character scalar")
  expect_error(.hzr_formula_update(~ 1, "add", c("x", "y")),
               "non-empty character scalar")
})

# .hzr_phase_update_formula -------------------------------------------------

test_that("phase_update_formula adds to an existing phase formula", {
  ph <- hzr_phase("cdf", t_half = 1, nu = 1, m = 1, formula = ~ age)
  out <- .hzr_phase_update_formula(ph, "add", "shock")
  expect_s3_class(out, "hzr_phase")
  expect_identical(.hzr_formula_rhs_terms(out$formula), c("age", "shock"))
})

test_that("phase_update_formula creates a formula when the phase had none", {
  ph <- hzr_phase("constant")
  expect_null(ph$formula)
  out <- .hzr_phase_update_formula(ph, "add", "age")
  expect_s3_class(out, "hzr_phase")
  expect_identical(.hzr_formula_rhs_terms(out$formula), "age")
})

test_that("phase_update_formula drops the only term, leaving `~ 1`", {
  # Not NULL: a NULL phase formula inherits the global design, so dropping a
  # phase's last covariate would bring the global ones back (#284).
  ph <- hzr_phase("cdf", t_half = 1, nu = 1, m = 1, formula = ~ age)
  out <- .hzr_phase_update_formula(ph, "drop", "age")
  expect_identical(deparse(out$formula), "~1")
})

test_that("phase_update_formula on NULL-formula drop is a no-op", {
  ph <- hzr_phase("constant")
  out <- .hzr_phase_update_formula(ph, "drop", "age")
  expect_null(out$formula)
})

test_that("phase_update_formula preserves other phase attributes", {
  ph <- hzr_phase("cdf", t_half = 0.5, nu = 2, m = 1,
                  fixed = "shapes", formula = ~ age)
  out <- .hzr_phase_update_formula(ph, "add", "shock")
  expect_identical(out$type, "cdf")
  expect_identical(out$t_half, 0.5)
  expect_identical(out$nu, 2)
  expect_identical(out$m, 1)
  expect_identical(out$fixed, c("t_half", "nu", "m"))
})

test_that("phase_update_formula rejects non-phase input", {
  expect_error(.hzr_phase_update_formula(list(), "add", "x"),
               "must be an `hzr_phase` object")
})

# .hzr_scope_current_vars ---------------------------------------------------

test_that("scope_current_vars lists vars for a single-distribution formula fit", {
  set.seed(1)
  n <- 100
  df <- data.frame(time = rexp(n), status = rep(1L, n),
                   age = rnorm(n), nyha = rnorm(n))
  fit <- hazard(
    Surv(time, status) ~ age + nyha,
    data = df,
    theta = c(0.5, 1.0, 0, 0),
    dist = "weibull",
    fit = TRUE
  )
  expect_identical(.hzr_scope_current_vars(fit), c("age", "nyha"))
})

test_that("scope_current_vars falls back to colnames(x) for non-formula fits", {
  fit <- .fit_weibull(n = 100L, betas = 0.3, seed = 1L)
  expect_identical(.hzr_scope_current_vars(fit), "x1")
})

test_that("scope_current_vars returns per-phase list for multiphase", {
  data(avc)
  set.seed(1)
  fit <- hazard(
    Surv(int_dead, dead) ~ 1,
    data = avc,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes", formula = ~ age),
      constant = hzr_phase("constant", formula = ~ age + mal)
    ),
    fit = TRUE,
    control = list(n_starts = 2L, maxit = 500L)
  )

  per_phase <- .hzr_scope_current_vars(fit)
  expect_named(per_phase, c("early", "constant"))
  expect_identical(per_phase$early, "age")
  expect_identical(per_phase$constant, c("age", "mal"))

  expect_identical(.hzr_scope_current_vars(fit, phase = "early"), "age")
  expect_identical(.hzr_scope_current_vars(fit, phase = "constant"),
                   c("age", "mal"))
})

test_that("scope_current_vars rejects a bad phase name for multiphase", {
  data(avc)
  set.seed(1)
  fit <- hazard(
    Surv(int_dead, dead) ~ 1, data = avc,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE,
    control = list(n_starts = 2L, maxit = 500L)
  )
  expect_error(.hzr_scope_current_vars(fit, phase = "late"),
               "Unknown phase")
})

test_that("scope_current_vars rejects `phase` on single-distribution fits", {
  fit <- .fit_weibull(n = 100L, seed = 1L)
  expect_error(.hzr_scope_current_vars(fit, phase = "early"),
               "only meaningful for multiphase")
})

test_that("scope_current_vars rejects non-hazard input", {
  expect_error(.hzr_scope_current_vars(list()),
               "must be a `hazard` object")
})
