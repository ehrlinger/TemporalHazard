# predict(newdata = ) on a multiphase fit whose phase carries its own formula.
#
# The phase design at newdata must be the fit's design: the same factor
# levels, contrasts and columns. Rebuilding it with a bare model.matrix()
# failed on a factor given as one label, and silently recoded a factor whose
# levels were a subset or in another order.
#
# Reference: the structural identity H_j(t | x) = exp(x beta_j) H0_j(t), with
# the baseline H0_j from a time-only newdata and beta_j from the fit's theta
# by name. It does not go through the design rebuild under test.

phase_formula_fit <- function(sum_contrasts = FALSE) {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  d$grp <- factor(ifelse(d$age > 100, "old", "young"))
  if (sum_contrasts) stats::contrasts(d$grp) <- stats::contr.sum(2)
  hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ grp),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  )
}

# Expected total cumulative hazard for rows with the given `young` indicator.
reference_cumhaz <- function(fit, time, young) {
  base <- predict(fit, newdata = data.frame(time = time),
                  type = "cumulative_hazard", decompose = TRUE)
  beta <- fit$fit$theta[["early.grpyoung"]]
  exp(beta * young) * base$early + base$constant
}

test_that("the reference can fail: the grp effect is not negligible", {
  fit <- phase_formula_fit()
  expect_true(fit$fit$converged)
  # A beta near 0 would let an ignored covariate pass every test below.
  expect_gt(abs(fit$fit$theta[["early.grpyoung"]]), 0.5)
  tt <- c(0.5, 2)
  expect_gt(min(abs(reference_cumhaz(fit, tt, 1) /
                      reference_cumhaz(fit, tt, 0) - 1)), 0.05)
})

test_that("a factor given as a single label codes to the fit's levels", {
  fit <- phase_formula_fit()
  tt <- c(0.5, 2, 6)
  young <- predict(fit, newdata = data.frame(time = tt, grp = "young"),
                   type = "cumulative_hazard")
  old <- predict(fit, newdata = data.frame(time = tt, grp = "old"),
                 type = "cumulative_hazard")
  expect_equal(young / reference_cumhaz(fit, tt, 1), rep(1, 3),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(old / reference_cumhaz(fit, tt, 0), rep(1, 3),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("a full-level factor, in any level order, codes to the fit's levels", {
  fit <- phase_formula_fit()
  tt <- c(1, 1)
  want <- reference_cumhaz(fit, tt, c(0, 1))
  same <- data.frame(time = tt,
                     grp = factor(c("old", "young"), levels = c("old", "young")))
  reversed <- data.frame(time = tt,
                         grp = factor(c("old", "young"),
                                      levels = c("young", "old")))
  subset <- data.frame(time = 1, grp = factor("young"))
  expect_equal(predict(fit, newdata = same, type = "cumulative_hazard") / want,
               c(1, 1), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(predict(fit, newdata = reversed, type = "cumulative_hazard") /
                 want, c(1, 1), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(predict(fit, newdata = subset, type = "cumulative_hazard") /
                 want[2], 1, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("extra and reordered newdata columns are ignored", {
  fit <- phase_formula_fit()
  tt <- c(0.5, 2)
  nd <- data.frame(age = c(5, 500), grp = c("young", "old"),
                   junk = c("a", "b"), time = tt)
  expect_equal(predict(fit, newdata = nd, type = "cumulative_hazard") /
                 reference_cumhaz(fit, tt, c(1, 0)),
               c(1, 1), tolerance = 1e-10, ignore_attr = TRUE)
  # decompose = TRUE and se.fit = TRUE read the same per-phase design.
  dec <- predict(fit, newdata = nd, type = "cumulative_hazard",
                 decompose = TRUE)
  base <- predict(fit, newdata = data.frame(time = tt),
                  type = "cumulative_hazard", decompose = TRUE)
  beta <- fit$fit$theta[["early.grpyoung"]]
  expect_equal(dec$early / (exp(beta * c(1, 0)) * base$early), c(1, 1),
               tolerance = 1e-10, ignore_attr = TRUE)
  se <- predict(fit, newdata = nd, type = "cumulative_hazard", se.fit = TRUE)
  expect_equal(se$fit / reference_cumhaz(fit, tt, c(1, 0)), c(1, 1),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("newdata carrying the fit's design columns by name is accepted", {
  fit <- phase_formula_fit()
  tt <- c(0.5, 2)
  nd <- data.frame(time = tt, grpyoung = c(1, 0))
  expect_equal(predict(fit, newdata = nd, type = "cumulative_hazard") /
                 reference_cumhaz(fit, tt, c(1, 0)),
               c(1, 1), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("an unseen level or a missing covariate is an error", {
  fit <- phase_formula_fit()
  expect_error(
    predict(fit, newdata = data.frame(time = 1, grp = "middle"),
            type = "cumulative_hazard"),
    "new level"
  )
  expect_error(
    predict(fit, newdata = data.frame(time = 1, age = 5),
            type = "cumulative_hazard"),
    "lacks the covariate column\\(s\\) 'grp'"
  )
})

test_that("a factor's non-default contrasts are the fit's, not the default", {
  # Under contr.sum the fit's column is grp1: old = +1, young = -1. Rebuilt
  # with default contrasts, a single label would code grpyoung instead.
  fit <- phase_formula_fit(sum_contrasts = TRUE)
  expect_true(fit$fit$converged)
  expect_identical(colnames(fit$fit$x_list$early), "grp1")
  tt <- c(0.5, 2)
  base <- predict(fit, newdata = data.frame(time = tt),
                  type = "cumulative_hazard", decompose = TRUE)
  beta <- fit$fit$theta[["early.grp1"]]
  expect_gt(abs(beta), 0.5)
  for (lab in c("old", "young")) {
    x <- if (lab == "old") 1 else -1
    got <- predict(fit, newdata = data.frame(time = tt, grp = lab),
                   type = "cumulative_hazard")
    expect_equal(got / (exp(beta * x) * base$early + base$constant),
                 c(1, 1), tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("a missing phase covariate is not filled from the workspace (#268)", {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  # A same-named object in the phase formula's own environment, which is
  # where a rebuild without the guard would look it up.
  mal <- 1
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ age, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes"),
      constant = hzr_phase("constant", formula = ~ mal)
    ),
    fit = TRUE
  )
  expect_error(
    predict(fit, newdata = data.frame(time = 2, age = 60),
            type = "cumulative_hazard"),
    "lacks the covariate column\\(s\\) 'mal' that phase 'constant' uses"
  )
})

test_that("a phase formula's environment constant still reaches newdata", {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  cutoff <- 100
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ I(age > cutoff)),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  )
  expect_true(fit$fit$converged)
  tt <- c(0.5, 2)
  base <- predict(fit, newdata = data.frame(time = tt),
                  type = "cumulative_hazard", decompose = TRUE)
  beta <- fit$fit$theta[["early.I(age > cutoff)TRUE"]]
  expect_gt(abs(beta), 0.5)
  got <- predict(fit, newdata = data.frame(time = tt, age = c(150, 50)),
                 type = "cumulative_hazard")
  expect_equal(got / (exp(beta * c(1, 0)) * base$early + base$constant),
               c(1, 1), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("all columns present: predictions unchanged from main 8a26c0e", {
  # Two numeric phase formulas, every covariate supplied: a case the old
  # rebuild got right. Values computed on main at 8a26c0e, before this fix.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ mal),
      constant = hzr_phase("constant", formula = ~ age)
    ),
    fit = TRUE
  )
  expect_true(fit$fit$converged)
  nd <- data.frame(time = c(0.5, 2, 6), mal = c(0, 1, 1), age = c(30, 60, 120))
  # Optimizer results, so a tolerance well above the optimizer's noise.
  expect_equal(predict(fit, newdata = nd, type = "cumulative_hazard"),
               c(0.08822505274, 0.42642619050, 0.49399211094),
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(predict(fit, newdata = nd, type = "survival"),
               c(0.9155548054, 0.6528380494, 0.6101855969),
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(predict(fit, newdata = nd, type = "hazard"),
               c(0.088515938428, 0.043098312832, 0.006837292565),
               tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("a fit without the stored phase design still predicts", {
  fit <- phase_formula_fit()
  fit$fit$x_design <- NULL
  tt <- c(1, 1)
  nd <- data.frame(time = tt,
                   grp = factor(c("old", "young"), levels = c("old", "young")))
  expect_equal(predict(fit, newdata = nd, type = "cumulative_hazard") /
                 reference_cumhaz(fit, tt, c(0, 1)),
               c(1, 1), tolerance = 1e-10, ignore_attr = TRUE)
})
