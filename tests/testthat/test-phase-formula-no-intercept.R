# A phase formula without an intercept (`~ 0 + x`, `~ x - 1`) builds the
# same design as `~ x` (#303). A phase has no free intercept of its own (its
# scale `mu` plays that role), so the design always drops one. Before the fix
# the design dropped its first column by position, which without an
# intercept is the first covariate: the phase was fitted without it, silently.

test_that("an intercept-free phase formula keeps its first term (#303)", {
  d <- data.frame(age = c(50, 61, 72, 43, 58, 66),
                  mal = c(0, 1, 0, 1, 1, 0),
                  grp = factor(c("a", "b", "c", "a", "b", "c")))
  control <- .hzr_formula_design(~ age + mal, d)$x
  expect_identical(colnames(control), c("age", "mal"))

  for (f in list(~ 0 + age + mal, ~ age + mal - 1, ~ -1 + age + mal)) {
    x <- .hzr_formula_design(f, d)$x
    expect_identical(x, control, info = deparse(f))
  }

  # A single term: the whole design, not an empty one.
  expect_identical(colnames(.hzr_formula_design(~ 0 + age, d)$x), "age")

  # A factor codes against a reference level, exactly as with an intercept.
  # Keeping all its levels would be collinear with the phase scale.
  expect_identical(.hzr_formula_design(~ 0 + grp, d)$x,
                   .hzr_formula_design(~ grp, d)$x)
  expect_identical(colnames(.hzr_formula_design(~ 0 + grp, d)$x),
                   c("grpb", "grpc"))
  # A factor after the first term lost the term and kept every level.
  expect_identical(colnames(.hzr_formula_design(~ 0 + age + grp, d)$x),
                   c("age", "grpb", "grpc"))

  # `.` expands against the data first.
  expect_identical(.hzr_formula_design(~ 0 + ., d)$x,
                   .hzr_formula_design(~ ., d)$x)
})

test_that("an intercept-free phase fits and predicts as its control (#303)", {
  skip_on_cran()
  d <- na.omit(avc[, c("int_dead", "dead", "age", "mal")])
  fitp <- function(fe) {
    hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
           phases = list(
             early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                               fixed = "m", formula = fe),
             constant = hzr_phase("constant")
           ),
           fit = TRUE, control = list(n_starts = 1L, conserve = FALSE))
  }
  control <- fitp(~ age)
  expect_identical(colnames(control$fit$x_list$early), "age")
  nd <- data.frame(time = c(0.5, 2, 5), age = c(40, 60, 80))
  p_control <- predict(control, newdata = nd, type = "cumulative_hazard")

  for (fe in list(~ 0 + age, ~ age - 1)) {
    fit <- fitp(fe)
    expect_identical(colnames(fit$fit$x_list$early), "age")
    expect_identical(fit$fit$theta, control$fit$theta)
    expect_identical(fit$fit$objective, control$fit$objective)
    expect_identical(
      predict(fit, newdata = nd, type = "cumulative_hazard"), p_control
    )
  }

  # The control is not the covariate-free model, so the equalities above
  # compare a fitted covariate and not two empty designs.
  none <- fitp(NULL)
  expect_gt(control$fit$objective - none$fit$objective, 1)
})
