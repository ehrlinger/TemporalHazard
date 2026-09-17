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
  early <- function(fe) {
    hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m", formula = fe)
  }
  fitp <- function(early) {
    hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
           phases = list(early = early, constant = hzr_phase("constant")),
           fit = TRUE, control = list(n_starts = 1L, conserve = FALSE))
  }
  control <- fitp(early(~ age))
  expect_identical(colnames(control$fit$x_list$early), "age")
  nd <- data.frame(time = c(0.5, 2, 5), age = c(40, 60, 80))
  p_control <- predict(control, newdata = nd, type = "cumulative_hazard")

  for (fe in list(~ 0 + age, ~ age - 1)) {
    expect_warning(ph <- early(fe), "removes the intercept")
    fit <- fitp(ph)
    expect_identical(colnames(fit$fit$x_list$early), "age")
    expect_identical(fit$fit$theta, control$fit$theta)
    expect_identical(fit$fit$objective, control$fit$objective)
    expect_identical(
      predict(fit, newdata = nd, type = "cumulative_hazard"), p_control
    )
  }

  # The control is not the covariate-free model, so the equalities above
  # compare a fitted covariate and not two empty designs.
  none <- fitp(early(NULL))
  expect_gt(control$fit$objective - none$fit$objective, 1)
})

test_that("hzr_phase() warns once that an intercept removal is ignored (#303)", {
  for (f in list(~ 0 + age, ~ age - 1, ~ -1 + age, ~ 0 + .)) {
    expect_warning(hzr_phase("constant", formula = f), "removes the intercept",
                   info = deparse(f))
  }
  for (f in list(~ age, ~ 1, ~ age + mal, ~ ., ~ 0)) {
    expect_no_warning(hzr_phase("constant", formula = f))
  }
})

test_that("the score test builds an intercept-free phase as the fit did (#303)", {
  skip_on_cran()
  # The score test rebuilds every phase but the candidate's from its stored
  # formula, so it must build that formula as the fit did.
  d <- na.omit(avc[, c("int_dead", "dead", "age", "mal", "com_iv", "inc_surg")])
  d$grp <- factor(ifelse(d$inc_surg > 2, "high", "low"),
                  levels = c("low", "high"))
  fitc <- function(constant) {
    hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
           phases = list(
             early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                               fixed = "m", formula = ~ age),
             constant = constant
           ),
           fit = TRUE, control = list(n_starts = 1L, conserve = FALSE))
  }
  control <- fitc(hzr_phase("constant", formula = ~ mal + grp))
  expect_warning(ph <- hzr_phase("constant", formula = ~ 0 + mal + grp),
                 "removes the intercept")
  fit <- fitc(ph)

  x <- .hzr_score_expand(fit, "com_iv", "early", d)$x_list$constant
  expect_identical(colnames(x), c("mal", "grphigh"))
  q <- .hzr_score_q(fit, "com_iv", "early", d)
  q_control <- .hzr_score_q(control, "com_iv", "early", d)
  expect_true(is.finite(q_control$stat))
  expect_identical(q$stat, q_control$stat)

  # A single-term phase lost its only column, so candidates in every other
  # phase could not be scored.
  expect_warning(ph <- hzr_phase("constant", formula = ~ 0 + age),
                 "removes the intercept")
  q <- .hzr_score_q(fitc(ph), "mal", "early", d)
  expect_true(is.na(q$reason))
  expect_true(is.finite(q$stat))
})

test_that("the score test declines a phase saved with the old design (#303)", {
  skip_on_cran()
  # A fit saved before #303 holds `~ 0 + o`, for an ordered `o`, as the
  # dummies `om`, `oh`. It now rebuilds as `o.L`, `o.Q`: the same count, so
  # without a name check the old coefficients were scored against the new
  # columns.
  d <- na.omit(avc[, c("int_dead", "dead", "age", "com_iv", "inc_surg")])
  d$o <- factor(cut(d$inc_surg, c(-Inf, 2, 4, Inf),
                    labels = c("l", "m", "h")), ordered = TRUE)
  expect_warning(ph <- hzr_phase("constant", formula = ~ 0 + o),
                 "removes the intercept")
  fit <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                dist = "multiphase",
                phases = list(
                  early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                    fixed = "m", formula = ~ age),
                  constant = ph
                ),
                fit = TRUE, control = list(n_starts = 1L, conserve = FALSE))
  expect_identical(colnames(fit$fit$x_list$constant), c("o.L", "o.Q"))
  expect_true(is.na(.hzr_score_q(fit, "com_iv", "early", d)$reason))

  old <- fit
  old$fit$x_list$constant <-
    stats::model.matrix(~ 0 + o, data = d)[, -1L, drop = FALSE]
  expect_identical(colnames(old$fit$x_list$constant), c("om", "oh"))
  expect_null(.hzr_score_expand(old, "com_iv", "early", d))
  expect_false(is.finite(.hzr_score_q(old, "com_iv", "early", d)$stat))
})
