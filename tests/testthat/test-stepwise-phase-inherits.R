# tests/testthat/test-stepwise-phase-inherits.R
# A multiphase phase with no formula of its own inherits the global design.
# A stepwise step that entered a variable into such a phase built `~ var`
# from scratch, so it REPLACED the inherited covariates: `early.age` became
# `early.mal`, the log-likelihood fell, and the step's p-value belonged to
# that smaller model (#284). Every refit here uses the same `control` as the
# explicit reference fit, so the two are the same `hazard()` call.

avc_284 <- na.omit(avc[, c("int_dead", "dead", "age", "mal")])
ctl_284 <- list(n_starts = 1L, conserve = FALSE)

ph_284 <- function(fe = NULL, fc = NULL) {
  list(
    early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m",
                         formula = fe),
    constant = hzr_phase("constant", formula = fc)
  )
}

fit_284 <- function(formula, phases) {
  suppressWarnings(hazard(formula, data = avc_284, dist = "multiphase",
                          phases = phases, fit = TRUE, control = ctl_284))
}

# The covariate coefficients of one phase, without its shape parameters.
phase_covs <- function(fit, phase) {
  nm <- names(coef(fit))
  shapes <- paste0(phase, ".", c("log_mu", "log_t_half", "nu", "m"))
  setdiff(grep(paste0("^", phase, "\\."), nm, value = TRUE), shapes)
}

wald_p <- function(fit, name) {
  i <- which(names(coef(fit)) == name)
  z <- coef(fit)[[i]] / sqrt(diag(fit$fit$vcov))[[i]]
  2 * stats::pnorm(-abs(z))
}

step_284 <- function(base, scope, direction = "forward", ...) {
  suppressWarnings(hzr_stepwise(base, data = avc_284, scope = scope,
                                direction = direction, trace = FALSE,
                                control = ctl_284, ...))
}

test_that("the phase update starts from the terms a phase inherits", {
  inh <- ~ age
  p0 <- hzr_phase("constant")
  added <- .hzr_phase_update_formula(p0, "add", "mal", inherited = inh)
  expect_identical(all.vars(added$formula), c("age", "mal"))
  # Emptying an inheriting phase leaves `~ 1`: NULL would inherit again.
  dropped <- .hzr_phase_update_formula(p0, "drop", "age", inherited = inh)
  expect_identical(deparse(dropped$formula), "~1")
  # A phase with its own formula ignores what it would inherit.
  own <- hzr_phase("constant", formula = ~ zz)
  kept <- .hzr_phase_update_formula(own, "add", "mal", inherited = inh)
  expect_identical(all.vars(kept$formula), c("zz", "mal"))
})

test_that("a Wald step adds to the inherited terms and reports their p-value", {
  base <- fit_284(survival::Surv(int_dead, dead) ~ age, ph_284())
  expect_identical(phase_covs(base, "early"), "early.age")
  sw <- step_284(base, list(early = ~ mal), criterion = "wald",
                 slentry = 0.99)
  explicit <- fit_284(survival::Surv(int_dead, dead) ~ age,
                      ph_284(~ age + mal))

  expect_identical(sw$steps$variable, "mal")
  expect_identical(sw$steps$phase, "early")
  expect_identical(phase_covs(sw, "early"), c("early.age", "early.mal"))
  expect_identical(phase_covs(sw, "constant"), "constant.age")
  expect_equal(sw$fit$objective, explicit$fit$objective, tolerance = 1e-8)
  # The entered model contains the base one, so its log-likelihood cannot
  # fall; under this control it did, from -196.44 to -204.43, when `age`
  # was replaced.
  expect_gt(sw$fit$objective, base$fit$objective)
  expect_equal(sw$steps$p_value, wald_p(explicit, "early.mal"),
               tolerance = 1e-6)
})

test_that("the score screen enters the variable instead of stopping", {
  base <- fit_284(survival::Surv(int_dead, dead) ~ age, ph_284())
  sw <- step_284(base, list(early = ~ mal), slentry = 0.99)
  explicit <- fit_284(survival::Surv(int_dead, dead) ~ age,
                      ph_284(~ age + mal))

  entered <- sw$steps[sw$steps$action == "enter", ]
  expect_identical(entered$variable, "mal")
  expect_identical(entered$phase, "early")
  expect_identical(phase_covs(sw, "early"), c("early.age", "early.mal"))
  expect_equal(sw$fit$objective, explicit$fit$objective, tolerance = 1e-8)
  # The same base model written with `early ~ age` as the phase's own formula
  # goes through the path that already worked, so its score p-value is the
  # reference. The old code could not score the inheriting candidate at all.
  own <- step_284(fit_284(survival::Surv(int_dead, dead) ~ age,
                          ph_284(~ age)),
                  list(early = ~ mal), slentry = 0.99)
  p_own <- own$steps$p_value[own$steps$action == "enter"]
  expect_true(is.finite(entered$p_value))
  expect_equal(entered$p_value, p_own, tolerance = 1e-8)
})

test_that("a drop starts from the terms the phase inherits", {
  base <- fit_284(survival::Surv(int_dead, dead) ~ age + mal, ph_284())
  pairs <- vapply(.hzr_stepwise_drop_candidates(base),
                  function(c) paste0(c$var, "@", c$phase), character(1))
  expect_setequal(pairs, c("age@early", "mal@early",
                           "age@constant", "mal@constant"))

  refit <- suppressWarnings(.hzr_refit_with_scope(
    base, "drop", "mal", phase = "early", data = avc_284, control = ctl_284
  ))
  explicit <- fit_284(survival::Surv(int_dead, dead) ~ age + mal,
                      ph_284(~ age))
  expect_identical(phase_covs(refit, "early"), "early.age")
  expect_identical(phase_covs(refit, "constant"),
                   c("constant.age", "constant.mal"))
  expect_equal(refit$fit$objective, explicit$fit$objective, tolerance = 1e-8)
})

test_that("a phase with its own formula is unchanged", {
  base <- fit_284(survival::Surv(int_dead, dead) ~ 1, ph_284(~ age))
  sw <- step_284(base, list(early = ~ mal), criterion = "wald",
                 slentry = 0.99)
  explicit <- fit_284(survival::Surv(int_dead, dead) ~ 1,
                      ph_284(~ age + mal))
  expect_identical(phase_covs(sw, "early"), c("early.age", "early.mal"))
  expect_identical(phase_covs(sw, "constant"), character(0))
  expect_equal(sw$fit$objective, explicit$fit$objective, tolerance = 1e-8)
})

test_that("a vector-interface fit whose phases inherit a direct `x` is refused", {
  # A refit has no terms to rebuild `x` from, so candidate refits dropped it
  # from every phase that inherited it.
  mp <- suppressWarnings(hazard(
    data = avc_284, time = int_dead, status = dead,
    x = as.matrix(avc_284["age"]), dist = "multiphase", phases = ph_284(),
    fit = TRUE, control = ctl_284
  ))
  expect_error(step_284(mp, list(early = ~ mal), criterion = "wald",
                        slentry = 0.99),
               "passed directly as `x`")

  # When every phase has its own formula, no phase reads `x`, and the
  # screen runs.
  mp2 <- suppressWarnings(hazard(
    data = avc_284, time = int_dead, status = dead,
    x = as.matrix(avc_284["age"]), dist = "multiphase",
    phases = ph_284(~ age, ~ age), fit = TRUE, control = ctl_284
  ))
  sw <- step_284(mp2, list(early = ~ mal), criterion = "wald",
                 slentry = 0.99)
  expect_identical(phase_covs(sw, "early"), c("early.age", "early.mal"))
})

test_that("an inherit refusal names its own remedy only", {
  mp <- suppressWarnings(hazard(
    data = avc_284, time = int_dead, status = dead,
    x = as.matrix(avc_284["age"]), dist = "multiphase", phases = ph_284(),
    fit = TRUE, control = ctl_284
  ))
  err <- tryCatch(step_284(mp, list(early = ~ mal), criterion = "wald",
                           slentry = 0.99),
                  error = conditionMessage)
  expect_match(err, "hzr_phase(formula = ~ ...)", fixed = TRUE)
  # The generic multiphase remedy asks for `time` and `status`, which this
  # caller already supplied.
  expect_false(grepl("supplying `time`", err, fixed = TRUE))
})

test_that("the score test pins the candidate where its column is", {
  # model.matrix() puts main effects before interactions, so a phase with
  # `age * mal` gains `opmos` BEFORE `age:mal`, not last. Pinning the zero in
  # the last slot scored `age:mal` as if it were `opmos` (p = 0.433 against
  # 0.251). The same model with a pre-built `age * mal` column is the check.
  d5 <- na.omit(avc[, c("int_dead", "dead", "age", "mal", "opmos")])
  d5$am <- d5$age * d5$mal
  score_p <- function(formula, phases) {
    base <- suppressWarnings(hazard(formula, data = d5, dist = "multiphase",
                                    phases = phases, fit = TRUE,
                                    control = ctl_284))
    sw <- suppressWarnings(hzr_stepwise(base, data = d5,
                                        scope = list(early = ~ opmos),
                                        direction = "forward",
                                        slentry = 0.99, trace = FALSE,
                                        control = ctl_284))
    sw$steps$p_value[sw$steps$action == "enter"]
  }
  # Inherited from the global formula.
  expect_equal(score_p(survival::Surv(int_dead, dead) ~ age * mal, ph_284()),
               score_p(survival::Surv(int_dead, dead) ~ age + mal + am,
                       ph_284()),
               tolerance = 1e-6)
  # The phase's own formula.
  expect_equal(score_p(survival::Surv(int_dead, dead) ~ 1,
                       ph_284(~ age * mal)),
               score_p(survival::Surv(int_dead, dead) ~ 1,
                       ph_284(~ age + mal + am)),
               tolerance = 1e-6)
})

test_that("an inherited time-varying design is refused", {
  # Rebuilt from its term labels, the phase would lose its windows: `age_w1`
  # and `age_w2` would become one constant `age` effect.
  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = avc_284,
    dist = "multiphase", phases = ph_284(), time_windows = 1, fit = TRUE,
    control = ctl_284
  ))
  expect_error(step_284(base, list(early = ~ mal), criterion = "wald",
                        slentry = 0.99),
               "time_windows")
})

test_that("an inherited term with several columns is refused before output", {
  d6 <- avc_284
  d6$grp <- cut(d6$age, 3)
  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ grp, data = d6, dist = "multiphase",
    phases = ph_284(), fit = TRUE, control = ctl_284
  ))
  err <- NULL
  out <- utils::capture.output(
    err <- tryCatch(
      hzr_stepwise(base, data = d6, scope = list(early = ~ mal),
                   criterion = "wald", slentry = 0.99, trace = TRUE,
                   control = ctl_284),
      error = conditionMessage
    )
  )
  expect_match(err, "more than one column")
  expect_identical(out, character(0))
})

test_that("emptying a phase's own formula does not make it inherit again", {
  base <- fit_284(survival::Surv(int_dead, dead) ~ age, ph_284())
  refit <- function(fit, action, var) {
    suppressWarnings(.hzr_refit_with_scope(fit, action, var, phase = "early",
                                           data = avc_284,
                                           control = ctl_284))
  }
  r1 <- refit(base, "drop", "age")
  r2 <- refit(r1, "add", "mal")
  r3 <- refit(r2, "drop", "mal")
  expect_identical(deparse(r3$spec$phases$early$formula), "~1")
  expect_identical(phase_covs(r3, "early"), character(0))
  expect_identical(phase_covs(r3, "constant"), "constant.age")
})

test_that("a scope for an inheriting phase is checked where its refit reads it", {
  # The refit builds the phase formula in the stored formula's environment,
  # so `zz` here reaches it even though neither the scope's frame nor the
  # search path can see it (#278's check, kept in step with #284).
  local_base <- function() {
    zz <- avc_284$age + sin(seq_len(nrow(avc_284))) # nolint: object_usage_linter.
    fit_284(survival::Surv(int_dead, dead) ~ age, ph_284())
  }
  expect_error(
    hzr_bootstrap(local_base(), n_boot = 2, seed = 1,
                  scope = list(early = ~ zz), criterion = "wald",
                  slentry = 0.99),
    "uses 'zz', which is not a column"
  )
})

test_that("a forward screen that never steps the inheriting phase still runs", {
  # The refit hands the inheriting phase the same global design back, so the
  # inherit refusals apply only to a phase the screen can step.
  d7 <- avc_284
  d7$grp <- cut(d7$age, 3)
  run <- function(formula, ..., criterion = "wald") {
    base <- suppressWarnings(hazard(
      formula, data = d7, dist = "multiphase",
      phases = ph_284(fc = ~ age), fit = TRUE, control = ctl_284, ...
    ))
    suppressWarnings(hzr_stepwise(base, data = d7,
                                  scope = list(constant = ~ mal),
                                  direction = "forward", criterion = criterion,
                                  slentry = 0.99, trace = FALSE,
                                  control = ctl_284))
  }
  sw <- run(survival::Surv(int_dead, dead) ~ grp)
  expect_identical(sw$steps$variable, "mal")
  expect_identical(phase_covs(sw, "constant"),
                   c("constant.age", "constant.mal"))
  expect_length(phase_covs(sw, "early"), 2L)

  # The inheriting phase keeps its windowed design under either criterion.
  # The score test rebuilt it from the plain global matrix instead, so the
  # screen stopped, blaming the candidate.
  for (crit in c("wald", "score")) {
    tw <- run(survival::Surv(int_dead, dead) ~ age, time_windows = 1,
              criterion = crit)
    expect_identical(tw$steps$variable, "mal")
    expect_true("constant.mal" %in% names(coef(tw)))
    expect_true(all(c("early.age_w1", "early.age_w2") %in% names(coef(tw))))
  }
})

test_that("the score expansion refuses rather than guess a column's slot", {
  # With no column names to match, the candidate's slot cannot be found; the
  # last slot is wrong after an interaction, so the expansion declines.
  base <- fit_284(survival::Surv(int_dead, dead) ~ age, ph_284())
  base$fit$x_list$early <- unname(base$fit$x_list$early)
  expect_null(.hzr_score_expand(base, "mal", "early", avc_284))
})

test_that("the refit itself refuses to step an unrebuildable inheriting phase", {
  # hzr_stepwise() refuses first, so this drives the refit's own check.
  d8 <- avc_284
  d8$grp <- cut(d8$age, 3)
  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ grp, data = d8, dist = "multiphase",
    phases = ph_284(), fit = TRUE, control = ctl_284
  ))
  expect_error(
    .hzr_refit_with_scope(base, "add", "mal", phase = "early", data = d8,
                          control = ctl_284),
    "more than one column"
  )
})
