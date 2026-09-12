# tests/testthat/test-phase-specific-covariates.R
#
# Phase 7d: correctness tests for phase-specific covariates (hzr_phase(formula=)).
#
# hzr_phase(formula = ~ ...) lets a covariate modulate the scale of a single
# phase, so the same covariate can carry different effects in different phases.
# The shipped smoke test only checks that such a fit runs; it does not check
# that the coefficients are correct or that a phase-specific covariate stays
# confined to its phase.
#
# These are simulation-based recovery tests -- the honest substitute for a SAS
# parity fixture when no reference output exists. They guard against this
# package's recurring "accepts the formal but never applies it" failure mode
# (seen with weights, counting-process times, and CoE weights) by confirming a
# covariate entered into a specific phase (a) recovers its true coefficient and
# (b) affects only that phase.

# ---------------------------------------------------------------------------
# Simulation helper: additive two-phase model (CDF early + constant), with
# optional phase-specific covariate effects on each phase's scale.
#   H(t | x) = mu_e(x) * Phi_early(t) + mu_c(x) * t
# Event times are drawn by inverting the additive cumulative hazard on a grid.
# ---------------------------------------------------------------------------
.sim_phase_cov <- function(n, beta_early = NULL, beta_const = NULL,
                           log_mu_e = -1.5, log_mu_c = -2.5,
                           t_half = 0.5, tmax = 80, seed = 1) {
  set.seed(seed)
  x      <- rnorm(n)
  tgrid  <- seq(0.001, tmax, length.out = 4000)
  Phi_e  <- hzr_phase_cumhaz(tgrid, t_half = t_half, nu = 1, m = 1, type = "cdf")
  U      <- -log(runif(n))
  sim_t  <- numeric(n)
  status <- integer(n)
  for (i in seq_len(n)) {
    mu_e <- exp(log_mu_e + if (!is.null(beta_early)) beta_early * x[i] else 0)
    mu_c <- exp(log_mu_c + if (!is.null(beta_const)) beta_const * x[i] else 0)
    H    <- mu_e * Phi_e + mu_c * tgrid
    # First grid index where the accumulated hazard reaches the target U[i];
    # NA only when the event never occurs within the grid (right-censored).
    j    <- match(TRUE, H >= U[i])
    if (is.na(j)) {
      sim_t[i]  <- max(tgrid)
      status[i] <- 0L            # right-censored: event beyond the grid
    } else {
      sim_t[i]  <- tgrid[j]
      status[i] <- 1L            # event observed (even if on the last point)
    }
  }
  data.frame(
    int_dead = sim_t,
    dead     = status,
    x        = x
  )
}

# ---------------------------------------------------------------------------
# Structural tests (fast, always run): a phase-specific covariate produces a
# coefficient for its own phase only.
# ---------------------------------------------------------------------------

test_that("a covariate in one phase produces a beta for that phase only", {
  df <- .sim_phase_cov(n = 200, beta_early = 0.8, seed = 3)

  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = df,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes", formula = ~ x),
      constant = hzr_phase("constant")
    ),
    fit     = TRUE,
    control = list(n_starts = 1L, maxit = 50L)
  )

  nm <- names(coef(fit))
  expect_true("early.x" %in% nm)
  expect_false("constant.x" %in% nm)
  # covariate_counts reflects the per-phase covariate placement
  expect_equal(unname(fit$fit$covariate_counts[["early"]]), 1L)
  expect_equal(unname(fit$fit$covariate_counts[["constant"]]), 0L)
})

test_that("the same covariate can enter two phases independently", {
  df <- .sim_phase_cov(n = 200, beta_early = 0.8, beta_const = -0.4, seed = 4)

  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = df,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes", formula = ~ x),
      constant = hzr_phase("constant", formula = ~ x)
    ),
    fit     = TRUE,
    control = list(n_starts = 1L, maxit = 50L)
  )

  nm <- names(coef(fit))
  expect_true(all(c("early.x", "constant.x") %in% nm))
  expect_equal(unname(fit$fit$covariate_counts[["early"]]), 1L)
  expect_equal(unname(fit$fit$covariate_counts[["constant"]]), 1L)
})

# ---------------------------------------------------------------------------
# Recovery tests (heavier simulations; skip on CRAN): the fitted coefficients
# match the data-generating truth.
# ---------------------------------------------------------------------------

test_that("phase-specific covariate in the early phase is recovered", {
  skip_on_cran()
  beta_early <- 0.8
  log_mu_e   <- -1.5
  log_mu_c   <- -3.0
  df <- .sim_phase_cov(n = 3000, beta_early = beta_early,
                       log_mu_e = log_mu_e, log_mu_c = log_mu_c, seed = 7)

  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = df,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes", formula = ~ x),
      constant = hzr_phase("constant")
    ),
    fit     = TRUE,
    control = list(n_starts = 3L, maxit = 1000L)
  )
  co <- coef(fit)

  expect_equal(unname(co[["early.x"]]),        beta_early, tolerance = 0.12)
  expect_equal(unname(co[["early.log_mu"]]),   log_mu_e,   tolerance = 0.15)
  expect_equal(unname(co[["constant.log_mu"]]), log_mu_c,  tolerance = 0.20)
})

test_that("the same covariate recovers opposite-sign effects across phases", {
  skip_on_cran()
  beta_early <- 1.0
  beta_const <- -0.6
  log_mu_e   <- -1.5
  log_mu_c   <- -2.5
  df <- .sim_phase_cov(n = 5000, beta_early = beta_early, beta_const = beta_const,
                       log_mu_e = log_mu_e, log_mu_c = log_mu_c, seed = 11)

  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = df,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes", formula = ~ x),
      constant = hzr_phase("constant", formula = ~ x)
    ),
    fit     = TRUE,
    control = list(n_starts = 3L, maxit = 1000L)
  )
  co <- coef(fit)

  # The defining feature: one covariate, two different (here opposite-sign)
  # phase-specific effects, both recovered.
  expect_equal(unname(co[["early.x"]]),    beta_early, tolerance = 0.12)
  expect_equal(unname(co[["constant.x"]]), beta_const, tolerance = 0.12)
  expect_true(co[["early.x"]] > 0)
  expect_true(co[["constant.x"]] < 0)
})

test_that("a covariate confined to one phase does not leak into the other", {
  skip_on_cran()
  # Truth: x affects ONLY the early phase. A model that also offers x to the
  # constant phase should estimate the constant-phase effect near zero.
  df <- .sim_phase_cov(n = 5000, beta_early = 1.0, beta_const = NULL,
                       log_mu_e = -1.5, log_mu_c = -2.5, seed = 13)

  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = df,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes", formula = ~ x),
      constant = hzr_phase("constant", formula = ~ x)
    ),
    fit     = TRUE,
    control = list(n_starts = 3L, maxit = 1000L)
  )
  co <- coef(fit)

  expect_equal(unname(co[["early.x"]]),    1.0, tolerance = 0.15)
  expect_equal(unname(co[["constant.x"]]), 0.0, tolerance = 0.15)
})

# ---------------------------------------------------------------------------
# Regression test: formula preprocessor parse-tree walk (not string regex)
# ---------------------------------------------------------------------------

test_that("formula with base-R function name matching phase name is not wrongly stripped", {
  # Regression for: the old string-regex `grepl("\\blog\\s*\\(", rhs_txt)`
  # would fire on `~ log(age)` when a phase happened to be named "log".
  # The parse-tree walk in .hzr_formula_has_phase_scope() must NOT strip a
  # formula whose only use of `log` is as a standard R function, not as a
  # phase-scoped call — even when the phase list contains no phase named "log".
  # Here we test the helper directly with a contrived phase name.
  df <- .sim_phase_cov(n = 200, beta_early = 0.5, seed = 99)
  df$log_age <- log(df$x + 2)  # a predictor that uses log()

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                          fixed = "shapes", formula = ~ log_age),
    constant = hzr_phase("constant")
  )
  # The global formula uses log() as a base-R function on the RHS.
  # Phase names are "early" and "constant" — neither is "log".
  # The preprocessor must leave the formula intact (no stripping to ~ 1).
  # If it wrongly strips, model.matrix() gets `~ 1` and the fit will lack
  # the log_age covariate, causing covariate_counts[["early"]] == 0.
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data = df, dist = "multiphase", phases = phases,
    fit = TRUE, control = list(n_starts = 1L)
  )
  # Covariate was routed via hzr_phase(formula = ~ log_age), not the global
  # formula, so covariate_counts should still be 1 regardless.
  expect_equal(unname(fit$fit$covariate_counts[["early"]]), 1L,
               label = "log_age covariate routed correctly via phase formula")

  # Now test the helper directly: a formula ~ log(age) with phases named
  # "early"/"constant" must NOT trigger phase-scope detection.
  rhs_no_scope <- quote(log(age) + offset(z))
  expect_false(
    .hzr_formula_has_phase_scope(rhs_no_scope, c("early", "constant")),
    label = "log() call not mistaken for phase scope"
  )

  # And a genuine phase-scoped formula MUST trigger it.
  rhs_with_scope <- quote(early(age + x) + constant(z))
  expect_true(
    .hzr_formula_has_phase_scope(rhs_with_scope, c("early", "constant")),
    label = "early() call correctly detected as phase scope"
  )
})

# ---------------------------------------------------------------------------
# Phase-scoped terms in the global formula are refused (#275)
# ---------------------------------------------------------------------------
# `constant(age)` in the global formula was stripped with the rest of the
# right-hand side and never reached its phase, so the fit converged without
# it. hazard() now stops instead of routing the term.

avc_275 <- na.omit(avc[, c("int_dead", "dead", "age", "mal")])

fit_275 <- function(formula, phase_formula = NULL) {
  hazard(formula, data = avc_275, dist = "multiphase",
         phases = list(
           early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                fixed = "m", formula = phase_formula),
           constant = hzr_phase("constant", formula = phase_formula)
         ),
         fit = TRUE, control = list(n_starts = 1L, conserve = FALSE))
}

test_that("a phase-scoped term in the global formula stops, naming the phase", {
  expect_error(fit_275(survival::Surv(int_dead, dead) ~ constant(age)),
               "phase 'constant'.*hzr_phase\\(.*formula = ~")
  expect_error(
    fit_275(survival::Surv(int_dead, dead) ~ early(age) + constant(mal)),
    "phases 'early', 'constant'"
  )
})

test_that("the refusal does not fire when `phases` is ignored", {
  # Outside dist = "multiphase" the phases are dropped with a warning, so a
  # phase named `log` must not stop `log(age)` from being a covariate.
  expect_warning(
    fit <- hazard(survival::Surv(int_dead, dead) ~ log(age),
                  data = avc_275, dist = "weibull",
                  phases = list(log = hzr_phase("constant")),
                  theta = c(0.3, 1, 0), fit = TRUE),
    "'phases' is ignored"
  )
  expect_identical(colnames(fit$data$x), "log(age)")
})

# The next two guard the fits the refusal must leave alone. Each global-formula
# fit is compared with the same model written through hzr_phase(formula = ),
# an independent route to the same design, not with a recorded number: from
# these starts BFGS stops short of the optimum, so a pinned log-likelihood
# would record where it stopped rather than a property of the model.

test_that("`~ log(age)`, with no phase named `log`, fits as a covariate", {
  fit <- fit_275(survival::Surv(int_dead, dead) ~ log(age))
  ref <- fit_275(survival::Surv(int_dead, dead) ~ 1,
                 phase_formula = ~ log(age))
  expect_true(fit$fit$converged)
  expect_identical(
    names(coef(fit)),
    c("early.log_mu", "early.log_t_half", "early.nu", "early.m",
      "early.log(age)", "constant.log_mu", "constant.log(age)")
  )
  expect_equal(fit$fit$objective, ref$fit$objective, tolerance = 1e-8)
  expect_equal(unname(coef(fit) / coef(ref)), rep(1, 7), tolerance = 1e-6)
})

test_that("a plain global covariate in a multiphase fit is unchanged", {
  fit <- fit_275(survival::Surv(int_dead, dead) ~ age)
  ref <- fit_275(survival::Surv(int_dead, dead) ~ 1, phase_formula = ~ age)
  expect_true(fit$fit$converged)
  expect_identical(
    names(coef(fit)),
    c("early.log_mu", "early.log_t_half", "early.nu", "early.m",
      "early.age", "constant.log_mu", "constant.age")
  )
  expect_equal(fit$fit$objective, ref$fit$objective, tolerance = 1e-8)
  expect_equal(unname(coef(fit) / coef(ref)), rep(1, 7), tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# `.` in a phase formula (#277)
# ---------------------------------------------------------------------------
# hzr_phase(formula = ~ .) was expanded against every column of `data`, so the
# response entered the phase design and the fit converged on it. hazard() now
# writes `.` out once, before fitting, against `data` without the Surv()
# variables, so every later reader of the phase formula sees explicit terms.

fit_277 <- function(phase_formula, data = avc_275) {
  hazard(survival::Surv(int_dead, dead) ~ 1, data = data, dist = "multiphase",
         phases = list(
           early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                fixed = "m", formula = phase_formula),
           constant = hzr_phase("constant")
         ),
         fit = TRUE, control = list(n_starts = 1L, conserve = FALSE))
}

test_that("`.` in a phase formula excludes the Surv() response columns", {
  fit <- fit_277(~ .)
  expect_identical(colnames(fit$fit$x_list$early), c("age", "mal"))
  expect_identical(all.vars(fit$spec$phases$early$formula), c("age", "mal"))
})

test_that("a `~ .` phase fit is the explicit phase fit", {
  fit_dot <- fit_277(~ .)
  fit_exp <- fit_277(~ age + mal)
  expect_true(fit_exp$fit$converged)
  expect_equal(fit_dot$fit$objective, fit_exp$fit$objective,
               tolerance = 1e-8)
  expect_equal(unname(coef(fit_dot) / coef(fit_exp)),
               rep(1, length(coef(fit_exp))), tolerance = 1e-6)
})

test_that("`. - mal` in a phase formula gives the `~ age` phase design", {
  fit <- fit_277(~ . - mal)
  expect_identical(colnames(fit$fit$x_list$early), "age")
})

test_that("predict() on a `~ .` phase fit needs no response columns", {
  fit_dot <- fit_277(~ .)
  fit_exp <- fit_277(~ age + mal)
  nd <- data.frame(time = 1, age = 60, mal = 1)
  s_dot <- predict(fit_dot, newdata = nd, type = "survival")
  expect_length(s_dot, 1L)
  expect_equal(s_dot, predict(fit_exp, newdata = nd, type = "survival"),
               tolerance = 1e-8)
})

test_that("`.` in a phase formula over response-only data", {
  d0 <- avc_275[, c("int_dead", "dead")]
  fit <- fit_277(~ ., data = d0)
  expect_identical(unname(fit$fit$covariate_counts[["early"]]), 0L)
  age <- avc_275$age
  expect_error(suppressWarnings(fit_277(~ . + age, data = d0)),
               "stands for no column")
})

test_that("`.` in a phase formula is refused on the vector interface", {
  # With `time =` and `status =` there is no Surv() term to say which columns
  # of `data` hold the response, so `.` cannot be expanded safely.
  expect_error(
    hazard(time = int_dead, status = dead, data = avc_275,
           dist = "multiphase",
           phases = list(
             early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                  fixed = "m", formula = ~ .),
             constant = hzr_phase("constant")
           ),
           fit = TRUE, control = list(n_starts = 1L, conserve = FALSE)),
    "formula interface"
  )
})
