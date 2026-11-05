# Too slow for the CRAN check time budget; runs on CI and locally.
# See "Where the check's time actually goes" in AGENTS.md.
skip_on_cran()

# `.fit_overfitted()` lives in helper-stepwise.R (auto-sourced by
# testthat).

# .hzr_stepwise_drop_candidates --------------------------------------------

test_that("drop_candidates enumerates all in-model vars for single-dist", {
  obj <- .fit_overfitted()
  cands <- .hzr_stepwise_drop_candidates(obj$fit)
  vars  <- vapply(cands, `[[`, character(1L), "var")
  expect_setequal(vars, c("x1", "x2", "x3"))
  # Every candidate has phase = NULL
  expect_true(all(vapply(cands, function(c) is.null(c$phase), logical(1L))))
})

test_that("drop_candidates emits (var, phase) pairs for multiphase", {
  data(avc)
  set.seed(1)
  fit <- hazard(
    Surv(int_dead, dead) ~ 1, data = avc,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes", formula = ~ age),
      constant = hzr_phase("constant", formula = ~ age + mal)
    ),
    fit = TRUE,
    control = list(n_starts = 2L, maxit = 500L)
  )
  cands <- .hzr_stepwise_drop_candidates(fit)
  pairs <- vapply(cands, function(c) paste0(c$var, "@", c$phase),
                  character(1L))
  expect_setequal(pairs, c("age@early", "age@constant", "mal@constant"))
})


# .hzr_stepwise_backward_step ----------------------------------------------

test_that("backward step drops the weakest noise variable under Wald", {
  obj <- .fit_overfitted()
  step <- .hzr_stepwise_backward_step(
    obj$fit, data = obj$data,
    criterion = "wald", slstay = 0.20
  )
  expect_true(step$accepted)
  # The dropped var must be one of the noise covariates
  expect_true(step$variable %in% c("x2", "x3"))
  # Refit should lack that coefficient
  expect_false(step$variable %in% colnames(step$fit$data$x))
  # all_scores must cover every in-model var, including the strong one
  expect_setequal(step$all_scores$variable, c("x1", "x2", "x3"))
})

test_that("backward step keeps every variable when they're all strong", {
  # High-strength signal on every covariate -> nothing clears drop
  # threshold.
  set.seed(77L)
  n <- 500L
  df <- data.frame(
    time   = rexp(n, rate = 1),
    status = rep(1L, n),
    x1     = rnorm(n),
    x2     = rnorm(n)
  )
  df$time <- df$time * exp(-0.8 * df$x1 - 0.7 * df$x2)

  fit <- hazard(
    Surv(time, status) ~ x1 + x2, data = df,
    theta = c(0.5, 1.0, 0, 0),
    dist = "weibull", fit = TRUE
  )

  step <- .hzr_stepwise_backward_step(
    fit, data = df, criterion = "wald", slstay = 0.20
  )
  expect_false(step$accepted)
  expect_identical(step$fit, fit)
  expect_true(all(step$all_scores$p_value < 0.05))
})

test_that("backward step under AIC uses ΔAIC rule and the W - 2*df shortcut", {
  obj <- .fit_overfitted()
  step <- .hzr_stepwise_backward_step(
    obj$fit, data = obj$data, criterion = "aic"
  )
  # Noise vars have z^2 ~ 1 on average, so z^2 - 2 < 0 -> drop.
  expect_true(step$accepted)
  expect_lt(step$score, 0)
  expect_equal(step$score, step$delta_aic)
})

test_that("force_in variables are scored but not dropped", {
  obj <- .fit_overfitted()
  step <- .hzr_stepwise_backward_step(
    obj$fit, data = obj$data, criterion = "wald", slstay = 0.20,
    force_in = c("x2", "x3")
  )
  # Only x1 is droppable, but x1 is the strong signal -> no drop.
  expect_false(step$accepted)
  # Every in-model var appears in all_scores with force_in flagged
  # correctly.
  scores <- step$all_scores
  expect_true(scores$force_in[scores$variable == "x2"])
  expect_true(scores$force_in[scores$variable == "x3"])
  expect_false(scores$force_in[scores$variable == "x1"])
})

test_that("force_in can lock out every candidate entirely", {
  obj <- .fit_overfitted()
  step <- .hzr_stepwise_backward_step(
    obj$fit, data = obj$data, criterion = "wald", slstay = 0.20,
    force_in = c("x1", "x2", "x3")
  )
  expect_false(step$accepted)
  # All rows flagged as force_in
  expect_true(all(step$all_scores$force_in))
})

test_that("backward step returns unchanged when there are no covariates", {
  set.seed(1)
  n <- 100L
  df <- data.frame(time = rexp(n), status = rep(1L, n))
  fit <- hazard(
    Surv(time, status) ~ 1, data = df,
    theta = c(0.5, 1.0), dist = "weibull", fit = TRUE
  )
  step <- .hzr_stepwise_backward_step(
    fit, data = df, criterion = "wald", slstay = 0.20
  )
  expect_false(step$accepted)
  expect_identical(step$fit, fit)
  expect_identical(nrow(step$all_scores), 0L)
})

test_that("backward step works on multiphase with phase-specific drop", {
  # Note: avoid 3+ covariates per phase — there's a pre-existing
  # recycling warning in .hzr_multiphase_cumhaz() triggered at that
  # threshold (flagged as a separate spawned task).
  data(avc)
  set.seed(1)
  fit <- hazard(
    Surv(int_dead, dead) ~ 1, data = avc,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes",
                            formula = ~ age + mal),
      constant = hzr_phase("constant")
    ),
    fit = TRUE,
    control = list(n_starts = 3L, maxit = 1000L)
  )

  step <- .hzr_stepwise_backward_step(
    fit, data = avc, criterion = "wald", slstay = 0.20,
    control = list(n_starts = 2L, maxit = 500L)
  )

  # If a drop happens, it's from the early phase (the only phase with
  # covariates), and the refit should be multiphase and converged.
  if (step$accepted) {
    expect_identical(step$phase, "early")
    expect_s3_class(step$fit, "hazard")
    expect_identical(step$fit$spec$dist, "multiphase")
    expect_true(step$fit$fit$converged)
  }
  # Regardless, every (var, phase) pair in the early phase must be
  # represented in all_scores.
  early_vars <- step$all_scores$variable[step$all_scores$phase == "early"]
  expect_setequal(early_vars, c("age", "mal"))
})
