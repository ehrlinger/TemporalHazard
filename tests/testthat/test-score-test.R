# Score-test (Q-statistic) unit and parity tests.
#
# The multiphase test is an ACCEPTANCE GATE: it compares R's per-step Q against
# SAS's own Q values, captured in inst/fixtures/stepwise-avc-forward-wald.rds
# from a PROC HAZARD run on the bundled avc data. If it fails, the statistic is
# wrong -- do not widen the tolerance.
#
# SAS's Q is deliberately approximate: its listing states "variances are
# approximate because shaping parameter covariances are ignored". R matches that
# approximation on purpose; the efficient score would NOT match.

# The early phase's shapes must be SAS's own fixed values (the unconditional fit
# values under CONDITION=14), not round numbers: they enter Phi_j(t) and so move
# Q directly. These are the same values test-stepwise-parity.R uses.
score_test_base_fit <- function(data, constant_phase = hzr_phase("constant")) {
  suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = data, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                           fixed = "shapes"),
      constant = constant_phase
    ),
    fit = TRUE
  ))
}

score_test_fixture_step1 <- function() {
  fx <- .hzr_load_stepwise_fixture("avc-forward-wald")
  testthat::skip_if(is.null(fx), "stepwise-avc-forward-wald fixture not installed")
  # Fixture step 1 is the first ENTER: SAS scored every candidate against the
  # intercept-only two-phase model and entered the largest Q.
  fx$steps[fx$steps$step_num == 1L, ]
}

test_that(".hzr_score_q reproduces SAS's first-step Q for the multiphase fixture", {
  step1 <- score_test_fixture_step1()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  base <- score_test_base_fit(avc)
  q <- .hzr_score_q(base, var = step1$variable, phase = step1$phase, data = avc)

  expect_equal(q$df, 1L)

  # Tolerance is 1e-3, not 1e-2, and it is load-bearing. `expect_equal()`'s
  # tolerance is RELATIVE here, and the three measured values are:
  #
  #   shipped (SAS-approximate)  34.33998  vs SAS 34.3396 -> rel 1.1e-5
  #   efficient (shapes in V)    34.6072   vs SAS 34.3396 -> rel 7.8e-3
  #
  # so at 1e-2 the efficient score PASSES and the gate fails to defend the
  # approximation it exists to lock in. 1e-3 is the loosest tolerance that
  # still rejects the efficient score while accepting the shipped statistic.
  # This is a tightening; do not loosen it back.
  expect_equal(q$stat, step1$stat, tolerance = 1e-3)

  # The fixture's p_value is 0.0001 -- SAS's truncated "<.0001" PRINT FLOOR, not
  # a value (the true p is ~4.6e-9). Comparing against it is vacuous: all.equal
  # switches to an absolute scale when mean(abs(target)) < tolerance, so any Q
  # above ~11 would pass. `stat` is the load-bearing comparison; assert only
  # that p lands under the same floor SAS printed.
  expect_lt(q$p_value, 1e-4)
})

test_that("the nuisance partition excludes shapes -- the efficient score does not match SAS", {
  # Guard for a deliberate approximation. SAS ignores shaping-parameter
  # covariances during selection, so .hzr_score_free_idx() must exclude shapes
  # from the nuisance block. Computing the EFFICIENT score (all parameters in
  # the nuisance block) is a defect, not an improvement: it is a different
  # statistic and it does not reproduce SAS. Pinned here so that widening the
  # partition can never regress silently past the parity gate.
  step1 <- score_test_fixture_step1()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- score_test_base_fit(avc)

  # Approximate partition, via the real code path.
  q_approx <- .hzr_score_q(base, var = step1$variable, phase = step1$phase,
                           data = avc)
  expect_equal(q_approx$stat, step1$stat, tolerance = 1e-3)

  # Efficient score, computed inline: same formula, but the nuisance index set
  # is ALL parameters rather than .hzr_score_free_idx()'s shape-free subset.
  exp_ <- .hzr_score_expand(base, step1$variable, step1$phase, avc)
  grad <- .hzr_score_gradient(base, exp_)
  info_exp <- .hzr_score_information_expanded(base, exp_)
  info_cur <- .hzr_score_information(base, theta = base$fit$theta)
  idx_all <- seq_along(base$fit$theta)
  inv_all <- solve(info_cur[idx_all, idx_all, drop = FALSE])

  b <- exp_$beta_idx
  i_bt <- info_exp[b, exp_$theta_idx[idx_all], drop = FALSE]
  v_eff <- info_exp[b, b] - as.numeric(i_bt %*% inv_all %*% t(i_bt))
  q_eff <- (grad[b]^2) / v_eff

  # Measured: q_eff = 34.6072 against SAS's 34.3396 -- rel 7.8e-3.
  expect_gt(abs(q_eff - step1$stat) / abs(step1$stat), 1e-3)
})

test_that(".hzr_score_free_idx keeps a covariate named like a shape parameter", {
  # The partition must be structural, not name-based. A covariate literally
  # named `m` yields the theta name `constant.m`, which a regex on shape names
  # drops from the nuisance block -- silently under-adjusting V_beta for every
  # other candidate in the step. A `constant` phase has no shape parameters at
  # all, so `constant.m` can only ever be a covariate.
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  avc$m <- avc$age

  fit <- score_test_base_fit(avc, constant_phase = hzr_phase("constant",
                                                             formula = ~ m))
  nm <- names(fit$fit$theta)
  expect_identical(nm, c("early.log_mu", "early.log_t_half", "early.nu",
                         "early.m", "constant.log_mu", "constant.m"))

  idx <- .hzr_score_free_idx(fit)

  # Keeps both log_mu positions and the constant phase's covariate; drops
  # exactly the early phase's three shape slots (positions 2:4).
  expect_identical(as.integer(idx), c(1L, 5L, 6L))
  expect_true(match("constant.m", nm) %in% idx)
  expect_false(match("early.m", nm) %in% idx)
})

test_that(".hzr_score_q returns NA when the nuisance block cannot be inverted", {
  # A block that exists but is not invertible must yield NA, not a wrongly
  # unadjusted (too large) Q. Distinct from the legitimately unadjusted case
  # where no nuisance parameters exist at all.
  step1 <- score_test_fixture_step1()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- score_test_base_fit(avc)

  singular <- list(inv = NULL, idx = .hzr_score_free_idx(base), ok = FALSE)
  q <- .hzr_score_q(base, var = step1$variable, phase = step1$phase,
                    data = avc, nuisance = singular)

  expect_true(is.na(q$stat))
  expect_true(is.na(q$p_value))
  expect_equal(q$df, 1L)
})

test_that(".hzr_score_q leaves v_beta unadjusted when no nuisance block exists", {
  # The mirror of the singular case: `ok = TRUE` with `inv = NULL` means there
  # are legitimately no nuisance parameters (e.g. an intercept-only model), and
  # the UNADJUSTED v_beta is the correct answer -- not NA. The two states must
  # stay distinguishable, so `inv = NULL` alone can never decide this.
  step1 <- score_test_fixture_step1()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- score_test_base_fit(avc)

  empty <- list(inv = NULL, idx = integer(0), ok = TRUE)
  q <- .hzr_score_q(base, var = step1$variable, phase = step1$phase,
                    data = avc, nuisance = empty)

  expect_false(is.na(q$stat))
  expect_true(is.finite(q$stat))
  expect_gt(q$stat, 0)
  # Unadjusted v_beta is larger than the adjusted one, so this Q is SMALLER
  # than the shipped statistic -- it is a different number, not SAS's.
  expect_lt(q$stat, step1$stat)
})

test_that(".hzr_score_q returns NA for a collinear candidate", {
  # A candidate that duplicates a column already in the model drives
  # I_bt %*% solve(I_tt) %*% I_tb -> I_bb, so v_beta collapses towards zero from
  # ABOVE. The old absolute `v_beta <= 0` guard never fired: v_beta landed at
  # ~1e-14, Q blew up to ~1e15, and the collinear candidate WON the step. The
  # relative floor is what makes the NA contract real.
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  avc$age_copy <- avc$age

  fit <- score_test_base_fit(avc, constant_phase = hzr_phase("constant",
                                                             formula = ~ age))
  q <- .hzr_score_q(fit, var = "age_copy", phase = "constant", data = avc)

  expect_true(is.na(q$stat))
  expect_true(is.na(q$p_value))
  expect_equal(q$df, 1L)
})

test_that(".hzr_score_q errors when data is not row-aligned with the fit", {
  # Row misalignment would otherwise make every candidate return NA, and
  # stepwise would report nothing significant -- a plausible-looking wrong
  # answer rather than a failure.
  step1 <- score_test_fixture_step1()
  data(avc, package = "TemporalHazard")
  avc_clean <- na.omit(avc)
  base <- score_test_base_fit(avc_clean)

  expect_error(
    .hzr_score_q(base, var = step1$variable, phase = step1$phase,
                 data = avc_clean[-1L, , drop = FALSE]),
    "row-aligned"
  )
})

# --- Defect guards: .hzr_score_free_idx() ---------------------------------

# NOTE: this is the ONLY test that pins the SAS-approximation partition (shapes
# excluded from the nuisance block) and the intercept fix (the defect where
# exponential's index set returned `2`, dropping `log_lambda`, its only
# baseline parameter) across all four single-distribution families. The
# "score Q agrees with the refit Wald chi-square" test above cannot substitute
# for it: mutation testing measured Q = 14.01 (correct baseline), 14.10
# (efficient score, shapes included -- a defect), and 8.44 (intercept dropped
# -- the pre-fix bug) against a Wald chi^2 of 13.04, and that test's
# `tolerance = 0.5` passes all three. Do not weaken or delete this test.
test_that(".hzr_score_free_idx keeps the intercept for single-distribution fits", {
  # The partition excludes SHAPES, not the intercept. Weibull's theta is
  # [mu, nu, betas...]: `mu` is the scale/intercept and must stay in the
  # nuisance block; only `nu` is a shape. Returning just the trailing beta
  # slots drops `mu`, under-adjusts V_beta, and makes Q too SMALL -- candidates
  # then silently fail to enter. Exponential is the sharpest case: its theta is
  # [log_lambda, betas...] with NO shape parameter at all, so nothing is
  # dropped. This mirrors the multiphase rule in the same function, where a
  # `constant` phase has no shapes and its `log_mu` is kept.
  #
  # log-logistic and log-normal are covered too: their theta names
  # ([log_alpha, log_beta, ...] and [mu, log_sigma, ...]) are the most
  # misleading of the four -- neither pair reads as "intercept, shape" on
  # sight -- which is exactly why pinning their index sets matters.
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  w <- hazard(survival::Surv(int_dead, dead) ~ age, data = avc,
              dist = "weibull", theta = c(0.01, 0.5, 0), fit = TRUE)
  expect_identical(as.integer(.hzr_score_free_idx(w)), c(1L, 3L))

  e <- hazard(survival::Surv(int_dead, dead) ~ age, data = avc,
              dist = "exponential", theta = c(0.01, 0), fit = TRUE)
  expect_identical(as.integer(.hzr_score_free_idx(e)), c(1L, 2L))

  ll <- hazard(survival::Surv(int_dead, dead) ~ age, data = avc,
               dist = "loglogistic", theta = c(0.01, 0.5, 0), fit = TRUE)
  expect_identical(as.integer(.hzr_score_free_idx(ll)), c(1L, 3L))

  ln <- hazard(survival::Surv(int_dead, dead) ~ age, data = avc,
               dist = "lognormal", theta = c(0.01, 0.5, 0), fit = TRUE)
  expect_identical(as.integer(.hzr_score_free_idx(ln)), c(1L, 3L))
})

test_that(".hzr_score_free_idx errors on an empty theta rather than mis-indexing", {
  # `seq.int(length(theta) - p_cov + 1L, length(theta))` with an empty theta
  # evaluates to seq.int(-1, 0) -> c(-1L, 0L). R does not error on that:
  # info[c(-1, 0), c(-1, 0)] silently DROPS row/column 1 and returns a
  # plausible-looking matrix. Fail loudly instead.
  fake <- structure(
    list(
      spec = list(dist = "weibull", control = list()),
      data = list(x = matrix(1, nrow = 2L, ncol = 1L)),
      fit  = list(theta = numeric(0))
    ),
    class = "hazard"
  )
  expect_error(.hzr_score_free_idx(fake), "theta")
})

test_that(".hzr_score_free_idx errors when multiphase covariate_counts are absent", {
  # Zero counts feed .hzr_log_mu_positions(), so every phase start after the
  # first would be wrong and the index set silently corrupt. A fitted
  # multiphase object always carries counts, so the old all-zero fallback only
  # converted a should-be-impossible state into a wrong V_beta.
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  avc$m <- avc$age

  fit <- score_test_base_fit(avc, constant_phase = hzr_phase("constant",
                                                             formula = ~ m))
  expect_identical(as.integer(.hzr_score_free_idx(fit)), c(1L, 5L, 6L))

  fit$fit$covariate_counts <- NULL
  expect_error(.hzr_score_free_idx(fit), "covariate_counts")
})

# --- Tier 2: single-distribution families, numeric oracle ------------------
# No SAS reference exists for these yet (PROC HAZARD sources in this repo are
# all multiphase). They are validated against numDeriv and against the existing
# refit path, which is the standard the analytic Hessians already meet.

.score_oracle_fit <- function(dist) {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  theta0 <- switch(
    dist,
    exponential = c(mu = 0.01),
    weibull     = c(mu = 0.01, nu = 0.5),
    loglogistic = c(mu = 0.01, nu = 0.5),
    lognormal   = c(mu = 0.01, nu = 0.5)
  )
  list(
    fit = hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = dist, theta = theta0, fit = TRUE),
    data = avc
  )
}

test_that("score U_beta matches numDeriv for each single-distribution family", {
  skip_if_not_installed("numDeriv")
  for (d in c("exponential", "weibull", "loglogistic", "lognormal")) {
    o <- .score_oracle_fit(d)
    exp_ <- .hzr_score_expand(o$fit, "age", phase = NULL, data = o$data)
    expect_false(is.null(exp_), label = paste0("expand for dist = ", d))

    # The analytic score at (theta_hat, beta = 0), against numDeriv on the same
    # expanded negative log-likelihood. This is the oracle: the analytic
    # gradient is what makes the score test cheaper than a refit, so it is the
    # thing that has to be right.
    grad <- .hzr_score_gradient(o$fit, exp_)
    nll <- function(par) {
      .hzr_score_single_nll(o$fit, exp_$x, par)
    }
    # SIGN: the .hzr_gradient_* functions return the score of the POSITIVE
    # log-likelihood, while .hzr_score_single_nll() is the negated (objective)
    # scale the information is a Hessian of. Q uses U^2, so the sign cannot
    # show up there -- which is exactly why it is asserted here instead.
    num <- -numDeriv::grad(nll, exp_$theta)
    expect_equal(as.numeric(grad), as.numeric(num), tolerance = 1e-5,
                 label = paste0("analytic vs numDeriv score for dist = ", d))

    # The reduced-model score is zero at the current MLE, so only the
    # candidate's own component survives. The tolerance here is the
    # optimizer's convergence floor (control$reltol = 1e-5 on the objective),
    # not the score's accuracy: theta_hat is a numerical MLE, so its score is
    # near zero rather than exactly zero. Compare against the candidate's own
    # component, which is ~1e3 -- four orders of magnitude larger.
    expect_lt(max(abs(grad[exp_$theta_idx])), 1e-2)

    q <- .hzr_score_q(o$fit, var = "age", phase = NULL, data = o$data)
    expect_true(is.finite(q$stat),
                label = paste0("finite Q for dist = ", d))
    expect_equal(q$df, 1L, label = paste0("df for dist = ", d))
    expect_true(q$p_value >= 0 && q$p_value <= 1,
                label = paste0("p in [0,1] for dist = ", d))
  }
})

test_that("score Q agrees with the refit Wald chi-square for a weak effect", {
  skip_if_not_installed("numDeriv")
  # Both are asymptotically chi^2(1) for the same hypothesis, so on a candidate
  # with a small true effect they must agree to within sampling tolerance. The
  # refit path is the trusted oracle here.
  #
  # NOTE: this only checks order-of-magnitude agreement, not discrimination.
  # Mutation testing showed `tolerance = 0.5` on a Wald chi^2 of ~13.04 admits
  # Q anywhere in ~[6.5, 19.6] -- wide enough that both the efficient score
  # (14.10, the shape-inclusion defect) and the intercept-dropped bug (8.44)
  # pass alongside the correct baseline (14.01). The
  # `.hzr_score_free_idx() keeps the intercept` test below is what actually
  # pins those cases; do not rely on this test to catch either regression.
  o <- .score_oracle_fit("weibull")
  q <- .hzr_score_q(o$fit, var = "age", phase = NULL, data = o$data)
  refit <- .hzr_refit_with_scope(o$fit, action = "add", var = "age",
                                 phase = NULL, data = o$data)
  # Single-distribution fits store UNNAMED theta; .hzr_wald_p() regenerates
  # canonical names via .hzr_parameter_names(), so the sole covariate is
  # `beta1`, never `age`.
  w <- .hzr_wald_p(refit, "beta1")
  expect_equal(q$stat, w$stat^2, tolerance = 0.5)
})

test_that("score returns NA for degenerate candidates rather than selecting them", {
  o <- .score_oracle_fit("weibull")
  d <- o$data
  d$const_col <- 1
  q <- .hzr_score_q(o$fit, var = "const_col", phase = NULL, data = d)
  expect_true(is.na(q$stat))
  expect_true(is.na(q$p_value))
})

test_that("score returns NA for a collinear single-distribution candidate", {
  skip_if_not_installed("numDeriv")
  # The collinear guard must hold on this path too: v_beta collapses towards
  # zero from ABOVE, so without the relative floor the duplicate column would
  # win the step with Q ~ 1e15.
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  avc$age_copy <- avc$age

  fit <- hazard(survival::Surv(int_dead, dead) ~ age, data = avc,
                dist = "weibull", theta = c(0.01, 0.5, 0), fit = TRUE)
  q <- .hzr_score_q(fit, var = "age_copy", phase = NULL, data = avc)
  expect_true(is.na(q$stat))
  expect_true(is.na(q$p_value))
})

test_that("score rejects a candidate already in the single-distribution model", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  fit <- hazard(survival::Surv(int_dead, dead) ~ age, data = avc,
                dist = "weibull", theta = c(0.01, 0.5, 0), fit = TRUE)
  expect_null(.hzr_score_expand(fit, "age", phase = NULL, data = avc))
})
