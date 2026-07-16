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
  fx <- readRDS(test_path("..", "..", "inst", "fixtures",
                          "stepwise-avc-forward-wald.rds"))
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
