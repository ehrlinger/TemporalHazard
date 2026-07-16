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

test_that(".hzr_score_q reproduces SAS's first-step Q for the multiphase fixture", {
  fx <- readRDS(test_path("..", "..", "inst", "fixtures",
                          "stepwise-avc-forward-wald.rds"))
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  # Fixture step 1 is the first ENTER: SAS scored every candidate against the
  # intercept-only two-phase model and entered the largest Q.
  step1 <- fx$steps[fx$steps$step_num == 1L, ]

  # The early phase's shapes must be SAS's own fixed values (the unconditional
  # fit values under CONDITION=14), not round numbers: they enter Phi_j(t) and
  # so move Q directly.  These are the same values test-stepwise-parity.R uses.
  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))

  q <- .hzr_score_q(base, var = step1$variable, phase = step1$phase, data = avc)

  expect_equal(q$df, 1L)
  expect_equal(q$stat, step1$stat, tolerance = 1e-2)
  expect_equal(q$p_value, step1$p_value, tolerance = 1e-3)
})
