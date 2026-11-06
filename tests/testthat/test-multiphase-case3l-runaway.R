# Issue #600: on avc, the M-fixed M=0 NU=1 model from t_half 0.5, nu 1,
# mu 1 with the default five starts selected a start at log mu 130,
# log t_half 289, nu -7.66, reporting objective +6380.86 and converged = TRUE.
# Three of the other four reached -214.875438, which a single start from the
# same point also reaches. The selected point was not an optimum
# but an underflow: Case 3L of hzr_decompos() rounded G to 0 there while g
# stayed positive, so the cumulative hazard vanished and the hazard did not.

test_that("a multi-start does not select a Case 3L underflow point (#600)", {
  skip_on_cran()
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0,
                                    fixed = "m")),
    theta = c(early.log_mu = 0, early.log_t_half = log(0.5), early.nu = 1,
              early.m = 0),
    fit = TRUE,
    control = list(n_starts = 5L, conserve = FALSE, start_seed = 3L)
  ))
  st <- fit$fit$starts
  # Five starts were scored, so the comparison below is over all of them.
  expect_identical(nrow(st), 5L)
  expect_true(all(is.finite(st$objective)))
  # Every start, not only the selected one: the runaway start must not report
  # a log-likelihood above the optimum the others reach.
  expect_true(all(st$objective <= -214.875438 + 1e-6))
  expect_equal(fit$fit$objective, -214.875438, tolerance = 1e-8)
  expect_equal(unname(fit$fit$theta[c("early.log_mu", "early.nu")]),
               c(-1.1089213, 3.6436668), tolerance = 1e-5)
})
