# A finite, enormous, meaningless Q ------------------------------------------
#
# Defect A of issue #134: a production screen accepted a candidate with
# `stat` = 92,211 on 1 df and `p = 0.000`, and the refit that followed made
# the model worse. Neither existing guard reaches that case -- `v_beta` stays
# positive and well above the collinearity floor -- so the statistic is
# reported as evidence.
#
# The cause is not near-collinearity on its own. Measured on `avc`, as a
# candidate approaches collinearity `u_beta` shrinks with `v_beta` and
# Q stays around 1-2 before the collinear guard takes over. Q only explodes
# when the model being scored against is NOT at its optimum: the reduced-model
# score is then no longer zero, so `u_beta` is inflated while `v_beta` stays
# small. That is exactly the state a failed refit leaves behind, which is why
# #134 saw it at step 9 after steps 8 and 9 had lowered the log-likelihood.
#
# The reference guards this and TemporalHazard did not. `dqstat.c` forms the
# implied coefficient and rejects it beyond +/-50:
#
#     *qbeta = (*qz)*(*qse);          /* = u_beta / v_beta */
#     if (fabs(*qbeta) > 50.0e0) { *qflag = 4; return; }

.diverge_case <- function(noise_sd = 0.01) {
  data("avc", package = "TemporalHazard", envir = environment())
  set.seed(1L)
  fit0 <- hazard(
    Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1L, maxit = 500L))
  # `age` in the model, then a near-copy of it offered as a candidate.
  with_age <- .hzr_refit_with_scope(
    fit0, action = "add", var = "age", phase = "early",
    data = avc, control = list(n_starts = 1L, maxit = 500L))
  set.seed(99L)
  d <- avc
  d$agee <- d$age + stats::rnorm(nrow(avc), sd = noise_sd)
  list(fit = with_age, data = d)
}

.moved_off_optimum <- function(fit, shift) {
  fit$fit$theta[1] <- fit$fit$theta[1] + shift
  fit
}

test_that("a candidate scored against a non-optimal model is declined", {
  skip_if_not_installed("numDeriv")
  o <- .diverge_case()

  # At the optimum the same candidate scores normally. This is the control:
  # without it the test could not tell the guard from the fixture.
  at_opt <- .hzr_score_q(o$fit, "agee", phase = "early", data = o$data)
  expect_false(is.na(at_opt$stat))
  expect_true(is.na(at_opt$reason))

  # Moved off the optimum -- what a failed refit leaves behind -- the score
  # statistic runs to millions and is now refused rather than reported.
  broken <- .moved_off_optimum(o$fit, 0.25)
  q <- .hzr_score_q(broken, "agee", phase = "early", data = o$data)
  expect_true(is.na(q$stat))
  expect_identical(q$reason, "coefficient_diverging")
})

test_that("the guard fires on the implied coefficient, not on Q being large", {
  skip_if_not_installed("numDeriv")
  o <- .diverge_case()
  broken <- .moved_off_optimum(o$fit, 0.25)

  nui <- .hzr_score_nuisance(broken)
  ex <- .hzr_score_expand(broken, "agee", phase = "early", data = o$data)
  info <- .hzr_score_information_expanded(broken, ex)
  grad <- .hzr_score_gradient(broken, ex)
  b <- ex$beta_idx
  t_idx <- ex$theta_idx[nui$idx]
  i_bt <- info[b, t_idx, drop = FALSE]
  v_beta <- info[b, b] - as.numeric(i_bt %*% nui$inv %*% t(i_bt))

  # Neither existing guard reaches this: v_beta is positive and far above the
  # collinearity floor. That is why the statistic was reportable.
  expect_gt(v_beta, 0)
  expect_gt(v_beta, info[b, b] * sqrt(.Machine$double.eps))
  # It is the implied coefficient that is absurd.
  expect_gt(abs(grad[b] / v_beta), 50)
})

test_that("a legitimate candidate is well inside the threshold", {
  skip_if_not_installed("numDeriv")
  # A guard that clipped ordinary candidates would be worse than none. With
  # a mildly correlated candidate the implied coefficient is single digits,
  # against a threshold of 50.
  o <- .diverge_case(noise_sd = 0.1)
  q <- .hzr_score_q(o$fit, "agee", phase = "early", data = o$data)
  expect_false(is.na(q$stat))
  expect_true(is.na(q$reason))
})

test_that("the diverging reason is distinct and points at the failed refit", {
  txt <- .hzr_score_reason_text("coefficient_diverging")
  expect_match(txt, "infinity")
  # The actionable half: this is often a symptom of an earlier bad step.
  expect_match(txt, "NOT at its optimum")
  expect_false(identical(txt, .hzr_score_reason_text("collinear")))
  expect_false(identical(txt, .hzr_score_reason_text("information_indefinite")))
})
