# Wald fallback for candidates the score statistic cannot test (#130).
#
# The score criterion computes SAS HAZARD's Q exactly: Q = grad^2 * I22, the
# reciprocal Schur complement of the OBSERVED information at beta = 0
# (src/vars/q1.c).  When a candidate's true effect is far from zero the
# log-likelihood is convex there, the Schur complement turns negative, and the
# statistic is undefined -- so the criterion declines candidates in proportion
# to how predictive they are.
#
# SAS has the same defect.  q1.c documents it ("IT IS POSSIBLE THAT THE
# PROGRAM WILL RETURN A NEGATIVE Q VALUE ... THE USER SHOULD USE THE MORE
# EXPENSIVE Q2 AS AN ALTERNATIVE") and dqstat.c declines the candidate with
# p = 1.  Q2 is referenced once in the C tree and never implemented.
#
# So the statistic stays bit-faithful and only the *handling* diverges: a
# candidate the score cannot test is refit and tested by Wald, which is what
# the unbuilt Q2 was for.  These tests pin down both halves -- that the
# statistic is unchanged, and that the screen no longer misses strong
# candidates.

planted <- function(seed = 11, n = 400, beta = 0.9) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  x3 <- stats::rnorm(n)
  data.frame(
    tt = stats::rexp(n, rate = 0.3 * exp(beta * x1)) + 0.01,
    ev = stats::rbinom(n, 1, 0.85),
    x1 = x1, x2 = x2, x3 = x3
  )
}

# Formula interface deliberately: hzr_stepwise() refuses a vector-interface
# base fit up front (.hzr_refit_blocker(), #159), because every candidate
# refit would fail.  Nothing about the Wald fallback depends on the
# interface -- the score criterion sees the same fitted MLE either way.
planted_fit <- function(D) {
  suppressWarnings(hazard(
    survival::Surv(tt, ev) ~ 1, data = D, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                  const = hzr_phase("constant")),
    fit = TRUE))
}

test_that("the score statistic itself is unchanged for an untestable candidate", {
  skip_on_cran()
  # Parity guard. The divergence from SAS is in the handling, not in Q, so
  # .hzr_score_q() must still decline x1 with the same reason as before.
  D <- planted()
  q <- .hzr_score_q(planted_fit(D), var = "x1", phase = "const", data = D)
  expect_true(is.na(q$stat))
  expect_identical(q$reason, "information_indefinite")
})

test_that("the planted strong effect is selected, not the noise variables", {
  skip_on_cran()
  # x1 has beta_true = 0.9 (LR ~ 178, p ~ 0); x2 and x3 are pure noise. Before
  # the fallback the screen entered x2 and x3 and never tested x1 -- a
  # confident, plausible, wrong selection.
  D <- planted()
  sw <- suppressWarnings(hzr_stepwise(
    fit = planted_fit(D), scope = list(const = ~ x1 + x2 + x3),
    data = D, direction = "both", slentry = 0.05))

  steps <- as.data.frame(sw)
  entered <- steps$variable[toupper(steps$action) == "ENTER"]
  expect_true("x1" %in% entered)
  # Asserting the noise stays out is the half that can fail for the right
  # reason: a fallback that simply refit everything would enter all three.
  expect_false("x2" %in% entered)
  expect_false("x3" %in% entered)
})

test_that("the Wald fallback is reported rather than silent", {
  skip_on_cran()
  D <- planted()
  sw <- suppressWarnings(hzr_stepwise(
    fit = planted_fit(D), scope = list(const = ~ x1 + x2 + x3),
    data = D, direction = "both", slentry = 0.05))
  # A criterion that quietly switched itself for some candidates would be the
  # same class of defect as the one being fixed.
  expect_true(sw$criteria$n_wald_fallbacks >= 1L)
})

test_that("only untestable-because-strong candidates fall back", {
  skip_on_cran()
  # A constant column is declined for a reason a refit cannot rescue, so it
  # must not cost one. Paying a refit per degenerate candidate would give back
  # the whole point of the score criterion.
  D <- planted()
  D$flat <- 1
  sw <- suppressWarnings(hzr_stepwise(
    fit = planted_fit(D), scope = list(const = ~ x1 + flat),
    data = D, direction = "both", slentry = 0.05))
  entered <- as.data.frame(sw)$variable[toupper(as.data.frame(sw)$action) == "ENTER"]
  expect_false("flat" %in% entered)
  expect_true("x1" %in% entered)
})

test_that("fallback reasons are classified, not lumped together", {
  # Every reason eligible for a refit must be one .hzr_score_reason_text()
  # knows, or the warning would print a bare code.
  expect_false(any(
    .hzr_score_reason_text(.hzr_score_fallback_reasons) %in%
      .hzr_score_fallback_reasons
  ))
  # Degenerate causes must stay out of the fallback set: a refit cannot make a
  # constant or collinear column testable.
  expect_false("collinear" %in% .hzr_score_fallback_reasons)
  expect_false("constant" %in% .hzr_score_fallback_reasons)
  expect_false("non_numeric" %in% .hzr_score_fallback_reasons)
  expect_true("information_indefinite" %in% .hzr_score_fallback_reasons)
})
