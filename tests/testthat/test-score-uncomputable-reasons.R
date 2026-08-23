# Why a candidate could not be scored, not just how many ----------------------
#
# .hzr_score_q() returned a bare NA for every failure, and hzr_stepwise()
# reported a bare count. Two of those failures mean opposite things:
#
#   collinear               drop the column
#   information_indefinite  the observed information at beta = 0 is not
#                           positive, which happens when the candidate's
#                           effect is LARGE -- it is probably the best
#                           variable on offer, and it was passed over
#
# Reporting the second as "a degenerate or collinear candidate column", which
# is what the old warning said, tells a user to discard their strongest
# variable. See issue #130.

.reason_fixture <- function(beta_true) {
  set.seed(4242)
  n <- 260
  d <- data.frame(z = rnorm(n), w = rnorm(n))
  tt <- rexp(n, 0.25 * exp(beta_true * d$z))
  d$t <- pmin(tt, 12)
  d$code <- ifelse(tt >= 12, 0L, 1L)
  ph <- list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 1,
                               fixed = "shapes"),
             const = hzr_phase("constant"))
  fit <- suppressWarnings(hazard(
    survival::Surv(t, code) ~ 1, data = d, dist = "multiphase", phases = ph,
    fit = TRUE, control = list(n_starts = 1L, maxit = 800L)))
  list(fit = fit, data = d)
}

test_that("a scorable candidate reports no reason", {
  fx <- .reason_fixture(0.2)
  q <- .hzr_score_q(fx$fit, "z", phase = "early", data = fx$data)
  expect_false(is.na(q$stat))
  expect_true(is.na(q$reason))
})

test_that("a strong candidate reports an indefinite information matrix", {
  # beta_true = 0.8 puts the candidate far enough from 0 that the
  # log-likelihood curves upward there and v_beta goes negative.
  fx <- .reason_fixture(0.8)
  q <- .hzr_score_q(fx$fit, "z", phase = "early", data = fx$data)
  expect_true(is.na(q$stat))
  expect_identical(q$reason, "information_indefinite")

  # Not a property of the data being unusable: wald tests the same candidate.
  w <- suppressWarnings(hzr_stepwise(
    fx$fit, scope = list(early = ~ z, const = ~ z), data = fx$data,
    direction = "forward", criterion = "wald", slentry = 0.2, trace = FALSE,
    control = list(n_starts = 1L)))
  expect_gt(length(setdiff(names(coef(w)), names(coef(fx$fit)))), 0L)
})

test_that("degenerate and collinear candidates keep their own reasons", {
  fx <- .reason_fixture(0.2)
  d <- fx$data
  d$k <- 1
  expect_identical(
    .hzr_score_q(fx$fit, "k", phase = "early", data = d)$reason, "constant")

  # Enter z, then offer an exact copy of it: adjusted variance collapses
  # towards zero from ABOVE, which is the opposite sign to the case above.
  entered <- suppressWarnings(hzr_stepwise(
    fx$fit, scope = list(early = ~ z, const = ~ z), data = d,
    direction = "forward", criterion = "score", slentry = 0.5, trace = FALSE,
    control = list(n_starts = 1L)))
  # Two exactly-collinear candidates, not one. For a perfect duplicate the
  # true v_beta is exactly 0, so its computed SIGN is decided by rounding --
  # this assertion passed on macOS and failed on Linux with
  # "information_indefinite" until the magnitude test was put ahead of the
  # sign test. A scaled copy rounds differently from an identical one, so the
  # pair samples both sides of the boundary.
  d$zdup <- d$z
  d$zscaled <- 2 * d$z
  for (v in c("zdup", "zscaled")) {
    expect_identical(
      .hzr_score_q(entered, v, phase = "early", data = d)$reason,
      "collinear",
      info = v)
  }
})

test_that("hzr_stepwise now tests the strong candidate instead of passing it over", {
  # Was: this run warned "passed over rather than tested", left
  # n_uncomputable_scores > 0, and selected nothing -- the candidate with the
  # largest effect on offer was the one the criterion could not test (#130).
  # The Wald fallback now refits exactly those candidates, so the screen
  # reaches them.  The reporting path below has not gone away; it is now only
  # for candidates whose fallback refit itself fails.
  fx <- .reason_fixture(0.8)
  r <- suppressWarnings(hzr_stepwise(
    fx$fit, scope = list(early = ~ z + w, const = ~ z + w), data = fx$data,
    direction = "forward", criterion = "score", slentry = 0.2,
    trace = FALSE, control = list(n_starts = 1L)))

  expect_false(isTRUE(r$criteria$stopped_uncomputable))
  # Rescued, not merely counted: nothing is left uncomputable, and the
  # fallback is reported rather than silently switching criterion.
  expect_identical(r$criteria$n_uncomputable_scores, 0L)
  expect_gt(r$criteria$n_wald_fallbacks, 0L)
  # z carries the planted effect; asserting it ENTERED is what distinguishes
  # "the fallback ran" from "the fallback got the right answer".
  steps <- as.data.frame(r)
  expect_true("z" %in% steps$variable[toupper(steps$action) == "ENTER"])
})

test_that("the reason text tells a user to keep the candidate, not drop it", {
  txt <- .hzr_score_reason_text("information_indefinite")
  expect_match(txt, "STRONG")
  expect_match(txt, "wald")
  # And it must not be reachable from the collinear code, which is the
  # opposite advice.
  expect_no_match(.hzr_score_reason_text("collinear"), "STRONG")

  # An unknown code reports as itself rather than becoming NA or "".
  expect_identical(.hzr_score_reason_text("not_a_real_code"), "not_a_real_code")
})

test_that("a degenerate-only screen does not claim a strong candidate", {
  fx <- .reason_fixture(0.2)
  d <- fx$data
  d$k <- 1
  # `k` is constant, `w` is noise: nothing enters, nothing is indefinite.
  r <- suppressWarnings(hzr_stepwise(
    fx$fit, scope = list(early = ~ k, const = ~ k), data = d,
    direction = "forward", criterion = "score", slentry = 0.05,
    trace = FALSE, control = list(n_starts = 1L)))
  expect_false("information_indefinite" %in%
                 names(r$criteria$uncomputable_reasons))
  expect_identical(names(r$criteria$uncomputable_reasons), "constant")
})


# The candidate's own curvature, separately from the model's ------------------
#
# SAS tests I_bb <= 0 before the tolerance test and reports it as a distinct
# diagnostic (q1.c err 2 vs err 3). The two say different things: err 2 is
# "this candidate's own observed information is not positive", err 3 is "the
# candidate is unusable given what is already in the model". Folding them
# together blames the current model for a fault in the candidate.
#
# The branch is REACHABLE on real data, and censoring rather than effect size
# is what reaches it: an interval-censored row contributes a difference of
# survival terms rather than a log density, so its second derivative is not
# sign-constrained. Measured on a 2-phase fit over a grid of sample size,
# interval share, interval width and true beta, 16 of 36 configurations gave
# i_bb < 0, the most negative being -0.27 (n = 300, 90% interval-censored,
# width 3, beta 1).
#
# It is NOT tested through such a fixture, and that is deliberate. i_bb is
# evaluated at the fitted theta, so its value carries the optimizer's
# cross-platform variation -- which is basin-scale, not rounding-scale. A
# fixture at i_bb = -6e-4 passed on macOS and failed on Linux and Windows
# because the fit landed slightly elsewhere and i_bb came out positive. A
# knife-edge numeric is the wrong instrument for a branch test; the branch is
# stubbed instead, following test-hessian-unavailable.R.

test_that("a candidate whose own information is not positive says so", {
  skip_if_not_installed("numDeriv")
  fx <- .reason_fixture(0.2)

  # Flip only the candidate's own diagonal, leaving everything else as the
  # real fit produced it, so the test isolates exactly the new branch.
  # Captured BEFORE mocking: afterwards this name resolves to the stub. Left
  # unqualified so it binds in the same namespace local_mocked_bindings()
  # writes to -- the package under test -- whether the suite runs under
  # load_all() or against the installed package.
  orig <- .hzr_score_information_expanded
  local_mocked_bindings(
    .hzr_score_information_expanded = function(current, exp_) {
      info <- orig(current, exp_)
      if (!is.null(info)) {
        b <- exp_$beta_idx
        info[b, b] <- -1e-3
      }
      info
    }
  )

  q <- .hzr_score_q(fx$fit, "z", phase = "early", data = fx$data)
  expect_true(is.na(q$stat))
  expect_identical(q$reason, "information_nonpositive")
})

test_that("the same candidate scores normally without the flipped diagonal", {
  # The control for the test above: without the stub this candidate is
  # scorable, so the stub is what produced the reason and not the fixture.
  skip_if_not_installed("numDeriv")
  fx <- .reason_fixture(0.2)
  q <- .hzr_score_q(fx$fit, "z", phase = "early", data = fx$data)
  expect_false(is.na(q$stat))
  expect_true(is.na(q$reason))
})

test_that("the non-positive and collinear diagnoses stay distinct", {
  txt <- .hzr_score_reason_text("information_nonpositive")
  expect_match(txt, "own observed information")
  # It must not read as collinearity, which carries the opposite remedy, nor
  # as the strong-candidate case.
  expect_no_match(txt, "collinear with")
  expect_no_match(txt, "STRONG")
  expect_false(identical(txt, .hzr_score_reason_text("collinear")))
  expect_false(identical(txt, .hzr_score_reason_text("information_indefinite")))
})
