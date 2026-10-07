# A candidate tested by NEITHER criterion is reported as such.
#
# criterion = "score" declines a candidate whose observed information is
# indefinite at beta = 0 -- which happens when the effect is LARGE (#130). The
# Wald fallback then refits it. But that refit can converge, producing a
# perfectly good point estimate, while its Hessian is singular: no standard
# error, so .hzr_wald_p() returns NA and the Wald test cannot run either.
#
# The fallback loop used to `next` there silently. The row kept the SCORE's
# reason, describing the first of two independent failures and saying nothing
# about the second, and no counter moved. A strongly predictive variable
# vanished and the screen rendered as an honest "nothing met slentry" -- the
# indistinguishability #159 and #130 were both about.

# An IDENTIFIED base (#565): the early phase's shapes are fixed. With them
# free the base had no maximum (exponential data), and whether x1's score
# could be computed depended on where the fit stopped on that ridge: on one
# CI platform it could, and the premise of this file failed. Both x1 and x2
# carry a strong effect here, because on an identified base a noise
# candidate is always scored normally and never reaches the fallback; x2 is
# the candidate the fallback rescues. The premises are asserted in
# `expect_fnv_premises()`. A rescued NOISE candidate is tested separately,
# with its score reason mocked.
fnv_fixture <- function(beta = 1.2, beta2 = 1.2, n = 400, seed = 7) {
  withr::local_seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  data.frame(tt = stats::rexp(n, 0.3 * exp(beta * x1 + beta2 * x2)) + 0.01,
             ev = stats::rbinom(n, 1, 0.75),
             x1 = x1, x2 = x2)
}

fnv_phases <- function() {
  list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0,
                         fixed = "shapes"),
       const = hzr_phase("constant"))
}

# The base is a maximum, its nuisance block inverts, and each named candidate
# is untestable by the score because of its own information.
expect_fnv_premises <- function(fit, D, indefinite = c("x1", "x2")) {
  expect_true(isTRUE(fit$fit$converged))
  expect_lte(fit$fit$rel_gradient, .Machine$double.eps^(1 / 3))
  expect_true(isTRUE(.hzr_score_nuisance(fit)$ok))
  for (v in indefinite) {
    expect_identical(.hzr_score_q(fit, var = v, phase = "const",
                                  data = D)$reason,
                     "information_indefinite", label = v)
  }
}

fnv_fit <- function(D) {
  suppressWarnings(hazard(survival::Surv(tt, ev) ~ 1, data = D,
                          dist = "multiphase", phases = fnv_phases(),
                          fit = TRUE))
}

# x1's score information is indefinite, so the score criterion sends it to the
# Wald fallback. Its rescuing refit used to converge with a singular Hessian,
# but that was the cold-started refit: since #551 the refit keeps the better
# of two starts and has a variance. The Hessian failure is therefore made
# here -- the x1 refit converges, keeps its estimate, and carries no variance
# matrix -- and the premise block asserts each of those.
fnv_no_variance <- function(env = parent.frame()) {
  orig <- .hzr_refit_with_scope
  local_mocked_bindings(
    .hzr_refit_with_scope = function(current, action = c("add", "drop"), var,
                                     phase = NULL, data, ...) {
      r <- orig(current, action = action, var = var, phase = phase,
                data = data, ...)
      if (identical(action, "add") && identical(var, "x1")) {
        r$fit$vcov <- NA
        r$fit$se <- NA
      }
      r
    },
    .env = env
  )
}

test_that("the fallback reaches x1, whose refit then has no variance", {
  skip_on_cran()
  D <- fnv_fixture()
  expect_fnv_premises(fnv_fit(D), D)
  # Unmocked: the score cannot test x1, so the fallback refits it, and that
  # refit (from the better of two starts) has a variance and tests it.
  out <- suppressWarnings(.hzr_stepwise_forward_step(
    current = fnv_fit(D), scope = list(const = ~ x1 + x2), data = D,
    criterion = "score", slentry = 0.05))
  x1_row <- out$all_scores[out$all_scores$variable == "x1", ]
  expect_true(x1_row$fallback)
  expect_identical(x1_row$stat_type, "wald_z")
  # With the variance removed, the same refit is the case below.
  fnv_no_variance()
  refit <- suppressWarnings(.hzr_refit_with_scope(
    fnv_fit(D), action = "add", var = "x1", phase = "const", data = D))
  expect_true(isTRUE(refit$fit$converged))
  expect_false(is.matrix(refit$fit$vcov))
})

test_that("the rescuing refit converges but has no variance to test with", {
  skip_on_cran()
  fnv_no_variance()
  # Pin the mechanism, not just the label. If the refit ever starts failing
  # outright, or starts producing an SE, this fixture stops exercising the
  # path the tests below describe and they would pass for a different reason.
  D <- fnv_fixture()
  refit <- .hzr_refit_with_scope(fnv_fit(D), action = "add", var = "x1",
                                 phase = "const", data = D)
  expect_true(isTRUE(refit$fit$converged))
  cname <- .hzr_candidate_coef_name(refit, "x1", "const")
  # A real estimate: the candidate is not degenerate, it is strong.
  expect_true(is.finite(refit$fit$par[[cname]]))
  # But no usable variance, so the Wald test cannot be computed.
  w <- .hzr_candidate_score(criterion = "wald", mode = "entry",
                            current = fnv_fit(D), candidate = refit,
                            names = cname)
  expect_true(is.na(w$score))
})

test_that("the row says the rescue failed, not just that the score did", {
  skip_on_cran()
  fnv_no_variance()
  D <- fnv_fixture()
  out <- suppressWarnings(.hzr_stepwise_forward_step(
    current = fnv_fit(D), scope = list(const = ~ x1 + x2), data = D,
    criterion = "score", slentry = 0.05))

  x1_row <- out$all_scores[out$all_scores$variable == "x1", ]
  expect_equal(nrow(x1_row), 1L)
  expect_true(is.na(x1_row$score))
  # THE FIX. Previously "information_indefinite", which is the score's reason
  # and describes only the first failure.
  expect_identical(x1_row$reason, "fallback_no_variance")

  # And x2 IS rescued and tested -- so the fallback still works and this is
  # not "the rescue stopped running". x2 carries a real effect on this
  # fixture, so the Wald test passes it; a rescued NOISE candidate declined
  # on its merits is the mocked case below.
  x2_row <- out$all_scores[out$all_scores$variable == "x2", ]
  expect_true(x2_row$fallback)
  expect_false(is.na(x2_row$score))
  expect_lt(x2_row$p_value, 0.05)
})

test_that("n_wald_fallbacks counts rescues, not attempts", {
  skip_on_cran()
  # x1 is attempted and yields no test; x2 is attempted and does. Only x2 is a
  # fallback. Miscounting here would overstate how much of the selection was
  # decided by the substituted criterion.
  fnv_no_variance()
  D <- fnv_fixture()
  out <- suppressWarnings(.hzr_stepwise_forward_step(
    current = fnv_fit(D), scope = list(const = ~ x1 + x2), data = D,
    criterion = "score", slentry = 0.05))
  expect_identical(out$n_wald_fallbacks, 1L)
  expect_identical(out$n_uncomputable, 1L)
  expect_false(out$all_scores$fallback[out$all_scores$variable == "x1"])
})

test_that("a rescued noise candidate is declined on its merits, and counted once", {
  skip_on_cran()
  # Plumbing, by construction: on an identified base a noise candidate is
  # always scored, so its score reason is mocked as information_indefinite to
  # send it to the fallback. x1 is untestable by the score for real and its
  # refit has no variance, as above. So x1 is an attempt that yields no test
  # and x2 a rescue that does: one fallback, one uncomputable.
  D <- fnv_fixture(beta2 = 0)
  base <- fnv_fit(D)
  expect_fnv_premises(base, D, indefinite = "x1")
  # The premise of the mock: x2 really is scored normally.
  expect_false(is.na(.hzr_score_q(base, var = "x2", phase = "const",
                                  data = D)$stat))
  fnv_no_variance()
  orig_q <- .hzr_score_q
  local_mocked_bindings(.hzr_score_q = function(current, var, ...) {
    r <- orig_q(current, var, ...)
    if (identical(var, "x2")) {
      r$stat <- NA_real_
      r$p_value <- NA_real_
      r$reason <- "information_indefinite"
    }
    r
  })
  out <- suppressWarnings(.hzr_stepwise_forward_step(
    current = base, scope = list(const = ~ x1 + x2), data = D,
    criterion = "score", slentry = 0.05))
  x2_row <- out$all_scores[out$all_scores$variable == "x2", ]
  expect_true(x2_row$fallback)
  expect_identical(x2_row$stat_type, "wald_z")
  expect_false(is.na(x2_row$score))
  expect_gt(x2_row$p_value, 0.05)
  expect_identical(out$n_wald_fallbacks, 1L)
  expect_identical(out$n_uncomputable, 1L)
})

test_that("a screen that tested nothing says so rather than looking clean", {
  skip_on_cran()
  fnv_no_variance()
  D <- fnv_fixture()
  # x1 alone in scope: x2 carries a real effect on this fixture and would
  # enter, and this test is about a screen with nothing it could test.
  w <- testthat::capture_warnings(
    sw <- hzr_stepwise(fit = fnv_fit(D), scope = list(const = ~ x1),
                       data = D, direction = "both", criterion = "score",
                       slentry = 0.05, trace = FALSE))
  # Zero steps is the honest outcome here -- neither criterion could test x1.
  expect_equal(nrow(as.data.frame(sw)), 0L)
  # But it must not be SILENT about it, or it is indistinguishable from
  # "nothing met slentry".
  expect_true(any(grepl("NEITHER criterion", w)))
  expect_identical(
    unname(sw$criteria$uncomputable_reasons["fallback_no_variance"]), 1L)
  # refit_failures stays empty -- the refit did not fail -- so a reader sent
  # there by the old wording would have found nothing.
  expect_length(sw$criteria$refit_failures %||% character(0), 0L)
})

test_that("the new reason has prose, and it names the right mechanism", {
  txt <- .hzr_score_reason_text("fallback_no_variance")
  # An unmapped code passes through as itself; that would be a bare token in a
  # user-facing warning.
  expect_false(identical(txt, "fallback_no_variance"))
  expect_match(txt, "no usable variance")
  expect_match(txt, "NEITHER")
  # And information_indefinite no longer claims the refit failed, which is
  # what it said before this case was distinguished.
  expect_false(grepl("that refit also failed",
                     .hzr_score_reason_text("information_indefinite")))
})

test_that("a weaker effect is still tested normally", {
  skip_on_cran()
  # The contrast that makes the above meaningful: same fixture, smaller beta,
  # and the score tests x1 directly. Without this, every assertion here could
  # pass for a screen that had simply stopped working.
  D <- fnv_fixture(beta = 0.4, beta2 = 0)
  sw <- suppressWarnings(hzr_stepwise(
    fit = fnv_fit(D), scope = list(const = ~ x1 + x2), data = D,
    direction = "both", criterion = "score", slentry = 0.05, trace = FALSE))
  st <- as.data.frame(sw)
  expect_true("x1" %in% st$variable[toupper(st$action) == "ENTER"])
  expect_identical(st$stat_type[1], "score_q")
})
