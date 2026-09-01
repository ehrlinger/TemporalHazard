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

fnv_fixture <- function(beta = 0.5, n = 400, seed = 7) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  data.frame(tt = stats::rexp(n, 0.3 * exp(beta * x1)) + 0.01,
             ev = stats::rbinom(n, 1, 0.75),
             x1 = x1, x2 = stats::rnorm(n))
}

fnv_phases <- function() {
  list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
       const = hzr_phase("constant"))
}

fnv_fit <- function(D) {
  suppressWarnings(hazard(survival::Surv(tt, ev) ~ 1, data = D,
                          dist = "multiphase", phases = fnv_phases(),
                          fit = TRUE))
}

test_that("the rescuing refit converges but has no variance to test with", {
  skip_on_cran()
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

  # And the noise variable IS rescued, scored and declined on its merits --
  # so the fallback still works and this is not "the rescue stopped running".
  x2_row <- out$all_scores[out$all_scores$variable == "x2", ]
  expect_true(x2_row$fallback)
  expect_false(is.na(x2_row$score))
  expect_gt(x2_row$p_value, 0.05)
})

test_that("n_wald_fallbacks counts rescues, not attempts", {
  skip_on_cran()
  # x1 is attempted and yields no test; x2 is attempted and does. Only x2 is a
  # fallback. Miscounting here would overstate how much of the selection was
  # decided by the substituted criterion.
  D <- fnv_fixture()
  out <- suppressWarnings(.hzr_stepwise_forward_step(
    current = fnv_fit(D), scope = list(const = ~ x1 + x2), data = D,
    criterion = "score", slentry = 0.05))
  expect_identical(out$n_wald_fallbacks, 1L)
  expect_identical(out$n_uncomputable, 1L)
  expect_false(out$all_scores$fallback[out$all_scores$variable == "x1"])
})

test_that("a screen that tested nothing says so rather than looking clean", {
  skip_on_cran()
  D <- fnv_fixture()
  w <- testthat::capture_warnings(
    sw <- hzr_stepwise(fit = fnv_fit(D), scope = list(const = ~ x1 + x2),
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
  D <- fnv_fixture(beta = 0.4)
  sw <- suppressWarnings(hzr_stepwise(
    fit = fnv_fit(D), scope = list(const = ~ x1 + x2), data = D,
    direction = "both", criterion = "score", slentry = 0.05, trace = FALSE))
  st <- as.data.frame(sw)
  expect_true("x1" %in% st$variable[toupper(st$action) == "ENTER"])
  expect_identical(st$stat_type[1], "score_q")
})
