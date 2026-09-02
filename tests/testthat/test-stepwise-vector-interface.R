# Stepwise over a vector-interface base fit (#160).
#
# hzr_stepwise()'s refit path required a formula-interface base fit for every
# distribution. For multiphase models that restriction was incidental: a scope
# change rewrites the *phase* formula, and the global formula only ever
# carried the response. Refusing the vector interface shut out every
# translated SAS SELECTION job, because SAS's censoring statements map onto
# this package's -1/0/1/2 status coding, which survival::Surv() does not
# share -- and 25 of the 26 such jobs also use ICENSOR, so requiring a formula
# would have forced exactly the status round-trip AGENTS.md records as having
# shipped a wrong answer.
#
# The single-distribution path genuinely does mutate the global formula, so it
# still requires one. That difference is asserted here so it stays deliberate.

vec_phases <- function() {
  list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
       const = hzr_phase("constant"))
}

vec_data <- function(n = 300, seed = 3, beta = 0) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  data.frame(tt = stats::rexp(n, 0.3 * exp(beta * x1)) + 0.01,
             ev = stats::rbinom(n, 1, 0.75),
             x1 = x1, x2 = stats::rnorm(n))
}

vec_fit <- function(D, ...) {
  suppressWarnings(hazard(time = D$tt, status = D$ev, dist = "multiphase",
                          phases = vec_phases(), fit = TRUE, ...))
}

test_that("the blocker refuses by distribution, not by interface alone", {
  skip_on_cran()
  D <- vec_data()
  fv <- vec_fit(D)
  fw <- suppressWarnings(hazard(time = D$tt, status = D$ev, dist = "weibull",
                                theta = c(0.3, 1), fit = TRUE))
  # Guard the guard: both really are vector-interface fits, or the contrast
  # below is between two formula fits and says nothing.
  expect_null(fv$call$formula)
  expect_null(fw$call$formula)

  # Multiphase: scope lives in the phase formulas, so no global one is needed.
  expect_null(.hzr_refit_blocker(fv))
  # Single-distribution: scope lives in the global formula, so one is.
  expect_match(.hzr_refit_blocker(fw), "vector interface")
})

test_that("a multiphase vector fit with no stored response is still refused", {
  # hazard() always records time/status, so this is a hand-assembled object --
  # but without the check an absent `time` reaches hazard() as a missing
  # argument several frames down, where the message names nothing useful.
  hollow <- structure(
    list(call = list(), spec = list(dist = "multiphase"), data = list()),
    class = "hazard")
  expect_match(.hzr_refit_blocker(hollow), "no `time` / `status` vectors")
})

test_that("a vector-interface multiphase fit refits with an added covariate", {
  skip_on_cran()
  D <- vec_data()
  fv <- vec_fit(D)
  out <- suppressWarnings(
    .hzr_refit_with_scope(fv, action = "add", var = "x1",
                          data = D, phase = "early"))
  expect_s3_class(out, "hazard")
  expect_true("early.x1" %in% names(out$fit$par))
  # It must stay on the vector interface rather than acquiring a formula it
  # was never given.
  expect_null(out$call$formula)
})

test_that("a vector-interface refit can drop a covariate again", {
  skip_on_cran()
  D <- vec_data()
  ph <- vec_phases()
  ph$early <- hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0, formula = ~x1)
  fv <- suppressWarnings(hazard(time = D$tt, status = D$ev, data = D,
                                dist = "multiphase", phases = ph, fit = TRUE))
  expect_true("early.x1" %in% names(fv$fit$par))
  out <- suppressWarnings(
    .hzr_refit_with_scope(fv, action = "drop", var = "x1",
                          data = D, phase = "early"))
  expect_false("early.x1" %in% names(out$fit$par))
})

test_that("a vector-interface refit carries weights and time_lower through", {
  skip_on_cran()
  D <- vec_data()
  D$w <- rep(c(1, 2), length.out = nrow(D))
  D$entry <- 0
  fv <- suppressWarnings(hazard(time = D$tt, status = D$ev,
                                time_lower = D$entry, weights = D$w,
                                dist = "multiphase", phases = vec_phases(),
                                fit = TRUE))
  out <- suppressWarnings(
    .hzr_refit_with_scope(fv, action = "add", var = "x1",
                          data = D, phase = "early"))
  expect_s3_class(out, "hazard")
  # Dropping either silently changes the likelihood -- a left-truncated cohort
  # refit without time_lower is at risk from 0. Assert they survived rather
  # than that the refit merely returned.
  expect_equal(out$data$weights, D$w)
  expect_equal(out$data$time_lower, D$entry)
})

test_that("a vector-interface refit carries BOTH censoring bounds through", {
  skip_on_cran()
  # The test above uses status 0/1 and entry 0, so `time_upper` is NULL there
  # and asserting on it would compare NULL with NULL -- an assertion that
  # cannot fail. This fixture makes both bounds load bearing: interval rows
  # with genuinely wide intervals, and a non-zero entry time on the rest.
  D <- vec_data()
  n <- nrow(D)
  D$st <- rep(c(1, 0, 2), length.out = n)
  D$lo <- ifelse(D$st == 2, D$tt, 0.05 + seq_len(n) / (20 * n))
  D$up <- ifelse(D$st == 2, D$tt + 0.75, D$tt)

  # Guard the guard: if the fixture degenerated to zero entry times or
  # zero-width intervals, dropping either bound would not change the
  # likelihood and the assertions below would prove nothing.
  expect_true(all(D$lo[D$st != 2] > 0))
  expect_true(all(D$up[D$st == 2] > D$lo[D$st == 2]))
  expect_gt(sum(D$st == 2), 0L)

  fv <- suppressWarnings(hazard(time = D$tt, status = D$st,
                                time_lower = D$lo, time_upper = D$up,
                                dist = "multiphase", phases = vec_phases(),
                                fit = TRUE))
  out <- suppressWarnings(
    .hzr_refit_with_scope(fv, action = "add", var = "x1",
                          data = D, phase = "early"))

  expect_s3_class(out, "hazard")
  # Dropping `time_lower` refits a left-truncated cohort as at risk from 0;
  # dropping `time_upper` turns every interval row into something else. Both
  # are silent -- the refit still returns a populated fit.
  expect_equal(out$data$time_lower, D$lo)
  expect_equal(out$data$time_upper, D$up)
  expect_equal(out$data$status, D$st)
})

test_that("a single-distribution vector fit still cannot be refit", {
  skip_on_cran()
  D <- vec_data()
  fw <- suppressWarnings(hazard(time = D$tt, status = D$ev, dist = "weibull",
                                theta = c(0.3, 1), fit = TRUE))
  # Deliberate, not incidental: the single-distribution path rewrites the
  # global formula, so there is nothing to mutate.
  expect_error(
    .hzr_refit_with_scope(fw, action = "add", var = "x1", data = D),
    "vector interface")
})

test_that("a stepwise screen over a vector base fit TAKES a step", {
  skip_on_cran()
  # THE ACCEPTANCE CRITERION from #160. A call that merely executes is the
  # defect the issue exists to avoid: before this, every candidate refit
  # failed, the forward step downgraded those failures to warnings, and the
  # screen returned zero steps -- rendering clean and selecting nothing,
  # indistinguishable from an honest "nothing met slentry".
  # beta = 0.4, not larger: at 0.5+ the score statistic goes indefinite on x1
  # (the #130 case -- the criterion declines candidates in proportion to how
  # predictive they are) and the run depends on the Wald fallback instead.
  # That machinery is tested in test-score-wald-fallback.R; this test is about
  # the refit INTERFACE, so it uses an effect the score can test directly.
  # Measured here: p = 1.8e-08 against slentry = 0.05.
  D <- vec_data(n = 400, seed = 7, beta = 0.4)   # x1 planted, x2 noise
  sw <- suppressWarnings(hzr_stepwise(
    fit = vec_fit(D), scope = list(const = ~ x1 + x2), data = D,
    direction = "both", criterion = "score", slentry = 0.05, trace = FALSE))

  steps <- as.data.frame(sw)
  entered <- steps$variable[toupper(steps$action) == "ENTER"]
  expect_gte(nrow(steps), 1L)
  expect_true("x1" %in% entered)
  # And it entered the phase the scope named, not somewhere else.
  expect_identical(steps$phase[steps$variable == "x1"][1], "const")
  # The noise stays out, so this is selection rather than "refit everything".
  expect_false("x2" %in% entered)
  # No candidate was silently skipped: a refit failure here would look like a
  # clean screen, which is the failure mode being fixed.
  expect_length(sw$criteria$refit_failures %||% character(0), 0L)
  # And the score criterion actually scored it, so this exercises the ordinary
  # path rather than leaning on the fallback.
  expect_identical(sw$criteria$n_wald_fallbacks, 0L)
  expect_identical(steps$stat_type[1], "score_q")
})

test_that("hzr_stepwise()'s pre-flight names a remedy that fits the model", {
  skip_on_cran()
  D <- vec_data()
  fw <- suppressWarnings(hazard(time = D$tt, status = D$ev, dist = "weibull",
                                theta = c(0.3, 1), fit = TRUE))
  # A single-distribution fit is told to rebuild with a formula...
  expect_error(
    suppressWarnings(hzr_stepwise(fit = fw, scope = ~ x1, data = D,
                                  direction = "forward")),
    "Surv\\(time, status\\)")
  # ...and a multiphase one is never told that, because for a translated SAS
  # job it would force the -1/0/1/2 -> Surv() status round-trip.
  expect_null(.hzr_refit_blocker(vec_fit(D)))
})
