# Broader integration tests for hzr_stepwise() — Step 8.7 /
# STEPWISE-DESIGN.md §6.1.  These exercise end-to-end scenarios that
# individual-helper tests do not cover: phase-specific selection,
# warm-start refit efficiency, and the candidate-failure warning path.


# Phase-specific entry -----------------------------------------------------

test_that("phase-specific candidates are enumerated distinctly per phase", {
  # The core requirement the design doc calls out is that a single
  # variable can be considered for entry into each phase independently.
  # Test the mechanism directly: after adding `x` to one phase, the
  # remaining candidates for the other phase still include `x`.
  data(avc)
  set.seed(1L)
  fit <- hazard(
    Surv(int_dead, dead) ~ 1, data = avc,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE,
    control = list(n_starts = 2L, maxit = 500L)
  )

  # Scope advertises `age` for both phases.
  scope <- list(early = ~ age, constant = ~ age)
  before <- .hzr_stepwise_candidates(fit, scope = scope, data = avc)
  pairs_before <- vapply(before,
                         function(c) paste0(c$var, "@", c$phase),
                         character(1L))
  expect_setequal(pairs_before, c("age@early", "age@constant"))

  # Add age into the early phase and re-enumerate.
  with_early <- .hzr_refit_with_scope(
    fit, action = "add", var = "age", phase = "early",
    data = avc, control = list(n_starts = 2L, maxit = 500L)
  )
  after <- .hzr_stepwise_candidates(with_early, scope = scope,
                                      data = avc)
  pairs_after <- vapply(after,
                        function(c) paste0(c$var, "@", c$phase),
                        character(1L))
  # age@constant should still be a candidate; age@early must not be.
  expect_true("age@constant" %in% pairs_after)
  expect_false("age@early" %in% pairs_after)
})


# Warm-start refit efficiency ----------------------------------------------

test_that("warm-started refits converge in far fewer iterations than cold starts", {
  # Fit a baseline model with two covariates, then add a third via the
  # refit wrapper. Compare optim counts to a cold hazard() call at the
  # same target.
  set.seed(17L)
  n <- 300L
  df <- data.frame(
    time   = rexp(n, 1),
    status = rep(1L, n),
    x1     = rnorm(n),
    x2     = rnorm(n),
    x3     = rnorm(n)
  )
  df$time <- df$time * exp(-0.8 * df$x1 - 0.4 * df$x2)

  base <- hazard(
    Surv(time, status) ~ x1 + x2, data = df,
    theta = c(0.5, 1.0, 0, 0), dist = "weibull", fit = TRUE
  )

  warm <- .hzr_refit_with_scope(
    base, action = "add", var = "x3", data = df
  )

  # Cold fit to the same target — starts every beta at 0 and shape
  # params at (0.5, 1.0) rather than the MLE.
  cold <- hazard(
    Surv(time, status) ~ x1 + x2 + x3, data = df,
    theta = c(0.5, 1.0, 0, 0, 0), dist = "weibull", fit = TRUE
  )

  expect_true(warm$fit$converged)
  # Warm and cold both reach the same MLE (essential correctness
  # guarantee).
  expect_equal(warm$fit$objective, cold$fit$objective, tolerance = 1e-3)

  # Warm should use no more function evaluations than cold.  The
  # design-doc promise is "typically one BFGS pass (~0.1 s)"; on
  # trivial single-distribution problems the gap can be small because
  # cold-start is already near the MLE (weibull starts at mu = 0.5,
  # nu = 1 which fits rexp(·, 1) data well).  Regressions in warm-
  # start correctness would show up as a strict increase here.
  warm_iters <- warm$fit$counts["function"]
  cold_iters <- cold$fit$counts["function"]
  expect_lte(warm_iters, cold_iters)
})


# Candidate-failure warning -------------------------------------------------

test_that("multi-level factor candidate errors with a main-effects note", {
  # A 3-level factor expands to 2 columns in the design matrix — v1
  # scope is main effects only, so stepwise should error cleanly
  # rather than silently score one contrast.
  set.seed(223L)
  n <- 200L
  df <- data.frame(
    time   = rexp(n),
    status = rep(1L, n),
    x1     = rnorm(n),
    f      = factor(sample(letters[1:3], n, replace = TRUE))
  )
  df$time <- df$time * exp(-0.8 * df$x1)

  # Fit with the factor already in the model so we can test the
  # single-dist guard path directly via .hzr_candidate_coef_name().
  fit <- hazard(
    Surv(time, status) ~ x1 + f, data = df,
    theta = c(0.5, 1.0, 0, 0, 0),
    dist = "weibull", fit = TRUE
  )

  expect_error(
    TemporalHazard:::.hzr_candidate_coef_name(fit, "f", NULL),
    "expands to multiple coefficients"
  )
})

test_that("one bad candidate does not prevent others from entering", {
  # Set up a scope where one candidate is broken (via a nonsense name
  # triggering the `var not in data` branch, which surfaces as an
  # error inside tryCatch) while a real one succeeds. The forward step
  # should warn about the broken candidate and still accept the good
  # one.
  obj <- .fit_driver_base(seed = 444L)
  df2 <- obj$data
  # Introduce a candidate that refers to a column *not* in `data` —
  # hazard() will error inside the candidate refit's tryCatch.
  expect_warning(
    step <- .hzr_stepwise_forward_step(
      obj$fit,
      scope = c("x1", "nonexistent"),
      data = df2,
      criterion = "wald", slentry = 0.30
    ),
    "candidate refit failed for nonexistent"
  )
  expect_true(step$accepted)
  expect_identical(step$variable, "x1")
  expect_identical(step$refit_failures, "nonexistent")
})


# Base fit whose formula was passed by symbol -------------------------------
#
# `hazard()` stores its call via match.call(), so a formula written into a
# variable first is stored as a *symbol*, not a formula.  Everything that
# rebuilds the call for a refit has to resolve it.

test_that("scope refit works when the base fit's formula is a symbol", {
  data(avc)
  set.seed(1L)
  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )

  f   <- Surv(int_dead, dead) ~ 1
  fit <- hazard(f, data = avc, dist = "multiphase", phases = phases,
                fit = TRUE, control = list(n_starts = 1L, maxit = 500L))

  # The premise of the test: the call really does hold a symbol.
  expect_true(is.symbol(fit$call$formula))

  refit <- .hzr_refit_with_scope(fit, action = "add", var = "age",
                                 phase = "early", data = avc,
                                 control = list(n_starts = 1L, maxit = 500L))

  expect_s3_class(refit, "hazard")
  expect_true("age" %in% names(coef(refit)) ||
                any(grepl("age", names(coef(refit)))))
})

test_that("stepwise selects identically whether the formula was literal", {
  data(avc)
  phases <- function() {
    list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    )
  }
  ctl   <- list(n_starts = 1L, maxit = 500L)
  scope <- list(early = ~ age, constant = ~ age)

  set.seed(2L)
  lit <- hazard(Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
                phases = phases(), fit = TRUE, control = ctl)
  sw_lit <- hzr_stepwise(lit, scope = scope, data = avc, direction = "forward",
                         slentry = 0.2, trace = FALSE)

  set.seed(2L)
  f   <- Surv(int_dead, dead) ~ 1
  sym <- hazard(f, data = avc, dist = "multiphase",
                phases = phases(), fit = TRUE, control = ctl)
  sw_sym <- hzr_stepwise(sym, scope = scope, data = avc, direction = "forward",
                         slentry = 0.2, trace = FALSE)

  # Guard against both screens being empty, which would compare equal while
  # testing nothing.
  expect_gt(length(names(coef(sw_lit))), length(names(coef(lit))))
  expect_equal(names(coef(sw_sym)), names(coef(sw_lit)))
})

test_that("an unresolvable stored formula errors clearly, not obscurely", {
  # A fit saved with saveRDS() and reloaded in a fresh session loses the
  # binding its `formula` symbol pointed at, so evaluating the stored call
  # raises "object 'f' not found" -- which says nothing about which fit is
  # at fault or what to do.  Inside hzr_stepwise()/hzr_bootstrap() the
  # per-candidate tryCatch() then swallows it into a generic "candidate
  # refit failed" and an empty screen.  Blanking `call_env` reproduces that
  # state without a second R session.
  data(avc)
  set.seed(3L)
  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  f   <- Surv(int_dead, dead) ~ 1
  fit <- hazard(f, data = avc, dist = "multiphase", phases = phases,
                fit = TRUE, control = list(n_starts = 1L, maxit = 400L))

  fit$call_env <- new.env(parent = emptyenv())

  expect_error(
    .hzr_refit_with_scope(fit, action = "add", var = "age", phase = "early",
                          data = avc,
                          control = list(n_starts = 1L, maxit = 400L)),
    "could not be resolved"
  )
})


# Uncomputable scores are not the same as "nothing met slentry" -------------
#
# Under criterion = "score", a candidate whose Q statistic cannot be computed
# (degenerate or collinear column, uninvertible nuisance block) yields NA and
# is dropped from `valid`.  When that happens to every candidate the step
# returns the same shape as a legitimate "no candidate cleared slentry" stop,
# so a broken screen and a finished one are indistinguishable.

.uncomputable_fixture <- function() {
  data(avc)
  avc$constcol <- 1.0        # degenerate: Q is NA, but the column is numeric
  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  fit <- hazard(Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
                phases = phases, fit = TRUE,
                control = list(n_starts = 1L, maxit = 400L))
  list(fit = fit, data = avc)
}

test_that("a step whose scores are all uncomputable says so", {
  fx <- .uncomputable_fixture()

  step <- .hzr_stepwise_forward_step(
    fx$fit, scope = list(early = ~ constcol), data = fx$data,
    criterion = "score", slentry = 0.5
  )

  expect_false(step$accepted)
  expect_identical(step$stop_reason, "scores_uncomputable")
  expect_equal(step$n_uncomputable, 1L)
})

test_that("a step that simply found nothing good enough is distinguishable", {
  fx <- .uncomputable_fixture()

  # `age` scores fine (p ~ 1e-4) but cannot clear an impossible slentry.
  step <- .hzr_stepwise_forward_step(
    fx$fit, scope = list(early = ~ age), data = fx$data,
    criterion = "score", slentry = 1e-12
  )

  expect_false(step$accepted)
  expect_identical(step$stop_reason, "no_candidate_met_slentry")
  expect_equal(step$n_uncomputable, 0L)
})

test_that("hzr_stepwise warns when a screen stops on uncomputable scores", {
  fx <- .uncomputable_fixture()

  expect_warning(
    sw <- hzr_stepwise(fx$fit, scope = list(early = ~ constcol),
                       data = fx$data, direction = "forward",
                       criterion = "score", slentry = 0.5, trace = FALSE),
    "could not be computed"
  )
  expect_gt(sw$criteria$n_uncomputable_scores, 0L)
})


# A formula passed by variable must work everywhere, not just in the refit ----
#
# #117 fixed .hzr_refit_with_scope(), and NEWS said the by-symbol defect
# "also affected hzr_stepwise() directly" and was fixed. Two sibling sites in
# .hzr_stepwise_candidates() still deparsed the stored call, so the default
# scope = NULL path raised `invalid formula "f": not a call` -- verbatim the
# string NEWS claims no longer occurs. Both distributions reach it.

test_that("scope = NULL works when the base formula was passed by variable", {
  data(avc)
  avc <- na.omit(avc)

  f    <- Surv(int_dead, dead) ~ age
  base <- hazard(f, data = avc, dist = "weibull", fit = TRUE,
                 theta = c(mu = 0.01, nu = 0.5, 0))
  expect_true(is.symbol(base$call$formula))

  sw <- hzr_stepwise(base, scope = NULL, data = avc, direction = "forward",
                     slentry = 0.05, trace = FALSE)
  expect_s3_class(sw, "hazard")
  # The Surv() response columns must not be offered as candidates.
  expect_false(any(c("int_dead", "dead") %in% sw$steps$variable))
})

test_that("scope = NULL works by variable for multiphase too", {
  data(avc)
  avc <- na.omit(avc)

  g   <- Surv(int_dead, dead) ~ 1
  ph  <- list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                                fixed = "shapes"),
              constant = hzr_phase("constant"))
  fit <- hazard(g, data = avc, dist = "multiphase", phases = ph, fit = TRUE,
                control = list(n_starts = 1L, maxit = 400L))
  expect_true(is.symbol(fit$call$formula))

  cands <- .hzr_stepwise_candidates(fit, scope = NULL, data = avc)
  expect_gt(length(cands), 0L)
  expect_false(any(vapply(cands, function(c) c$var, character(1L)) %in%
                     c("int_dead", "dead")))
})

test_that("scope = NULL keeps logical columns and drops unmodellable ones", {
  # Whether a 0/1 field arrives logical or numeric depends on the reader that
  # built the frame, not on the variable. Dropping logicals would make the
  # default scope silently narrower for one of two equivalent read paths.
  d <- data.frame(
    t    = c(2, 4, 6, 8, 10, 12, 14, 16),
    ev   = c(1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L),
    num  = c(1, 2, 3, 4, 5, 6, 7, 8),
    flag = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE),
    chr  = letters[1:8],
    stringsAsFactors = FALSE
  )
  expect_equal(.hzr_modellable_vars(d, c("num", "flag", "chr")),
               c("num", "flag"))
})


# A zero-step screen must say WHY it is empty (#159) -------------------------
#
# The refit-failure counterpart of the uncomputable-score tests above.  Every
# accepted step goes through .hzr_refit_with_scope(); when the refit fails the
# forward/backward step downgrades it to a warning and returns nothing.  The
# step functions have always reported `refit_failures`, but hzr_stepwise()
# never read it, so a screen that could not fit a single candidate returned
# the same zero-row `steps` as one that fit them all and liked none.

test_that("hzr_stepwise rejects a vector-interface base fit up front", {
  # The reproduction from #159: `.hzr_refit_with_scope()` cannot rebuild a
  # call that has no stored formula, so EVERY candidate refit fails.  Before
  # the fix this ran to completion and returned a zero-step result with
  # nothing on the object to distinguish it from an honest empty screen.
  data(avc)
  avc <- na.omit(avc)
  vecfit <- hazard(time = avc$int_dead, status = avc$dead, dist = "weibull",
                   fit = TRUE, theta = c(mu = 0.01, nu = 0.5))
  # Premise of the test: the call really does lack a formula.
  expect_null(vecfit$call$formula)

  expect_error(
    hzr_stepwise(vecfit, scope = ~ age + mal, data = avc,
                 direction = "forward", trace = FALSE),
    "vector interface"
  )
})

test_that(".hzr_refit_blocker is the single answer both callers read", {
  data(avc)
  avc <- na.omit(avc)
  vecfit <- hazard(time = avc$int_dead, status = avc$dead, dist = "weibull",
                   fit = TRUE, theta = c(mu = 0.01, nu = 0.5))
  frmfit <- hazard(Surv(int_dead, dead) ~ 1, data = avc, dist = "weibull",
                   fit = TRUE, theta = c(mu = 0.01, nu = 0.5))

  expect_null(.hzr_refit_blocker(frmfit))
  expect_match(.hzr_refit_blocker(vecfit), "vector interface")

  # The predicate is distribution-aware (#160): a multiphase fit's scope lives
  # in the phase formulas, so it needs no global one and is NOT blocked. Both
  # vector fits above and below are the same interface; only `dist` differs,
  # which is what makes this a contrast rather than two spellings of one case.
  vecmp <- suppressWarnings(hazard(
    time = avc$int_dead, status = avc$dead, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                  const = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1)))
  expect_null(vecmp$call$formula)
  expect_null(.hzr_refit_blocker(vecmp))
})

# Fixture: an intercept-only formula fit whose only candidate names a column
# that is not in `data`, so hazard() errors inside the candidate tryCatch.
.refit_failure_fixture <- function() {
  data(avc)
  avc <- na.omit(avc)
  fit <- hazard(Surv(int_dead, dead) ~ 1, data = avc, dist = "weibull",
                fit = TRUE, theta = c(mu = 0.01, nu = 0.5))
  list(fit = fit, data = avc)
}

test_that("a screen emptied by refit failures records that on the object", {
  fx <- .refit_failure_fixture()

  warns <- character()
  sw <- withCallingHandlers(
    hzr_stepwise(fx$fit, scope = c("nonexistent"), data = fx$data,
                 direction = "forward", criterion = "wald",
                 trace = FALSE, control = list(n_starts = 1L)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # The run must say out loud that the empty screen is not a null result.
  expect_true(any(grepl("never tested", warns, fixed = TRUE)))

  # Assert on the OBJECT, not on printed output: this is what a caller has.
  expect_equal(nrow(sw$steps), 0L)
  expect_true(sw$criteria$stopped_refit_failed)
  expect_equal(sw$criteria$n_refit_failures, 1L)
  expect_identical(sw$criteria$refit_failures, "nonexistent")
})

test_that("an honest empty screen records no refit failure", {
  # The discriminating negative control: same zero-row `steps`, but nothing
  # failed -- the candidate was tested and simply did not clear slentry.
  fx <- .refit_failure_fixture()

  sw <- hzr_stepwise(fx$fit, scope = ~ age, data = fx$data,
                     direction = "forward", criterion = "score",
                     slentry = 1e-12, trace = FALSE,
                     control = list(n_starts = 1L))

  expect_equal(nrow(sw$steps), 0L)
  expect_false(sw$criteria$stopped_refit_failed)
  expect_equal(sw$criteria$n_refit_failures, 0L)
  expect_identical(sw$criteria$refit_failures, character())
})

test_that("the trace does not claim 'no further action' when nothing was fit", {
  fx <- .refit_failure_fixture()

  sw <- suppressWarnings(
    hzr_stepwise(fx$fit, scope = c("nonexistent"), data = fx$data,
                 direction = "forward", criterion = "wald",
                 trace = FALSE, control = list(n_starts = 1L))
  )
  trace <- stepwise_trace(sw)

  expect_false(any(grepl("no further action", trace, fixed = TRUE)))
  expect_true(any(grepl("refit FAILED", trace, fixed = TRUE)))
  expect_true(any(grepl("nonexistent", trace, fixed = TRUE)))

  # print() shows exactly the trace, so the reader sees the same thing.
  expect_output(print(sw), "refit FAILED")
})

test_that("a score screen stopped on uncomputable scores says so in the trace", {
  fx <- .uncomputable_fixture()

  sw <- suppressWarnings(
    hzr_stepwise(fx$fit, scope = list(early = ~ constcol), data = fx$data,
                 direction = "forward", criterion = "score", slentry = 0.5,
                 trace = FALSE)
  )
  trace <- stepwise_trace(sw)

  expect_false(any(grepl("no further action", trace, fixed = TRUE)))
  expect_true(any(grepl("could be COMPUTED", trace, fixed = TRUE)))
})

test_that("partial refit failure keeps the steps it did accept", {
  # Partial failure is not total failure: one candidate is unfittable, the
  # other enters.  Erroring out here would discard a real accepted step, so
  # the run continues and records the failure instead.
  fx <- .refit_failure_fixture()

  sw <- suppressWarnings(
    hzr_stepwise(fx$fit, scope = c("age", "nonexistent"), data = fx$data,
                 direction = "forward", criterion = "wald", trace = FALSE,
                 control = list(n_starts = 1L))
  )

  expect_equal(nrow(sw$steps), 1L)
  expect_identical(sw$steps$variable, "age")
  expect_gt(sw$criteria$n_refit_failures, 0L)
  expect_true(all(sw$criteria$refit_failures == "nonexistent"))
})
