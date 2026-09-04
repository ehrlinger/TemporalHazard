# objective = "sas" -- PROC HAZARD's interval-censored contribution.
#
# Design: inst/dev/SAS-INTERVAL-OBJECTIVE-DESIGN.md, §5 tests 2-5.  Tests 1 and
# 6 (end-to-end parity and iteration-trace evaluation against SAS's printed
# log-likelihood) live with the `uslife2023` fixture.
#
# Everything here is self-contained: no SAS reference, no study mount, no RNG.
# The data are built from deterministic sequences so a change in R's random
# number stream cannot move an assertion, and the identities are algebraic, so
# they pin the formula rather than a stored number.

# --- Deterministic interval data -------------------------------------------
# 24 rows, 24 distinct widths, lower bounds strictly positive.  Distinct widths
# matter: with one shared width the -log(u - l) term is an additive constant
# and a per-row divisor is indistinguishable from a global one.
sas_iv_lo <- seq(0.10, 2.40, length.out = 24)
sas_iv_up <- sas_iv_lo + seq(0.05, 6.05, length.out = 24)
sas_iv_Ll <- seq(0.02, 1.20, length.out = 24)
sas_iv_Lu <- sas_iv_Ll + seq(1e-04, 0.85, length.out = 24)
sas_iv_w  <- seq(0.75, 40.00, length.out = 24)

sas_iv <- function(objective) {
  TemporalHazard:::.hzr_logl_interval(
    cumhaz_lower = sas_iv_Ll, cumhaz_upper = sas_iv_Lu,
    lower = sas_iv_lo, upper = sas_iv_up, weights = sas_iv_w,
    objective = objective)
}

# --- Shared multiphase spec for the threading / gradient tests --------------
sas_mp_phases <- function() {
  list(early    = hzr_phase("cdf", t_half = 0.8, nu = 0.4, m = -0.6),
       constant = hzr_phase("constant"))
}
sas_mp_cc    <- c(early = 0L, constant = 0L)
sas_mp_xl    <- list(early = NULL, constant = NULL)
sas_mp_theta <- c(log(0.06), log(0.8), 0.4, -0.6, -2.7)

# Mixed status: exact events, right-censored and interval-censored together,
# so the tests also confirm the other two branches are left alone.
sas_mp_n      <- 24L
sas_mp_status <- rep(c(1, 0, 2), length.out = sas_mp_n)
sas_mp_lo     <- seq(0.10, 2.40, length.out = sas_mp_n)
sas_mp_up     <- sas_mp_lo + seq(0.20, 4.00, length.out = sas_mp_n)
sas_mp_time   <- ifelse(sas_mp_status == 2, sas_mp_up, sas_mp_lo + 0.4)
sas_mp_w      <- seq(0.50, 5.00, length.out = sas_mp_n)

sas_mp_ll <- function(theta = sas_mp_theta, objective = "likelihood") {
  TemporalHazard:::.hzr_logl_multiphase(
    theta = theta, time = sas_mp_time, status = sas_mp_status,
    time_lower = sas_mp_lo, time_upper = sas_mp_up, weights = sas_mp_w,
    phases = sas_mp_phases(), covariate_counts = sas_mp_cc,
    x_list = sas_mp_xl, objective = objective)
}


# ---------------------------------------------------------------------------
# Test 2 -- analytic identity.  Pins the algebra with no SAS dependency.
# ---------------------------------------------------------------------------

test_that('objective = "sas" differs from the likelihood by the closed-form term', {
  dL <- sas_iv_Lu - sas_iv_Ll
  expect_equal(
    sas_iv("sas") - sas_iv("likelihood"),
    sum(sas_iv_w * log(dL / (exp(dL) - 1))) - sum(sas_iv_w * log(sas_iv_up - sas_iv_lo)),
    tolerance = 1e-12)
})

test_that('objective = "sas" is the interval-mean-hazard density term', {
  # The definitional form, written as SAS writes it: log[S(u) dLambda / (u-l)].
  # The implementation expands this in log space to avoid underflowing S(u);
  # at these magnitudes the two must agree to machine precision.
  expect_equal(
    sas_iv("sas"),
    sum(sas_iv_w * log(exp(-sas_iv_Lu) * (sas_iv_Lu - sas_iv_Ll) /
                         (sas_iv_up - sas_iv_lo))),
    tolerance = 1e-12)
})

test_that("the log-space form survives an S(u) that underflows to zero", {
  # Lambda(u) > ~745 sends exp(-Lambda(u)) to exactly 0, so the definitional
  # product form above returns -Inf.  The implementation must not.
  big <- TemporalHazard:::.hzr_logl_interval(
    cumhaz_lower = 799.5, cumhaz_upper = 800, lower = 0, upper = 1,
    weights = 1, objective = "sas")
  expect_true(is.finite(big))
  expect_equal(big, -800 + log(0.5), tolerance = 1e-12)
  expect_identical(log(exp(-800) * 0.5 / 1), -Inf)   # the form it replaces
})

test_that("the width term switches off when every interval has width one", {
  # This is what makes the US life table the clean anchor: log(u - l) = 0, so
  # the S(u)dLambda core is isolated.
  up1 <- sas_iv_lo + 1
  expect_equal(
    TemporalHazard:::.hzr_logl_interval(
      cumhaz_lower = sas_iv_Ll, cumhaz_upper = sas_iv_Lu,
      lower = sas_iv_lo, upper = up1, weights = sas_iv_w, objective = "sas"),
    sum(sas_iv_w * (-sas_iv_Lu + log(sas_iv_Lu - sas_iv_Ll))),
    tolerance = 1e-12)
})


# ---------------------------------------------------------------------------
# Test 3 -- gradient/objective consistency.  This is the check that catches the
# duplicated-formula hazard the .hzr_logl_interval() extraction removed: the
# objective and the finite-difference closure inside the gradient used to carry
# separate copies of the interval contribution.
# ---------------------------------------------------------------------------

test_that("the analytic gradient matches numDeriv under both objectives", {
  skip_if_not_installed("numDeriv")
  for (obj in c("likelihood", "sas")) {
    ana <- TemporalHazard:::.hzr_gradient_multiphase(
      theta = sas_mp_theta, time = sas_mp_time, status = sas_mp_status,
      time_lower = sas_mp_lo, time_upper = sas_mp_up, weights = sas_mp_w,
      phases = sas_mp_phases(), covariate_counts = sas_mp_cc,
      x_list = sas_mp_xl, objective = obj)
    num <- numDeriv::grad(function(th) sas_mp_ll(th, obj), sas_mp_theta)
    # The interval contribution reaches the gradient by one-sided finite
    # difference with a sqrt(.Machine$double.eps) step, so ~1e-6 absolute is
    # the accuracy floor of the method, not slack.
    expect_equal(ana, num, tolerance = 1e-5,
                 info = paste("objective =", obj))
  }
})

test_that("the two objectives disagree, so a passing gradient is not vacuous", {
  # Without this, every assertion above would still pass if `objective` were
  # threaded but ignored -- the "sas" tests would silently be re-testing the
  # likelihood against itself.
  expect_false(isTRUE(all.equal(sas_mp_ll(objective = "sas"),
                                sas_mp_ll(objective = "likelihood"))))
  expect_false(isTRUE(all.equal(sas_iv("sas"), sas_iv("likelihood"))))
})


# ---------------------------------------------------------------------------
# Test 4 -- guards.  The stop()/-Inf split is deliberate: -Inf is a value the
# optimizer is expected to walk away from, whereas a data defect must not be
# survivable.
# ---------------------------------------------------------------------------

test_that('objective = "sas" stops on a non-positive interval width', {
  expect_error(
    TemporalHazard:::.hzr_logl_interval(
      cumhaz_lower = c(0.1, 0.2), cumhaz_upper = c(0.4, 0.5),
      lower = c(1, 2), upper = c(1, 2), weights = c(1, 1), objective = "sas"),
    "requires upper > lower")
})

test_that("the width guard is checked before the feasibility guard", {
  # A row with u <= l almost always also yields dLambda <= 0.  If the order
  # were reversed the feasibility guard would return -Inf first and silently
  # convert a data defect into a region the optimizer merely avoids.
  expect_error(
    TemporalHazard:::.hzr_logl_interval(
      cumhaz_lower = 0.5, cumhaz_upper = 0.5,   # dLambda == 0 as well
      lower = 2, upper = 1, weights = 1, objective = "sas"),
    "requires upper > lower")
})

test_that("non-positive dLambda is parameter infeasibility, not an error", {
  for (obj in c("likelihood", "sas")) {
    expect_identical(
      TemporalHazard:::.hzr_logl_interval(
        cumhaz_lower = c(0.5, 0.5), cumhaz_upper = c(0.5, 0.9),
        lower = c(0, 0), upper = c(1, 1), weights = c(1, 1), objective = obj),
      -Inf, info = paste("objective =", obj))
  }
})

test_that("the width-guard message counts offenders against rows checked", {
  # The count after "of" is the number of interval rows examined, not the
  # number that failed -- reporting the failure count twice reads as though
  # every row failed.
  err <- tryCatch(
    TemporalHazard:::.hzr_logl_interval(
      cumhaz_lower = c(0.1, 0.2, 0.3), cumhaz_upper = c(0.4, 0.5, 0.6),
      lower = c(1, 2, 3), upper = c(1, 5, 3), weights = c(1, 1, 1),
      objective = "sas"),
    error = conditionMessage)
  expect_match(err, "2 of 3 interval row\\(s\\) fail this")
  expect_match(err, "at index/indices 1, 3")
})

test_that('objective = "sas" rejects left-censored rows', {
  # PROC HAZARD has no left-censoring statement, so no SAS run corresponds to
  # the result.
  expect_error(
    TemporalHazard:::.hzr_logl_multiphase(
      theta = sas_mp_theta, time = c(1, 2), status = c(-1, 2),
      time_lower = c(0, 1), time_upper = c(1, 2), weights = c(1, 1),
      phases = sas_mp_phases(), covariate_counts = sas_mp_cc,
      x_list = sas_mp_xl, objective = "sas"),
    "does not support left-censored rows")
})

test_that("the gradient refuses exactly the data the objective refuses", {
  # Guarding only the objective would leave the gradient computing happily for
  # data the objective rejects. The gradient is reachable on its own (the score
  # test calls it directly), so the objective's refusal is not guaranteed to
  # come first.
  args <- list(theta = sas_mp_theta, time = c(1, 2), status = c(-1, 2),
               time_lower = c(0, 1), time_upper = c(1, 2), weights = c(1, 1),
               phases = sas_mp_phases(), covariate_counts = sas_mp_cc,
               x_list = sas_mp_xl, objective = "sas")
  expect_error(do.call(TemporalHazard:::.hzr_logl_multiphase, args),
               "does not support left-censored rows")
  expect_error(do.call(TemporalHazard:::.hzr_gradient_multiphase, args),
               "does not support left-censored rows")

  # ...and under the default, both still accept it.
  args$objective <- "likelihood"
  expect_no_error(do.call(TemporalHazard:::.hzr_logl_multiphase, args))
  expect_no_error(do.call(TemporalHazard:::.hzr_gradient_multiphase, args))
})

test_that('objective = "sas" applies only to dist = "multiphase"', {
  expect_error(
    hazard(time = 1:5, status = rep(1, 5), dist = "weibull",
           objective = "sas", fit = FALSE),
    'applies only to dist = "multiphase"')
})

test_that("an unknown objective is rejected", {
  expect_error(
    TemporalHazard:::.hzr_logl_interval(
      cumhaz_lower = 0.1, cumhaz_upper = 0.2, lower = 0, upper = 1,
      weights = 1, objective = "probability"),
    "should be one of")
})


# ---------------------------------------------------------------------------
# Test 5 -- the default is unchanged.  The extraction of .hzr_logl_interval()
# had to be behavior-preserving; these pin the default against the formula it
# replaced rather than against a stored number.
# ---------------------------------------------------------------------------

test_that("the default objective is the interval probability", {
  expect_equal(
    sas_iv("likelihood"),
    sum(sas_iv_w * (-sas_iv_Ll +
                      TemporalHazard::hzr_log1mexp(sas_iv_Lu - sas_iv_Ll))),
    tolerance = 1e-12)
  # log(S(l) - S(u)), the same quantity written directly.
  expect_equal(
    sas_iv("likelihood"),
    sum(sas_iv_w * log(exp(-sas_iv_Ll) - exp(-sas_iv_Lu))),
    tolerance = 1e-10)
})

test_that("omitting objective is identical to naming the default", {
  expect_identical(
    TemporalHazard:::.hzr_logl_interval(
      cumhaz_lower = sas_iv_Ll, cumhaz_upper = sas_iv_Lu,
      lower = sas_iv_lo, upper = sas_iv_up, weights = sas_iv_w),
    sas_iv("likelihood"))
  expect_identical(sas_mp_ll(), sas_mp_ll(objective = "likelihood"))
})


# ---------------------------------------------------------------------------
# Tests 1 and 6 -- end-to-end parity against SAS's printed log-likelihood, at
# the optimum and at three off-optimum iterates.
#
# Reference: /studies/general/uslife/table2023/distributions/hz.icall.lst,
# reproduced here through the shipped `uslife2023` fixture so the assertions
# outlive the SAS license. Every row of that fit is interval-censored, so the
# whole log-likelihood IS the interval contribution -- there is no exact-event
# or right-censored term to net out.
#
# SAS prints the log-likelihood to 6 significant figures, so the admissible
# error at this magnitude is a half-ulp: +-0.5. A tighter tolerance would be
# asserting more precision than SAS reported.
# ---------------------------------------------------------------------------

# Cumulative hazard of the three-phase model. With TAU = GAMMA = ALPHA = 1 the
# late phase is a pure Weibull, G3 = t^eta. Built from the package's own shape
# functions rather than re-derived.
sas_uslife_Lambda <- function(t, mue, thalf, nu, muc, mul, eta) {
  mue * hzr_decompos(t, t_half = thalf, nu = nu, m = 0)$G +
    muc * t +
    mul * hzr_decompos_g3(t, tau = 1, gamma = 1, alpha = 1, eta = eta)$G3
}

test_that("the uslife2023 fixture gates on the SAS job's own printed figures", {
  # Protocol: gate the cohort before trusting any parity number. These are the
  # figures the .lst prints, not values recomputed from the fixture.
  expect_identical(nrow(uslife2023), 124L)
  expect_true(all(uslife2023$age_u - uslife2023$age_l == 1))
  expect_equal(sum(uslife2023$d_all), 100000.0125, tolerance = 1e-9)
  expect_equal(min(uslife2023$d_all), 0.2352, tolerance = 1e-4)
  expect_equal(max(uslife2023$d_all), 3620.335, tolerance = 1e-4)
})

test_that('objective = "sas" reproduces SAS hz.icall at its final estimates', {
  lo <- uslife2023$age_l
  up <- uslife2023$age_u
  w  <- uslife2023$d_all
  # Natural-scale estimates as printed in the .lst parameter table.
  lam <- function(t) {
    sas_uslife_Lambda(t, mue = 4.504471e-03, thalf = 1.533746e-03,
                      nu = -2.159600, muc = 1.082926e-03,
                      mul = 5.309996e-17, eta = 8.386895)
  }
  ll <- TemporalHazard:::.hzr_logl_interval(
    lam(lo), lam(up), lo, up, w, objective = "sas")

  expect_lt(abs(ll - (-410414)), 0.5)             # SAS print precision
  expect_lt(abs(ll - (-410414)), 0.01)            # regression canary: -410414.0025

  # Conservation of events: SAS reports sum_i d_i Lambda(u_i) = 100000, which
  # validates the whole cumulative hazard at once rather than one point of it.
  expect_equal(sum(w * lam(up)), 100000, tolerance = 1e-6)

  # The default objective must NOT match. Without this, the test would also
  # pass if `objective` were ignored and both branches computed the same thing.
  ll_lik <- TemporalHazard:::.hzr_logl_interval(
    lam(lo), lam(up), lo, up, w, objective = "likelihood")
  expect_lt(abs(ll_lik - (-406268)), 0.5)
  expect_gt(abs(ll_lik - (-410414)), 4000)
})

test_that('objective = "sas" tracks SAS along its own iteration trace', {
  # At the optimum a wrong objective and a different optimum are
  # indistinguishable. Off-optimum iterates separate them, and cost four
  # objective evaluations with no refitting.
  #
  # Trace columns are the internal scale: E2 = log(t_half), E3 with
  # nu = -exp(E3), E0/C0/L0 = log(mu) per phase, L4 = log(eta). The nu mapping
  # is read off the .lst directly: -exp(0.7699208) = -2.15960 against a printed
  # NU of -2.159600.
  trace <- list(
    list(it =  0, sas = -410612,
         p = c(-5.15857, 1.0677340, 2.102790, -5.21282, -7.18509, -36.5550)),
    list(it =  4, sas = -410419,
         p = c(-5.17612, 1.0637730, 2.124433, -5.38322, -6.83621, -37.3892)),
    list(it =  8, sas = -410414,
         p = c(-5.83365, 0.9155983, 2.126492, -5.39579, -6.82636, -37.4678)),
    list(it = 13, sas = -410414,
         p = c(-6.48004, 0.7699208, 2.126670, -5.40268, -6.82809, -37.4744))
  )
  lo <- uslife2023$age_l
  up <- uslife2023$age_u
  w  <- uslife2023$d_all

  for (tr in trace) {
    p <- tr$p
    lam <- function(t) {
      sas_uslife_Lambda(t, mue = exp(p[4]), thalf = exp(p[1]),
                        nu = -exp(p[2]), muc = exp(p[5]),
                        mul = exp(p[6]), eta = exp(p[3]))
    }
    ll <- TemporalHazard:::.hzr_logl_interval(
      lam(lo), lam(up), lo, up, w, objective = "sas")
    expect_lt(abs(ll - tr$sas), 0.5)
  }

  # The trace spans 198 log-likelihood units, so agreement across it is not
  # an artefact of every iterate sitting near the same value.
  expect_gt(abs(trace[[1]]$sas - trace[[4]]$sas), 190)
})


# ---------------------------------------------------------------------------
# Objective is recorded on the fit and survives every refit path.
#
# `objective` was accepted, match.arg'd and threaded to the optimizer, but
# never stored -- so `.hzr_refit_with_scope()`, which rebuilds the hazard()
# call from `$spec`, refit every stepwise candidate under the likelihood while
# the base fit's own objective was the SAS density. That differences
# `delta_logLik` and `aic` across two estimands and reports a full `$steps`
# table with no warning. The tests below assert the recording, the default for
# objects that predate the argument, and the forwarding.
# ---------------------------------------------------------------------------

sas_obj_data <- function(n = 60, seed = 11) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  data.frame(
    tt = stats::rexp(n, 0.4) + 0.05,
    ev = rep(c(1, 0), length.out = n),
    x1 = x1,
    x2 = stats::rnorm(n)
  )
}

sas_obj_phases <- function() {
  list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
       const = hzr_phase("constant"))
}

test_that("hazard() records the objective it was fitted under", {
  skip_on_cran()
  D <- sas_obj_data()
  f_lik <- suppressWarnings(hazard(time = D$tt, status = D$ev,
                                   dist = "multiphase",
                                   phases = sas_obj_phases(), fit = TRUE))
  f_sas <- suppressWarnings(hazard(time = D$tt, status = D$ev,
                                   dist = "multiphase",
                                   phases = sas_obj_phases(),
                                   objective = "sas", fit = TRUE))
  # Distinct values, not merely non-NULL: a spec that hardcoded either one
  # would pass a `!is.null()` check.
  expect_identical(f_lik$spec$objective, "likelihood")
  expect_identical(f_sas$spec$objective, "sas")
})

test_that(".hzr_fit_objective() defaults only when nothing was recorded", {
  skip_on_cran()
  D <- sas_obj_data()
  f <- suppressWarnings(hazard(time = D$tt, status = D$ev, dist = "multiphase",
                               phases = sas_obj_phases(), objective = "sas",
                               fit = TRUE))
  expect_identical(TemporalHazard:::.hzr_fit_objective(f), "sas")

  # A fit from before the argument existed carries no `objective`; it was
  # necessarily estimated under the likelihood.
  legacy <- f
  legacy$spec$objective <- NULL
  expect_identical(TemporalHazard:::.hzr_fit_objective(legacy), "likelihood")
})

test_that("a scope refit keeps the base fit's objective", {
  skip_on_cran()
  D <- sas_obj_data()
  f_sas <- suppressWarnings(hazard(Surv(tt, ev) ~ 1, data = D,
                                   dist = "multiphase",
                                   phases = sas_obj_phases(),
                                   objective = "sas", fit = TRUE))
  expect_identical(f_sas$spec$objective, "sas")

  out <- suppressWarnings(
    TemporalHazard:::.hzr_refit_with_scope(f_sas, action = "add", var = "x1",
                                           data = D, phase = "early"))
  # The refit is where the estimand used to change silently.
  expect_identical(out$spec$objective, "sas")

  f_lik <- suppressWarnings(hazard(Surv(tt, ev) ~ 1, data = D,
                                   dist = "multiphase",
                                   phases = sas_obj_phases(), fit = TRUE))
  out_lik <- suppressWarnings(
    TemporalHazard:::.hzr_refit_with_scope(f_lik, action = "add", var = "x1",
                                           data = D, phase = "early"))
  expect_identical(out_lik$spec$objective, "likelihood")
})

test_that("the two objectives are far enough apart for the mix-up to matter", {
  # Guards the tests above against becoming cosmetic: if the two objectives
  # ever agreed numerically, forwarding the label would protect nothing and a
  # test asserting only the label would still pass.
  expect_gt(abs(sas_mp_ll(objective = "sas") -
                  sas_mp_ll(objective = "likelihood")), 1)
})

test_that("a vector-interface sas fit keeps its objective across a scope refit", {
  skip_on_cran()
  # The intersection of #160 (multiphase models refit from a vector-interface
  # base fit) and the objective fix. Neither change tested it on its own, and
  # it is the combination SAS SELECTION jobs actually land on: 25 of 26 use
  # ICENSOR, so they arrive interval-censored, through the vector interface,
  # and are the reason objective = "sas" exists.
  set.seed(7)
  n <- 60L
  D <- data.frame(tt = stats::rexp(n, 0.4) + 0.05,
                  x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  D$st <- rep(c(1, 0, 2), length.out = n)
  D$lo <- ifelse(D$st == 2, D$tt, 0.02)
  D$up <- ifelse(D$st == 2, D$tt + 0.6, D$tt)

  # Guard the guard: interval rows must be present or objective = "sas" is
  # indistinguishable from the likelihood and this test proves nothing.
  expect_gt(sum(D$st == 2), 0L)

  f <- suppressWarnings(hazard(time = D$tt, status = D$st, time_lower = D$lo,
                               time_upper = D$up, dist = "multiphase",
                               phases = sas_obj_phases(), objective = "sas",
                               fit = TRUE))
  expect_null(f$call$formula)
  expect_identical(f$spec$objective, "sas")

  out <- suppressWarnings(
    TemporalHazard:::.hzr_refit_with_scope(f, action = "add", var = "x1",
                                           data = D, phase = "early"))
  # Both halves have to survive together: the response vectors rebuilt from
  # `$data` (#160) and the estimand read from `$spec` (this fix).
  expect_identical(out$spec$objective, "sas")
  expect_equal(out$data$time_upper, D$up)
  expect_equal(out$data$time_lower, D$lo)
})

test_that("a refit refuses a conflicting objective and tolerates a redundant one", {
  skip_on_cran()
  # `objective` reaches .hzr_refit_with_scope() through `...` -- hzr_stepwise()
  # documents `...` as "Passed to the underlying hazard() refits" -- while the
  # refit also supplies it from the base fit. Before the guard this matched the
  # same formal twice and died inside do.call() with a message about argument
  # matching that says nothing about estimands.
  D <- sas_obj_data()
  f <- suppressWarnings(hazard(Surv(tt, ev) ~ 1, data = D, dist = "multiphase",
                               phases = sas_obj_phases(), objective = "sas",
                               fit = TRUE))
  refit <- function(...) {
    TemporalHazard:::.hzr_refit_with_scope(f, action = "add", var = "x1",
                                           data = D, phase = "early", ...)
  }

  # Conflicting: refused, and the message has to say WHY a refit cannot change
  # estimand or the next reader deletes the guard as pedantry.
  expect_error(suppressWarnings(refit(objective = "likelihood")),
               "cannot be changed in a refit")
  expect_error(suppressWarnings(refit(objective = "likelihood")),
               "between two estimands rather than between two models")
  # Specifically NOT the bare R argument-matching error it used to raise: the
  # regression would still "error", so asserting only that it errors would pass
  # against the broken code.
  msg <- tryCatch(suppressWarnings(refit(objective = "likelihood")),
                  error = conditionMessage)
  expect_false(grepl("matched by multiple actual arguments", msg, fixed = TRUE))

  # Redundant: dropped, not refused -- restating the base fit's own objective
  # asks for nothing the refit was not already going to do.
  out <- suppressWarnings(refit(objective = "sas"))
  expect_identical(out$spec$objective, "sas")

  # Absent: unchanged.
  expect_identical(suppressWarnings(refit())$spec$objective, "sas")
})

# ---------------------------------------------------------------------------
# Test 5 -- the SAS data guards are entry guards (#213).
#
# Both conditions are pure functions of the data: identical at every start, so
# reporting them through the optimizer's per-start tryCatch framed a data
# defect as a convergence problem and invited raising `n_starts`, a remedy that
# cannot work. They are checked in hazard() before any optimization.
# ---------------------------------------------------------------------------

test_that("left-censored rows under sas are reported as a data defect, not a fit failure", {
  ph <- list(early = hzr_phase("cdf"), late = hzr_phase("constant"))
  err <- tryCatch(
    suppressWarnings(hazard(
      time = c(1, 2, 3, 4, 5, 6), status = c(1, 0, -1, 1, 0, 1),
      dist = "multiphase", phases = ph, fit = TRUE, objective = "sas")),
    error = conditionMessage)

  expect_match(err, "does not support left-censored rows")
  # The defect: the message used to arrive wrapped in "produced no usable fit
  # from N starts: N errored", which reads as a convergence problem.
  expect_no_match(err, "no usable fit")
  expect_no_match(err, "start")
})

test_that("the sas guards fire with fit = FALSE as well", {
  # The combination can never produce a fit, so it is rejected when it is
  # specified rather than when it is first evaluated. Before this, fit = FALSE
  # returned a normal hazard object carrying an objective it could never use.
  ph <- list(early = hzr_phase("cdf"), late = hzr_phase("constant"))
  expect_error(
    suppressWarnings(hazard(
      time = c(1, 2, 3, 4, 5, 6), status = c(1, 0, -1, 1, 0, 1),
      dist = "multiphase", phases = ph, fit = FALSE, objective = "sas")),
    "does not support left-censored rows")
})

test_that("the interval-width guard reports the data row, not the interval-subset row", {
  # The inner guard sees only the interval rows, so it numbered the offender 1.
  # Checked on the full data, the same row is 3 -- which is the number the user
  # can look up.
  ph <- list(early = hzr_phase("cdf"), late = hzr_phase("constant"))
  err <- tryCatch(
    suppressWarnings(hazard(
      time = c(1, 2, 3, 4, 5, 6), status = c(1, 0, 2, 1, 0, 1),
      time_lower = c(1, 2, 3, 4, 5, 6), time_upper = c(1, 2, 3, 4, 5, 6),
      dist = "multiphase", phases = ph, fit = TRUE, objective = "sas")),
    error = conditionMessage)

  expect_match(err, "requires upper > lower")
  expect_match(err, "at index/indices 3")
  expect_no_match(err, "no usable fit")
})

test_that("the entry guards are specific to objective = 'sas'", {
  # Negative control: the same data are legitimate under the default
  # objective, which is the statistically consistent one. A guard that fired
  # here would be rejecting valid work.
  ph <- list(early = hzr_phase("cdf"), late = hzr_phase("constant"))
  expect_no_error(suppressWarnings(hazard(
    time = c(1, 2, 3, 4, 5, 6), status = c(1, 0, -1, 1, 0, 1),
    dist = "multiphase", phases = ph, fit = FALSE, objective = "likelihood")))
})

test_that("the entry guards accept data that sas does support", {
  # Negative control against over-firing: interval rows with a positive width
  # and no left-censoring are exactly what objective = "sas" is for. fit = TRUE
  # so the objective is actually evaluated on the data the guard declared
  # valid -- at fit = FALSE the guard's verdict is never tested against the
  # objective's.
  ph <- list(early = hzr_phase("cdf"), late = hzr_phase("constant"))
  expect_no_error(suppressWarnings(hazard(
    time = c(1, 2, 3, 4, 5, 6), status = c(1, 0, 2, 1, 0, 1),
    time_lower = c(1, 2, 2.0, 4, 5, 6), time_upper = c(1, 2, 3.5, 4, 5, 6),
    dist = "multiphase", phases = ph, fit = TRUE, objective = "sas")))
})

test_that("interval rows with no bounds supplied are caught at entry", {
  # Both bounds default to `time`, so every interval row has zero width. This
  # is the commonest way to hit the guard and the earlier tests missed it:
  # they passed time_lower/time_upper explicitly, so an entry check that gave
  # up whenever a bound was NULL passed the whole file.
  ph <- list(early = hzr_phase("cdf"), late = hzr_phase("constant"))
  err <- tryCatch(
    suppressWarnings(hazard(
      time = c(1, 2, 3, 4, 5, 6), status = c(1, 0, 2, 1, 0, 1),
      dist = "multiphase", phases = ph, fit = TRUE, objective = "sas")),
    error = conditionMessage)
  expect_match(err, "requires upper > lower")
  expect_match(err, "at index/indices 3")
  expect_no_match(err, "no usable fit")
})

test_that("the entry check reads the bounds, not `time`", {
  # Pins the parity claim in .hzr_check_sas_data()'s @details. Row 3 has a
  # valid interval (5, 7) but a `time` of 1, so an implementation that
  # normalised from `time` instead of the supplied bounds would reject it;
  # row 6 is the mirror -- bounds are degenerate while `time` looks fine.
  ph <- list(early = hzr_phase("cdf"), late = hzr_phase("constant"))
  expect_no_error(suppressWarnings(hazard(
    time        = c(1, 2, 1, 4, 5, 6), status = c(1, 0, 2, 1, 0, 1),
    time_lower  = c(1, 2, 5, 4, 5, 6),
    time_upper  = c(1, 2, 7, 4, 5, 6),
    dist = "multiphase", phases = ph, fit = FALSE, objective = "sas")))

  err <- tryCatch(
    suppressWarnings(hazard(
      time        = c(1, 2, 3, 4, 5, 6), status = c(1, 0, 1, 1, 0, 2),
      time_lower  = c(1, 2, 3, 4, 5, 9),
      time_upper  = c(1, 2, 3, 4, 5, 9),
      dist = "multiphase", phases = ph, fit = FALSE, objective = "sas")),
    error = conditionMessage)
  expect_match(err, "requires upper > lower")
  expect_match(err, "at index/indices 6")
})

test_that("the entry check and the objective agree on which rows offend", {
  # The @details claims the two cannot disagree. This must call the objective
  # DIRECTLY: hazard() now runs .hzr_check_sas_data() before optimization, so
  # comparing the entry check against hazard() would compare the entry check
  # with itself and pass no matter what the objective did.
  tt <- c(1, 2, 3, 4, 5, 6)
  st <- c(1, 0, 2, 1, 0, 1)
  objective_errors <- function(lo, up) {
    tryCatch({
      TemporalHazard:::.hzr_logl_multiphase(
        theta = sas_mp_theta, time = tt, status = st,
        time_lower = lo, time_upper = up, weights = rep(1, length(tt)),
        phases = sas_mp_phases(), covariate_counts = sas_mp_cc,
        x_list = sas_mp_xl, objective = "sas")
      FALSE
    }, error = function(e) TRUE)
  }
  entry_errors <- function(lo, up) {
    tryCatch({
      TemporalHazard:::.hzr_check_sas_data(st, tt, lo, up, "sas")
      FALSE
    }, error = function(e) TRUE)
  }

  cases <- list(
    # Offenders: every interval row has zero width once bounds default to time.
    list(lo = NULL,                    up = NULL,                    bad = TRUE),
    list(lo = c(1, 2, 3, 4, 5, 6),     up = NULL,                    bad = TRUE),
    list(lo = NULL,                    up = c(1, 2, 3, 4, 5, 6),     bad = TRUE),
    list(lo = c(1, 2, 3, 4, 5, 6),     up = c(1, 2, 3, 4, 5, 6),     bad = TRUE),
    # Accepted: row 3 is the only interval row and has a positive width.
    list(lo = c(1, 2, 2.0, 4, 5, 6),   up = c(1, 2, 3.5, 4, 5, 6),   bad = FALSE),
    # Accepted, and `time` disagrees with the bounds on the interval row --
    # an implementation normalising from `time` would reject this.
    list(lo = c(1, 2, 5.0, 4, 5, 6),   up = c(1, 2, 7.0, 4, 5, 6),   bad = FALSE)
  )

  for (i in seq_along(cases)) {
    cs <- cases[[i]]
    entry <- entry_errors(cs$lo, cs$up)
    inner <- objective_errors(cs$lo, cs$up)
    expect_identical(entry, inner, info = paste("case", i))
    expect_identical(entry, cs$bad, info = paste("case", i))
  }
})

test_that("an NA status under sas names the argument and the row", {
  # any(status == -1) is NA-poisoned, so without this the inner guard's `if`
  # threw a bare "missing value where TRUE/FALSE needed".
  ph <- list(early = hzr_phase("cdf"), late = hzr_phase("constant"))
  err <- tryCatch(
    suppressWarnings(hazard(
      time = c(1, 2, 3, 4, 5, 6), status = c(1, 0, NA, 1, 0, 1),
      dist = "multiphase", phases = ph, fit = FALSE, objective = "sas")),
    error = conditionMessage)
  expect_match(err, "complete 'status' vector")
  expect_match(err, "at index/indices 3")
  expect_no_match(err, "missing value where TRUE/FALSE needed")
})
