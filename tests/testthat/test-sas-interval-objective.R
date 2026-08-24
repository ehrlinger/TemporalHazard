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
# optimiser is expected to walk away from, whereas a data defect must not be
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
  # convert a data defect into a region the optimiser merely avoids.
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
# had to be behaviour-preserving; these pin the default against the formula it
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
# outlive the SAS licence. Every row of that fit is interval-censored, so the
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
