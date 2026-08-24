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
