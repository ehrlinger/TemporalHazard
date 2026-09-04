# A phase can leave the model without leaving the output. These tests pin the
# two ways it happens, which have DIFFERENT consequences:
#
#   absent    -- the phase never starts; mu and shape are both unidentified.
#   saturated -- the phase finished before the first observation; it is then a
#                constant offset, so mu stays well identified and only the
#                shape parameters go flat.
#
# The second is the one that misled a parity investigation (temporal_hazard#143,
# where the report had it as "mu stops being identified" -- mu is in fact the
# one parameter that survives).

ident_ns <- function(f) get(f, asNamespace("TemporalHazard"))

ident_setup <- function(t_half, nu = 0, m = -0.4) {
  ph <- ident_ns(".hzr_validate_phases")(
    list(early    = hzr_phase("cdf", t_half = t_half, nu = nu, m = m),
         constant = hzr_phase("constant")))
  list(phases = ph,
       cc     = c(early = 0L, constant = 0L),
       xl     = list(early = NULL, constant = NULL),
       theta  = c(log(0.045), log(t_half), nu, m, log(0.036)))
}

ident_time <- function() {
  data("avc", package = "TemporalHazard", envir = environment())
  t <- avc$int_dead
  t[t > 0]
}

ident_check <- function(t_half, nu = 0, m = -0.4, tol = 1e-8) {
  s <- ident_setup(t_half, nu, m)
  ident_ns(".hzr_check_phase_identifiability")(
    s$theta, ident_time(), s$phases, s$cc, s$xl, tol = tol)
}


test_that("a healthy two-phase fit raises nothing", {
  # t_half = 0.003 is temporal_hazard#143's own reproducer. The early phase is
  # still climbing across the observed range, so nothing is unidentified -- the
  # guard must not fire here, or it would confirm the report it was written to
  # correct.
  expect_no_warning(ident_check(0.003))
  sh <- suppressWarnings(ident_check(0.003))
  expect_gt(sh$variation[sh$phase == "early"], 0.5)
})

test_that("a saturated phase warns, and says mu is the survivor", {
  w <- capture_warnings(ident_check(1e-6))
  expect_length(w, 1L)
  expect_match(w, "constant across the observed times")
  expect_match(w, "'mu' remains identified but the shape parameters do not")
  expect_match(w, "'early'")
})

test_that("saturation is detected as exactly zero variation", {
  sh <- suppressWarnings(ident_check(1e-6))
  expect_identical(sh$variation[sh$phase == "early"], 0)
  # ...and mu is NOT flagged absent: the phase still supplies most of Lambda.
  expect_gt(sh$share[sh$phase == "early"], 0.5)
})

test_that("an absent phase warns about mu AND shape", {
  # Half-life far beyond follow-up: G(t) never lifts off, so the phase
  # contributes essentially none of Lambda anywhere.
  w <- capture_warnings(ident_check(1e8))
  expect_true(any(grepl("has not started by the end of follow-up", w)))
  expect_true(any(grepl("neither its 'mu' nor its shape", w)))
})

test_that("an absent phase is reported once, as absent rather than as flat", {
  # A phase that never starts is trivially also flat. Reporting both would be
  # two warnings for one condition, and "absent" is the more informative.
  w <- capture_warnings(ident_check(1e8))
  expect_length(w, 1L)
})

test_that("the tolerance is honoured in both directions", {
  expect_no_warning(ident_check(1e-6, tol = 0))       # nothing can fall below 0
  expect_warning(ident_check(0.003, tol = 0.9),       # 0.843 < 0.9
                 "constant across the observed times")
})

test_that("variation is withheld, not guessed, for a phase carrying covariates", {
  # With covariates mu differs by row, so the spread of the contribution mixes
  # covariate variation with the shape's. Reporting it would invite reading
  # covariate spread as shape identifiability.
  data("avc", package = "TemporalHazard", envir = environment())
  ph <- ident_ns(".hzr_validate_phases")(
    list(early    = hzr_phase("cdf", t_half = 1e-6, nu = 0, m = -0.4),
         constant = hzr_phase("constant")))
  xm <- matrix(as.numeric(seq_len(nrow(avc)) %% 2), ncol = 1,
               dimnames = list(NULL, "grp"))
  sh <- ident_ns(".hzr_phase_shares")(
    c(log(0.045), log(1e-6), 0, -0.4, 0.3, log(0.036)),
    avc$int_dead,
    ph, c(early = 1L, constant = 0L), list(early = xm, constant = NULL))
  expect_true(is.na(sh$variation[sh$phase == "early"]))
  expect_false(is.na(sh$share[sh$phase == "early"]))
})

test_that("the fit carries the shares so the warning is checkable", {
  skip_if_not_installed("survival")
  data("avc", package = "TemporalHazard", envir = environment())
  ph <- list(early = hzr_phase("cdf", t_half = 0.003, nu = 0, m = -0.4,
                               fixed = c("nu", "m", "t_half")),
             constant = hzr_phase("constant"))
  f <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
              dist = "multiphase", phases = ph,
              theta = c(log(0.045), log(0.003), 0, -0.4, log(0.036)),
              fit = TRUE, control = list(n_starts = 1, conserve = TRUE))
  expect_s3_class(f$fit$phase_share, "data.frame")
  expect_setequal(f$fit$phase_share$phase, c("early", "constant"))
})


# --- control$phase_share_tol validation -------------------------------------
# NA is the dangerous input, not the obviously-wrong one: `shares < NA` is NA,
# which() drops it, and the guard would pass every phase while looking like it
# ran. A check written against silent failure must not fail silently itself.

ident_fit_tol <- function(tol) {
  data("avc", package = "TemporalHazard", envir = environment())
  ph <- list(early = hzr_phase("cdf", t_half = 1e-6, nu = 0, m = -0.4,
                               fixed = c("nu", "m", "t_half")),
             constant = hzr_phase("constant"))
  hazard(survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
         phases = ph, theta = c(log(0.045), log(1e-6), 0, -0.4, log(0.036)),
         fit = TRUE,
         control = list(n_starts = 1, conserve = TRUE, phase_share_tol = tol))
}

test_that("a bad phase_share_tol is rejected, naming the option", {
  skip_if_not_installed("survival")
  for (bad in list(NA_real_, NA, "1e-8", -1, c(1e-8, 1e-9), numeric(0), Inf)) {
    # fixed = TRUE: the option name contains a `$`, and escaping it for a
    # regex is exactly the kind of quiet mismatch this test exists to catch.
    expect_error(ident_fit_tol(bad), "'control$phase_share_tol'",
                 fixed = TRUE, info = deparse(bad))
  }
})

test_that("phase_share_tol = 0 silences the check rather than erroring", {
  skip_if_not_installed("survival")
  expect_no_warning(ident_fit_tol(0))
})

test_that("a valid phase_share_tol still warns", {
  skip_if_not_installed("survival")
  expect_warning(ident_fit_tol(1e-8), "constant across the observed times")
})


# --- return type is stable --------------------------------------------------

test_that("shares are a data.frame even when nothing can be measured", {
  # No observed time carries a positive total, so there is nothing to measure.
  # The frame must still have the same shape -- fit$phase_share is one type for
  # a caller to handle, not two -- and NA must mean "not measured" rather than
  # being confusable with a measured zero.
  s <- ident_setup(0.003)
  sh <- ident_ns(".hzr_phase_shares")(s$theta, numeric(0), s$phases, s$cc, s$xl)
  expect_s3_class(sh, "data.frame")
  expect_identical(sh$phase, c("early", "constant"))
  expect_true(all(is.na(sh$share)))
  expect_true(all(is.na(sh$variation)))
})

test_that("unmeasurable shares warn about nothing", {
  s <- ident_setup(0.003)
  expect_no_warning(
    ident_ns(".hzr_check_phase_identifiability")(
      s$theta, numeric(0), s$phases, s$cc, s$xl))
})


# ---------------------------------------------------------------------------
# One distinct observation time (#211).
#
# `variation` is (max - min) / max of a phase's contribution, which is 0 for
# EVERY phase when there is only one distinct time -- by construction, not
# because any phase saturated. The scan then fired and described a cause that
# had not occurred, including telling a `constant` phase that its shape
# parameters were unidentified when it has none. The real condition is worse
# than either of the two being reported and went unnamed.
# ---------------------------------------------------------------------------

ident_tied <- function(times) {
  s <- ident_setup(t_half = 0.5)
  warns <- character(0)
  withCallingHandlers(
    ident_ns(".hzr_check_phase_identifiability")(
      s$theta, times, s$phases, s$cc, s$xl),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  warns
}

test_that("tied observation times are reported as their own condition", {
  w <- ident_tied(rep(2, 20))

  expect_length(w, 1L)
  expect_match(w, "one distinct observation time")

  # The three claims that were wrong for this input.
  expect_no_match(w, "shape parameters")
  expect_no_match(w, "finished before the first observation")
  expect_no_match(w, "half-life")
})

test_that("a single-row fit reports the same condition", {
  w <- ident_tied(2)
  expect_length(w, 1L)
  expect_match(w, "one distinct observation time")
})

test_that("the tied-time condition names both phases, not one", {
  # Every phase is unidentified here, so reporting one of them would imply the
  # other is fine.
  w <- ident_tied(rep(2, 20))
  expect_match(w, "early")
  expect_match(w, "constant")
})

test_that("the tied-time guard still returns the shares", {
  # summary() and fit$phase_share read this; an early return must not drop it.
  s <- ident_setup(t_half = 0.5)
  sh <- suppressWarnings(ident_ns(".hzr_check_phase_identifiability")(
    s$theta, rep(2, 20), s$phases, s$cc, s$xl))
  expect_s3_class(sh, "data.frame")
  expect_identical(sh$phase, c("early", "constant"))
})

test_that("two distinct times are enough to run the ordinary scan", {
  # Boundary: the guard must not swallow the real diagnostics one time point
  # too eagerly. With t_half = 0.5 and times 2 and 3 the early phase has long
  # since saturated, which is the ordinary saturated message.
  w <- ident_tied(c(2, 3))
  expect_false(any(grepl("one distinct observation time", w)))
})

test_that("no measurable rows is not the same condition as one distinct time", {
  # Both leave `variation` uninformative, but they are different problems and
  # only one of them is worth a warning. The distinct-time count is NA when
  # nothing was measured, so it cannot be confused with a measured 0 -- the
  # same NA-versus-zero rule the shares themselves follow.
  s <- ident_setup(0.003)
  none <- ident_ns(".hzr_phase_shares")(s$theta, numeric(0), s$phases, s$cc, s$xl)
  tied <- ident_ns(".hzr_phase_shares")(s$theta, rep(2, 20), s$phases, s$cc, s$xl)

  expect_identical(attr(none, "n_distinct_time"), NA_integer_)
  expect_identical(attr(tied, "n_distinct_time"), 1L)
})

test_that("the distinct-time count is taken over measured rows, not all rows", {
  # Rows at t = 0 carry a zero total, so `ok` drops them and no measure is
  # taken there. `time` then has two distinct values while the measurements
  # rest on one -- and it is the measurements that make `variation` degenerate.
  # Counting over all of `time` would call this fit measurable and hand it back
  # to the saturation scan, which is the misdiagnosis being fixed.
  s <- ident_setup(t_half = 0.5)
  tt <- c(0, 0, 2, 2)
  sh <- ident_ns(".hzr_phase_shares")(s$theta, tt, s$phases, s$cc, s$xl)

  expect_identical(length(unique(tt)), 2L)
  expect_identical(attr(sh, "n_distinct_time"), 1L)

  w <- character(0)
  withCallingHandlers(
    ident_ns(".hzr_check_phase_identifiability")(
      s$theta, tt, s$phases, s$cc, s$xl),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_length(w, 1L)
  expect_match(w, "one distinct observation time")
})
