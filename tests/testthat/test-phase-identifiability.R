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
