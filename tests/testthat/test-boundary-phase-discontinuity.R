# #448 -- a cdf/hazard phase whose fitted nu collapses toward zero becomes a
# step function at t_half.  The decomposition already refuses nu == 0 exactly
# (decomposition.R, "the nu -> 0 limit is degenerate"), but a fit can sit at
# nu = -1.4e-16, take the ordinary m > 0 && nu < 0 branch, and return
# converged = TRUE with no warning.  The guard was at the exact point; the
# pathology is a neighbourhood.
#
# These tests state a FACT (the phase is indistinguishable from a step at this
# data's resolution) and carry the MAGNITUDE in $detail, rather than tuning a
# threshold.  See also the family record contract on fit$fit$boundary.

avc_fit_448 <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  tt <- avc$int_dead
  st <- avc$dead
  k <- is.finite(tt) & is.finite(st) & tt > 0
  hazard(
    time = tt[k], status = st[k], dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
      late  = hzr_phase("g3", tau = 5, gamma = 1, alpha = 1, eta = 1,
                        constraint = "alpha_gamma_eta")
    ),
    fit = TRUE, control = list(n_starts = 1L)
  )
}

test_that("a cdf phase driven to nu ~ 0 warns, and the warning is classed", {
  skip_on_cran()
  # the assignment is INSIDE expect_warning: outside, `fit` would be bound to
  # the condition object and every later assertion would silently test that.
  expect_warning(
    fit <- avc_fit_448(),
    class = "hzr_phase_discontinuity"
  )
  expect_s3_class(
    tryCatch(avc_fit_448(), hzr_boundary = function(e) e),
    "hzr_boundary"
  )
  # the fit still returns: we warn, we do not refuse (SAS completes too)
  expect_true(is.finite(fit$fit$objective))
})

test_that("the discontinuity is recorded on fit$fit$boundary with its magnitude", {
  skip_on_cran()
  suppressWarnings(fit <- avc_fit_448())

  expect_true(is.list(fit$fit$boundary))
  rec <- Filter(function(r) identical(r$mechanism, "phase_discontinuity"),
                fit$fit$boundary)
  expect_length(rec, 1L)
  rec <- rec[[1L]]

  expect_identical(rec$phase, "early")
  expect_identical(rec$parameter, "nu")
  expect_true(is.character(rec$detail) && nzchar(rec$detail))

  # the magnitude does the discriminating, so it must be IN the message:
  # this fit sits at |nu| = 1.4e-16 with 5 event times at t_half.
  expect_match(rec$detail, "nu", fixed = TRUE)
  expect_match(rec$detail, "[0-9]")
})

test_that("a healthy cdf fit is examined and found clean (NULL, not NA)", {
  skip_on_cran()
  data(avc, package = "TemporalHazard", envir = environment())
  tt <- avc$int_dead
  st <- avc$dead
  k <- is.finite(tt) & is.finite(st) & tt > 0
  fit <- suppressWarnings(hazard(
    time = tt[k], status = st[k], dist = "multiphase",
    phases = list(early = hzr_phase("cdf"), late = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1L)
  ))
  # the census measured |nu| >= 0.8175 on every healthy fit, against 1.4e-16
  # here, so this must NOT trip.
  expect_true(abs(fit$fit$theta[["early.nu"]]) > 1e-3)
  # `expect_null()` alone cannot tell "examined, found nothing" from "never
  # implemented".  A presence check on names() CANNOT close that gap in R:
  # assigning NULL to a list element DELETES it, so the clean state is
  # genuinely absent from names() -- `$weak` behaves the same way, by
  # convention (verified: is.null(f$fit$weak) TRUE, "weak" %in% names FALSE).
  # The tri-state is therefore pinned on its OBSERVABLE CONSEQUENCE in the
  # degraded record: NULL is examined-and-clean and is NOT a lost capability.
  expect_null(fit$fit$boundary)
  expect_false("boundary_check" %in% fit$degraded)
})

test_that("an unfitted object did NOT look, and says so", {
  skip_on_cran()
  data(avc, package = "TemporalHazard", envir = environment())
  tt <- avc$int_dead
  st <- avc$dead
  k <- is.finite(tt) & is.finite(st) & tt > 0
  unfit <- hazard(
    time = tt[k], status = st[k], dist = "multiphase",
    phases = list(early = hzr_phase("cdf"), late = hzr_phase("constant")),
    fit = FALSE
  )
  # This is the assertion that makes the pair above non-hollow: it FAILS until
  # the capability is registered, so "clean fits are not degraded" cannot pass
  # merely because nothing ever writes the record.
  expect_true(is.na(unfit$fit$boundary))
  expect_true("boundary_check" %in% unfit$degraded)
})
