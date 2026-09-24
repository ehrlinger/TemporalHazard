# A "hazard" phase is -log(1 - G(t)), which diverges as G -> 1: the phase is
# UNBOUNDED BY CONSTRUCTION. Nothing constrains t_half to stay inside the
# observed support, so the optimizer can walk it below the data, the whole
# observed range lands where G is essentially 1, and the objective runs away.
# On the shipped `cabgkul` data that produced a log-likelihood of +290082 with
# converged = TRUE and no warning naming the cause (#444).
#
# The fit is NOT blocked and no bound is imposed -- that would move existing
# estimates. It warns, and records what it found in $fit$boundary.

ub_data <- function(n = 200, seed = 1) {
  set.seed(seed)
  data.frame(t = rweibull(n, 1.4, 5), s = rbinom(n, 1, 0.8))
}

ub_fit <- function(d, phases, ...) {
  suppressWarnings(hazard(time = d$t, status = d$s, dist = "multiphase",
                          phases = phases, fit = TRUE,
                          control = list(n_starts = 1L), ...))
}

test_that("a hazard phase fitted below the observed support is recorded and warned", {
  d <- ub_data()
  ph <- list(early = hzr_phase("hazard", t_half = 0.5, nu = 1, m = 0),
             late  = hzr_phase("constant"))
  # Start the early phase far below the data so the optimizer keeps it there.
  # The WARNING is what a user sees; the RECORD is what tooling reads.
  # Caught by CLASS ONLY. A bare `warning =` handler catches whichever
  # warning fires FIRST -- this fit also warns about a non-positive-definite
  # Hessian -- and would report simpleWarning, which looks like the class
  # having been stripped when it has not.
  w <- tryCatch(
    hazard(time = d$t, status = d$s, dist = "multiphase", phases = ph,
           fit = TRUE, control = list(n_starts = 1L),
           theta = c(log(0.1), log(min(d$t) / 1000), 1, 0, log(0.05))),
    hzr_unbounded_phase = function(w) w
  )
  expect_s3_class(w, "hzr_unbounded_phase")
  expect_s3_class(w, "hzr_boundary")
  expect_equal(class(w), c("hzr_unbounded_phase", "hzr_boundary", "warning",
                           "condition"))

  f <- ub_fit(d, ph, theta = c(log(0.1), log(min(d$t) / 1000), 1, 0, log(0.05)))
  b <- f$fit$boundary
  expect_true(is.list(b))
  expect_length(b, 1L)
  expect_equal(b[[1L]]$mechanism, "unbounded_phase")
  expect_equal(b[[1L]]$phase, "early")
  expect_equal(b[[1L]]$parameter, "t_half")
  # The message states a FACT; the magnitude in $detail discriminates, so no
  # threshold is tuned. 1 - G(t_min) is the mechanism itself: measured 0.053
  # for a fit hugging the edge of its data, and 3.8e-12 for the cabgkul case
  # that produced +290082.
  expect_match(b[[1L]]$detail, "t_half", fixed = TRUE)
  expect_match(b[[1L]]$detail, "1 - G", fixed = TRUE)
})

test_that("a fit that was examined and found nothing reads NULL, not NA", {
  # THE DISTINCTION THE TRI-STATE EXISTS FOR. A two-state field cannot tell
  # "looked and found nothing" from "never looked", and NULL is the common
  # case, so the confusion would be invisible.
  d <- ub_data()
  f <- ub_fit(d, list(early = hzr_phase("hazard", t_half = 3, nu = 1, m = 0),
                      late  = hzr_phase("constant")))
  # expect_null() alone would pass on a field that is ABSENT, so it cannot
  # tell "examined, found nothing" from "never implemented". Presence in
  # names() CANNOT close that here: `x$a <- NULL` DELETES the element in R,
  # so the NULL state is genuinely absent from names() -- and `$weak` does
  # the same, so this is the package's convention, not a defect of this
  # field. The distinction is pinned on its OBSERVABLE CONSEQUENCE instead:
  # NULL is "examined, nothing found" and is NOT a degraded capability,
  # while NA is "did not run" and IS one, with a cause.
  expect_null(f$fit$boundary)
  expect_false("boundary_check" %in% f$degraded)
  expect_false(.hzr_is_na_scalar(f$fit$boundary))
})

test_that("a fit that was never examined reads NA, and says why", {
  d <- ub_data()
  f <- hazard(time = d$t, status = d$s, dist = "multiphase",
              phases = list(early = hzr_phase("hazard", t_half = 3, nu = 1, m = 0),
                            late  = hzr_phase("constant")),
              fit = FALSE)
  expect_true(.hzr_is_na_scalar(f$fit$boundary))
  expect_true("boundary_check" %in% f$degraded)
  expect_match(f$degraded_causes[["boundary_check"]], "not fitted", fixed = TRUE)
})

test_that("only an unbounded phase type is checked", {
  # A cdf phase below the first observation is NOT this defect: cdf is
  # bounded, so a small t_half is an ordinary estimate, not a runaway. The
  # predicate is keyed on the TYPE, not on the parameter name.
  expect_true(.hzr_phase_type_unbounded("hazard"))
  expect_false(.hzr_phase_type_unbounded("cdf"))
  expect_false(.hzr_phase_type_unbounded("constant"))
  expect_false(.hzr_phase_type_unbounded("g3"))
  d <- ub_data()
  f <- ub_fit(d, list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0),
                      late  = hzr_phase("constant")),
              theta = c(log(0.1), log(min(d$t) / 1000), 1, 0, log(0.05)))
  # expect_null() alone would pass on a field that is ABSENT, so it cannot
  # tell "examined, found nothing" from "never implemented". Presence in
  # names() CANNOT close that here: `x$a <- NULL` DELETES the element in R,
  # so the NULL state is genuinely absent from names() -- and `$weak` does
  # the same, so this is the package's convention, not a defect of this
  # field. The distinction is pinned on its OBSERVABLE CONSEQUENCE instead:
  # NULL is "examined, nothing found" and is NOT a degraded capability,
  # while NA is "did not run" and IS one, with a cause.
  expect_null(f$fit$boundary)
  expect_false("boundary_check" %in% f$degraded)
})

test_that("the check keys on the FITTED t_half, not the starting value", {
  # A guard keyed on where a parameter STARTS is a different guard from one
  # keyed on where it CONVERGED. A fit that starts low and converges inside
  # the data is not the defect.
  d <- ub_data()
  f <- ub_fit(d, list(early = hzr_phase("hazard", t_half = 3, nu = 1, m = 0),
                      late  = hzr_phase("constant")))
  th <- exp(unname(f$fit$theta[["early.log_t_half"]]))
  expect_gt(th, min(d$t))          # converged inside the data
  # expect_null() alone would pass on a field that is ABSENT, so it cannot
  # tell "examined, found nothing" from "never implemented". Presence in
  # names() CANNOT close that here: `x$a <- NULL` DELETES the element in R,
  # so the NULL state is genuinely absent from names() -- and `$weak` does
  # the same, so this is the package's convention, not a defect of this
  # field. The distinction is pinned on its OBSERVABLE CONSEQUENCE instead:
  # NULL is "examined, nothing found" and is NOT a degraded capability,
  # while NA is "did not run" and IS one, with a cause.
  expect_null(f$fit$boundary)      # so nothing is recorded
  expect_false("boundary_check" %in% f$degraded)
})

test_that("summary reports the finding, and says nothing when there is none", {
  d <- ub_data()
  ph <- list(early = hzr_phase("hazard", t_half = 0.5, nu = 1, m = 0),
             late  = hzr_phase("constant"))
  bad <- ub_fit(d, ph, theta = c(log(0.1), log(min(d$t) / 1000), 1, 0, log(0.05)))
  out <- paste(utils::capture.output(print(summary(bad))), collapse = " ")
  expect_match(out, "unbounded", fixed = TRUE)
  ok <- ub_fit(d, list(early = hzr_phase("hazard", t_half = 3, nu = 1, m = 0),
                       late  = hzr_phase("constant")))
  out2 <- paste(utils::capture.output(print(summary(ok))), collapse = " ")
  expect_false(grepl("unbounded", out2, fixed = TRUE))
})

test_that("the suite's own marginal case is pinned, and it is marginal", {
  # THE ONE FIT IN THIS PACKAGE'S SUITE THAT TRIPS. It lives in
  # test-theta-names.R, inside a blanket suppressWarnings() covering a loop
  # over three phase specs of which only this one warns, so that file needs no
  # change and gets none -- the behaviour is pinned HERE instead, where it can
  # be asserted without making an unrelated test conditional.
  #
  # It is also the case that shows why the magnitude belongs in the record.
  # Its t_half sits a factor of 1.18 below the first observation with 5.3% of
  # its mass still beyond it -- a phase hugging the edge of its data. The
  # cabgkul fit that returned +290082 sits at 93x with 3.8e-12. Both trip the
  # same factual test; only one is degenerate, and the record says which
  # without a threshold being tuned.
  set.seed(1)
  n <- 120
  tt <- stats::rexp(n, 0.3) + 0.05
  st <- rep(1L, n)
  ph <- list(a = hzr_phase("hazard"), b = hzr_phase("constant"),
             c = hzr_phase("cdf"))
  fit <- suppressWarnings(hazard(time = tt, status = st, dist = "multiphase",
                                 phases = ph, fit = TRUE,
                                 control = list(n_starts = 1)))
  b <- fit$fit$boundary
  expect_true(is.list(b))
  expect_equal(b[[1L]]$phase, "a")
  expect_equal(b[[1L]]$mechanism, "unbounded_phase")
  t_half <- exp(unname(fit$fit$theta[["a.log_t_half"]]))
  ratio <- min(tt) / t_half
  expect_gt(ratio, 1)            # it does trip
  expect_lt(ratio, 2)            # and it is marginal, unlike cabgkul's 93x
})
