# A "hazard" phase is -log(1 - G(t)), which diverges as G -> 1: the phase is
# UNBOUNDED BY CONSTRUCTION. Nothing constrains t_half to stay inside the
# observed support, so the optimizer can walk it below the data, the whole
# observed range lands where G is essentially 1, and the objective runs away.
# On the shipped `cabgkul` data that produced a log-likelihood of +290082 with
# converged = TRUE and no warning naming the cause (#444).
#
# The fit is NOT blocked and no bound is imposed -- that would move existing
# estimates. It warns, and records what it found in $fit$boundary.

# ON THE TRI-STATE ASSERTIONS BELOW. expect_null() passes on a field that is
# ABSENT, so it cannot tell "examined, found nothing" from "never
# implemented". Presence in names() cannot close that here either: `x$a <-
# NULL` DELETES the element in R, and `$weak` behaves the same, so the NULL
# state is genuinely absent from names() by the package's own convention.
# The distinction is therefore pinned on its OBSERVABLE CONSEQUENCE --
# NULL is not a degraded capability, NA is one and carries a cause.
#
# Numbers in $detail are asserted ANCHORED to their phrase. On one fixture
# 1 - G underflows to exactly 0, so format(., digits = 4) is "0", which is a
# substring of nearly every other figure in the message: an unanchored match
# passed while a mutation reporting G instead of 1 - G survived.

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
  # ASSERT THE NUMBERS, not the template. Matching "t_half" and "1 - G"
  # only matched string literals in the paste0() that builds this, so
  # inverting the ratio or reporting G instead of 1 - G both SURVIVED as
  # mutations while these passed. The figures are the whole payload of the
  # record, so they are compared against an independent calculation.
  t_half <- exp(unname(f$fit$theta[["early.log_t_half"]]))
  t_min <- min(d$t[d$t > 0])
  dec <- hzr_decompos(t_min, t_half = t_half,
                      nu = unname(f$fit$theta[["early.nu"]]),
                      m = unname(f$fit$theta[["early.m"]]))
  want_mass <- exp(dec$log_surv)
  # The mass is reported from log(1 - G), not from 1 - G (#578). On this
  # fixture it is about 2e-13: formed as 1 - G it kept three digits, and
  # before #578 the fit sat further out, where it rounded to exactly 0 and
  # the record printed "0". So the figure must be positive, and 1 - G formed
  # the old way must agree with it to the digits it has.
  expect_gt(want_mass, 0)
  # As a ratio: both are about 2e-13, below any tolerance, so an absolute
  # comparison would accept 1 - G rounding to 0 (Copilot, #583).
  expect_equal((1 - dec$G) / want_mass, 1, tolerance = 1e-2)
  # ANCHOR the number to its phrase. An unanchored match passed while a
  # mutation reporting G instead of 1 - G survived.
  expect_gt(nchar(format(t_min / t_half, digits = 3)), 2L)
  expect_match(b[[1L]]$detail,
               paste0("a factor of ", format(t_min / t_half, digits = 3), "."),
               fixed = TRUE)
  expect_match(b[[1L]]$detail,
               paste0("1 - G(t_min), is ", format(want_mass, digits = 4), "."),
               fixed = TRUE)
  expect_gt(t_min / t_half, 1)          # the ratio is stated the right way up
})

test_that("a fit that was examined and found nothing reads NULL, not NA", {
  # THE DISTINCTION THE TRI-STATE EXISTS FOR. A two-state field cannot tell
  # "looked and found nothing" from "never looked", and NULL is the common
  # case, so the confusion would be invisible.
  d <- ub_data()
  f <- ub_fit(d, list(early = hzr_phase("hazard", t_half = 3, nu = 1, m = 0),
                      late  = hzr_phase("constant")))
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

test_that("the suite's one former trip was the clamp, and no longer trips (#578)", {
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
  # UNTIL #578 THIS FIT TRIPPED, and the trip was the defect's own work. The
  # hazard phase's cumulative hazard was formed from 1 - G and clamped once
  # G rounded to 1: the hazard came back near 1e290, each event added about
  # +668 to the log-likelihood, and the optimizer converged on that at
  # +60479.26 with 120 events, t_half at 0.133 below a first time of 0.174.
  # With log(1 - G) carried (#578) that point no longer exists. The fit ends
  # with the phase inside its data and an ordinary log-likelihood (-260.495
  # when measured), so there is nothing below the support to report.
  expect_lt(fit$fit$objective, 0)
  expect_true(all(is.finite(predict(fit, type = "hazard"))))
  b <- fit$fit$boundary
  # By mechanism, not position: $boundary can hold several records in no
  # guaranteed order (#261).
  ub <- Filter(function(r) identical(r$mechanism, "unbounded_phase"),
               if (is.list(b)) b else list())
  expect_length(ub, 0L)
  t_half <- exp(unname(fit$fit$theta[["a.log_t_half"]]))
  expect_gt(t_half, min(tt))
  # A known positive for the detector would be a fixture whose hazard phase
  # genuinely fits below its support. None was found: nine fits of late
  # follow-up and left-shifted Weibull data (three seeds each, hazard phase
  # with m = 0 and m = 1, five starts) all put t_half inside the data. The
  # fixture above ("a hazard phase fitted below the observed support") still
  # trips, but from a start at t_min / 1000 that the optimizer does not
  # leave: its log-likelihood is -430.69 there against -410.13 from a start
  # inside the data. Whether the detector has a genuine case to catch is the
  # review of #444, done separately.
})

test_that("every observed time counts, not just `time`", {
  # t_min was min(time) alone. On a left-truncated fit time_lower is the
  # ENTRY time and can lie below min(time), so a t_half between the two was
  # reported as "below the first observed time" while observations existed
  # below it -- a false claim in the record's own words (#444).
  th <- c(a.log_mu = log(0.1), a.log_t_half = log(0.5), a.nu = 1, a.m = 0,
          c.log_mu = log(0.05))
  ph <- list(a = hzr_phase("hazard"), c = hzr_phase("constant"))
  tt <- c(1, 2, 3)
  # With `time` alone, t_half = 0.5 is below min(time) = 1 and trips.
  r1 <- .hzr_boundary_check_impl(th, ph, tt, fitted = TRUE)
  expect_true(is.list(r1$boundary))
  # An entry time of 0.2 means observation began BELOW t_half, so it does not.
  r2 <- .hzr_boundary_check_impl(th, ph, tt, fitted = TRUE,
                                 time_lower = c(0.2, 0.2, 0.2))
  expect_null(r2$boundary)
})

test_that("a duplicated phase name yields one record, not two", {
  # The loop ran over names(phases), and `phases[[nm]]` resolves a duplicated
  # name to the FIRST element -- so c("a", "a") iterated twice and emitted two
  # identical records for one phase. It runs by index now.
  # The TYPES must differ, or the two implementations agree by accident: with
  # both duplicates the same type, looping by name emits two records because
  # it reads the first element twice, and looping by index emits two because
  # there really are two -- the same count for opposite reasons. An earlier
  # version of this test asserted that count and did NOT kill the revert.
  th <- c(a.log_mu = log(0.1), a.log_t_half = log(0.001), a.nu = 1, a.m = 0,
          c.log_mu = log(0.05))
  ph <- list(a = hzr_phase("hazard"), a = hzr_phase("cdf"),
             c = hzr_phase("constant"))
  r <- .hzr_boundary_check_impl(th, ph, c(1, 2, 3), fitted = TRUE)
  # One record: the "hazard" duplicate. Looping by name reads phases[["a"]]
  # -- the FIRST element -- for both, sees "hazard" twice, and emits two.
  expect_length(r$boundary, 1L)
  expect_identical(r$boundary[[1L]]$phase, "a")
  # An UNNAMED phase is skipped rather than producing a "<NA>.log_t_half" key.
  ph2 <- list(a = hzr_phase("hazard"), hzr_phase("hazard"))
  names(ph2) <- c("a", "")
  r2 <- .hzr_boundary_check_impl(th, ph2, c(1, 2, 3), fitted = TRUE)
  expect_length(r2$boundary, 1L)
})

test_that("a weight-0 row does not supply the first observed time", {
  # The likelihood drops a weight-0 row, so its time is not an observed time
  # of this fit. Before, a weight-0 row placed below the fitted t_half moved
  # t_min below it and silently removed the record.
  d <- ub_data()
  ph <- list(early = hzr_phase("hazard", t_half = 0.5, nu = 1, m = 0),
             late  = hzr_phase("constant"))
  th <- c(log(0.1), log(min(d$t) / 1000), 1, 0, log(0.05))
  base <- ub_fit(d, ph, theta = th)
  padded <- ub_fit(list(t = c(d$t, 1e-12), s = c(d$s, 0)), ph, theta = th,
                   weights = c(rep(1, length(d$t)), 0))
  # The premise: the ghost row sits below the fitted t_half.
  expect_lt(1e-12, exp(unname(padded$fit$theta[["early.log_t_half"]])))
  expect_true(is.list(padded$fit$boundary))
  expect_identical(padded$fit$boundary[[1L]]$detail,
                   base$fit$boundary[[1L]]$detail)
})

test_that("a row the designs drop for an NA covariate is not an observed time", {
  # A phase formula is evaluated in `data`, and the multiphase designs drop a
  # row with an NA there before the likelihood sees it, while hazard() keeps
  # it in `time` (the global formula has no covariate to lose it on).
  d <- ub_data()
  n <- length(d$t)
  withr::local_seed(2)
  df <- data.frame(t = c(d$t, 1e-12), s = c(d$s, 0),
                   z = c(stats::rnorm(n), NA))
  ph <- list(early = hzr_phase("hazard", t_half = 0.5, nu = 1, m = 0),
             late  = hzr_phase("constant", formula = ~ z))
  fit <- suppressWarnings(hazard(
    survival::Surv(t, s) ~ 1, data = df, dist = "multiphase",
    phases = ph, fit = TRUE, control = list(n_starts = 1L),
    theta = c(log(0.1), log(min(d$t) / 1000), 1, 0, log(0.05), 0)))
  t_half <- exp(unname(fit$fit$theta[["early.log_t_half"]]))
  expect_lt(1e-12, t_half)  # the premise: the dropped row is below t_half
  expect_true(is.list(fit$fit$boundary))
  expect_match(fit$fit$boundary[[1L]]$detail,
               paste0("below the first observed time (",
                      format(min(d$t), digits = 4), ")"), fixed = TRUE)
})

test_that("a bound the likelihood never reads is not an observed time", {
  # `time_lower` is ignored on a left-censored row, so a value there is not a
  # point the fit evaluates at, and must not become the first observed time.
  d <- ub_data()
  n <- length(d$t)
  tt <- c(d$t, stats::median(d$t))
  st <- c(d$s, -1)
  lo <- c(rep(0, n), 1e-12)
  ph <- list(early = hzr_phase("hazard", t_half = 0.5, nu = 1, m = 0),
             late  = hzr_phase("constant"))
  fit <- suppressWarnings(hazard(
    time = tt, status = st, time_lower = lo, dist = "multiphase",
    phases = ph, fit = TRUE, control = list(n_starts = 1L),
    theta = c(log(0.1), log(min(d$t) / 1000), 1, 0, log(0.05))))
  t_half <- exp(unname(fit$fit$theta[["early.log_t_half"]]))
  expect_lt(1e-12, t_half)  # the premise: the unused bound is below t_half
  expect_true(is.list(fit$fit$boundary))
  expect_match(fit$fit$boundary[[1L]]$detail,
               paste0("below the first observed time (",
                      format(min(d$t), digits = 4), ")"), fixed = TRUE)
})

test_that("a `time` an explicit bound replaces is not an observed time", {
  # A left-censored row is evaluated at its `time_upper`; its `time` is not
  # read, so a tiny value there must not become the first observed time.
  d <- ub_data()
  n <- length(d$t)
  tt <- c(d$t, 1e-12)
  st <- c(d$s, -1)
  up <- c(d$t, stats::median(d$t))
  ph <- list(early = hzr_phase("hazard", t_half = 0.5, nu = 1, m = 0),
             late  = hzr_phase("constant"))
  fit <- suppressWarnings(hazard(
    time = tt, status = st, time_upper = up, dist = "multiphase",
    phases = ph, fit = TRUE, control = list(n_starts = 1L),
    theta = c(log(0.1), log(min(d$t) / 1000), 1, 0, log(0.05))))
  t_half <- exp(unname(fit$fit$theta[["early.log_t_half"]]))
  expect_lt(1e-12, t_half)
  expect_true(is.list(fit$fit$boundary))
  expect_match(fit$fit$boundary[[1L]]$detail,
               paste0("below the first observed time (",
                      format(min(d$t), digits = 4), ")"), fixed = TRUE)
})
