# The record on real fits (#242). Every test asserts a guard first -- that the
# fit really is in the state the record is about -- and that no entry fell
# back to "cause not recorded", which exists for paths the reasons miss and
# which no path this package knows about should reach.

weib_data <- function(n = 80) {
  set.seed(7)
  data.frame(t = stats::rweibull(n, 1.4, 3), d = rep(1L, n))
}

# Copied from test-coe-applied.R's coe_fit(): explicit starts, n_starts = 1.
mp_fit <- function(status_pattern, conserve = TRUE) {
  set.seed(42)
  n <- 60
  time <- stats::rexp(n, 0.2) + 0.05
  status <- rep(status_pattern, length.out = n)
  args <- list(time = time, status = status)
  if (any(status == 2)) {
    args$time_lower <- time
    args$time_upper <- ifelse(status == 2, time * 1.5, time)
  }
  suppressWarnings(do.call(hazard, c(args, list(
    phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                  const = hzr_phase("constant")),
    dist = "multiphase", fit = TRUE,
    control = list(conserve = conserve, n_starts = 1)
  ))))
}

test_that("a clean fit records nothing", {
  fit <- hazard(survival::Surv(t, d) ~ 1, data = weib_data(), dist = "weibull",
                fit = TRUE, theta = c(mu = 1, nu = 1))
  expect_true(is.matrix(fit$fit$vcov))             # guard: SEs exist
  expect_null(fit$fit$weak)                        # guard: examined, clean
  expect_identical(fit$degraded, character(0))
  expect_length(fit$degraded_causes, 0L)
})

test_that("an unfitted object records that it was not fitted, and why", {
  dat <- weib_data()
  fit0 <- hazard(time = dat$t, status = dat$d, dist = "weibull", fit = FALSE)
  expect_identical(fit0$degraded,
                   c("fitting", "standard_errors", "weak_direction_check"))
  expect_identical(fit0$degraded_causes[["fitting"]],
                   "not requested (fit = FALSE)")
  expect_false("cause not recorded" %in% fit0$degraded_causes)
})

test_that("fit = TRUE without theta says the optimizer did not run", {
  # hazard() runs the single-distribution optimizer only when theta is given,
  # and has never said so. If this call errors or fits, the premise has
  # changed: stop and report, do not adjust the test.
  dat <- weib_data()
  fit1 <- hazard(time = dat$t, status = dat$d, dist = "weibull", fit = TRUE)
  expect_null(fit1$fit$counts)                     # guard: no optimizer ran
  expect_identical(fit1$degraded_causes[["fitting"]],
                   "no starting values (theta = NULL); the optimizer did not run")
})

test_that("CoE switched off by the user is listed as not requested", {
  skip_on_cran()
  fit <- mp_fit(c(1, 0), conserve = FALSE)
  expect_false(fit$spec$control$conserve_applied)  # guard
  expect_identical(fit$degraded_causes[["conservation_of_events"]],
                   "not requested (conserve = FALSE)")
  expect_false("cause not recorded" %in% fit$degraded_causes)
})

test_that("interval rows disable CoE, and the record says why", {
  skip_on_cran()
  fit <- mp_fit(c(1, 0, 2))
  expect_false(is.null(fit$fit$theta))             # guard: the fit ran
  expect_identical(fit$degraded_causes[["conservation_of_events"]],
                   "status outside {0, 1} (left or interval censoring)")
  expect_false("cause not recorded" %in% fit$degraded_causes)
})

test_that("no analytic Hessian and no numDeriv: SEs and the weak check are listed", {
  skip_on_cran()
  # Interval rows make the analytic multiphase Hessian decline by design, so
  # the fit needs numDeriv; mocked absent, no Hessian exists at all.
  local_mocked_bindings(.hzr_numderiv_available = function() FALSE)
  fit <- mp_fit(c(1, 0, 2))
  expect_false(is.matrix(fit$fit$vcov))            # guard: SEs really were lost
  expect_identical(
    fit$degraded_causes[c("standard_errors", "weak_direction_check")],
    c(standard_errors = "numDeriv not installed and no analytic Hessian",
      weak_direction_check = "standard errors unavailable"))
  expect_false("cause not recorded" %in% fit$degraded_causes)
})

test_that("summary() carries the record", {
  dat <- weib_data()
  fit0 <- hazard(time = dat$t, status = dat$d, dist = "weibull", fit = FALSE)
  s <- summary(fit0)
  expect_identical(s$degraded, fit0$degraded)
  expect_identical(s$degraded_causes, fit0$degraded_causes)
})

test_that("a failed full-information recompute under CoE is recorded", {
  skip_on_cran()
  # With the analytic Hessian declining and numDeriv absent, both the fit's
  # own vcov and the CoE recompute have no Hessian to work from.
  local_mocked_bindings(
    .hzr_hessian_multiphase = function(...) NULL,
    .hzr_numderiv_available = function() FALSE
  )
  fit <- mp_fit(c(1, 0))
  expect_true(fit$spec$control$conserve_applied)   # guard: the recompute ran
  expect_identical(fit$degraded_causes[["conserved_phase_variance"]],
                   "numDeriv not installed")
  expect_false("cause not recorded" %in% fit$degraded_causes)
})

outhaz_fixture <- function() {
  system.file("extdata", "outhaz-fixture.rds", package = "TemporalHazard")
}

test_that("an imported SAS fit records that R did not fit or examine it", {
  obj <- .hzr_outhaz_to_spec(hzr_read_outhaz(outhaz_fixture()), need_vcov = TRUE)
  expect_true(is.matrix(obj$fit$vcov))             # guard: covariance was read
  expect_identical(
    obj$degraded_causes,
    c(fitting = "imported from SAS output; not fitted in R",
      weak_direction_check = "imported from SAS output; no R Hessian"))
})

test_that("an import read without its covariance also lists standard errors", {
  obj <- .hzr_outhaz_to_spec(hzr_read_outhaz(outhaz_fixture()), need_vcov = FALSE)
  expect_false(is.matrix(obj$fit$vcov))            # guard
  expect_identical(obj$degraded,
                   c("fitting", "standard_errors", "weak_direction_check"))
  expect_identical(obj$degraded_causes[["standard_errors"]],
                   "covariance not imported")
})
