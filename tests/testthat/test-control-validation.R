# hazard() checks `control` against the elements the fitter reads (#376).
# An element it would ignore used to be accepted silently: a typo such as
# `n_startz` left the default in force, and `control$fix` returned a fit
# identical to the unconstrained one, the "fixed" parameter having moved.

cv_data <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc)[1:150, ]
}

cv_weibull <- function(control, fit = FALSE) {
  d <- cv_data()
  hazard(survival::Surv(int_dead, dead) ~ age, data = d, dist = "weibull",
         theta = c(mu = 0.01, nu = 0.5, 0.004), fit = fit,
         control = control)
}

cv_multiphase <- function(control, fit = FALSE) {
  d <- cv_data()
  suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = fit, control = control
  ))
}

test_that("an unknown control element is refused and named", {
  expect_error(cv_weibull(list(n_startz = 99)),
               "'n_startz'.*maxit, reltol, abstol, shape_param_count")
  expect_error(cv_multiphase(list(n_startz = 99, nonsense = "x")),
               "'n_startz', 'nonsense'")
  # The multiphase list names the multiphase elements too.
  expect_error(cv_multiphase(list(n_startz = 99)), "n_starts, conserve")
})

test_that("control$fix is refused, pointing at hzr_phase(fixed = )", {
  for (make in list(cv_weibull, cv_multiphase)) {
    expect_error(make(list(fix = 2L)), "hzr_phase\\(fixed = \\)")
  }
  # It says what happened before, so a user of an old fit knows it was
  # unconstrained, and it is not reported as a mere unknown element.
  expect_error(cv_weibull(list(fix = 2L)), "never read")
  expect_error(cv_weibull(list(fix = 2L, maxit = 50)), "control\\$fix")
})

test_that("control$quasi is refused as never read", {
  expect_error(cv_multiphase(list(quasi = TRUE)), "control\\$quasi.*never read")
})

test_that("a multiphase-only element on a single distribution says so", {
  expect_error(cv_weibull(list(n_starts = 3)),
               "'n_starts'.*only for dist = \"multiphase\"")
})

test_that("control must be a named list", {
  expect_error(cv_weibull(list(200)), "named")
})

test_that("the elements the fitter reads are accepted", {
  # Known positive: legitimate lists still construct and fit.
  expect_s3_class(cv_weibull(list(maxit = 200, reltol = 1e-8, abstol = 1e-6)),
                  "hazard")
  expect_s3_class(cv_weibull(list(shape_param_count = 2L)), "hazard")
  expect_s3_class(
    cv_multiphase(list(n_starts = 1L, conserve = FALSE, phase_share_tol = 0,
                       start_seed = 3L, maxit = 200)),
    "hazard"
  )
  expect_s3_class(cv_weibull(list()), "hazard")
  fit <- cv_weibull(list(maxit = 500), fit = TRUE)
  expect_true(is.finite(fit$fit$objective))
})
