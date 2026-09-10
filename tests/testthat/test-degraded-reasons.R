# Reasons carried up to the record (#242): the Hessian path's reason for losing
# standard errors, and the weak-direction check's reason for not looking.
# Each cause string is the spec's, verbatim.

# A one-parameter exponential model, log-rate scale, small enough to fit in
# milliseconds and with a Hessian known in closed form: exp(theta) * sum(time).
exp_logl <- function(theta, time, status, time_lower = NULL, time_upper = NULL,
                     x = NULL, weights = NULL, return_gradient = FALSE) {
  sum(status * theta[1] - exp(theta[1]) * time)
}
exp_grad <- function(theta, time, status, ...) {
  sum(status) - exp(theta[1]) * sum(time)
}
fit_exp <- function(hess) {
  set.seed(3)
  time <- stats::rexp(40, 0.5)
  suppressWarnings(.hzr_optim_generic(
    logl_fn = exp_logl, gradient_fn = exp_grad,
    time = time, status = rep(1, 40), theta_start = 0,
    hessian_fn = function(p) hess(p, time)
  ))
}

test_that(".hzr_safe_solve() says why it returned no covariance", {
  expect_identical(suppressWarnings(.hzr_safe_solve(matrix(NaN, 2, 2)))$reason,
                   "Hessian has non-finite entries")
  expect_identical(suppressWarnings(.hzr_safe_solve(matrix(0, 2, 2)))$reason,
                   "Hessian not invertible")
  expect_identical(.hzr_safe_solve(diag(2))$reason, NA_character_)
})

test_that("an analytic Hessian that inverts leaves no reason", {
  res <- fit_exp(function(p, time) matrix(exp(p) * sum(time), 1, 1))
  expect_true(is.matrix(res$vcov))                 # guard: SEs exist
  expect_identical(res$se_unavailable_reason, NA_character_)
})

test_that("each way the Hessian path loses standard errors names its cause", {
  res <- fit_exp(function(p, time) matrix(NaN, 1, 1))
  expect_false(is.matrix(res$vcov))
  expect_identical(res$se_unavailable_reason, "Hessian has non-finite entries")

  res <- fit_exp(function(p, time) matrix(0, 1, 1))
  expect_false(is.matrix(res$vcov))
  expect_identical(res$se_unavailable_reason, "Hessian not invertible")

  local_mocked_bindings(.hzr_numderiv_available = function() FALSE)
  res <- fit_exp(function(p, time) NULL)
  expect_false(is.matrix(res$vcov))
  expect_identical(res$se_unavailable_reason,
                   "numDeriv not installed and no analytic Hessian")
})

test_that("a numDeriv failure is named, not folded into 'not installed'", {
  skip_if_not_installed("numDeriv")
  local_mocked_bindings(hessian = function(...) stop("boom"),
                        .package = "numDeriv")
  res <- fit_exp(function(p, time) NULL)
  expect_false(is.matrix(res$vcov))
  expect_identical(res$se_unavailable_reason, "numDeriv::hessian() failed")
})

test_that("every 'could not look' exit of the weak-direction check names why", {
  why <- function(...) .hzr_weak_direction_impl(...)$reason
  expect_identical(why(NULL, NA_real_), "standard errors unavailable")
  expect_identical(why(diag(2), NA_real_), "Hessian condition number unavailable")
  expect_identical(why(NA, 1e-10), "standard errors unavailable")
  # Finite, positive diagonal, but a non-finite entry among estimated params.
  expect_identical(why(matrix(c(1, Inf, Inf, 1), 2), 1e-10),
                   "covariance has non-finite entries")
  # Finite covariance whose correlation overflows: 1e300 / (1e-150)^2.
  expect_identical(why(matrix(c(1e-300, 1e300, 1e300, 1e-300), 2), 1e-10),
                   "covariance has non-finite entries")
})

test_that("the split leaves .hzr_weak_direction() returning what it did", {
  expect_identical(.hzr_weak_direction(NULL, NA_real_), NA)
  expect_null(.hzr_weak_direction(diag(2), 0.5))
  res <- .hzr_weak_direction_impl(diag(2), 0.5)
  expect_true("weak" %in% names(res))     # a NULL result keeps its slot
  expect_null(res$weak)
  expect_identical(res$reason, NA_character_)
})
