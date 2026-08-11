# Standard errors going missing must say why.
#
# Three layers combine into a silent failure: the analytic Hessian declines
# for left/interval-censored rows (by design), the optimizer falls back to
# numDeriv, and numDeriv is a Suggests so it may simply be absent. The result
# was rcond = NA, pd = NA and vcov() returning a bare logical, with nothing
# naming the cause -- the user sees diag(vcov(fit)) complain about 'nrow'.

test_that("a failing numDeriv fallback warns instead of returning NA mutely", {
  # The analytic Hessian declining is routine; numDeriv is the documented
  # fallback. When that fallback also fails there is no third option, and the
  # fit used to return rcond = NA / pd = NA / vcov = bare logical in silence.
  # Two warnings are expected: the numDeriv failure, then the no-Hessian
  # consequence.
  skip_if_not_installed("numDeriv")

  logl <- function(theta, time, status, ...) -sum((theta - 1)^2)
  grad <- function(theta, time, status, ...) -2 * (theta - 1)

  local_mocked_bindings(
    hessian = function(...) stop("synthetic numDeriv failure"),
    .package = "numDeriv"
  )

  expect_warning(
    res <- TemporalHazard:::.hzr_optim_generic(
      logl_fn     = logl,
      gradient_fn = grad,
      time        = c(1, 2, 3),
      status      = c(1, 0, 1),
      theta_start = c(0.5, 0.5),
      hessian_fn  = function(...) NULL
    ),
    "numDeriv::hessian\\(\\) failed"
  )
})

test_that("the no-Hessian branch reports NA diagnostics loudly", {
  skip_if_not_installed("numDeriv")

  logl <- function(theta, time, status, ...) -sum((theta - 1)^2)
  grad <- function(theta, time, status, ...) -2 * (theta - 1)

  local_mocked_bindings(
    hessian = function(...) stop("synthetic numDeriv failure"),
    .package = "numDeriv"
  )

  res <- suppressWarnings(TemporalHazard:::.hzr_optim_generic(
    logl_fn     = logl,
    gradient_fn = grad,
    time        = c(1, 2, 3),
    status      = c(1, 0, 1),
    theta_start = c(0.5, 0.5),
    hessian_fn  = function(...) NULL
  ))

  # Behaviour is unchanged -- still NA diagnostics -- only now it is announced.
  expect_true(is.na(res$rcond))
  expect_true(is.na(res$pd))
})

test_that("numDeriv present means a Hessian is computed even for interval-censored fits", {
  # The distinction that matters, and the one the server hit:
  #   rcond = NA        -> no Hessian was produced AT ALL (numDeriv missing)
  #   rcond = <number>  -> a Hessian was produced; whether it inverts is a
  #                        separate question about conditioning
  # Only the first is the numDeriv gap. Asserting vcov() is a matrix would
  # conflate the two, since a well-formed but singular Hessian still yields
  # vcov = NA by design.
  skip_if_not_installed("numDeriv")

  e <- new.env()
  utils::data("cabgkul", package = "TemporalHazard", envir = e)
  d <- e$cabgkul

  # Mark a few rows interval-censored so the analytic multiphase Hessian
  # declines (hessian-multiphase.R coverage contract) and numDeriv is the
  # only remaining route.
  status <- as.integer(d$dead)
  idx <- which(status == 0L)[1:5]
  status[idx] <- 2L
  lower <- d$int_dead
  lower[idx] <- pmax(d$int_dead[idx] - 0.1, 0)

  fit <- suppressWarnings(hazard(
    time = d$int_dead, status = status,
    time_lower = lower, time_upper = d$int_dead,
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.19, nu = 1.4, m = 1),
      late  = hzr_phase("g3", tau = 1, gamma = 1, alpha = 1, eta = 1.7,
                        fixed = c("tau", "gamma", "alpha"))
    ),
    fit = TRUE, control = list(n_starts = 1, maxit = 300)))

  expect_true(fit$fit$converged)
  expect_false(is.na(fit$fit$rcond))
})
