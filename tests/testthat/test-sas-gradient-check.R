# test-sas-gradient-check.R -- SAS/C's acceptance test in .hzr_optim_generic()
#
# SAS/C HAZARD accepts an optimum only when the relative gradient
#   max_i |g_i| * max(|x_i|, 1) / max(|f|, 1)
# is at most gradtl = eps^(1/3) (src/optim/umstop.c). optim()'s BFGS stops on
# the relative change in the objective instead, which reports convergence
# short of the optimum on a flat valley. .hzr_optim_generic() now polishes
# such a stop with stats::nlm() at SAS's tolerances.
#
# The objective is a Rosenbrock valley shifted by 1e4. The offset makes
# optim's relative criterion, reltol * |f|, about 0.1, so BFGS stops well
# down the valley and still reports convergence 0 -- the shape of the defect
# on real fits, in two parameters with a known optimum at (1, 1).

gradtl <- .Machine$double.eps^(1 / 3)

rosen_logl <- function(theta, ...) {
  -(1e4 + 100 * (theta[2] - theta[1]^2)^2 + (1 - theta[1])^2)
}
rosen_score <- function(theta, ...) {
  c(400 * theta[1] * (theta[2] - theta[1]^2) + 2 * (1 - theta[1]),
    -200 * (theta[2] - theta[1]^2))
}
# Hessian of the negative log-likelihood, so no numDeriv is needed.
rosen_hessian <- function(theta) {
  matrix(c(1200 * theta[1]^2 - 400 * theta[2] + 2, -400 * theta[1],
           -400 * theta[1], 200), 2, 2)
}
rel_gradient <- function(theta) {
  max(abs(rosen_score(theta)) * pmax(abs(theta), 1)) /
    max(abs(rosen_logl(theta)), 1)
}
start <- c(-1.2, 1)
plain_bfgs <- function() {
  stats::optim(start, function(th) -rosen_logl(th),
               function(th) -rosen_score(th), method = "BFGS",
               control = list(maxit = 1000, reltol = 1e-5))
}
# Log-likelihood below the optimum, whose value is exactly -1e4.
ll_gap <- function(theta) -1e4 - rosen_logl(theta)

test_that("plain BFGS at the default reltol stops short of SAS's test here", {
  # The premise of the next test. If this ever passes gradtl, the fixture no
  # longer exercises the polish and must be made harder.
  bfgs <- plain_bfgs()
  expect_equal(bfgs$convergence, 0L)
  expect_gt(rel_gradient(bfgs$par), 10 * gradtl)
})

test_that("a converged BFGS stop that fails SAS's test is polished to pass it", {
  fit <- .hzr_optim_generic(
    logl_fn = rosen_logl, gradient_fn = rosen_score,
    time = 1, status = 1, theta_start = start, hessian_fn = rosen_hessian
  )
  expect_equal(fit$convergence, 0L)
  expect_lte(fit$rel_gradient, gradtl)
  expect_equal(fit$polish_code, 1L)
  # SAS's test is relative to |f|: at |f| = 1e4 it accepts |g_i * x_i| up to
  # about 0.06, so a passing point is within SAS's tolerance of the optimum,
  # not at it. What it buys is the log-likelihood, so assert on that: the
  # polished point is within 1e-4 of the maximum and at least 100 times
  # closer than plain BFGS stopped.
  expect_lt(ll_gap(fit$par), 1e-4)
  expect_lt(ll_gap(fit$par), ll_gap(plain_bfgs()$par) / 100)
  expect_true(all(abs(fit$par - c(1, 1)) < 1e-2),
              label = paste("par =", paste(signif(fit$par, 8), collapse = ", ")))
  # The recorded statistic is the one at the returned point.
  expect_equal(fit$rel_gradient, rel_gradient(fit$par), tolerance = 1e-8)
})

test_that("hazard() warns only on the polish's hard failures, and records every result", {
  # The optimizer's return is the real one with the two fields overridden, so
  # everything else hazard() reads from it is genuine.
  set.seed(7)
  df <- data.frame(time = stats::rexp(80, 0.4), status = rep(c(1, 1, 0), length.out = 80),
                   z = stats::rnorm(80))
  real <- .hzr_optim_exponential
  fit_with <- function(code, rel = 1e-2) {
    testthat::local_mocked_bindings(.hzr_optim_exponential = function(...) {
      r <- real(...)
      r$polish_code <- code
      r$rel_gradient <- rel
      r
    })
    hazard(survival::Surv(time, status) ~ z, data = df, dist = "exponential",
           theta = c(log_rate = 0, z = 0), fit = TRUE)
  }
  # nlm code 5, SAS/C's "unbounded": warns.
  expect_warning(f5 <- fit_with(5L), "kept rising")
  expect_identical(f5$fit$polish_code, 5L)
  # nlm code 4, SAS/C's "reached no convergence": warns and says how to go on.
  expect_warning(fit_with(4L), "iteration limit.*control\\$maxit")
  # nlm code 3, where SAS/C prints a caution and retries: recorded, not warned.
  expect_no_warning(f3 <- fit_with(3L), message = "relative gradient")
  expect_identical(f3$fit$polish_code, 3L)
  expect_equal(f3$fit$rel_gradient, 1e-2)
  expect_output(print(f3), "relative 0.01 .*not met, nlm code 3")
  expect_output(print(summary(f3)), "relative 0.01 .*not met, nlm code 3")
  # A passing fit says so.
  f_ok <- fit_with(NA_integer_, rel = 1e-9)
  expect_output(print(f_ok), "relative 1e-09 .*; met\\)")
})

test_that("the bounded (L-BFGS-B) path is not polished", {
  fit <- .hzr_optim_generic(
    logl_fn = rosen_logl, gradient_fn = rosen_score,
    time = 1, status = 1, theta_start = start, hessian_fn = rosen_hessian,
    use_bounds = TRUE, lower_bounds = c(-Inf, -Inf)
  )
  expect_true(is.na(fit$rel_gradient))
  expect_true(is.na(fit$polish_code))
})
