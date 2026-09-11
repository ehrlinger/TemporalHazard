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
  fit_with <- function(code, rel = 1e-2, conv = 0L) {
    testthat::local_mocked_bindings(.hzr_optim_exponential = function(...) {
      r <- real(...)
      r$polish_code <- code
      r$rel_gradient <- rel
      r$convergence <- conv
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
  # A fit that did not converge is not "converged over a failing point": no
  # gradient warning even with a hard-failure code attached.
  expect_no_warning(fit_with(4L, conv = 1L), message = "relative-gradient test")
  # nlm code 3, where SAS/C prints a caution and retries: recorded, not warned.
  expect_no_warning(f3 <- fit_with(3L), message = "relative-gradient test")
  expect_identical(f3$fit$polish_code, 3L)
  expect_equal(f3$fit$rel_gradient, 1e-2)
  expect_output(print(f3), "relative 0.01 .*not met, nlm code 3")
  expect_output(print(summary(f3)), "relative 0.01 .*not met, nlm code 3")
  # A passing fit says so, and shows the nlm() code when there is one.
  f_ok <- fit_with(NA_integer_, rel = 1e-9)
  expect_output(print(f_ok), "relative 1e-09 .*; met\\)")
  expect_output(print(fit_with(1L, rel = 1e-9)), "; met, nlm code 1\\)")
  # Code 4 whose statistic nonetheless meets the test -- possible when
  # nlm() stops on a partial score, as under CoE -- does not warn.
  expect_no_warning(fit_with(4L, rel = 1e-9), message = "relative-gradient test")
})

test_that("a polished fit keeps the caller's parameter names and says it went on", {
  # nlm() returns an unnamed estimate; optim() keeps names, and the
  # single-distribution fits return par as is, so a lost name surfaced as an
  # unnamed coef() only on fits that happened to be polished.
  fit <- .hzr_optim_generic(
    logl_fn = rosen_logl, gradient_fn = rosen_score, time = 1, status = 1,
    theta_start = c(a = -1.2, b = 1), hessian_fn = rosen_hessian
  )
  expect_equal(fit$polish_code, 1L)
  expect_named(fit$par, c("a", "b"))
  expect_match(fit$message, "continued with nlm\\(\\)")
})

test_that("rel_gradient is NA, never a pass, where the gradient cannot be trusted", {
  # A likelihood that is -Inf everywhere: the objective is clamped to 1e10 and
  # the wrapped gradient returns zeros, which read as a relative gradient of 0.
  flat <- .hzr_optim_generic(
    logl_fn = function(theta, ...) -Inf, gradient_fn = rosen_score,
    time = 1, status = 1, theta_start = start,
    hessian_fn = function(theta) diag(2)
  )
  expect_true(is.na(flat$rel_gradient))
  # A score with a NaN component: the wrapper zeroes it, and BFGS stops far
  # from the optimum with an apparently tiny gradient.
  nan_score <- suppressWarnings(.hzr_optim_generic(
    logl_fn = rosen_logl, gradient_fn = function(theta, ...) {
      c(NaN, rosen_score(theta)[2])
    },
    time = 1, status = 1, theta_start = start, hessian_fn = rosen_hessian
  ))
  expect_true(is.na(nan_score$rel_gradient))
  # Converged with no usable gradient: said, not left blank.
  expect_match(.hzr_format_gradient_test(nan_score$rel_gradient,
                                         nan_score$polish_code),
               "not evaluated")
  # Not converged, or no record at all: no line.
  expect_null(.hzr_format_gradient_test(NA_real_, NA_integer_,
                                        converged = FALSE))
  expect_null(.hzr_format_gradient_test(NULL, NULL))
})

test_that("with an inexact score the recorded test uses the objective's own gradient", {
  # Conservation of Events hands .hzr_optim_generic() the partial score at the
  # conserved theta, which omits how the conserved scale moves: a gradient
  # that is not the objective's. Model that with a score that is wrong in a
  # known way: a hundredth of the truth in its first component, which is the
  # one that sets SAS's maximum near this optimum, and right in its second.
  # The continuation still follows that score; the recorded test must not.
  biased <- function(theta, ...) rosen_score(theta) * c(0.01, 1)
  rel_biased <- function(theta) {
    max(abs(biased(theta)) * pmax(abs(theta), 1)) /
      max(abs(rosen_logl(theta)), 1)
  }
  fit <- .hzr_optim_generic(
    logl_fn = rosen_logl, gradient_fn = biased, time = 1, status = 1,
    theta_start = start, hessian_fn = rosen_hessian, gradient_exact = FALSE
  )
  expect_equal(fit$convergence, 0L)
  # Premise: at the returned point the two statistics differ, so the next
  # assertion can tell them apart. If this fails the fixture no longer
  # discriminates and must change.
  expect_gt(abs(rel_biased(fit$par) / rel_gradient(fit$par) - 1), 0.1)
  # The recorded statistic is the objective's own, to finite-difference
  # accuracy, not the biased score's.
  expect_lt(abs(fit$rel_gradient / rel_gradient(fit$par) - 1), 1e-3)
})

test_that("multiphase passes gradient_exact = FALSE exactly when CoE is applied", {
  # The flag is set at .hzr_optim_multiphase()'s call; spy on it there, so
  # reverting that one line fails this test.
  set.seed(3)
  n <- 300
  tt <- c(stats::rexp(n / 2, 3), stats::rexp(n / 2, 0.15))
  cens <- stats::runif(n, 1, 15)
  d <- data.frame(time = pmin(tt, cens), status = as.integer(tt <= cens))
  real <- .hzr_optim_generic
  seen <- logical(0)
  testthat::local_mocked_bindings(
    .hzr_optim_generic = function(..., gradient_exact = TRUE) {
      seen <<- c(seen, gradient_exact)
      real(..., gradient_exact = gradient_exact)
    }
  )
  fit_with <- function(conserve) {
    seen <<- logical(0)
    f <- suppressWarnings(hazard(
      survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
      phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0),
                    constant = hzr_phase("constant")),
      fit = TRUE, control = list(n_starts = 1, conserve = conserve)
    ))
    list(applied = isTRUE(f$spec$control$conserve_applied), seen = seen)
  }
  on <- fit_with(TRUE)
  off <- fit_with(FALSE)
  expect_true(on$applied)
  expect_false(off$applied)
  expect_true(length(on$seen) > 0 && !any(on$seen))
  expect_true(length(off$seen) > 0 && all(off$seen))
})

test_that("a finite difference that lands on the clamp is NA, not a fabricated gradient", {
  # The maximum sits on the edge of the region where the likelihood is
  # defined: beyond theta[1] = 1 it is -Inf, so the objective is clamped to
  # 1e10 there, and a central difference taken at the edge has one side on
  # the clamp. Without the guard that side makes a "gradient" of 1e10 / (2h).
  edge_logl <- function(theta, ...) {
    if (theta[1] > 1) -Inf else -(1e4 + (theta[1] - 2)^2 + theta[2]^2)
  }
  edge_score <- function(theta, ...) c(-2 * (theta[1] - 2), -2 * theta[2])
  fit <- suppressWarnings(.hzr_optim_generic(
    logl_fn = edge_logl, gradient_fn = edge_score, time = 1, status = 1,
    theta_start = c(0, 0.5), hessian_fn = function(theta) diag(2),
    gradient_exact = FALSE
  ))
  expect_equal(fit$convergence, 0L)
  # Premise: the fit stopped within one difference step of the edge.
  expect_true(fit$par[1] <= 1 &&
                fit$par[1] > 1 - .Machine$double.eps^(1 / 3))
  expect_true(is.na(fit$rel_gradient))
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
