# The G3 second derivatives are central differences, and a fixed additive
# step is too small for a large shape and too large for a small one (#332).
# At gamma = 220 the step moved gamma by 5e-7 of itself and rounding took
# over; at eta = 0.009 it moved eta by 1.3%. The Hessian entries were off by
# up to 6% and the standard errors of gamma and eta by 24%.
#
# The reference is the Jacobian of the analytic gradient. It divides by the
# step once, not twice, and agrees with numDeriv's Hessian of the
# log-likelihood at d = 1e-2 and with a nested numerical gradient to about
# 1e-5 at these points. numDeriv's default d = 0.1 does not: at gamma = 200
# it is itself off by a factor of 3.

g3_steps_data <- function() {
  withr::local_seed(20260916)
  th0 <- c(log(0.3), log(10), 4, 1.5, 0.5)
  n <- 1500
  grid <- seq(0.001, 60, length.out = 20000)
  cumhaz <- exp(th0[1]) *
    hzr_decompos_g3(grid, exp(th0[2]), th0[3], th0[4], th0[5])$G3
  t_event <- stats::approx(cumhaz, grid, xout = -log(stats::runif(n)),
                           rule = 2, ties = "ordered")$y
  cens <- stats::runif(n, 2, 25)
  list(time = pmin(t_event, cens), status = as.integer(t_event <= cens))
}

g3_steps_hessians <- function(d, th) {
  names(th) <- paste0("late.", c("log_mu", "log_tau", "gamma", "alpha", "eta"))
  phases <- list(late = hzr_phase("g3"))
  counts <- c(late = 0L)
  x_list <- list(late = NULL)
  analytic <- .hzr_hessian_multiphase(th, d$time, d$status, phases = phases,
                                      covariate_counts = counts,
                                      x_list = x_list)
  grad <- function(p) {
    .hzr_gradient_multiphase(p, d$time, d$status, phases = phases,
                             covariate_counts = counts, x_list = x_list,
                             sanitize = FALSE)
  }
  reference <- numDeriv::jacobian(grad, th,
                                  method.args = list(d = 1e-3, r = 6))
  list(analytic = analytic, reference = -(reference + t(reference)) / 2)
}

test_that("g3 standard errors match the reference at an extreme shape (#332)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  d <- g3_steps_data()
  # #327's regression point: gamma ~ 220, eta ~ 0.009.
  h <- g3_steps_hessians(d, c(-0.91648, 2.6238, 220.08, 1.2174, 0.0090875))
  se <- sqrt(diag(solve(h$analytic)))
  se_ref <- sqrt(diag(solve(h$reference)))
  expect_true(all(is.finite(se_ref)))
  # gamma and eta are the weakly identified pair the fixed step got wrong.
  expect_lt(max(abs(se / se_ref - 1)), 0.01)
  # A fixed eta step moved eta by 1.3% here and put this entry 1.4% off; it
  # is 2e-4 off with the relative step, below the reference's own 1.5e-3.
  rel <- abs(h$analytic - h$reference) / abs(h$reference)
  expect_lt(rel["late.log_tau", "late.eta"], 2e-3)
})

test_that("g3 Hessian entries match the reference across shapes (#332)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  d <- g3_steps_data()
  rel <- function(h) {
    max(abs(h$analytic - h$reference) / pmax(abs(h$reference), 1e-6))
  }
  # Large gamma, small eta: rounding of a step too small for gamma.
  expect_lt(rel(g3_steps_hessians(d, c(log(0.25), log(9), 100, 1.5, 0.02))),
            1e-4)
  expect_lt(rel(g3_steps_hessians(d, c(log(0.25), log(9), 200, 1, 0.01))),
            1e-4)
  # Small alpha: a step of 1.2e-4 is 0.6% of alpha = 0.02.
  expect_lt(rel(g3_steps_hessians(d, c(log(0.3), log(10), 4, 0.02, 0.5))),
            1e-3)
})

test_that("a fit at an extreme g3 shape reports the reference standard errors (#332)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  d <- g3_steps_data()
  fit <- suppressWarnings(hazard(
    time = d$time, status = d$status, dist = "multiphase",
    phases = list(late = hzr_phase("g3", tau = exp(2.6238), gamma = 220.08,
                                   alpha = 1.2174, eta = 0.0090875)),
    fit = TRUE, control = list(n_starts = 1L, conserve = FALSE)
  ))
  expect_true(isTRUE(fit$fit$converged))
  expect_gt(fit$fit$theta[["late.gamma"]], 100)
  h <- g3_steps_hessians(d, unname(fit$fit$theta))
  se <- sqrt(diag(fit$fit$vcov))
  se_ref <- sqrt(diag(solve(h$reference)))
  expect_true(all(is.finite(se_ref)))
  # The fixed step reported SE(gamma) = 52 and SE(eta) = 0.0028 here, where
  # the reference gives about 750 and 0.035: an order of magnitude too
  # confident. The reference itself moves about 5% with its step size.
  expect_lt(max(abs(se / se_ref - 1)), 0.2)
})

test_that("the g3 score is right for a small positive alpha (#332)", {
  # A forward step of 1e-5 at alpha <= 1e-5 moved alpha by 100% or more of
  # itself: dPhi/dalpha and the log-likelihood gradient in alpha were 50-70%
  # off, and the optimizer used that gradient. The reference is a Richardson
  # central difference on the closed form and on the log-likelihood, with a
  # step of 1% of alpha, so it never crosses alpha = 0; at 0.1% it agrees to
  # about 1e-8.
  withr::local_seed(20260917)
  tau <- 1e5
  gamma <- 1.5
  eta <- 0.8
  grid <- 10^seq(-3, 2, length.out = 20000)
  cumhaz <- 0.3 * hzr_decompos_g3(grid, tau, gamma, 1e-5, eta)$G3
  t_event <- stats::approx(cumhaz, grid, xout = -log(stats::runif(800)),
                           rule = 2, ties = "ordered")$y
  cens <- stats::runif(800, 2, 25)
  time <- pmin(t_event, cens)
  status <- as.integer(t_event <= cens)
  phases <- list(late = hzr_phase("g3"))
  counts <- c(late = 0L)
  x_list <- list(late = NULL)
  loglik <- function(p) {
    sum(log(.hzr_multiphase_hazard(time[status == 1], p, phases, counts,
                                   x_list))) -
      sum(.hzr_multiphase_cumhaz(time, p, phases, counts, x_list))
  }
  richardson <- function(f, x) {
    h <- 0.01 * x
    d1 <- (f(x + h) - f(x - h)) / (2 * h)
    d2 <- (f(x + h / 2) - f(x - h / 2)) / h
    (4 * d2 - d1) / 3
  }
  for (alpha in c(5e-6, 1e-5)) {
    phi <- .hzr_g3_phase_derivatives(time[1:5], tau, gamma, alpha, eta)
    ref <- vapply(time[1:5], function(t1) {
      richardson(function(a) hzr_decompos_g3(t1, tau, gamma, a, eta)$G3, alpha)
    }, numeric(1))
    expect_lt(max(abs(phi$dPhi_dalpha / ref - 1)), 1e-6)

    theta <- c(log(0.3), log(tau), gamma, alpha, eta)
    score <- .hzr_gradient_multiphase(theta, time, status, phases = phases,
                                      covariate_counts = counts,
                                      x_list = x_list, sanitize = FALSE)
    ref_score <- richardson(function(a) {
      loglik(c(log(0.3), log(tau), gamma, a, eta))
    }, alpha)
    expect_lt(abs(score[4] / ref_score - 1), 1e-6)
  }
})
