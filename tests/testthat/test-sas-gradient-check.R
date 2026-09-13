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

test_that("the acceptance test asks for the unsanitised score", {
  # A gradient that, like .hzr_gradient_multiphase(), zeroes a component it
  # cannot evaluate unless asked not to. The optimizer needs the zero; the
  # test must see the NaN and report NA, never a small gradient.
  guarded <- function(theta, ..., sanitize = TRUE) {
    g <- c(NaN, rosen_score(theta)[2])
    if (sanitize) g[!is.finite(g)] <- 0
    g
  }
  fit <- suppressWarnings(.hzr_optim_generic(
    logl_fn = rosen_logl, gradient_fn = guarded, time = 1, status = 1,
    theta_start = start, hessian_fn = rosen_hessian
  ))
  expect_equal(fit$convergence, 0L)
  expect_true(is.na(fit$rel_gradient))
})

test_that(".hzr_gradient_multiphase returns NA, not zeros, when asked for the raw score", {
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
                 const = hzr_phase("constant"))
  set.seed(1)
  tt <- stats::rexp(50) + 0.01
  # m < 0 with nu < 0 is undefined: the gradient cannot be evaluated there.
  args <- list(c(log(0.5), log(1), -1, -0.5, log(0.1)), tt, rep(1, 50),
               phases = phases, covariate_counts = c(early = 0L, const = 0L),
               x_list = list(early = NULL, const = NULL))
  expect_true(all(do.call(.hzr_gradient_multiphase, args) == 0))
  expect_true(all(is.na(do.call(.hzr_gradient_multiphase,
                                c(args, sanitize = FALSE)))))

  # A g3 phase with gamma <= 0 is infeasible too (the g3 early exit).
  g3_args <- list(c(log(0.1), log(1), -1, 1, 1),
                  tt, rep(1, 50), phases = list(late = hzr_phase("g3")),
                  covariate_counts = c(late = 0L), x_list = list(late = NULL))
  expect_true(all(do.call(.hzr_gradient_multiphase, g3_args) == 0))
  expect_true(all(is.na(do.call(.hzr_gradient_multiphase,
                                c(g3_args, sanitize = FALSE)))))

  # A cumulative hazard that overflows (the non-finite H/h guard).
  big <- args
  big[[1]] <- c(log(0.5), log(1), 1, 0, 800)
  expect_true(all(do.call(.hzr_gradient_multiphase, big) == 0))
  expect_true(all(is.na(do.call(.hzr_gradient_multiphase,
                                c(big, sanitize = FALSE)))))

  # One non-finite component at a feasible point (the final zeroing, the
  # site that matters most): a shape derivative that comes back NaN.
  ok <- args
  ok[[1]] <- c(log(0.5), log(1), 1, 0.5, log(0.1))
  real_pd <- .hzr_phase_derivatives
  testthat::local_mocked_bindings(.hzr_phase_derivatives = function(...) {
    pd <- real_pd(...)
    pd$dPhi_dm[1] <- NaN
    pd
  })
  g_sane <- do.call(.hzr_gradient_multiphase, ok)
  g_raw <- do.call(.hzr_gradient_multiphase, c(ok, sanitize = FALSE))
  expect_true(all(is.finite(g_sane)))
  # NA exactly, as documented: a NaN passes is.na() but is not NA_real_.
  expect_identical(g_raw[4], NA_real_)
  expect_true(all(is.finite(g_raw[-4])))
})

test_that("a multiphase fit's acceptance test reaches the gradient unsanitised", {
  # Through the base closure and the fixed-mask wrapper (m is fixed, so the
  # wrapper is on the path): the test asks for the raw score at least once,
  # and the optimizer still gets the sanitised one.
  set.seed(3)
  n <- 300
  tt <- c(stats::rexp(n / 2, 3), stats::rexp(n / 2, 0.15))
  cens <- stats::runif(n, 1, 15)
  d <- data.frame(time = pmin(tt, cens), status = as.integer(tt <= cens))
  real <- .hzr_gradient_multiphase
  seen <- logical(0)
  testthat::local_mocked_bindings(
    .hzr_gradient_multiphase = function(..., sanitize = TRUE) {
      seen <<- c(seen, sanitize)
      real(..., sanitize = sanitize)
    }
  )
  suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0,
                                    fixed = "m"),
                  constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1, conserve = FALSE)
  ))
  expect_true(any(!seen))
  expect_true(any(seen))
})

test_that("the finite-difference check does not straddle m = 0", {
  # A kink at theta[2] = 0 -- slope +5 above, -3 below -- standing in for the
  # cusp where the cdf and hazard families meet.
  kinked <- function(theta) {
    (theta[1] - 1)^2 + ifelse(theta[2] >= 0, 5 * theta[2], -3 * theta[2])
  }
  h <- .Machine$double.eps^(1 / 3)
  for (m in c(h / 3, 0, -h / 3)) {
    g <- .hzr_fd_gradient(kinked, c(0.5, m), sign_bounded = 2L)
    expect_equal(g[[2]], if (m >= 0) 5 else -3, tolerance = 1e-6,
                 label = sprintf("slope at m = %g", m))
  }
  # Without the bound the central stencil mixes the two slopes (2.33 here).
  expect_gt(abs(.hzr_fd_gradient(kinked, c(0.5, h / 3))[[2]] - 5), 1)

  # The 1% cap below 0: a left branch that varies on the scale of |m|, as the
  # cusp does. A step of eps^(1/3) would dwarf |m| and miss the slope.
  cusp <- function(theta) {
    (theta[1] - 1)^2 + ifelse(theta[2] >= 0, 5 * theta[2], sqrt(-theta[2]))
  }
  m <- -h / 3
  expect_equal(.hzr_fd_gradient(cusp, c(0.5, m), sign_bounded = 2L)[[2]],
               -1 / (2 * sqrt(-m)), tolerance = 1e-3)
})

test_that("multiphase marks every free shape m, and only m, as sign-bounded", {
  set.seed(3)
  n <- 300
  tt <- c(stats::rexp(n / 2, 3), stats::rexp(n / 2, 0.15))
  cens <- stats::runif(n, 1, 15)
  d <- data.frame(time = pmin(tt, cens), status = as.integer(tt <= cens),
                  z = stats::rnorm(n))
  real <- .hzr_optim_generic
  got <- list()
  testthat::local_mocked_bindings(
    .hzr_optim_generic = function(..., theta_start, sign_bounded = integer(0)) {
      got[[length(got) + 1L]] <<- list(free = names(theta_start),
                                       bounded = names(theta_start)[sign_bounded])
      real(..., theta_start = theta_start, sign_bounded = sign_bounded)
    }
  )
  fit_phases <- function(phases) {
    got <<- list()
    suppressWarnings(hazard(
      survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
      phases = phases, fit = TRUE, control = list(n_starts = 1)
    ))
    got
  }
  # A: a parameter ahead of a free m leaves the optimizer's vector (e1's nu is
  # fixed), so e1.m's reduced position (3) differs from its full one (4).
  runs <- fit_phases(list(e1 = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0,
                                         formula = ~ z, fixed = "nu"),
                          constant = hzr_phase("constant")))
  expect_true(length(runs) > 0)
  # Premise: without a dropped parameter before e1.m this could not tell full
  # positions from reduced ones.
  expect_false("e1.nu" %in% runs[[1]]$free)
  for (r in runs) expect_identical(r$bounded, "e1.m")
  # B: a fixed m is not in the optimizer's vector and must not be marked.
  runs <- fit_phases(list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0),
                          mid = hzr_phase("hazard", t_half = 5, nu = 1, m = 0,
                                          fixed = "m"),
                          constant = hzr_phase("constant")))
  expect_true(length(runs) > 0)
  for (r in runs) expect_identical(r$bounded, "early.m")
  # C: both families that carry m are marked, cdf and hazard.
  runs <- fit_phases(list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0),
                          mid = hzr_phase("hazard", t_half = 5, nu = 1, m = 0),
                          constant = hzr_phase("constant")))
  expect_true(length(runs) > 0)
  for (r in runs) expect_setequal(r$bounded, c("early.m", "mid.m"))
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
