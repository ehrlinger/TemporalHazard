# The log_t_half derivative of an early phase at small t_half (#574).
#
# .hzr_phase_derivatives() stepped t_half linearly with an absolute floor on
# the step. Below t_half = 1e-4 the step stopped shrinking with t_half, and
# below about 6e-10 it was larger than t_half itself, so the difference went
# one-sided over several times t_half. At log_t_half = -23.28, nu = 0.2931,
# m = 120 the derivative came back at 0.27 of its value, and the multiphase
# score for log_t_half at -1.34 where finite differences of the
# log-likelihood give -28.1.
#
# The oracle is a Richardson derivative of hzr_decompos() in log_t_half,
# taken twice from base steps a hundred-fold apart. A cell is compared only
# where the two agree, so a disagreement with the package is the package's.

thalf_times <- c(0.02, 0.1, 0.5, 1, 3, 8)

thalf_value <- function(log_t_half, nu, m, type, what, t) {
  d <- hzr_decompos(t, t_half = exp(log_t_half), nu = nu, m = m)
  if (type == "cdf") {
    if (what == "Phi") d$G else d$g
  } else if (what == "Phi") {
    -log(pmax(1 - d$G, .Machine$double.xmin))
  } else {
    d$h
  }
}

thalf_oracle <- function(log_t_half, nu, m, type, what) {
  at <- function(d) {
    vapply(thalf_times, function(t) {
      numDeriv::grad(function(l) thalf_value(l, nu, m, type, what, t),
                     log_t_half, method.args = list(d = d))
    }, numeric(1))
  }
  list(fine = at(1e-4), coarse = at(1e-2))
}

# Worst disagreement over the times, relative to the largest derivative.
thalf_error <- function(log_t_half, nu, m, type, what) {
  o <- thalf_oracle(log_t_half, nu, m, type, what)
  scale <- max(abs(o$fine))
  pd <- .hzr_phase_derivatives(thalf_times, t_half = exp(log_t_half), nu = nu,
                               m = m, type = type)
  got <- pd[[paste0("d", what, "_dlog_thalf")]]
  list(error = max(abs(got - o$fine)) / scale,
       oracle_spread = max(abs(o$fine - o$coarse)) / scale,
       scale = scale)
}

test_that("the log_t_half derivative holds down to very small t_half (#574)", {
  skip_if_not_installed("numDeriv")
  # The issue's shape. The first two rows agreed before the fix and are the
  # known positives; from -18 down they did not.
  for (lth in c(-1.2, -9.2, -16, -18, -20, -22, -23.28, -25)) {
    for (what in c("Phi", "phi")) {
      e <- thalf_error(lth, nu = 0.2931, m = 120, type = "cdf", what = what)
      label <- paste0("cdf ", what, " at log_t_half = ", lth)
      # Premises: there is a derivative to compare, and the oracle is settled.
      expect_gt(e$scale, 1e-6, label = paste(label, "(oracle size)"))
      expect_lt(e$oracle_spread, 1e-7, label = paste(label, "(oracle spread)"))
      expect_lt(e$error, 1e-5, label = label)
    }
  }
})

test_that("cells that agreed before the fix still agree (#574)", {
  skip_if_not_installed("numDeriv")
  # Ordinary scales, where the step is what it always was.
  for (lth in c(-1.2, -9.2, 3)) {
    for (type in c("cdf", "hazard")) {
      e <- thalf_error(lth, nu = 1, m = 1, type = type, what = "Phi")
      expect_lt(e$error, 1e-6, label = paste(type, "Phi at log_t_half =", lth))
    }
  }
})

test_that("from t_half = 1e-4 down to the cap the step is the one it was (#574)", {
  # Near saturation the phase's values carry rounding noise, and a step held
  # at its ordinary size below t_half = 1e-4 divides that noise by too
  # little: against a Richardson oracle it lost digits in 18 of 1096 cells
  # that had agreed. So the step grows as 1e-4 / t_half there, as the old
  # linear step did, and the result must be the old one. The reference is
  # that old difference, written out.
  eps_rel <- .Machine$double.eps^(1 / 3)
  Phi <- function(th) {
    -log(pmax(1 - hzr_decompos(thalf_times, t_half = th, nu = 1, m = 1)$G,
              .Machine$double.xmin))
  }
  cells <- list(list(log_t_half = -12, tol = 5e-7),
                list(log_t_half = -14, tol = 4e-6))
  for (cell in cells) {
    th <- exp(cell$log_t_half)
    h <- eps_rel * 1e-4
    old <- th * (Phi(th + h) - Phi(th - h)) / (2 * h)
    got <- .hzr_phase_derivatives(thalf_times, t_half = th, nu = 1, m = 1,
                                  type = "hazard")$dPhi_dlog_thalf
    # The derivative is about -1 at every time here.
    expect_gt(max(abs(old)), 0.5)
    expect_lt(max(abs(got - old)) / max(abs(old)), cell$tol,
              label = paste("hazard Phi at log_t_half =", cell$log_t_half))
  }
})

test_that("a subnormal t_half is differenced over its real spacing (#574)", {
  # For nu = 1, m = 1 and t_half far below t, the cdf density is t_half / t^2
  # to leading order, and so is its derivative in log_t_half. Subnormal
  # points are rounded to a coarse grid, so dividing by the nominal step put
  # this 32% off at t_half = 1e-320. Times small enough that the density
  # itself is still well above the smallest subnormal.
  tt <- c(0.1, 0.5)
  for (th in c(1e-300, 1e-318, 3e-320, 1e-320)) {
    got <- .hzr_phase_derivatives(tt, t_half = th, nu = 1, m = 1,
                                  type = "cdf")$dphi_dlog_thalf
    expect_equal(got / (th / tt^2), c(1, 1), tolerance = 0.02,
                 label = paste("t_half =", th))
  }
})

test_that("a t_half with no room to step returns NaN, not a clean zero (#574)", {
  # In the denormal range the two points of the difference coincide. Zero
  # would read as a derivative; NaN cannot be mistaken for one.
  pd <- .hzr_phase_derivatives(c(0.5, 2), t_half = 5e-324, nu = 1, m = 1,
                               type = "cdf")
  expect_true(all(is.nan(pd$dPhi_dlog_thalf)))
  expect_true(all(is.nan(pd$dphi_dlog_thalf)))
  # The other two derivatives do not step t_half and are unaffected.
  expect_true(all(is.finite(pd$dPhi_dnu)))
  expect_true(all(is.finite(pd$dPhi_dm)))
  # At the largest double the upper point overflows and the spacing is
  # infinite: a finite difference over it is 0, where just below that the
  # derivative is as large as -0.05 at these times.
  big <- c(1e300, 1e305, 1e307)
  pd <- .hzr_phase_derivatives(big, t_half = .Machine$double.xmax, nu = 1,
                               m = 1, type = "cdf")
  expect_true(all(is.nan(pd$dPhi_dlog_thalf)))
  below <- .hzr_phase_derivatives(big, t_half = .Machine$double.xmax * (1 - 1e-5),
                                  nu = 1, m = 1, type = "cdf")
  expect_lt(min(below$dPhi_dlog_thalf), -0.01)
})

test_that("the multiphase score for log_t_half is right at the issue's point (#574)", {
  skip_if_not_installed("numDeriv")
  set.seed(11)
  n <- 300
  z <- rbinom(n, 1, 0.4)
  age <- round(rnorm(n), 2)
  t_event <- ifelse(runif(n) < 0.3, rexp(n, 4 * exp(0.5 * z)),
                    rexp(n, 0.15 * exp(0.3 * age)))
  cens <- runif(n, 0.5, 8)
  d <- data.frame(stop = pmin(t_event, cens) + 0.01,
                  event = as.integer(t_event <= cens), z = z, age = age)
  # What the optimizer is handed, captured without fitting.
  seen <- NULL
  local_mocked_bindings(.hzr_optim_generic = function(...) {
    seen <<- list(...)
    stop("captured")
  })
  theta <- c(2.524, -23.28, 0.2931, 120, 0.58, -1.5, 0.40)
  suppressWarnings(try(hazard(
    survival::Surv(stop, event) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, formula = ~ z),
      constant = hzr_phase("constant", formula = ~ age)
    ),
    theta = theta, fit = TRUE, control = list(n_starts = 1, conserve = FALSE)
  ), silent = TRUE))
  expect_false(is.null(seen))
  expect_equal(unname(seen$theta_start), theta)
  obj <- function(p) {
    seen$logl_fn(p, seen$time, seen$status, seen$time_lower, seen$time_upper,
                 seen$x, weights = seen$weights)
  }
  score <- seen$gradient_fn(theta, seen$time, seen$status, seen$time_lower,
                            seen$time_upper, seen$x, weights = seen$weights)
  oracle <- numDeriv::grad(obj, theta)
  # The component #574 is about is far from zero here.
  expect_gt(abs(oracle[2]), 1)
  expect_equal(unname(score)[2], oracle[2], tolerance = 1e-5)
  expect_equal(unname(score), oracle, tolerance = 1e-5)
})
