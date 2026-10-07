# Numerical Hessians at a covariate offset (#598).
#
# Every numDeriv::hessian() site used numDeriv's default step, 0.1 * |par|.
# A covariate offset by c moves the intercept by about -beta * c, so that
# step spans many standard errors of the intercept and the second difference
# reads curvature far from the maximum. On main 440c0ebf, a Weibull fit with
# left- and interval-censored rows and x offset by 100 reached the centred
# fit's log-likelihood to 3e-8 and reported SE(beta) 0.035 of the centred
# fit's, with no warning; an exponential fit at offset 1000 reported 1e-28.
#
# The oracle is the centred fit (the same model, so the same SEs for every
# parameter but the intercept), or an analytic Hessian where one exists.

.n598_data <- function(seed = 598) {
  set.seed(seed)
  n <- 300
  x <- stats::rnorm(n)
  t <- stats::rweibull(n, shape = 1.3, scale = 4 * exp(-0.5 * x / 1.3))
  cc <- stats::runif(n, 1, 12)
  time <- pmin(t, cc)
  status <- as.integer(t <= cc)
  tl <- time
  tu <- time
  ev <- which(status == 1)
  ic <- ev[seq(1, length(ev), by = 4)]
  rest <- setdiff(ev, ic)
  lc <- rest[seq(1, length(rest), by = 7)]
  status[ic] <- 2L
  tl[ic] <- 0.7 * time[ic]
  tu[ic] <- 1.2 * time[ic]
  status[lc] <- -1L
  tu[lc] <- 1.1 * time[lc]
  data.frame(time = time, status = status, tl = tl, tu = tu, x = x)
}

.n598_starts <- list(weibull = c(0.2, 1, 0), exponential = c(log(0.3), 0),
                     loglogistic = c(0, 0.2, 0), lognormal = c(1, 0, 0))

.n598_fit <- function(d, dist, off) {
  suppressWarnings(hazard(
    time = d$time, status = d$status, time_lower = d$tl, time_upper = d$tu,
    x = matrix(d$x + off, dimnames = list(NULL, "x")), dist = dist,
    theta = .n598_starts[[dist]], fit = TRUE))
}

# Standard errors from a Hessian of the negative log-likelihood.
.n598_se <- function(h) sqrt(diag(solve((h + t(h)) / 2)))

test_that("left- and interval-censored SEs do not move with a covariate offset (#598)", {
  d <- .n598_data()
  # Premise: the numerical Hessian is reached because left- and
  # interval-censored rows are present (the analytic forms decline them).
  expect_equal(as.vector(table(factor(d$status, levels = -1:2))),
               c(24, 80, 141, 55))
  for (dist in names(.n598_starts)) {
    c0 <- .n598_fit(d, dist, 0)
    p <- length(.n598_starts[[dist]])
    se0 <- sqrt(diag(c0$fit$vcov))
    for (off in c(100, 1000)) {
      f <- .n598_fit(d, dist, off)
      lab <- paste(dist, "offset", off)
      # Premise: the same maximum, so any SE difference is the Hessian's.
      expect_lt(abs(f$fit$objective - c0$fit$objective), 1e-6, label = lab)
      expect_true(is.matrix(f$fit$vcov), label = lab)
      if (!is.matrix(f$fit$vcov)) next
      se <- sqrt(diag(f$fit$vcov))
      # Every parameter but the intercept (index 1) is offset-invariant.
      expect_equal(unname(se[-1] / se0[-1]), rep(1, p - 1),
                   tolerance = 1e-3, label = lab)
    }
  }
})

test_that("a multiphase fit with interval-censored rows keeps its SEs at an offset (#598)", {
  d <- .n598_data()
  mp <- function(off) {
    dd <- transform(d, x = x + off)
    suppressWarnings(hazard(
      time = dd$time, status = dd$status, time_lower = dd$tl,
      time_upper = dd$tu, data = dd, dist = "multiphase",
      phases = list(hzr_phase("constant", formula = ~ x)),
      theta = c(log(0.2), 0), fit = TRUE))
  }
  c0 <- mp(0)
  f <- mp(100)
  expect_lt(abs(f$fit$objective - c0$fit$objective), 1e-6)
  expect_true(is.matrix(f$fit$vcov))
  se0 <- sqrt(diag(c0$fit$vcov))
  se <- sqrt(diag(f$fit$vcov))
  expect_equal(unname(se[2] / se0[2]), 1, tolerance = 1e-3)
})

test_that(".hzr_numeric_hessian() matches the exact Weibull Hessian at offsets (#598)", {
  set.seed(1)
  n <- 300
  x <- stats::rnorm(n)
  t <- stats::rweibull(n, shape = 1.4, scale = 3 * exp(-0.5 * x / 1.4))
  cc <- stats::runif(n, 1, 10)
  time <- pmin(t, cc)
  st <- as.integer(t <= cc)
  # The centred MLE on the internal scale (alpha, psi, beta), from the fit.
  f0 <- suppressWarnings(hazard(time = time, status = st, x = matrix(x),
                                dist = "weibull", theta = c(0.2, 1, 0),
                                fit = TRUE))
  th <- unname(f0$fit$theta)
  phi0 <- c(th[2] * log(th[1]), log(th[2]), th[3])
  # A singular Hessian (numDeriv's default at offset 1000) is infinitely off.
  err <- function(h, h_true) {
    se <- tryCatch(.n598_se(h), error = function(e) NA_real_)
    e <- max(abs(se / .n598_se(h_true) - 1))
    if (is.na(e)) Inf else e
  }
  default_err <- numeric(0)
  for (off in c(0, 10, 100, 1000)) {
    xo <- x + off
    # The same maximum at every offset: alpha absorbs -beta * offset.
    phi <- phi0 - c(phi0[3] * off, 0, 0)
    nll <- function(ph) {
      eta <- ph[1] + ph[3] * xo
      g <- exp(ph[2])
      -sum(st * (eta + log(g) + (g - 1) * log(time)) - exp(eta) * time^g)
    }
    h_true <- .hzr_hessian_weibull_internal(phi, time, st, x = matrix(xo))
    # Premise, at offset 0: this nll is the objective h_true differentiates.
    if (off == 0) {
      expect_lt(err(numDeriv::hessian(nll, phi, method.args = list(d = 1e-4)),
                    h_true), 1e-3)
    }
    default_err[as.character(off)] <- err(numDeriv::hessian(nll, phi), h_true)
    expect_lt(err(.hzr_numeric_hessian(nll, phi), h_true), 1e-4,
              label = paste("offset", off))
  }
  # Premise: the comparison can fail. numDeriv's default step is far off at
  # offset 100, so a helper that returned it would not pass the above.
  expect_gt(default_err[["100"]], 0.1)
})

test_that(".hzr_numeric_hessian() reports a saddle as a saddle (#598)", {
  fn <- function(p) 0.5 * (p[1]^2 - 2 * p[2]^2) + 0.3 * p[1] * p[2]
  h_true <- matrix(c(1, 0.3, 0.3, -2), 2)
  h <- .hzr_numeric_hessian(fn, c(0.4, -0.2))
  expect_equal(h, h_true, tolerance = 1e-6)
  expect_lt(min(eigen((h + t(h)) / 2, only.values = TRUE)$values), 0)
})

test_that(".hzr_numeric_hessian() caps the step on a direction with almost no information (#598)", {
  # Curvature 1e-4 along p[2]: the whitened step would be 0.1 / sqrt(1e-4)
  # = 10, where exp() is nothing like its quadratic at the point. numDeriv's
  # own first step there, 0.1 * |p[2]| = 0.1, reads it accurately. (A
  # polynomial would not show this: Richardson extrapolation is exact for a
  # quartic at any step.)
  fn <- function(p) 0.5 * p[1]^2 + 1e-4 * exp(p[2] - 1)
  h <- .hzr_numeric_hessian(fn, c(0.3, 1))
  expect_equal(h[2, 2] / 1e-4, 1, tolerance = 1e-6)
  expect_equal(h[1, 1], 1, tolerance = 1e-6)
  # Premise: uncapped, the step would read that far out.
  expect_gt(0.1 / sqrt(h[2, 2]) / (0.1 * 1), 50)
})

test_that("the score paths' numerical Hessians match the analytic ones at an offset (#598)", {
  set.seed(2)
  n <- 300
  x <- stats::rnorm(n) + 100
  t <- stats::rexp(n, 0.2 * exp(0.5 * (x - 100)))
  cc <- stats::runif(n, 1, 10)
  time <- pmin(t, cc)
  st <- as.integer(t <= cc)
  # Single distribution: right-censored exponential takes the analytic
  # Hessian in the fit, so its vcov is the oracle for the score path's
  # numerical one (on the same, stored, scale).
  fe <- suppressWarnings(hazard(time = time, status = st,
                                x = matrix(x, dimnames = list(NULL, "x")),
                                dist = "exponential", theta = c(-50, 0.5),
                                fit = TRUE))
  h <- .hzr_score_single_hessian(fe, fe$data$x, fe$fit$theta)
  expect_equal(unname(.n598_se(h) / sqrt(diag(fe$fit$vcov))), c(1, 1),
               tolerance = 1e-4)

  # Multiphase: the numerical fallback against the analytic Hessian it
  # stands in for, on right-censored data where both exist.
  dd <- data.frame(time = time, status = st, x = x)
  fm <- suppressWarnings(hazard(
    time = dd$time, status = dd$status, data = dd, dist = "multiphase",
    phases = list(hzr_phase("constant", formula = ~ x)),
    theta = c(-50, 0.5), fit = TRUE))
  th <- fm$fit$theta
  phases <- .hzr_score_phases(fm)
  h_num <- .hzr_score_multiphase_hessian(fm, th, phases,
                                         fm$fit$covariate_counts,
                                         fm$fit$x_list)
  h_an <- .hzr_hessian_multiphase(
    th, time = dd$time, status = dd$status, x = fm$data$x,
    weights = fm$data$weights, phases = phases,
    covariate_counts = fm$fit$covariate_counts, x_list = fm$fit$x_list)
  expect_equal(unname(.n598_se(h_num) / .n598_se(h_an)), c(1, 1),
               tolerance = 1e-4)
})
