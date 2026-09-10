# test-weibull-hessian.R -- numDeriv cross-check for .hzr_hessian_weibull_internal()
#
# Reconstruct the optimizer's internal objective (negative log-likelihood on the
# (alpha, psi, beta) scale). Post-PR#58 the natural-scale .hzr_logl_weibull
# equals logl_internal for event+right data, so we can map phi -> (mu, nu, beta)
# and reuse it.
.wb_obj_internal <- function(phi, time, status, time_lower = NULL,
                             time_upper = NULL, x = NULL, weights = NULL) {
  alpha <- phi[1]
  g <- exp(phi[2])
  beta <- if (length(phi) > 2) phi[3:length(phi)] else numeric(0)
  mu <- exp(alpha / g)
  theta_nat <- c(mu, g, beta)
  -.hzr_logl_weibull(theta_nat, time, status, time_lower, time_upper,
                     x = x, weights = weights)
}

test_that(".hzr_hessian_weibull_internal matches numDeriv (no covariates)", {
  skip_if_not_installed("numDeriv")
  set.seed(21)
  n <- 300
  time <- rexp(n, 0.6) + 0.01
  status <- rbinom(n, 1, 0.75)
  phi <- c(alpha = 0.4, psi = log(1.5))
  obj <- function(p) .wb_obj_internal(p, time, status)
  H_an <- .hzr_hessian_weibull_internal(phi, time, status)
  H_nd <- numDeriv::hessian(obj, phi)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_weibull_internal matches numDeriv (covariates + weights)", {
  skip_if_not_installed("numDeriv")
  set.seed(22)
  n <- 400
  x <- cbind(a = rnorm(n), b = rbinom(n, 1, 0.4))
  eta <- as.numeric(x %*% c(0.3, -0.4))
  time <- rweibull(n, shape = 1.4, scale = exp(-eta / 1.4)) + 0.01
  status <- rbinom(n, 1, 0.8)
  w <- runif(n, 0.5, 2)
  phi <- c(alpha = 0.2, psi = log(1.4), 0.3, -0.4)
  obj <- function(p) .wb_obj_internal(p, time, status, x = x, weights = w)
  H_an <- .hzr_hessian_weibull_internal(phi, time, status, x = x, weights = w)
  H_nd <- numDeriv::hessian(obj, phi)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_weibull_internal matches numDeriv (counting-process start)", {
  skip_if_not_installed("numDeriv")
  set.seed(23)
  n <- 250
  start <- runif(n, 0, 0.5)
  time <- start + rexp(n, 0.7) + 0.05
  status <- rbinom(n, 1, 0.7)
  phi <- c(alpha = 0.3, psi = log(1.2))
  obj <- function(p) .wb_obj_internal(p, time, status, time_lower = start)
  H_an <- .hzr_hessian_weibull_internal(phi, time, status, time_lower = start)
  H_nd <- numDeriv::hessian(obj, phi)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_weibull_internal returns NULL for left/interval censoring", {
  set.seed(24)
  n <- 40
  time <- rexp(n, 0.5) + 0.01
  status <- rep(1L, n)
  status[1:5] <- 2L
  expect_null(.hzr_hessian_weibull_internal(
    c(alpha = 0, psi = 0), time, status,
    time_lower = time, time_upper = time + 0.5))
  status2 <- rep(1L, n)
  status2[1:5] <- -1L
  expect_null(.hzr_hessian_weibull_internal(
    c(alpha = 0, psi = 0), time, status2, time_upper = time + 0.5))
})

test_that(".hzr_hessian_weibull_internal handles right-censored time = 0 rows", {
  skip_if_not_installed("numDeriv")
  set.seed(27)
  n <- 60
  time <- c(0, rexp(n - 1, 0.6) + 0.01)      # row 1: right-censored at time 0
  status <- c(0L, rbinom(n - 1, 1, 0.7))
  phi <- c(alpha = 0.3, psi = log(1.3))
  H_an <- .hzr_hessian_weibull_internal(phi, time, status)
  expect_true(all(is.finite(H_an)))
  obj <- function(p) .wb_obj_internal(p, time, status)
  H_nd <- numDeriv::hessian(obj, phi)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that("weibull fit vcov uses the analytic Hessian (matches numDeriv+delta)", {
  skip_if_not_installed("numDeriv")
  set.seed(25)
  n <- 500
  z <- rnorm(n)
  time <- rweibull(n, shape = 1.5, scale = exp(-0.5 * z / 1.5)) + 0.01
  status <- rbinom(n, 1, 0.85)
  fit <- hazard(
    survival::Surv(time, status) ~ z,
    data = data.frame(time, status, z),
    dist = "weibull", theta = c(mu = 1, nu = 1, z = 0), fit = TRUE
  )
  # Reference: numDeriv internal-scale Hessian -> invert -> delta-method J.
  th <- fit$fit$theta              # natural (mu, nu, beta)
  mu <- th[1]
  nu <- th[2]
  beta <- th[-(1:2)]
  alpha <- unname(nu * log(mu))
  psi <- unname(log(nu))
  phi <- c(alpha, psi, unname(beta))
  obj <- function(p) .wb_obj_internal(p, time, status, x = cbind(z))
  v_int <- solve(numDeriv::hessian(obj, phi))
  pdim <- length(phi)
  jmat <- diag(pdim)
  jmat[1, 1] <- mu / nu
  jmat[1, 2] <- -mu * alpha / nu
  jmat[2, 2] <- nu
  v_ref <- jmat %*% v_int %*% t(jmat)
  expect_equal(unname(fit$fit$vcov), unname(v_ref), tolerance = 1e-3)
  expect_true(isTRUE(fit$fit$pd))
})

test_that("a masked internal variance keeps the weibull vcov on the (mu, nu) scale", {
  # .hzr_safe_solve() NAs the row and column of a non-positive variance. The
  # delta-method transform used to be skipped whenever the vcov held any NA,
  # leaving the surviving SEs on the (alpha, psi) scale beside (mu, nu)
  # estimates. Force the mask on a healthy fit and compare to the unmasked SEs.
  set.seed(25)
  n <- 500
  z <- rnorm(n)
  time <- rweibull(n, shape = 1.5, scale = exp(-0.5 * z / 1.5)) + 0.01
  status <- rbinom(n, 1, 0.85)
  dat <- data.frame(time, status, z)
  fit_weibull <- function() {
    hazard(survival::Surv(time, status) ~ z, data = dat, dist = "weibull",
           theta = c(mu = 1, nu = 1, z = 0), fit = TRUE)
  }
  v_ref <- unname(fit_weibull()$fit$vcov)
  real_solve <- .hzr_safe_solve
  vcov_masked <- function(k) {
    local_mocked_bindings(.hzr_safe_solve = function(H, ...) {
      res <- real_solve(H, ...)
      res$vcov[k, ] <- NA_real_
      res$vcov[, k] <- NA_real_
      res
    })
    unname(fit_weibull()$fit$vcov)
  }

  # alpha masked: nu = exp(psi) depends on psi alone, so the (nu, z) block
  # survives and must equal the unmasked one -- SE(nu), not SE(psi) = SE(nu) / nu.
  v1 <- vcov_masked(1)
  expect_equal(is.na(v1), outer(1:3, 1:3, function(i, j) i == 1 | j == 1))
  expect_equal(v1[2:3, 2:3], v_ref[2:3, 2:3])

  # psi masked: mu depends on psi through dmu/dpsi, so it is masked as well.
  v2 <- vcov_masked(2)
  expect_equal(is.na(v2), outer(1:3, 1:3, function(i, j) i <= 2 | j <= 2))
  expect_equal(v2[3, 3], v_ref[3, 3])
})

test_that("weibull SEs are invariant to covariate rescaling", {
  set.seed(26)
  n <- 600
  z <- rnorm(n)
  time <- rweibull(n, shape = 1.3, scale = exp(-0.4 * z / 1.3)) + 0.01
  status <- rbinom(n, 1, 0.9)
  f1 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z),
               dist = "weibull", theta = c(mu = 1, nu = 1, z = 0), fit = TRUE)
  f2 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z / 100),
               dist = "weibull", theta = c(mu = 1, nu = 1, z = 0), fit = TRUE)
  se1 <- sqrt(diag(f1$fit$vcov))
  se2 <- sqrt(diag(f2$fit$vcov))
  # mu and nu SEs are invariant to covariate rescaling.
  expect_equal(unname(se1[1]), unname(se2[1]), tolerance = 1e-2)
  expect_equal(unname(se1[2]), unname(se2[2]), tolerance = 1e-2)
})

test_that(".hzr_hessian_weibull_internal errors on inconsistent theta/x", {
  set.seed(62)
  n <- 30
  time <- rexp(n, 0.5) + 0.01
  status <- rbinom(n, 1, 0.7)
  expect_error(
    .hzr_hessian_weibull_internal(c(alpha = 0, psi = 0, z = 0.3), time, status, x = NULL),
    "ncol")
  x2 <- cbind(a = rnorm(n), b = rnorm(n))
  expect_error(
    .hzr_hessian_weibull_internal(c(alpha = 0, psi = 0, z = 0.3), time, status, x = x2),
    "ncol")
})
