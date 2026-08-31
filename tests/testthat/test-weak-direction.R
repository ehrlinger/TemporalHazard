# Ridge (weak-identification) detection.
#
# An ill-conditioned Hessian already warns that *standard errors* are
# unreliable. That understates the problem when the ill-conditioning is a
# ridge: there the individual point estimates along the flat direction are not
# determined by the data either, only their combination is. These tests pin
# down the detector that names which parameters are jointly identified.
#
# The detector works on the correlation of the estimates rather than on raw
# variances, so it reports trade-offs between parameters. A single imprecise
# parameter is not a ridge and is left to the existing rcond warning.

ridge_vcov <- function(r = 0.99999, scale = c(1, 1)) {
  V <- matrix(c(1, r, r, 1), 2, 2)
  D <- diag(scale)
  V <- D %*% V %*% D
  dimnames(V) <- list(c("a", "b"), c("a", "b"))
  V
}

test_that("a well-conditioned fit has no weak direction", {
  V <- solve(matrix(c(4, 1, 1, 3), 2, 2))
  expect_null(.hzr_weak_direction(V, rcond = 0.3, param_names = c("a", "b")))
})

test_that("rcond above tolerance suppresses detection even on a ridge", {
  # The gate is the package's existing ill-conditioning threshold, so the
  # ridge note never fires where the rcond warning stays silent.
  expect_null(.hzr_weak_direction(ridge_vcov(), rcond = 1e-3,
                                  param_names = c("a", "b")))
})

test_that("a two-parameter ridge names both parameters", {
  w <- .hzr_weak_direction(ridge_vcov(), rcond = 1e-10,
                           param_names = c("a", "b"))
  expect_false(is.null(w))
  expect_setequal(w$params, c("a", "b"))
  expect_gt(abs(w$correlation), 0.99)
})

test_that("the flat direction is scale-invariant", {
  # THE REGRESSION THIS FILE EXISTS FOR. The same ridge, but with the second
  # parameter measured on a 1000x larger scale -- as m (~27) is against nu
  # (~0.027) on the multiphase fixture. A direction taken from the raw
  # covariance loads almost entirely on the larger-scale parameter and reports
  # the ridge as a single unidentified parameter; only the correlation-scaled
  # form recovers the true two-parameter combination.
  V <- ridge_vcov(scale = c(1, 1000))

  raw <- eigen(V, symmetric = TRUE)
  raw_w <- raw$vectors[, which.max(raw$values)]^2
  expect_gt(max(raw_w), 0.99)          # raw units: one parameter dominates

  w <- .hzr_weak_direction(V, rcond = 1e-10, param_names = c("a", "b"))
  expect_false(is.null(w))
  expect_setequal(w$params, c("a", "b"))
})

test_that("an imprecise but uncorrelated parameter is not a ridge", {
  # Large variance on its own is not a trade-off. The existing rcond warning
  # and the parameter's own standard error already cover this case.
  V <- diag(c(1, 1, 1e6))
  expect_null(.hzr_weak_direction(V, rcond = 1e-10,
                                  param_names = c("a", "b", "c")))
})

test_that("a moderate correlation is not a ridge", {
  # Two parameters that do span the flat direction but trade off only weakly.
  # This is the case that exercises the correlation threshold itself: the
  # uncorrelated test above exits earlier, at the "direction spans a single
  # parameter" guard, and so never reaches it.
  V <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  expect_null(.hzr_weak_direction(V, rcond = 1e-10,
                                  param_names = c("a", "b")))
  # Just over the threshold, the same shape is reported.
  V2 <- matrix(c(1, 0.995, 0.995, 1), 2, 2)
  expect_setequal(
    .hzr_weak_direction(V2, rcond = 1e-10, param_names = c("a", "b"))$params,
    c("a", "b")
  )
})

test_that("degenerate covariances yield no weak direction rather than an error", {
  nms <- c("a", "b")
  expect_null(.hzr_weak_direction(NULL, rcond = 1e-10, param_names = nms))
  expect_null(.hzr_weak_direction(NA, rcond = 1e-10, param_names = nms))
  expect_null(.hzr_weak_direction(ridge_vcov(), rcond = NA_real_,
                                  param_names = nms))
  expect_null(.hzr_weak_direction(matrix(c(1, NA, NA, 1), 2, 2),
                                  rcond = 1e-10, param_names = nms))
  # A non-estimated (fixed) parameter carries an NA row/column.
  V <- ridge_vcov()
  V3 <- matrix(NA_real_, 3, 3)
  V3[1:2, 1:2] <- V
  expect_setequal(
    .hzr_weak_direction(V3, rcond = 1e-10,
                        param_names = c("a", "b", "fixed"))$params,
    c("a", "b")
  )
  # Zero variance makes the correlation undefined.
  expect_null(.hzr_weak_direction(diag(c(0, 1)), rcond = 1e-10,
                                  param_names = nms))
})

test_that("unnamed parameters fall back to positional labels", {
  w <- .hzr_weak_direction(ridge_vcov(), rcond = 1e-10, param_names = NULL)
  expect_false(is.null(w))
  expect_length(w$params, 2)
})

test_that("the multiphase ridge fixture warns and records nu and m", {
  skip_on_cran()
  set.seed(42)
  n <- 100
  time   <- rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                 const = hzr_phase("constant"))

  w <- testthat::capture_warnings(
    fit <- hazard(time = time, status = status, dist = "multiphase",
                  phases = phases, fit = TRUE,
                  control = list(n_starts = 5, maxit = 500))
  )
  # The perturbed start wins and lands on the m*nu ridge, where m (~27) and nu
  # (~0.027) are identified only through their product.
  expect_true(any(grepl("only in combination", w)))
  expect_false(is.null(fit$fit$weak))
  expect_setequal(fit$fit$weak$params, c("early.nu", "early.m"))
})

test_that("a well-identified fit stays silent and records no weak direction", {
  skip_on_cran()
  set.seed(7)
  n <- 300
  time   <- rweibull(n, shape = 1.5, scale = 2)
  status <- rep(1L, n)
  # `theta` is required: hazard()'s single-distribution fit branch is
  # `else if (fit && !is.null(theta))`, so without it nothing is optimised,
  # $vcov and $rcond come back NULL, and .hzr_weak_direction() returns at its
  # first guard.  Both assertions below then pass with the whole detector
  # deleted -- which is what this test used to do.
  w <- testthat::capture_warnings(
    fit <- hazard(time = time, status = status, dist = "weibull",
                  theta = c(0.5, 1.0), fit = TRUE)
  )
  # Guard the guard: assert the fit actually reached the detector, or the
  # two assertions below are vacuous again.
  expect_false(is.null(fit$fit$rcond))
  expect_false(is.null(fit$fit$vcov))
  expect_false(any(grepl("only in combination", w)))
  expect_null(fit$fit$weak)
})

test_that("a ridge is found even when a larger correlated block outranks it", {
  # Regression for the leading-eigenvector bug: gating only the top
  # eigenvector on cor_tol missed a near-perfect two-parameter ridge whenever
  # some block of k >= 3 moderately correlated parameters carried more
  # variance.  An equicorrelated k-block has eigenvalue 1 + (k - 1) * r, so it
  # passes the ridge's ceiling of 2 as soon as r > 1 / (k - 1) -- 0.5 at
  # k = 3.  The detector then returned NULL, which is the documented signal
  # for "well identified".
  ridge_with_block <- function(rho) {
    R <- diag(5)
    R[1, 2] <- R[2, 1] <- 0.99999
    for (i in 3:5) {
      for (j in 3:5) if (i != j) R[i, j] <- rho
    }
    .hzr_weak_direction(R, rcond = 1e-10,
                        param_names = paste0("p", seq_len(5)))
  }

  # Below the crossover the block loses on variance and the old code worked.
  expect_setequal(ridge_with_block(0.45)$params, c("p1", "p2"))

  # Above it the block wins on variance but is not a trade-off. These are the
  # values that returned NULL before the scan was introduced.
  for (rho in c(0.51, 0.60, 0.75)) {
    w <- ridge_with_block(rho)
    expect_false(is.null(w))
    expect_setequal(w$params, c("p1", "p2"))
    expect_gt(abs(w$correlation), 0.999)
  }
})

test_that("summary() prints the ridge note when one is recorded", {
  skip_on_cran()
  set.seed(42)
  n <- 100
  time   <- rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                 const = hzr_phase("constant"))
  fit <- suppressWarnings(
    hazard(time = time, status = status, dist = "multiphase", phases = phases,
           fit = TRUE, control = list(n_starts = 5, maxit = 500)))
  # The note is strwrap()'d for the console, so collapse whitespace before
  # matching rather than assuming where the line breaks land.
  out <- paste(capture.output(print(summary(fit))), collapse = " ")
  out <- gsub("\\s+", " ", out)
  expect_match(out, "only in combination")
  expect_match(out, "early\\.nu")
})

test_that("a single-distribution ridge is reported under the real parameter names", {
  skip_on_cran()
  # optim() returns an unnamed `par` for the single-distribution fits, so
  # reading names(fit$par) alone labelled the warning "par3"/"par4" while
  # summary() printed "beta1"/"beta2" against the same numbers.  The names
  # must come from the same source summary() uses.
  set.seed(3)
  n <- 250
  x1 <- stats::rnorm(n)
  x2 <- x1 + stats::rnorm(n, sd = 1e-4)   # near-collinear: a genuine ridge
  df <- data.frame(
    tt = stats::rweibull(n, 1.4, 2),
    st = stats::rbinom(n, 1, 0.9),
    x1 = x1, x2 = x2
  )
  fit <- suppressWarnings(hazard(
    survival::Surv(tt, st) ~ x1 + x2, data = df, dist = "weibull",
    theta = c(0.5, 1, 0, 0), fit = TRUE
  ))

  # Assert the fixture is still a ridge, or the naming assertion below is
  # vacuous.
  expect_false(is.null(fit$fit$weak))
  expect_null(names(fit$fit$par))          # the condition the fallback exists for
  expect_setequal(fit$fit$weak$params, c("beta1", "beta2"))
})
