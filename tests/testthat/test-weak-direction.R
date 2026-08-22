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
  w <- testthat::capture_warnings(
    fit <- hazard(time = time, status = status, dist = "weibull", fit = TRUE)
  )
  expect_false(any(grepl("only in combination", w)))
  expect_null(fit$fit$weak)
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
