test_that("hazard() builds a hazard object", {
  x <- matrix(c(1, 0, 0, 1), ncol = 2)
  fit <- hazard(
    time = c(1, 2),
    status = c(1, 0),
    x = x,
    theta = c(0.2, -0.1),
    dist = "weibull",
    maxit = 50
  )

  expect_s3_class(fit, "hazard")
  expect_equal(fit$spec$dist, "weibull")
  expect_equal(fit$legacy_args$maxit, 50)
})

test_that("hazard() validates core dimensions", {
  expect_error(
    hazard(time = c(1, 2), status = 1),
    "same length"
  )

  expect_error(
    hazard(time = c(1, 2), status = c(1, 0), x = matrix(1, nrow = 3, ncol = 1)),
    "rows must match"
  )
})

test_that("predict.hazard returns linear predictor and hazard scale", {
  x <- matrix(c(1, 0, 0, 1), ncol = 2)
  fit <- hazard(time = c(1, 2), status = c(1, 0), x = x, theta = c(0.3, -0.2))

  eta <- predict(fit, type = "linear_predictor")
  hz <- predict(fit, type = "hazard")

  expect_equal(eta, c(0.3, -0.2), tolerance = 1e-12)
  expect_equal(hz, exp(eta), tolerance = 1e-12)
})

test_that("predict.hazard accepts newdata", {
  fit <- hazard(
    time = c(1, 2),
    status = c(1, 0),
    x = matrix(c(1, 2, 3, 4), ncol = 2),
    theta = c(0.5, 0.25)
  )

  newdata <- matrix(c(2, 1, 0, 1), ncol = 2)
  eta <- predict(fit, newdata = newdata, type = "linear_predictor")

  expect_equal(eta, as.numeric(newdata %*% c(0.5, 0.25)), tolerance = 1e-12)
})

test_that("predict.hazard errors when theta is missing", {
  fit <- hazard(time = c(1, 2), status = c(1, 0), x = matrix(c(1, 0), ncol = 1))
  expect_error(predict(fit), "No coefficients")
})

test_that("summary.hazard returns model summary metadata", {
  fit <- hazard(
    time = c(1, 2, 3),
    status = c(1, 0, 1),
    x = matrix(c(1, 0, 0, 1, 1, 1), ncol = 2),
    theta = c(0.3, -0.2),
    dist = "weibull"
  )

  s <- summary(fit)

  expect_s3_class(s, "summary.hazard")
  expect_equal(s$n, 3)
  expect_equal(s$p, 2)
  expect_equal(s$dist, "weibull")
  expect_equal(rownames(s$coefficients), c("beta1", "beta2"))
  expect_equal(s$coefficients$estimate, c(0.3, -0.2), tolerance = 1e-12)
})

test_that("summary.hazard includes standard errors when vcov is available", {
  fit <- hazard(
    time = c(1, 2),
    status = c(1, 0),
    x = matrix(c(1, 0), ncol = 1),
    theta = c(0.4),
    dist = "exponential"
  )
  fit$fit$vcov <- matrix(0.25, nrow = 1, ncol = 1)

  s <- summary(fit)

  expect_true(s$has_vcov)
  expect_equal(s$coefficients$std_error, 0.5, tolerance = 1e-12)
  expect_true(is.finite(s$coefficients$z_stat))
  expect_true(is.finite(s$coefficients$p_value))
})

test_that("optimizer produces valid vcov and SEs for all distributions", {
  skip_if_not_installed("numDeriv")

  set.seed(42)
  n <- 200
  time <- rexp(n, 0.5)
  cens <- runif(n, 0, quantile(time, 0.8))
  status <- as.integer(time <= cens)
  time <- pmin(time, cens)

  # Weibull
  fit_w <- hazard(time = time, status = status,
                  theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
  expect_true(is.matrix(fit_w$fit$vcov))
  expect_true(all(diag(fit_w$fit$vcov) > 0))
  expect_true(all(is.finite(fit_w$fit$se)))
  expect_equal(fit_w$fit$se, sqrt(diag(fit_w$fit$vcov)))

  # Exponential
  fit_e <- hazard(time = time, status = status,
                  theta = c(log(0.3)), dist = "exponential", fit = TRUE)
  expect_true(is.matrix(fit_e$fit$vcov))
  expect_true(all(diag(fit_e$fit$vcov) > 0))
  expect_true(all(is.finite(fit_e$fit$se)))

  # Log-logistic
  fit_ll <- hazard(time = time, status = status,
                   theta = c(0, 0.2), dist = "loglogistic", fit = TRUE)
  expect_true(is.matrix(fit_ll$fit$vcov))
  expect_true(all(diag(fit_ll$fit$vcov) > 0))
  expect_true(all(is.finite(fit_ll$fit$se)))

  # Log-normal
  fit_ln <- hazard(time = time, status = status,
                   theta = c(1, 0), dist = "lognormal", fit = TRUE)
  expect_true(is.matrix(fit_ln$fit$vcov))
  expect_true(all(diag(fit_ln$fit$vcov) > 0))
  expect_true(all(is.finite(fit_ln$fit$se)))
})

test_that("fit$se is conformable with fit$par whether or not vcov has NA rows", {
  skip_if_not_installed("numDeriv")

  # `$se` is derived from the fit's vcov. A multiphase fit legitimately carries
  # NA variance rows for parameters held fixed, and the guard used to collapse
  # the whole vector to a scalar NA as soon as any cell was NA -- so `$se` came
  # back length 1 against a length-5 `$par`, and naming it errored. This is the
  # `$se` half of the contract vcov.hazard() already keeps: a scalar NA means
  # there is no variance matrix at all, not that some of it is missing.

  set.seed(101)
  n <- 150
  df <- data.frame(time = rexp(n, 0.6) + 0.01,
                   status = rbinom(n, 1, 0.7))

  # Finite vcov -- one SE per parameter.
  fit_w <- hazard(time = df$time, status = df$status,
                  theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
  expect_false(anyNA(fit_w$fit$vcov))
  expect_length(fit_w$fit$se, length(fit_w$fit$par))
  expect_true(all(is.finite(fit_w$fit$se)))

  # vcov with NA rows (fixed shapes) -- still one SE per parameter.
  fit_mp <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = df, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1, maxit = 300, conserve = TRUE)))
  expect_true(is.matrix(fit_mp$fit$vcov))
  expect_true(anyNA(fit_mp$fit$vcov))
  expect_length(fit_mp$fit$se, length(fit_mp$fit$par))
  # Naming the SEs against the parameters is what the old shape broke.
  expect_named(stats::setNames(fit_mp$fit$se, names(fit_mp$fit$par)),
               names(fit_mp$fit$par))

  # Length and naming alone cannot fail on a vector of NA, which is exactly
  # what the first fix returned: one NA per parameter, conformable and empty.
  # Assert first that this fixture reaches both branches, then that the SEs
  # that were computable actually survive.
  d_mp <- diag(fit_mp$fit$vcov)
  computable <- is.finite(d_mp) & d_mp >= 0
  expect_true(any(computable))
  expect_true(any(!computable))
  expect_equal(unname(fit_mp$fit$se[computable]), unname(sqrt(d_mp[computable])))
  expect_true(all(is.na(fit_mp$fit$se[!computable])))

  # `$se` carries the variance matrix's names when it has them. No fitted path
  # currently produces a named `fit$fit$vcov`, so this guards against the two
  # diverging later rather than reproducing a fault seen today.
  expect_equal(names(fit_mp$fit$se), names(diag(fit_mp$fit$vcov)))

  # `$se` and summary() derive from the same matrix, so they cannot disagree.
  # They did: summary() reported real standard errors where `$se` was all NA.
  expect_equal(unname(fit_mp$fit$se),
               unname(summary(fit_mp)$coefficients$std_error))
})
test_that("vcov.hazard returns a named matrix and preserves NA rows", {
  # A multiphase fit legitimately has NA variance rows: parameters held fixed
  # (e.g. early shapes) and the Conservation-of-Events-conserved phase log_mu
  # carry no Hessian-based variance. The old contract returned a scalar NA for
  # the whole matrix whenever any cell was NA, making vcov() unusable for
  # multiphase models and discarding the finite free-parameter block.
  theta <- c(early.log_mu = -1.5, early.log_t_half = log(0.5),
             early.nu = 1, early.m = 1, early.x = 0.8,
             constant.log_mu = -3.0, constant.x = -0.4)
  V <- matrix(NA_real_, 7, 7)
  free <- c(1L, 5L, 7L)               # early.log_mu, early.x, constant.x
  V[free, free] <- matrix(c(0.0101, -0.0065, 0.0002,
                            -0.0065, 0.0072, -0.0004,
                            0.0002, -0.0004, 0.0012), 3, 3)
  fit <- structure(list(fit = list(theta = theta, vcov = V)), class = "hazard")

  Vout <- vcov(fit)
  expect_true(is.matrix(Vout))                       # not scalar NA
  expect_equal(dim(Vout), c(7L, 7L))
  expect_equal(rownames(Vout), names(theta))         # named for alignment
  expect_equal(colnames(Vout), names(theta))
  expect_true(anyNA(Vout))                           # NA rows preserved
  expect_true(all(is.na(diag(Vout)[c(2, 3, 4, 6)]))) # fixed/conserved -> NA
  expect_true(all(is.finite(diag(Vout)[free])))      # free block intact
})

test_that("vcov.hazard distinguishes a covariate shared across phases by name", {
  theta <- c(early.log_mu = -1.5, early.x = 0.8,
             constant.log_mu = -3.0, constant.x = -0.4)
  V <- diag(c(0.01, 0.02, NA_real_, 0.03))
  fit <- structure(list(fit = list(theta = theta, vcov = V)), class = "hazard")

  Vout <- vcov(fit)
  expect_setequal(rownames(Vout), c("early.log_mu", "early.x",
                                    "constant.log_mu", "constant.x"))
  # The two phase-specific x coefficients must resolve to distinct slots.
  expect_equal(Vout["early.x", "early.x"], 0.02)
  expect_equal(Vout["constant.x", "constant.x"], 0.03)
})

test_that("vcov.hazard returns scalar NA only when the matrix is truly absent", {
  fit_null <- structure(list(fit = list(theta = c(a = 1), vcov = NULL)),
                        class = "hazard")
  expect_true(is.na(vcov(fit_null)) && !is.matrix(vcov(fit_null)))

  fit_scalar <- structure(list(fit = list(theta = c(a = 1), vcov = NA)),
                          class = "hazard")
  expect_true(is.na(vcov(fit_scalar)) && !is.matrix(vcov(fit_scalar)))
})

test_that("vcov.hazard names single-phase (weibull) coefficients", {
  skip_if_not_installed("numDeriv")
  set.seed(7)
  n <- 200
  time <- rexp(n, 0.5)
  status <- rep(1L, n)
  fit <- hazard(time = time, status = status,
                theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
  Vout <- vcov(fit)
  expect_true(is.matrix(Vout))
  expect_false(is.null(rownames(Vout)))
  expect_equal(rownames(Vout), colnames(Vout))
})

test_that("summary shows finite z-stats and p-values from fitted model", {
  skip_if_not_installed("numDeriv")

  set.seed(99)
  n <- 150
  x <- matrix(rnorm(n), ncol = 1)
  time <- rexp(n, exp(0.5 * x))
  status <- rep(1L, n)

  fit <- hazard(time = time, status = status, x = x,
                theta = c(0.3, 1.0, 0), dist = "weibull", fit = TRUE)
  s <- summary(fit)

  expect_true(s$has_vcov)
  expect_true(all(is.finite(s$coefficients$std_error)))
  expect_true(all(is.finite(s$coefficients$z_stat)))
  expect_true(all(is.finite(s$coefficients$p_value)))
  expect_true(all(s$coefficients$p_value >= 0 & s$coefficients$p_value <= 1))
})

test_that("print.summary.hazard prints without error", {
  fit <- hazard(
    time = c(1, 2),
    status = c(1, 0),
    x = matrix(c(1, 0), ncol = 1),
    theta = c(0.4),
    dist = "exponential"
  )

  expect_output(print(summary(fit)), "hazard model summary")
})

test_that("hazard() accepts formula interface with Surv()", {
  df <- data.frame(
    time = c(1, 2, 3),
    status = c(1, 0, 1),
    x1 = c(0, 1, 0),
    x2 = c(1, 1, 0)
  )

  fit <- hazard(
    Surv(time, status) ~ x1 + x2,
    data = df,
    theta = c(0.3, -0.2),
    dist = "weibull"
  )

  expect_s3_class(fit, "hazard")
  expect_equal(fit$spec$dist, "weibull")
  expect_equal(length(fit$fit$theta), 2)
})

test_that("formula interface extracts predictors correctly", {
  df <- data.frame(
    time = c(1, 2, 3),
    status = c(1, 0, 1),
    x1 = c(0.5, 1.5, 2.5),
    x2 = c(1, 2, 3)
  )

  fit <- hazard(
    Surv(time, status) ~ x1 + x2,
    data = df,
    theta = c(-0.5, 0.3, -0.1),  # log(lambda), beta1, beta2 for exponential
    dist = "exponential"
  )

  pred <- predict(fit, type = "linear_predictor")
  expected_eta <- c(0.5 * 0.3 + 1 * (-0.1), 1.5 * 0.3 + 2 * (-0.1), 2.5 * 0.3 + 3 * (-0.1))
  expect_equal(pred, expected_eta, tolerance = 1e-12)
})

test_that("formula interface requires data when formula is provided", {
  expect_error(
    hazard(Surv(c(1, 2), c(1, 0)) ~ x, data = NULL),
    "'data' is required"
  )
})

test_that("formula interface works without intercept", {
  df <- data.frame(
    time = c(1, 2, 3),
    status = c(1, 0, 1),
    x = c(1, 0, 1)
  )

  fit <- hazard(
    Surv(time, status) ~ x - 1,
    data = df,
    theta = 0.4,
    dist = "exponential"
  )

  expect_equal(ncol(fit$data$x), 1)
  expect_equal(as.numeric(fit$data$x[, 1]), c(1, 0, 1))
})
