# A resample with no events at one level of a binary covariate has its
# maximum at infinity. The refit stops where the slope falls below the
# optimizer's tolerance and reports converged = TRUE with a small gradient,
# so it passed every convergence test, and its estimate was pooled (#654).

.mono_data <- withr::with_seed(654, {
  n <- 80
  x <- c(rep(1, 6), rep(0, n - 6))
  status <- stats::rbinom(n, 1, 0.6)
  status[x == 1] <- c(1, 0, 0, 0, 0, 0)
  data.frame(time = stats::rexp(n, 0.3), status = status, x = x)
})

.mono_fit <- function(d) {
  hazard(survival::Surv(time, status) ~ x, data = d, dist = "exponential",
         theta = c(0.3, 0), fit = TRUE)
}

test_that("the monotone-likelihood check separates a rising fit from a maximum", {
  d <- .mono_data
  expect_false(.hzr_monotone_likelihood(.mono_fit(d)))
  # Drop the one event at x = 1: the coefficient now runs to -Inf.
  sep <- d[!(d$x == 1 & d$status == 1), ]
  fit_sep <- .mono_fit(sep)
  expect_true(isTRUE(fit_sep$fit$converged))
  expect_true(.hzr_monotone_likelihood(fit_sep))
})

test_that("the check follows a separation along a combination of covariates", {
  # No events at (x1, x2) = (1, 0) and no rows at (0, 1): beta1 -> -Inf and
  # beta2 -> +Inf together, with beta1 + beta2 finite. A step in either
  # coefficient alone unbalances the (1, 1) rows, so only a ray along the
  # combination stays flat.
  withr::local_seed(6542)
  n <- 120
  grp <- rep(c("00", "10", "11"), each = n / 3)
  d <- data.frame(
    time = stats::rexp(n, 0.3),
    status = ifelse(grp == "10", 0L, stats::rbinom(n, 1, 0.6)),
    x1 = as.numeric(grp != "00"), x2 = as.numeric(grp == "11")
  )
  fit <- hazard(survival::Surv(time, status) ~ x1 + x2, data = d,
                dist = "exponential", theta = c(0.3, 0, 0), fit = TRUE)
  expect_true(isTRUE(fit$fit$converged))
  b <- unname(utils::tail(fit$fit$theta, 2))
  expect_true(b[1] < -5 && b[2] > 5)
  expect_true(.hzr_monotone_likelihood(fit))
})

# No events at the reference level: the intercept goes to -Inf and both
# other levels' coefficients to +Inf. No single coefficient's ray is flat.
.mono_reference <- withr::with_seed(6543, {
  n <- 120
  g <- factor(rep(c("a", "b", "c"), each = n / 3))
  data.frame(time = stats::rexp(n, 0.3),
             status = ifelse(g == "a", 0L, stats::rbinom(n, 1, 0.6)),
             g = g)
})

test_that("the check follows a separation through the intercept, with or without vcov", {
  # The fit warns, rightly, that it is weakly identified.
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ g, data = .mono_reference,
    dist = "exponential", theta = c(0.3, 0, 0), fit = TRUE))
  expect_true(isTRUE(fit$fit$converged))
  expect_lt(fit$fit$theta[[1]], -5)
  expect_true(.hzr_monotone_likelihood(fit))
  # At separation inverting the Hessian is what usually fails; the check
  # must still see it.
  no_vcov <- fit
  no_vcov$fit$vcov <- NA
  expect_true(.hzr_monotone_likelihood(no_vcov))
  # The known negative takes the same path.
  ok <- .mono_fit(.mono_data)
  ok$fit$vcov <- NA
  expect_false(.hzr_monotone_likelihood(ok))
})

test_that("the Weibull intercept is followed on the log scale", {
  # Its mu is a positive rate, so the intercept's run to -Inf is mu -> 0,
  # which no linear step in mu can follow.
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ g, data = .mono_reference,
    dist = "weibull", theta = c(0.3, 1, 0, 0), fit = TRUE))
  expect_true(isTRUE(fit$fit$converged))
  expect_lt(fit$fit$theta[[1]], 1e-8)
  expect_true(.hzr_monotone_likelihood(fit))
  fit$fit$vcov <- NA
  expect_true(.hzr_monotone_likelihood(fit))
})

test_that("one extreme covariate value does not read as a flat likelihood", {
  # A row at x = 1e7 sets the step that moves the linear predictor by 10 at
  # 1e-6, where a genuine maximum falls by about 4e-5. The step is at least
  # three standard errors, where it falls by about 4.5.
  d <- withr::with_seed(2, {
    x <- stats::rnorm(200)
    data.frame(time = stats::rexp(200, 0.3 * exp(0.3 * x)),
               status = stats::rbinom(200, 1, 0.7), x = x)
  })
  d$x[1] <- 1e7
  d$status[1] <- 0L
  d$time[1] <- 5
  fit <- hazard(survival::Surv(time, status) ~ x, data = d,
                dist = "exponential", theta = c(-1, 0), fit = TRUE)
  expect_true(isTRUE(fit$fit$converged))
  expect_false(.hzr_monotone_likelihood(fit))
  # Without vcov the curvature comes from finite differences, whose step in
  # x must not overflow exp() on that row.
  fit$fit$vcov <- NA
  expect_false(.hzr_monotone_likelihood(fit))
})

test_that("the check sees an exactly flat direction in the design", {
  # A factor level absent from a resample leaves a column of zeros, whose
  # coefficient stays where it started; an exact multiple of a column is
  # flat along a direction rays find only to rounding.
  d <- withr::with_seed(5, {
    n <- 100
    g <- factor(c(rep("rare", 3), rep(c("a", "b"), length.out = n - 3)),
                levels = c("a", "b", "rare"))
    data.frame(time = stats::rexp(n, 0.3), status = stats::rbinom(n, 1, 0.7),
               g = g)
  })[-(1:3), ]
  stopifnot(!any(d$g == "rare"))
  for (start in c(0, 0.5)) {
    fit <- suppressWarnings(hazard(
      survival::Surv(time, status) ~ g, data = d, dist = "exponential",
      theta = c(-1, 0, start), fit = TRUE))
    expect_identical(utils::tail(unname(fit$fit$theta), 1), start)
    expect_true(.hzr_monotone_likelihood(fit))
    fit$fit$vcov <- NA
    expect_true(.hzr_monotone_likelihood(fit))
  }
  d <- withr::with_seed(3, {
    x1 <- stats::rnorm(200)
    data.frame(time = stats::rexp(200, 0.3 * exp(0.3 * x1)),
               status = stats::rbinom(200, 1, 0.7), x1 = x1, x2 = 2 * x1)
  })
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ x1 + x2, data = d, dist = "exponential",
    theta = c(-1, 0, 0), fit = TRUE))
  expect_true(.hzr_monotone_likelihood(fit))
  fit$fit$vcov <- NA
  expect_true(.hzr_monotone_likelihood(fit))
  # Nearly collinear but determined is not flat.
  d$x2 <- 2 * d$x1 + withr::with_seed(4, stats::rnorm(200, sd = 0.05))
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ x1 + x2, data = d, dist = "exponential",
    theta = c(-1, 0, 0), fit = TRUE))
  expect_false(.hzr_monotone_likelihood(fit))
})

test_that("an offset covariate does not read as a flat likelihood", {
  # Shifted by 20000, the covariate's move along the combined direction
  # nearly cancels the intercept's. Sized on the covariate alone the step
  # shrank to nothing, and every replicate of a genuine maximum was dropped.
  x <- rep(c(-1, 0, 1), 100)
  for (shift in c(0, 20000)) {
    d <- data.frame(time = exp(-0.001 * x) / 0.3, status = 1, x = x + shift)
    fit <- suppressWarnings(hazard(
      survival::Surv(time, status) ~ x, data = d, dist = "exponential",
      theta = c(log(0.3) - 0.001 * shift, 0.001), fit = TRUE))
    expect_false(.hzr_monotone_likelihood(fit), label = paste("shift", shift))
    fit$fit$vcov <- NA
    expect_false(.hzr_monotone_likelihood(fit), label = paste("shift", shift))
  }
  b <- suppressWarnings(hzr_bootstrap(fit, n_boot = 3, seed = 1))
  expect_identical(b$n_success, 3L)
})

test_that("rescaling every weight leaves the verdict as it was", {
  # Weights scale the log-likelihood and every fall along a ray, but not
  # the maximum. At weight 5e-8 a genuine maximum's fall over the capped
  # step was below an absolute tolerance, and every replicate was dropped.
  d <- data.frame(time = 1:20, status = 1, w = 5e-8)
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = d, dist = "exponential",
    weights = w, theta = -log(10.5), fit = TRUE))
  expect_false(.hzr_monotone_likelihood(fit))
  b <- suppressWarnings(hzr_bootstrap(fit, n_boot = 3, seed = 1))
  expect_identical(b$n_success, 3L)
  # The same estimates under rescaled weights: a separation stays flat and
  # a genuine maximum stays a maximum. (Fitted at weight 5e-8 the optimizer
  # barely leaves its start, so the rescaling is applied to fitted objects.)
  sep <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ g, data = .mono_reference,
    dist = "exponential", theta = c(0.3, 0, 0), fit = TRUE))
  ok <- .mono_fit(.mono_data)
  for (s in c(5e-8, 1e6)) {
    sep$data$weights <- rep(s, length(sep$data$time))
    ok$data$weights <- rep(s, length(ok$data$time))
    expect_true(.hzr_monotone_likelihood(sep), label = paste("separated", s))
    expect_false(.hzr_monotone_likelihood(ok), label = paste("maximum", s))
  }
})

test_that("the check sees a multiphase phase with no events of its own", {
  skip_on_cran() # multiphase fit
  # With no deaths after six months the constant phase is not needed: its
  # log_mu runs to -Inf on its own, and only its own ray is flat.
  e <- new.env()
  utils::data("avc", package = "TemporalHazard", envir = e)
  d <- stats::na.omit(e$avc[, c("int_dead", "dead")])
  d$dead[d$int_dead > 6] <- 0L
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = phases, fit = TRUE,
    control = list(n_starts = 1L, conserve = FALSE)))
  expect_true(isTRUE(fit$fit$converged))
  expect_lt(utils::tail(fit$fit$theta, 1), -15)
  expect_true(.hzr_monotone_likelihood(fit))
  # And without vcov, where the step is sized on the linear predictor that
  # log_mu moves.
  fit$fit$vcov <- NA
  expect_true(.hzr_monotone_likelihood(fit))
})

test_that("the rank check sees a column collinear with a conserved log_mu", {
  skip_on_cran() # multiphase fits
  # Conservation of Events re-solves the constant phase's log_mu, and
  # fixed_mask marks it; a constant covariate there is not identified, and
  # two starts end at different values with the same likelihood.
  e <- new.env()
  utils::data("avc", package = "TemporalHazard", envir = e)
  d <- stats::na.omit(e$avc[, c("int_dead", "dead")])
  d$one <- 1
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m",
                      formula = ~ 1),
    constant = hzr_phase("constant", formula = ~ one)
  )
  fits <- lapply(c(0, 0.7), function(st) {
    suppressWarnings(hazard(
      survival::Surv(int_dead, dead) ~ one, data = d, dist = "multiphase",
      phases = phases, fit = TRUE, theta = c(-2.7, -1.9, 1.4, 1, -4, st),
      control = list(n_starts = 1L)))
  })
  b <- vapply(fits, function(f) utils::tail(unname(f$fit$theta), 1), 0)
  expect_gt(abs(diff(b)), 1)
  for (f in fits) {
    # The case is the conserved log_mu (5th) left masked with no finite
    # vcov, as CoE leaves it when the vcov recompute fails. Whether it fails
    # varies by platform (it succeeds on R-devel), so the state is set here.
    f$fit$fixed_mask[[5]] <- TRUE
    f$fit$vcov <- NA
    expect_true(.hzr_monotone_likelihood(f))
  }
})

test_that("an error in the check fails its replicate, not the run", {
  fit <- .mono_fit(.mono_data)
  calls <- 0L
  local_mocked_bindings(.hzr_monotone_likelihood = function(...) {
    calls <<- calls + 1L
    if (calls == 1L) FALSE else stop("probe")  # the first call is the original
  })
  b <- suppressWarnings(hzr_bootstrap(fit, n_boot = 5, seed = 1))
  expect_identical(b$n_success, 0L)
  expect_identical(b$failure_reasons[["likelihood check failed: probe"]], 5L)
})

test_that("hzr_bootstrap() drops replicates on a monotone likelihood (#654)", {
  fit <- .mono_fit(.mono_data)
  reason <- "likelihood flat or still rising (no unique finite maximum)"
  expect_warning(
    b <- hzr_bootstrap(fit, n_boot = 40, seed = 1),
    "replicates were dropped because their likelihood was flat or still"
  )
  n_mono <- b$failure_reasons[[reason]]
  expect_gt(n_mono, 0L)
  expect_equal(b$n_success + b$n_failed, 40L)

  # The same run with the check off pools every replicate, as before.
  local_mocked_bindings(.hzr_monotone_likelihood = function(...) FALSE)
  pooled <- hzr_bootstrap(fit, n_boot = 40, seed = 1)
  expect_equal(pooled$n_success, 40L)
  bx <- function(r) r$estimate[r$parameter == "x"]
  # The dropped replicates are exactly the ones that ran out to infinity,
  # and every one kept is unchanged.
  dropped <- setdiff(pooled$replicates$replicate, b$replicates$replicate)
  expect_length(dropped, n_mono)
  far <- pooled$replicates$parameter == "x" & pooled$replicates$estimate < -5
  expect_setequal(dropped, pooled$replicates$replicate[far])
  kept <- pooled$replicates[pooled$replicates$replicate %in%
                              b$replicates$replicate, ]
  expect_equal(kept$estimate, b$replicates$estimate)
  sd_of <- function(bt) bt$summary$sd[bt$summary$parameter == "x"]
  expect_gt(sd_of(pooled), 5 * sd_of(b))
  expect_lt(max(abs(bx(b$replicates))), 5)
})

test_that("hzr_bootstrap() warns when the fit itself is on a monotone likelihood", {
  d <- .mono_data
  fit_sep <- .mono_fit(d[!(d$x == 1 & d$status == 1), ])
  w <- character()
  withCallingHandlers(
    hzr_bootstrap(fit_sep, n_boot = 3, seed = 1),
    warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("fit being bootstrapped is flat or still rising", w)))
})
