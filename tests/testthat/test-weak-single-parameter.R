# test-weak-single-parameter.R
# The weak-direction detector also names ONE parameter the data do not
# determine (#415). The pairwise reading cannot: a correlation matrix
# normalises that parameter's variance away.

nm <- c("late.log_mu", "late.log_tau", "late.gamma", "late.alpha")
th <- c(-0.5, 1, 4, 0.8)
gated <- 1e-12  # below .hzr_rcond_tol, as on every fit this reading targets
shapes <- c("late.gamma", "late.alpha", "late.eta")

# gamma's standard error, relative to gamma, is `rel`; everything else is
# determined and uncorrelated with it.
single_vcov <- function(rel, gamma = th[3]) {
  diag(c(0.01, 0.02, (rel * gamma)^2, 0.001))
}

test_that("a parameter the data do not determine is named on its own", {
  w <- .hzr_weak_direction(single_vcov(50), gated, nm, theta = th,
                           shape_names = shapes)
  expect_true(is.list(w))
  expect_true(isTRUE(w$single))
  expect_identical(w$params, "late.gamma")
  expect_equal(w$se_metric, 50, tolerance = 1e-8)
  expect_match(.hzr_weak_direction_message(w), "near-flat direction on its own",
               fixed = TRUE)
})

test_that("without the fitted values the reading does not run", {
  # Every pre-#415 caller passes no theta, so its answer is unchanged.
  expect_null(.hzr_weak_direction(single_vcov(50), gated, nm))
})

test_that("a parameter's own precision is not the criterion; the gate is", {
  # The same dominant direction, relative standard error 0.3: at gamma-hat =
  # 1e7 that is what a flat-to-infinity likelihood reports, so it is named
  # when the Hessian is ill-conditioned ...
  w <- .hzr_weak_direction(single_vcov(0.3, gamma = 1e7), gated, nm,
                           theta = replace(th, 3, 1e7),
                           shape_names = shapes)
  expect_identical(w$params, "late.gamma")
  expect_match(.hzr_weak_direction_message(w), "'late.gamma' = 1e+07",
               fixed = TRUE)
  # ... and not when it is not.
  expect_null(.hzr_weak_direction(single_vcov(0.3), 1e-4, nm, theta = th,
                           shape_names = shapes))
})

test_that("no single parameter dominating means nothing is named", {
  # The flattest direction is shared by log_mu and log_tau (loading 0.707
  # each, correlation 0.5, below the pairwise reading's bar too).
  v <- diag(c(1, 1, (0.1 * th[3])^2, 0.001))
  v[1, 2] <- v[2, 1] <- 0.5
  expect_null(.hzr_weak_direction(v, gated, nm, theta = th,
                           shape_names = shapes))
})

test_that("the rcond gate still decides whether anything is read", {
  expect_null(.hzr_weak_direction(single_vcov(50), 1e-4, nm, theta = th,
                           shape_names = shapes))
})

test_that("a log-scale estimate near 0 is not scaled by its own value", {
  # gamma is the flat one (relative standard error 50). log_tau = 1e-4 with a
  # standard error of 0.5 is determined to within a factor of e, but scaling
  # it by its own value, as #415's candidate did, reads a relative standard
  # error of 5000 and names log_tau instead.
  th0 <- th
  th0[2] <- 1e-4
  v <- diag(c(0.01, 0.25, (50 * th[3])^2, 0.001))
  w <- .hzr_weak_direction(v, gated, nm, theta = th0,
                           shape_names = shapes)
  expect_identical(w$params, "late.gamma")
  value_scaled <- v / outer(abs(th0), abs(th0))
  expect_identical(nm[which.max(abs(eigen(value_scaled)$vectors[, 1]))],
                   "late.log_tau")
})

test_that("a pair the correlation reading finds is still reported as the pair", {
  v <- matrix(c(1, 0.999, 0.999, 1), 2) * 100
  w <- .hzr_weak_direction(v, gated, c("early.m", "early.nu"),
                           theta = c(27, 0.027),
                           shape_names = shapes)
  expect_false(isTRUE(w$single))
  expect_setequal(w$params, c("early.m", "early.nu"))
})

test_that("a coefficient is never named on its own; a log-scale parameter only when listed", {
  # A covariate in units of 1e-4 opens the rcond gate by scaling alone, and
  # its well-identified coefficient dominated the flattest direction.
  logs <- c("late.log_mu", "late.log_tau")
  v <- diag(c(0.01, 0.02, (0.1 * th[3])^2, 0.001, 565^2))
  w <- .hzr_weak_direction(v, gated, c(nm, "late.x2"),
                           theta = c(th, 2641),
                           shape_names = shapes, log_names = logs)
  expect_null(w)
  # A log_mu with a standard error of 100 is named once hazard() lists it as
  # log-scale (#673), and not by a caller that does not.
  v2 <- diag(c(1e4, 0.02, (0.1 * th[3])^2, 0.001))
  expect_null(.hzr_weak_direction(v2, gated, nm, theta = th,
                           shape_names = shapes))
  w2 <- .hzr_weak_direction(v2, gated, nm, theta = th, shape_names = shapes,
                            log_names = logs)
  expect_identical(w2$params, "late.log_mu")
  expect_true(isTRUE(w2$single))
  # The same dominant direction with a log-scale standard error of 5, under
  # the bar of 10: a determined log_mu, so nothing is named.
  v3 <- diag(c(25, 1e-6, (1e-3 * th[3])^2, 1e-6))
  expect_null(.hzr_weak_direction(v3, gated, nm, theta = th,
                                  shape_names = shapes, log_names = logs))
})

# A constant-phase fit whose covariate `x` is z in units of `unit`, so a
# large unit opens the rcond gate by scaling alone (#673 review).
unit_fit <- function(unit = 1e4, covariate = "x") {
  withr::local_seed(2)
  n <- 400
  z <- stats::rnorm(n)
  te <- stats::rexp(n, 0.3 * exp(0.4 * z))
  ce <- stats::runif(n, 0, 8)
  d <- data.frame(time = pmin(te, ce), status = as.integer(te <= ce))
  d[[covariate]] <- unit * z
  suppressWarnings(hazard(
    stats::as.formula(paste("survival::Surv(time, status) ~", covariate)),
    data = d, dist = "multiphase",
    phases = list(constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1)))
}

test_that("a determined log_mu is not named when covariate units open the gate (#673)", {
  fit <- unit_fit(1e4)
  # The premise: the gate is open, and log_mu is well determined.
  expect_lt(fit$fit$rcond, .hzr_rcond_tol)
  expect_lt(sqrt(diag(fit$fit$vcov))[[1]], 1)
  expect_null(fit$fit$weak)
})

test_that("a covariate named like a log-scale parameter is never named (#673 review)", {
  # A column called `log_mu` makes a second `constant.log_mu`, the
  # coefficient (z about 7). The name cannot say which parameter it is.
  fit <- unit_fit(1e-4, covariate = "log_mu")
  expect_identical(names(fit$fit$theta), rep("constant.log_mu", 2L))
  expect_lt(fit$fit$rcond, .hzr_rcond_tol)
  expect_null(fit$fit$weak)
})

test_that("a conserved phase switched off at the boundary is not identified (#673)", {
  # Nobody dies late, so the constant phase has nothing to fit: its log_mu
  # runs to about -28 with a standard error near 3.5e4, and the Hessian's
  # rcond is about 1e-11. Before #673 the single-parameter reading looked
  # only at g3 shapes and fit$fit$weak was NULL, "checked, identified".
  withr::local_seed(1)
  n <- 300
  te <- stats::rexp(n, 2)
  dies <- stats::runif(n) < 0.25
  d <- data.frame(time = ifelse(dies, pmin(te, 1.5), stats::runif(n, 2, 10)),
                  status = as.integer(dies))
  fit <- NULL
  w <- testthat::capture_warnings(fit <- hazard(
    survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 2, m = 1,
                                    fixed = c("nu", "m")),
                  constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1)))
  # The premise: the gate is open and log_mu is undetermined.
  expect_lt(fit$fit$rcond, .hzr_rcond_tol)
  expect_lt(fit$fit$theta[["constant.log_mu"]], -20)
  j <- match("constant.log_mu", names(fit$fit$theta))
  expect_gt(sqrt(diag(fit$fit$vcov))[[j]], 1e3)
  expect_identical(fit$fit$weak$params, "constant.log_mu")
  expect_true(isTRUE(fit$fit$weak$single))
  expect_true(any(grepl("'constant.log_mu'", w, fixed = TRUE)))
})

test_that("nu is not a named shape: it is signed and 0 is legitimate", {
  v <- diag(c(0.01, 0.05^2))
  expect_null(.hzr_weak_direction(v, gated, c("early.log_t_half", "early.nu"),
                                  theta = c(0, 1e-3),
                           shape_names = shapes))
})

test_that("an identified fit with a small-unit covariate names nothing", {
  skip_if_not_installed("numDeriv")
  withr::local_seed(1)
  n <- 400
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n) * 1e-4
  ev <- (stats::rexp(n) / exp(0.5 * x1 + 3000 * x2))^(1 / 1.5)
  cen <- stats::runif(n, 0, 3)
  d <- data.frame(time = pmin(ev, cen), status = as.numeric(ev <= cen),
                  x1 = x1, x2 = x2)
  fit <- suppressWarnings(hazard(survival::Surv(time, status) ~ x1 + x2,
                                 data = d, dist = "weibull",
                                 theta = c(0.5, 1, 0, 0), fit = TRUE))
  # The premise: the gate is open, so only the naming rule keeps this quiet.
  expect_lt(fit$fit$rcond, .hzr_rcond_tol)
  expect_null(fit$fit$weak)
})

test_that("a covariate named like a shape is not a shape", {
  # A cdf phase with a covariate called `gamma` has a coefficient
  # `early.gamma`; only the phase specs say which names are g3 shapes.
  v <- diag(c(0.01, (50 * 4)^2))
  expect_null(.hzr_weak_direction(v, gated, c("early.log_mu", "early.gamma"),
                                  theta = c(-1, 4), shape_names = shapes))
  expect_identical(
    .hzr_weak_direction(v, gated, c("early.log_mu", "early.gamma"),
                        theta = c(-1, 4), shape_names = "early.gamma")$params,
    "early.gamma")
})
