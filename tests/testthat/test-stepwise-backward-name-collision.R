# A dropped variable is tested on its own coefficient, found by its term (#315).
#
# The backward step looked a variable's coefficient up by its bare name.
# model.matrix() names a logical `flag`'s column `flagTRUE`, while a factor
# `fla` with level `g` owns a column named `flag`, so in `~ fla + flag` the
# name `flag` found the factor's dummy. It was loud only because `fla` then
# errored ("expands to multiple coefficients") before any decision was made;
# resolving `fla` alone would have made `flag`'s wrong lookup silent.
#
# Every test compares each variable's drop statistic against a z computed by
# position from the fit's theta and vcov, and first checks that the two
# columns' z are far apart: matching one of two indistinguishable numbers
# would prove nothing.

bnc_data <- function(n = 400L, seed = 11L) {
  set.seed(seed)
  d <- data.frame(
    fla  = factor(sample(c("f", "g"), n, replace = TRUE)),
    flag = sample(c(TRUE, FALSE), n, replace = TRUE)
  )
  d$time <- stats::rexp(n) * exp(-1.2 * d$flag)
  d$status <- 1L
  d
}

# By position in theta, never by name.
bnc_z <- function(fit, k) {
  unname(fit$fit$theta[k] / sqrt(fit$fit$vcov[k, k]))
}

bnc_stat <- function(step, var) {
  row <- step$all_scores[step$all_scores$variable == var, ]
  expect_equal(nrow(row), 1L)
  row$stat
}

bnc_single <- function(rhs, d) {
  hazard(stats::update(survival::Surv(time, status) ~ 1, rhs), data = d,
         dist = "weibull", theta = c(0.5, 1, 0, 0), fit = TRUE)
}

bnc_multi <- function(rhs, d) {
  hazard(survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
         phases = list(constant = hzr_phase("constant", formula = rhs)),
         fit = TRUE)
}

test_that("each variable is tested on its own column, either order (single dist)", {
  d <- bnc_data()
  # Theta position of each variable's column: shapes mu, nu come first.
  orders <- list(
    list(rhs = ~ fla + flag, cols = c("flag", "flagTRUE"), fla = 3L, flag = 4L),
    list(rhs = ~ flag + fla, cols = c("flagTRUE", "flag"), fla = 4L, flag = 3L)
  )
  for (o in orders) {
    fit <- bnc_single(o$rhs, d)
    # The collision: the factor's dummy carries the logical's bare name.
    expect_identical(colnames(fit$data$x), o$cols)
    z_fla <- bnc_z(fit, o$fla)
    z_flag <- bnc_z(fit, o$flag)
    expect_gt(abs(z_flag - z_fla), 5)

    step <- .hzr_stepwise_backward_step(fit, data = d, criterion = "wald",
                                        slstay = 0.1)
    expect_equal(bnc_stat(step, "flag"), z_flag, tolerance = 1e-8)
    expect_equal(bnc_stat(step, "fla"), z_fla, tolerance = 1e-8)

    # fla (p about 0.16) is the one to go; flag (p below 1e-29) stays.
    expect_true(step$accepted)
    expect_identical(step$variable, "fla")
    expect_equal(step$stat, z_fla, tolerance = 1e-8)
    expect_equal(step$p_value, 2 * stats::pnorm(-abs(z_fla)), tolerance = 1e-8)

    # The post-drop refit is the model without fla, with flag's own column.
    direct <- hazard(survival::Surv(time, status) ~ flag, data = d,
                     dist = "weibull", theta = c(0.5, 1, 0), fit = TRUE)
    expect_identical(colnames(step$fit$data$x), "flagTRUE")
    expect_equal(step$fit$fit$theta, direct$fit$theta, tolerance = 1e-4)
  }
})

test_that("each variable is tested on its own column, either order (multiphase)", {
  d <- bnc_data()
  for (rhs in list(~ fla + flag, ~ flag + fla)) {
    fit <- bnc_multi(rhs, d)
    nms <- names(stats::coef(fit))
    # Position of each column in theta, located by the name model.matrix()
    # gave it, not the variable's.
    k_fla <- match("constant.flag", nms)
    k_flag <- match("constant.flagTRUE", nms)
    expect_false(anyNA(c(k_fla, k_flag)))
    z_fla <- bnc_z(fit, k_fla)
    z_flag <- bnc_z(fit, k_flag)
    expect_gt(abs(z_flag - z_fla), 5)

    step <- .hzr_stepwise_backward_step(fit, data = d, criterion = "wald",
                                        slstay = 0.1)
    expect_equal(bnc_stat(step, "flag"), z_flag, tolerance = 1e-8)
    expect_equal(bnc_stat(step, "fla"), z_fla, tolerance = 1e-8)
    expect_true(step$accepted)
    expect_identical(step$variable, "fla")
    expect_equal(step$stat, z_fla, tolerance = 1e-8)

    direct <- bnc_multi(~ flag, d)
    expect_identical(names(stats::coef(step$fit)),
                     c("constant.log_mu", "constant.flagTRUE"))
    expect_equal(unname(stats::coef(step$fit)), unname(stats::coef(direct)),
                 tolerance = 1e-4)
  }
})

test_that("a phase inheriting the global design is resolved by its term", {
  # A phase with no formula of its own reads its columns from the global
  # design (fit$data$x_design), not from a phase design.
  d <- bnc_data()
  fit <- hazard(survival::Surv(time, status) ~ fla + flag, data = d,
                dist = "multiphase",
                phases = list(constant = hzr_phase("constant")), fit = TRUE)
  expect_true(.hzr_phase_inherits_global(fit, "constant"))
  expect_identical(.hzr_candidate_coef_name(fit, "flag", "constant"),
                   "constant.flagTRUE")
  expect_identical(.hzr_candidate_coef_name(fit, "fla", "constant"),
                   "constant.flag")
})

test_that("a phase formula beside a global design reads its own design", {
  # The phase has a term (`z`) the global formula lacks, so reading the
  # global design for it cannot find its columns; and the two designs order
  # the colliding columns differently.
  d <- bnc_data()
  d$z <- stats::rnorm(nrow(d))
  # Only name resolution is tested here. An exponential sample under a cdf
  # plus constant model has a singular Hessian, and hazard() says so.
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ fla + flag, data = d,
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0,
                        formula = ~ z + flag + fla),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  expect_false(.hzr_phase_inherits_global(fit, "early"))
  expect_true(.hzr_phase_inherits_global(fit, "constant"))
  expect_identical(colnames(fit$fit$x_list$early),
                   c("z", "flagTRUE", "flag"))
  expect_identical(colnames(fit$fit$x_list$constant), c("flag", "flagTRUE"))
  expect_identical(.hzr_candidate_coef_name(fit, "z", "early"), "early.z")
  expect_identical(.hzr_candidate_coef_name(fit, "flag", "early"),
                   "early.flagTRUE")
  expect_identical(.hzr_candidate_coef_name(fit, "fla", "early"),
                   "early.flag")
  expect_identical(.hzr_candidate_coef_name(fit, "flag", "constant"),
                   "constant.flagTRUE")
  expect_identical(.hzr_candidate_coef_name(fit, "fla", "constant"),
                   "constant.flag")
})

test_that("the whole backward run drops fla and keeps flag", {
  d <- bnc_data()
  fit <- bnc_single(~ fla + flag, d)
  sw <- hzr_stepwise(fit, data = d, direction = "backward",
                     criterion = "wald", slstay = 0.1, trace = FALSE)
  drops <- sw$steps[sw$steps$action == "drop", ]
  expect_identical(drops$variable, "fla")
  expect_equal(drops$stat, bnc_z(fit, 3L), tolerance = 1e-8)
  # hzr_stepwise() returns the final fit itself.
  expect_identical(colnames(sw$data$x), "flagTRUE")
})

test_that("a logical with no collision is found by its term", {
  # Its column is `flagTRUE`, so the bare name matched nothing and the drop
  # test stopped with "not found in the design matrix".
  d <- bnc_data()
  d$z <- stats::rnorm(nrow(d))
  fit <- hazard(survival::Surv(time, status) ~ z + flag, data = d,
                dist = "weibull", theta = c(0.5, 1, 0, 0), fit = TRUE)
  expect_identical(colnames(fit$data$x), c("z", "flagTRUE"))
  step <- .hzr_stepwise_backward_step(fit, data = d, criterion = "wald",
                                      slstay = 0.1)
  expect_equal(bnc_stat(step, "flag"), bnc_z(fit, 4L), tolerance = 1e-8)
  expect_equal(bnc_stat(step, "z"), bnc_z(fit, 3L), tolerance = 1e-8)
})

test_that("plain numeric terms resolve to the positions the name lookup gave", {
  # No collision: the term route must name exactly what the name lookup did.
  set.seed(3L)
  n <- 300L
  d <- data.frame(a = stats::rnorm(n), b = stats::rnorm(n))
  d$time <- stats::rexp(n) * exp(-0.8 * d$a)
  d$status <- 1L
  fit <- hazard(survival::Surv(time, status) ~ a + b, data = d,
                dist = "weibull", theta = c(0.5, 1, 0, 0), fit = TRUE)
  expect_identical(.hzr_candidate_coef_name(fit, "a", NULL), "beta1")
  expect_identical(.hzr_candidate_coef_name(fit, "b", NULL), "beta2")
  step <- .hzr_stepwise_backward_step(fit, data = d, criterion = "wald",
                                      slstay = 0.2)
  expect_equal(bnc_stat(step, "a"), bnc_z(fit, 3L), tolerance = 1e-8)
  expect_equal(bnc_stat(step, "b"), bnc_z(fit, 4L), tolerance = 1e-8)
})

test_that("a multi-column term still errors in the drop path", {
  set.seed(223L)
  n <- 200L
  d <- data.frame(time = stats::rexp(n), status = 1L, x1 = stats::rnorm(n),
                  f = factor(sample(letters[1:3], n, replace = TRUE)))
  fit <- hazard(survival::Surv(time, status) ~ x1 + f, data = d,
                dist = "weibull", theta = c(0.5, 1, 0, 0, 0), fit = TRUE)
  expect_error(
    .hzr_stepwise_backward_step(fit, data = d, criterion = "wald"),
    "expands to multiple coefficients \\(.fb., .fc.\\)"
  )
})

test_that("a fit with no stored design keeps the name lookup", {
  # The vector interface stores no formula design, so the name is all there is.
  d <- bnc_data()
  x <- cbind(a = stats::rnorm(nrow(d)), b = as.numeric(d$flag))
  fit <- hazard(time = d$time, status = d$status, x = x, dist = "weibull",
                theta = c(0.5, 1, 0, 0), fit = TRUE)
  expect_null(fit$data$x_design)
  expect_null(.hzr_term_coef_name(fit, "b", NULL, d))
  expect_identical(.hzr_candidate_coef_name(fit, "b", NULL), "beta2")
})
