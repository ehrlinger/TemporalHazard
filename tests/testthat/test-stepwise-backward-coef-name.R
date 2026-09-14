# A single-distribution coefficient is tested by its position, whatever the
# user named `theta` (#304).
#
# .hzr_candidate_coef_name() names a single-distribution coefficient `beta<k>`
# by its design column. .hzr_wald_p() took the coefficient names from
# `names(theta)` whenever they were all non-empty, so a theta named after its
# covariates (`c(mu =, nu =, gb = 0, z = 0)`) offered `gb`, `z` and no `beta1`:
# the backward step stopped with "Unknown coefficient name(s)". Worse, a
# covariate named like a shape parameter (`nu`) would have matched the shape.
#
# Every test compares against a z computed by position from the fit's theta
# and vcov, and first checks that the competing coefficient's z is far away:
# matching one of two indistinguishable numbers would prove nothing.

bcn_data <- function(n = 300L, seed = 7L) {
  set.seed(seed)
  d <- data.frame(gb = stats::runif(n), z = stats::rnorm(n))
  d$time <- stats::rexp(n, 0.2) * exp(-2 * d$gb)
  d$status <- stats::rbinom(n, 1, 0.7)
  d
}

z_at <- function(fit, k) {
  unname(fit$fit$theta[k] / sqrt(fit$fit$vcov[k, k]))
}

bcn_fit <- function(formula, d, theta) {
  hazard(formula, data = d, dist = "weibull", theta = theta, fit = TRUE)
}

test_that("backward drop is the same whatever theta is named", {
  d <- bcn_data()
  f <- survival::Surv(time, status) ~ gb + z
  thetas <- list(
    covariate = c(mu = 0.2, nu = 1, gb = 0, z = 0),
    beta      = c(mu = 0.2, nu = 1, beta1 = 0, beta2 = 0),
    unnamed   = c(0.2, 1, 0, 0)
  )
  for (nm in names(thetas)) {
    fit <- bcn_fit(f, d, thetas[[nm]])
    z_z <- z_at(fit, 4L)
    # gb has a real effect and z none, so the two tests are far apart.
    expect_gt(abs(z_at(fit, 3L) - z_z), 5)

    sw <- hzr_stepwise(fit, data = d, direction = "backward",
                       criterion = "wald", slstay = 0.2, trace = FALSE)
    drops <- sw$steps[sw$steps$action == "drop", ]
    expect_identical(drops$variable, "z", label = nm)
    expect_equal(drops$stat, z_z, tolerance = 1e-6, label = nm)
    expect_equal(drops$p_value, 2 * stats::pnorm(-abs(z_z)),
                 tolerance = 1e-6, label = nm)
  }
})

test_that("a covariate named like a shape parameter is tested as the covariate", {
  d <- bcn_data()
  d$nu <- d$z
  fit <- bcn_fit(survival::Surv(time, status) ~ gb + nu, d,
                 c(mu = 0.2, nu = 1, gb = 0, nu = 0))
  expect_identical(names(stats::coef(fit)), c("mu", "nu", "gb", "nu"))
  z_cov <- z_at(fit, 4L)
  expect_gt(abs(z_at(fit, 2L) - z_cov), 5)

  step <- .hzr_stepwise_backward_step(fit, data = d, criterion = "wald",
                                      slstay = 0.2)
  row <- step$all_scores[step$all_scores$variable == "nu", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$stat, z_cov, tolerance = 1e-6)
  expect_true(step$accepted)
  expect_identical(step$variable, "nu")
})

test_that("forward Wald from a covariate-named theta tests the entered column", {
  d <- bcn_data()
  base <- bcn_fit(survival::Surv(time, status) ~ z, d,
                  c(mu = 0.2, nu = 1, z = 0))
  direct <- bcn_fit(survival::Surv(time, status) ~ z + gb, d,
                    c(mu = 0.2, nu = 1, z = 0, gb = 0))
  z_gb <- z_at(direct, 4L)
  expect_gt(abs(z_gb - z_at(direct, 3L)), 5)

  step <- .hzr_stepwise_forward_step(base, scope = ~ gb, data = d,
                                     criterion = "wald", slentry = 0.05)
  row <- step$all_scores[step$all_scores$variable == "gb", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$stat, z_gb, tolerance = 1e-3)
  expect_true(step$accepted)
})
