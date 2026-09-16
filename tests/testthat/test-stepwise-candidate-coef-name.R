# The entered candidate is tested on the column the refit added.
#
# .hzr_candidate_coef_name() used to find a candidate's coefficient by its bare
# variable name. model.matrix() names a logical's column `<var>TRUE`, so when a
# factor already in the model has a dummy that carries the bare name -- factor
# `fla` with level `g` gives the column `flag` -- the lookup found the factor's
# dummy instead of the logical `flag`. The Wald candidate loop and the score
# path's Wald fallback then reported the dummy's test as the candidate's, a
# populated p-value over the wrong coefficient.
#
# Every test here compares against a z computed from a direct fit's coef() and
# vcov(), and first checks that the candidate's z and the dummy's z are far
# apart: matching one of two indistinguishable numbers would prove nothing.

ccn_data <- function(effect = 1.2, n = 400L, seed = 11L) {
  set.seed(seed)
  d <- data.frame(
    fla  = factor(sample(c("f", "g"), n, replace = TRUE)),
    flag = sample(c(TRUE, FALSE), n, replace = TRUE)
  )
  d$time <- stats::rexp(n) * exp(-effect * d$flag)
  d$status <- 1L
  d
}

# By position: a single-distribution coef() carries no names, only vcov() does.
wald_z <- function(fit, nm) {
  i <- match(nm, rownames(stats::vcov(fit)))
  unname(stats::coef(fit)[i] / sqrt(stats::vcov(fit)[i, i]))
}

ccn_phases <- function(formula) {
  list(constant = hzr_phase("constant", formula = formula))
}

scored_row <- function(step, var) {
  row <- step$all_scores[step$all_scores$variable == var, ]
  expect_equal(nrow(row), 1L)
  row
}

test_that("a logical candidate is tested on its own column (single dist)", {
  d <- ccn_data()
  base <- hazard(survival::Surv(time, status) ~ fla, data = d,
                 dist = "weibull", theta = c(0.5, 1, 0), fit = TRUE)
  direct <- hazard(survival::Surv(time, status) ~ fla + flag, data = d,
                   dist = "weibull", theta = c(0.5, 1, 0, 0), fit = TRUE)

  # The collision: the factor's dummy carries the logical's bare name.
  expect_identical(colnames(direct$data$x), c("flag", "flagTRUE"))
  z_cand <- wald_z(direct, "beta2")
  expect_gt(abs(z_cand - wald_z(direct, "beta1")), 5)

  step <- .hzr_stepwise_forward_step(base, scope = ~ flag, data = d,
                                     criterion = "wald", slentry = 0.05)
  row <- scored_row(step, "flag")
  expect_equal(row$stat, z_cand, tolerance = 1e-3)
  expect_lt(row$p_value, 1e-10)
  expect_true(step$accepted)
})

test_that("a logical candidate enters a single-dist fit with no collision", {
  # With no factor in the model the bare name matched nothing, and the lookup
  # stopped with "not found in the design matrix" on a column the refit had.
  d <- ccn_data()
  base <- hazard(survival::Surv(time, status) ~ 1, data = d,
                 dist = "weibull", theta = c(0.5, 1), fit = TRUE)
  direct <- hazard(survival::Surv(time, status) ~ flag, data = d,
                   dist = "weibull", theta = c(0.5, 1, 0), fit = TRUE)
  expect_identical(colnames(direct$data$x), "flagTRUE")

  step <- .hzr_stepwise_forward_step(base, scope = ~ flag, data = d,
                                     criterion = "wald", slentry = 0.05)
  row <- scored_row(step, "flag")
  expect_equal(row$stat, wald_z(direct, "beta1"), tolerance = 1e-3)
  expect_true(step$accepted)
})

test_that("a logical candidate is tested on its own column (multiphase)", {
  skip_on_cran()
  d <- ccn_data()
  base <- hazard(survival::Surv(time, status) ~ 1, data = d,
                 dist = "multiphase", phases = ccn_phases(~ fla), fit = TRUE)
  direct <- hazard(survival::Surv(time, status) ~ 1, data = d,
                   dist = "multiphase", phases = ccn_phases(~ fla + flag),
                   fit = TRUE)

  expect_identical(colnames(direct$fit$x_list$constant), c("flag", "flagTRUE"))
  z_cand <- wald_z(direct, "constant.flagTRUE")
  expect_gt(abs(z_cand - wald_z(direct, "constant.flag")), 5)

  step <- .hzr_stepwise_forward_step(base, scope = list(constant = ~ flag),
                                     data = d, criterion = "wald",
                                     slentry = 0.05)
  row <- scored_row(step, "flag")
  expect_equal(row$stat, z_cand, tolerance = 1e-3)
  expect_lt(row$p_value, 1e-10)
  expect_true(step$accepted)
})

test_that("the score path's Wald fallback tests the candidate's column", {
  skip_on_cran()
  # A large effect makes the score's information indefinite at beta = 0, so
  # the candidate is declined and rescued by a Wald refit (#130). That rescue
  # reported the factor dummy's z -- about -1.3, p = 0.18 -- for a candidate
  # this strong, and the screen rejected its best variable as noise.
  d <- ccn_data(effect = 5)
  base <- hazard(survival::Surv(time, status) ~ 1, data = d,
                 dist = "multiphase", phases = ccn_phases(~ fla), fit = TRUE)
  direct <- hazard(survival::Surv(time, status) ~ 1, data = d,
                   dist = "multiphase", phases = ccn_phases(~ fla + flag),
                   fit = TRUE)
  z_cand <- wald_z(direct, "constant.flagTRUE")
  expect_gt(abs(z_cand - wald_z(direct, "constant.flag")), 5)

  step <- suppressWarnings(.hzr_stepwise_forward_step(
    base, scope = list(constant = ~ flag), data = d,
    criterion = "score", slentry = 0.05
  ))
  row <- scored_row(step, "flag")
  # The fixture must reach the fallback, or this tests the score path instead.
  expect_true(row$fallback)
  expect_identical(row$stat_type, "wald_z")
  expect_equal(row$stat, z_cand, tolerance = 1e-3)
  expect_true(step$accepted)
})

test_that("a candidate whose column is not the last is still found", {
  # model.matrix() puts main effects before interactions, so a candidate `c`
  # added to `~ a * b` lands third of four. Every fixture above has the
  # candidate last, which "take the last column" would also pass.
  set.seed(29L)
  n <- 400L
  d <- data.frame(a = stats::rnorm(n), b = stats::rnorm(n),
                  c = stats::rnorm(n))
  d$time <- stats::rexp(n) * exp(-0.5 * d$c)
  d$status <- 1L
  base <- hazard(survival::Surv(time, status) ~ a * b, data = d,
                 dist = "weibull", theta = c(0.5, 1, 0, 0, 0), fit = TRUE)
  direct <- hazard(survival::Surv(time, status) ~ a * b + c, data = d,
                   dist = "weibull", theta = c(0.5, 1, 0, 0, 0, 0), fit = TRUE)
  expect_identical(colnames(direct$data$x), c("a", "b", "c", "a:b"))
  z_cand <- wald_z(direct, "beta3")
  expect_gt(abs(z_cand - wald_z(direct, "beta4")), 5)

  step <- .hzr_stepwise_forward_step(base, scope = ~ c, data = d,
                                     criterion = "wald", slentry = 0.05)
  expect_equal(scored_row(step, "c")$stat, z_cand, tolerance = 1e-3)
})

test_that("a candidate that only reparameterises the phase is refused", {
  skip_on_cran()
  # Adding `z` to `~ z:f` turns the columns `z:fa, z:fb` into `z, z:fb`: the
  # same column space, the same likelihood. One name is new, so the column
  # difference alone "found" the candidate, and the step reported a tiny
  # p-value and accepted a variable that changed nothing. The score path
  # requires the column count to rise by exactly one; so must this.
  set.seed(5L)
  n <- 300L
  d <- data.frame(z = stats::rnorm(n),
                  f = factor(sample(c("a", "b"), n, replace = TRUE)))
  d$time <- stats::rexp(n) * exp(-0.6 * d$z)
  d$status <- 1L
  base <- hazard(survival::Surv(time, status) ~ 1, data = d,
                 dist = "multiphase", phases = ccn_phases(~ z:f), fit = TRUE)
  refit <- .hzr_refit_with_scope(base, "add", "z", phase = "constant",
                                 data = d)
  # The mechanism: one new name, no new column, no change in fit.
  expect_identical(colnames(base$fit$x_list$constant), c("z:fa", "z:fb"))
  expect_identical(colnames(refit$fit$x_list$constant), c("z", "z:fb"))
  expect_equal(refit$fit$objective, base$fit$objective, tolerance = 1e-6)

  expect_error(
    .hzr_stepwise_forward_step(base, scope = list(constant = ~ z), data = d,
                               criterion = "wald", slentry = 0.05),
    "does not add a column"
  )
})

test_that("a candidate that adds more than one column still errors", {
  set.seed(223L)
  n <- 200L
  d <- data.frame(
    time   = stats::rexp(n),
    status = 1L,
    x1     = stats::rnorm(n),
    f      = factor(sample(letters[1:3], n, replace = TRUE))
  )
  current <- hazard(survival::Surv(time, status) ~ x1, data = d,
                    dist = "weibull", theta = c(0.5, 1, 0), fit = TRUE)
  fit <- hazard(survival::Surv(time, status) ~ x1 + f, data = d,
                dist = "weibull", theta = c(0.5, 1, 0, 0, 0), fit = TRUE)
  expect_error(
    .hzr_candidate_coef_name(fit, "f", NULL, current = current),
    "expands to multiple coefficients"
  )
})
