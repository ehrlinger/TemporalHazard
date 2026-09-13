# tests/testthat/test-gof-per-subject.R
# hzr_gof() expected events are per subject (issue #254): each subject
# contributes its own cumulative hazard at exit, minus the cumulative hazard
# at its counting-process entry time.  The covariate-mean curve (par_cumhaz)
# is for plotting only and must not drive cum_expected.

.gof_avc <- local({
  data(avc, package = "TemporalHazard")
  na.omit(avc)
})

.gof_phases <- function() {
  list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
}

test_that("covariate fit: total_expected is the per-subject cumulative hazard sum", {
  avc <- .gof_avc
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ age + mal,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5, beta_age = 0, beta_mal = 0),
    fit   = TRUE
  )
  expect_true(isTRUE(fit$fit$converged))

  # Independent path: predict() on a newdata frame built from the raw columns.
  nd <- avc[, c("age", "mal")]
  nd$time <- avc$int_dead
  h_subject <- predict(fit, newdata = nd, type = "cumulative_hazard")

  s <- attr(hzr_gof(fit), "summary")
  expect_equal(s$total_expected, sum(h_subject), tolerance = 1e-8)
  expect_equal(s$total_observed, sum(avc$dead))
})

test_that("multiphase covariate fit with conserve = TRUE gives E/O = 1", {
  avc <- .gof_avc
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ age + status + mal + com_iv,
    data    = avc,
    dist    = "multiphase",
    phases  = .gof_phases(),
    fit     = TRUE,
    control = list(n_starts = 1, maxit = 1000, conserve = TRUE)
  )
  expect_true(isTRUE(fit$fit$converged))
  expect_true(isTRUE(fit$spec$control$conserve_applied))

  s <- attr(hzr_gof(fit), "summary")
  expect_equal(s$total_observed, sum(avc$dead))
  # Compare a ratio to 1, so the assertion is relative, not absolute.
  expect_equal(s$total_expected / s$total_observed, 1, tolerance = 1e-6)
  expect_equal(s$final_residual / s$total_observed, 0, tolerance = 1e-6)
})

test_that("intercept-only fit: cum_expected matches the covariate-free formula", {
  data(cabgkul, package = "TemporalHazard")
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = cabgkul,
    dist  = "weibull",
    theta = c(mu = 0.10, nu = 1.0),
    fit   = TRUE
  )
  expect_true(isTRUE(fit$fit$converged))
  gof <- hzr_gof(fit)

  # With no covariates every subject shares one curve, so the per-subject sum
  # at each exit time is (events + censored) * H(t): the pre-#254 formula.
  old <- cumsum((gof$n_event + gof$n_censor) * gof$par_cumhaz)
  expect_equal(gof$cum_expected, old, tolerance = 1e-10)
  expect_equal(gof$cum_observed, cumsum(gof$n_event))
  # And it is not trivially zero: the fit predicts events.
  expect_gt(max(gof$cum_expected), 0.5 * sum(cabgkul$dead))
})

test_that("left-truncated fit subtracts the entry-time cumulative hazard", {
  avc <- .gof_avc
  set.seed(2)
  late <- stats::runif(nrow(avc)) < 0.4
  entry <- ifelse(late, avc$int_dead * stats::runif(nrow(avc), 0.1, 0.9), 0)
  expect_gt(sum(entry > 0), 50)

  fit <- hazard(
    time       = avc$int_dead,
    status     = avc$dead,
    time_lower = entry,
    x          = as.matrix(avc[, c("age", "mal")]),
    dist       = "multiphase",
    phases     = .gof_phases(),
    fit        = TRUE,
    control    = list(n_starts = 1, maxit = 1000, conserve = TRUE)
  )
  expect_true(isTRUE(fit$fit$converged))
  expect_true(isTRUE(fit$spec$control$conserve_applied))

  nd_exit  <- data.frame(age = avc$age, mal = avc$mal, time = avc$int_dead)
  nd_entry <- data.frame(age = avc$age, mal = avc$mal, time = entry)
  h_exit  <- predict(fit, newdata = nd_exit,  type = "cumulative_hazard")
  h_entry <- predict(fit, newdata = nd_entry, type = "cumulative_hazard")
  # The entry term is not negligible, so omitting it would be caught.
  expect_gt(sum(h_entry), 10)

  s <- attr(hzr_gof(fit), "summary")
  expect_equal(s$total_expected, sum(h_exit - h_entry), tolerance = 1e-8)
  expect_equal(s$total_expected / s$total_observed, 1, tolerance = 1e-6)
})

test_that("vector-interface covariate fit uses the stored design matrix", {
  avc <- .gof_avc
  fit <- hazard(
    time   = avc$int_dead,
    status = avc$dead,
    x      = as.matrix(avc[, c("age", "mal")]),
    dist   = "weibull",
    theta  = c(0.01, 0.5, 0, 0),
    fit    = TRUE
  )
  expect_true(isTRUE(fit$fit$converged))
  expect_null(fit$data$frame)

  nd <- data.frame(age = avc$age, mal = avc$mal, time = avc$int_dead)
  s <- attr(hzr_gof(fit), "summary")
  expect_equal(s$total_expected,
               sum(predict(fit, newdata = nd, type = "cumulative_hazard")),
               tolerance = 1e-8)
})

test_that("weighted fit: observed and expected are both weighted", {
  avc <- .gof_avc
  set.seed(1)
  w <- stats::runif(nrow(avc), 0.5, 2)
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ age + mal,
    data    = avc,
    weights = w,
    dist    = "multiphase",
    phases  = .gof_phases(),
    fit     = TRUE,
    control = list(n_starts = 1, maxit = 1000, conserve = TRUE)
  )
  expect_true(isTRUE(fit$fit$converged))
  expect_true(isTRUE(fit$spec$control$conserve_applied))

  s <- attr(hzr_gof(fit), "summary")
  # CoE under weights is sum(w * H) = sum(w * d); the unweighted sums differ.
  expect_equal(s$total_observed, sum(w * avc$dead), tolerance = 1e-10)
  expect_equal(s$total_expected / s$total_observed, 1, tolerance = 1e-6)
})

test_that("custom time_grid: observed and expected cover the same subjects", {
  avc <- .gof_avc
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ age + mal,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5, beta_age = 0, beta_mal = 0),
    fit   = TRUE
  )
  # A grid holding only the first 20 distinct event times.
  grid <- sort(unique(avc$int_dead[avc$dead == 1]))[1:20]
  on_grid <- avc$int_dead %in% grid

  nd <- avc[on_grid, c("age", "mal")]
  nd$time <- avc$int_dead[on_grid]
  gof <- hzr_gof(fit, time_grid = grid)
  s <- attr(gof, "summary")
  expect_equal(s$total_observed, sum(avc$dead[on_grid]))
  expect_equal(s$total_expected,
               sum(predict(fit, newdata = nd, type = "cumulative_hazard")),
               tolerance = 1e-8)
})
