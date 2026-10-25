# hzr_deciles() on a multiphase time_windows fit (#629): the horizon ranking
# must use the coefficients of the window that contains the horizon, not the
# window each subject exited in. The reference is predict(newdata = ), which
# rebuilds and expands the design at the time it is given.

test_that("multiphase time_windows: horizon risk uses the horizon's window", {
  set.seed(629)
  n <- 400
  x1 <- stats::rnorm(n)
  tt <- stats::rweibull(n, shape = 1.2, scale = 6 * exp(-0.5 * x1))
  cc <- stats::runif(n, 2, 12)
  d <- data.frame(time = pmin(tt, cc), status = as.integer(tt <= cc),
                  x1 = x1)
  fm <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ x1, data = d,
    dist = "multiphase", time_windows = 2,
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  expect_true(isTRUE(fm$fit$converged))
  horizon <- 6

  # Premise: subjects exit in both windows and the horizon is in the second,
  # so the exit-time design read at the horizon gives a different answer.
  expect_true(any(d$time <= 2) && any(d$time > 2))
  ref_h <- predict(fm, newdata = data.frame(time = horizon, x1 = d$x1),
                   type = "cumulative_hazard")
  exit_design_h <- .hzr_multiphase_cumhaz(
    rep(horizon, n), fm$fit$theta, fm$fit$phases, fm$fit$covariate_counts,
    fm$fit$x_list
  )
  expect_gt(max(abs(exit_design_h / ref_h - 1)), 0.5)

  dc <- hzr_deciles(fm, time = horizon, groups = 4)

  ref_group <- as.integer(floor(4 * (rank(ref_h, ties.method = "first") - 1) /
                                  n)) + 1L
  expect_equal(dc$n, as.integer(tabulate(ref_group, 4)))
  expect_equal(dc$mean_survival,
               as.numeric(tapply(exp(-ref_h), ref_group, mean)),
               tolerance = 1e-12)

  # Expected events stay at each subject's own exit time, now summed over the
  # horizon's partition.
  exit_h <- predict(fm, type = "cumulative_hazard")
  expect_equal(dc$expected, as.numeric(tapply(exit_h, ref_group, sum)),
               tolerance = 1e-12)
  expect_equal(attr(dc, "overall")$total_expected, sum(exit_h),
               tolerance = 1e-12)
})

test_that("multiphase time_windows: a phase's own formula keeps its design", {
  # Only a phase that took the global x is window-expanded. Here the early
  # phase's own formula has as many columns as the expanded global design
  # (x2, x3 against x1_w1, x1_w2), so expanding it too would not fail on
  # shape: only the values can show it.
  set.seed(6291)
  n <- 400
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                  x3 = stats::rnorm(n))
  tt <- stats::rweibull(n, shape = 1.2,
                        scale = 6 * exp(-0.5 * d$x1 - 0.4 * d$x2))
  cc <- stats::runif(n, 2, 12)
  d$time <- pmin(tt, cc)
  d$status <- as.integer(tt <= cc)
  fm <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ x1, data = d,
    dist = "multiphase", time_windows = 2,
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ x2 + x3),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  expect_true(isTRUE(fm$fit$converged))
  expect_identical(colnames(fm$fit$x_list$early), c("x2", "x3"))
  expect_identical(colnames(fm$fit$x_list$constant), c("x1_w1", "x1_w2"))

  horizon <- 6
  ref_h <- predict(fm, newdata = data.frame(time = horizon, x1 = d$x1,
                                            x2 = d$x2, x3 = d$x3),
                   type = "cumulative_hazard")
  dc <- hzr_deciles(fm, time = horizon, groups = 4)
  ref_group <- as.integer(floor(4 * (rank(ref_h, ties.method = "first") - 1) /
                                  n)) + 1L
  expect_equal(dc$n, as.integer(tabulate(ref_group, 4)))
  expect_equal(dc$mean_survival,
               as.numeric(tapply(exp(-ref_h), ref_group, mean)),
               tolerance = 1e-12)
})

test_that("multiphase time_windows: a legacy fit is rebuilt only where it inherited", {
  # A fit saved before attr(x_list, "from_formula") is classified by column
  # names. Here the early phase's own formula uses columns NAMED x1_w1 and
  # x1_w2, independent of x1, which match the global ~ x1 expanded at one cut
  # (Codex on #657). The stored own-formula design is not windowed, so it is
  # also the horizon design and must be kept.
  set.seed(6292)
  n <- 400
  d <- data.frame(x1 = stats::rnorm(n), x1_w1 = stats::rnorm(n),
                  x1_w2 = stats::rnorm(n))
  tt <- stats::rweibull(n, shape = 1.2,
                        scale = 6 * exp(-0.5 * d$x1 - 0.4 * d$x1_w1))
  cc <- stats::runif(n, 2, 12)
  d$time <- pmin(tt, cc)
  d$status <- as.integer(tt <= cc)
  fm <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ x1, data = d,
    dist = "multiphase", time_windows = 2,
    phases = list(
      early    = hzr_phase("cdf", t_half = 2, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ x1_w1 + x1_w2),
      constant = hzr_phase("constant", formula = ~ 1)
    ),
    fit = TRUE
  ))
  expect_true(isTRUE(fm$fit$converged))
  legacy <- fm
  attr(legacy$fit$x_list, "from_formula") <- NULL
  # Premise: the name fallback misreads the phase as inheriting.
  expect_true(.hzr_phase_inherits_global(legacy, "early"))

  horizon <- 6
  ref_h <- .hzr_multiphase_cumhaz(
    rep(horizon, n), fm$fit$theta, fm$fit$phases, fm$fit$covariate_counts,
    fm$fit$x_list
  )
  ref_group <- as.integer(floor(4 * (rank(ref_h, ties.method = "first") - 1) /
                                  n)) + 1L
  ref_ms <- as.numeric(tapply(exp(-ref_h), ref_group, mean))
  expect_equal(hzr_deciles(fm, time = horizon, groups = 4)$mean_survival,
               ref_ms, tolerance = 1e-12)
  expect_equal(hzr_deciles(legacy, time = horizon, groups = 4)$mean_survival,
               ref_ms, tolerance = 1e-12)

  # A legacy fit whose phase really inherited is still rebuilt at the horizon.
  inh <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ x1, data = d,
    dist = "multiphase", time_windows = 2,
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  attr(inh$fit$x_list, "from_formula") <- NULL
  inh_ref <- predict(inh, newdata = data.frame(time = horizon, x1 = d$x1),
                     type = "cumulative_hazard")
  inh_group <- as.integer(floor(4 * (rank(inh_ref, ties.method = "first") - 1) /
                                  n)) + 1L
  expect_equal(hzr_deciles(inh, time = horizon, groups = 4)$mean_survival,
               as.numeric(tapply(exp(-inh_ref), inh_group, mean)),
               tolerance = 1e-12)
})
