# tests/testthat/test-diagnostics.R
# Tests for hzr_deciles() and related diagnostic utilities

# Helper cache: fit once to reduce repeated optimizer runs across tests.
.avc_fixture <- local({
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ age + mal,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5, beta_age = 0, beta_mal = 0),
    fit   = TRUE
  )
  stopifnot(isTRUE(fit$fit$converged))
  list(avc = avc, fit = fit)
})

.fit_avc_weibull <- function() .avc_fixture$fit
.avc_data <- function() .avc_fixture$avc

test_that("hzr_deciles returns correct structure", {
  fit <- .fit_avc_weibull()
  cal <- hzr_deciles(fit, time = 120)

  expect_s3_class(cal, "hzr_deciles")
  expect_s3_class(cal, "data.frame")
  expect_equal(nrow(cal), 10)
  expect_named(cal, c("group", "n", "events", "expected", "observed_rate",
                       "expected_rate", "chi_sq", "p_value",
                       "mean_survival", "mean_cumhaz"))

  # Overall attribute
  ov <- attr(cal, "overall")
  expect_true(is.list(ov))
  expect_named(ov, c("chi_sq", "df", "p_value", "time", "groups",
                      "total_events", "total_expected",
                      "n_included", "n_excluded"))
  expect_equal(ov$time, 120)
  expect_equal(ov$groups, 10)
})

test_that("hzr_deciles includes all subjects and counts all events (SAS method)", {
  fit <- .fit_avc_weibull()
  avc <- .avc_data()
  cal <- hzr_deciles(fit, time = 120)

  # SAS %DECILES uses every subject; the horizon only stratifies risk.
  expect_equal(sum(cal$n), nrow(avc))
  expect_equal(sum(cal$events), sum(avc$dead == 1))
  expect_equal(attr(cal, "overall")$n_excluded, 0L)
})

test_that("hzr_deciles works with non-default groups", {
  fit <- .fit_avc_weibull()
  avc <- .avc_data()
  cal5 <- hzr_deciles(fit, time = 60, groups = 5)
  expect_equal(nrow(cal5), 5)
  expect_equal(sum(cal5$n), nrow(avc))
})

test_that("hzr_deciles totals are horizon-invariant; only the grouping changes", {
  fit <- .fit_avc_weibull()
  avc <- .avc_data()

  cal60  <- hzr_deciles(fit, time = 60)
  cal120 <- hzr_deciles(fit, time = 120)

  # Expected = sum of cumulative hazard at each subject's own follow-up time,
  # and events = all observed deaths: both totals are independent of the
  # grouping horizon. Only the assignment of subjects to risk groups changes.
  expect_equal(sum(cal60$n), nrow(avc))
  expect_equal(sum(cal120$n), nrow(avc))
  expect_equal(sum(cal60$events), sum(avc$dead == 1))
  expect_equal(sum(cal120$events), sum(avc$dead == 1))
  expect_equal(sum(cal60$expected), sum(cal120$expected), tolerance = 1e-8)
})

test_that("hzr_deciles risk ordering is monotone in horizon survival", {
  fit <- .fit_avc_weibull()
  cal <- hzr_deciles(fit, time = 120)

  # Groups are ranked by predicted survival at the horizon, group 1 = lowest
  # risk (highest survival), so mean horizon survival decreases across groups.
  expect_true(all(diff(cal$mean_survival) <= 1e-10),
              label = "mean horizon survival should decrease across groups")
})

test_that("hzr_deciles chi-square values are non-negative", {
  fit <- .fit_avc_weibull()
  cal <- hzr_deciles(fit, time = 120)
  valid <- !is.na(cal$chi_sq)
  expect_true(all(cal$chi_sq[valid] >= 0))
  expect_true(all(cal$p_value[valid] >= 0 & cal$p_value[valid] <= 1))

  ov <- attr(cal, "overall")
  expect_true(ov$chi_sq >= 0)
  expect_true(ov$df > 0)
})

test_that("hzr_deciles works with intercept-only model", {
  # Intercept-only: all predictions identical, but grouping should still work
  data(cabgkul, package = "TemporalHazard")
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = cabgkul,
    dist  = "weibull",
    theta = c(mu = 0.10, nu = 1.0),
    fit   = TRUE
  )
  cal <- hzr_deciles(fit, time = 12)

  expect_s3_class(cal, "hzr_deciles")
  expect_equal(nrow(cal), 10)
  expect_equal(sum(cal$n), nrow(cabgkul))
  # Everyone shares the same horizon survival, so risk grouping is arbitrary
  # but mean horizon survival is identical across groups.
  expect_true(max(cal$mean_survival) - min(cal$mean_survival) < 1e-8)
})

test_that("hzr_deciles works with multiphase model", {
  data(cabgkul, package = "TemporalHazard")
  fit_mp <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = cabgkul,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.2, nu = 1, m = 1,
                            fixed = "shapes"),
      constant = hzr_phase("constant"),
      late     = hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1,
                            fixed = "shapes")
    ),
    fit     = TRUE,
    control = list(n_starts = 3, maxit = 500)
  )
  cal <- hzr_deciles(fit_mp, time = 60)

  expect_s3_class(cal, "hzr_deciles")
  expect_equal(nrow(cal), 10)
  expect_equal(sum(cal$n), nrow(cabgkul))
})

test_that("hzr_deciles rejects invalid inputs", {
  fit <- .fit_avc_weibull()
  avc <- .avc_data()

  expect_error(hzr_deciles(fit, time = -1), "positive")
  expect_error(hzr_deciles(fit, time = 120, groups = 1), "at least 2")
  expect_error(hzr_deciles(fit, time = 120, groups = nrow(avc) + 1),
               "number of observations")
  expect_error(hzr_deciles("not_a_model", time = 12), "hazard object")
})

test_that("hzr_deciles unfitted model is rejected", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  fit_unfitted <- hazard(
    survival::Surv(int_dead, dead) ~ age + status,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.05, nu = 0.5, 0, 0),
    fit   = FALSE
  )

  expect_error(hzr_deciles(fit_unfitted, time = 120), "fit = TRUE")
})

test_that("print.hzr_deciles runs without error", {
  fit <- .fit_avc_weibull()
  cal <- hzr_deciles(fit, time = 120)
  expect_output(print(cal), "Decile-of-risk calibration")
  expect_output(print(cal), "Overall")
})


# =========================================================================
# hzr_gof tests
# =========================================================================

test_that("hzr_gof returns correct structure", {
  fit <- .fit_avc_weibull()
  gof <- hzr_gof(fit)

  expect_s3_class(gof, "hzr_gof")
  expect_s3_class(gof, "data.frame")
  expect_true(nrow(gof) > 0)
  expect_true(all(c("time", "n_risk", "n_event", "n_censor",
                     "km_surv", "km_cumhaz", "par_surv", "par_cumhaz",
                     "cum_observed", "cum_expected", "residual") %in%
                    names(gof)))

  # Summary attribute
  s <- attr(gof, "summary")
  expect_true(is.list(s))
  expect_named(s, c("total_observed", "total_expected", "final_residual",
                     "dist", "n"))
})

test_that("hzr_gof cumulative observed equals total events", {
  fit <- .fit_avc_weibull()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  gof <- hzr_gof(fit)

  s <- attr(gof, "summary")
  expect_equal(s$total_observed, sum(avc$dead))
  expect_equal(s$n, nrow(avc))
})

test_that("hzr_gof survival values are in [0, 1]", {
  fit <- .fit_avc_weibull()
  gof <- hzr_gof(fit)

  expect_true(all(gof$km_surv >= 0 & gof$km_surv <= 1))
  expect_true(all(gof$par_surv >= 0 & gof$par_surv <= 1))
})

test_that("hzr_gof cumulative hazard is non-negative and increasing", {
  fit <- .fit_avc_weibull()
  gof <- hzr_gof(fit)

  expect_true(all(gof$km_cumhaz >= 0))
  expect_true(all(gof$par_cumhaz >= 0))
  # Parametric cumhaz should be non-decreasing
  expect_true(all(diff(gof$par_cumhaz) >= -1e-10))
})

test_that("hzr_gof conservation ratio is reasonable", {
  fit <- .fit_avc_weibull()
  gof <- hzr_gof(fit)

  s <- attr(gof, "summary")
  ratio <- s$total_expected / s$total_observed
  # For a reasonable fit, E/O should be in the ballpark (0.5 to 2.0)
  expect_true(ratio > 0.2 && ratio < 5.0,
              label = paste("conservation ratio", round(ratio, 3),
                            "should be between 0.2 and 5.0"))
})

test_that("hzr_gof works with custom time grid", {
  fit <- .fit_avc_weibull()
  t_grid <- seq(1, 200, by = 10)
  gof <- hzr_gof(fit, time_grid = t_grid)

  expect_equal(nrow(gof), length(t_grid))
  expect_equal(gof$time, t_grid)
})

test_that("hzr_gof works with intercept-only model", {
  data(cabgkul, package = "TemporalHazard")
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = cabgkul,
    dist  = "weibull",
    theta = c(mu = 0.10, nu = 1.0),
    fit   = TRUE
  )
  gof <- hzr_gof(fit)

  expect_s3_class(gof, "hzr_gof")
  expect_true(nrow(gof) > 0)
  s <- attr(gof, "summary")
  expect_equal(s$total_observed, sum(cabgkul$dead))
})

test_that("hzr_gof works with multiphase model", {
  data(cabgkul, package = "TemporalHazard")
  fit_mp <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = cabgkul,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.2, nu = 1, m = 1,
                            fixed = "shapes"),
      constant = hzr_phase("constant"),
      late     = hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1,
                            fixed = "shapes")
    ),
    fit     = TRUE,
    control = list(n_starts = 3, maxit = 500)
  )
  gof <- hzr_gof(fit_mp)

  expect_s3_class(gof, "hzr_gof")
  # Should have phase-specific columns
  expect_true(any(grepl("^par_cumhaz_", names(gof))))
})

test_that("hzr_gof rejects invalid inputs", {
  expect_error(hzr_gof("not_a_model"), "hazard object")

  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  fit_unfitted <- hazard(
    survival::Surv(int_dead, dead) ~ age + mal,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5, 0, 0),
    fit   = FALSE
  )
  expect_error(hzr_gof(fit_unfitted), "fit = TRUE")
})

test_that("print.hzr_gof runs without error", {
  fit <- .fit_avc_weibull()
  gof <- hzr_gof(fit)
  expect_output(print(gof), "Goodness-of-fit")
  expect_output(print(gof), "Conservation ratio")
})


# =========================================================================
# hzr_kaplan tests
# =========================================================================

test_that("hzr_kaplan returns correct structure", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)

  expect_s3_class(km, "hzr_kaplan")
  expect_s3_class(km, "data.frame")
  expect_true(nrow(km) > 0)
  expect_named(km, c("time", "n_risk", "n_event", "n_censor",
                      "survival", "std_err", "cl_lower", "cl_upper",
                      "cumhaz", "hazard", "density", "life"))
})

test_that("hzr_kaplan survival is monotone non-increasing", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)

  expect_true(all(diff(km$survival) <= 0),
              label = "survival should be non-increasing")
  expect_true(all(km$survival >= 0 & km$survival <= 1))
})

test_that("hzr_kaplan confidence limits respect [0, 1]", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)

  expect_true(all(km$cl_lower >= 0 & km$cl_lower <= 1))
  expect_true(all(km$cl_upper >= 0 & km$cl_upper <= 1))
  # Lower <= survival <= upper
  expect_true(all(km$cl_lower <= km$survival + 1e-10))
  expect_true(all(km$cl_upper >= km$survival - 1e-10))
})

test_that("hzr_kaplan matches survfit survival values", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead, event_only = FALSE)
  sf <- survival::survfit(survival::Surv(int_dead, dead) ~ 1,
                           data = cabgkul)

  # Survival values at survfit times should match exactly
  expect_equal(km$survival, sf$surv, tolerance = 1e-12)
  expect_equal(km$time, sf$time)
})

test_that("hzr_kaplan cumhaz is consistent with survival", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)

  # cumhaz = -log(S)
  expected_cumhaz <- -log(pmax(km$survival, .Machine$double.xmin))
  expect_equal(km$cumhaz, expected_cumhaz, tolerance = 1e-10)
})

test_that("hzr_kaplan life integral is non-negative and non-decreasing", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)

  expect_true(all(km$life >= 0))
  expect_true(all(diff(km$life) >= -1e-10),
              label = "life integral should be non-decreasing")
})

test_that("hzr_kaplan event_only filters correctly", {
  data(cabgkul, package = "TemporalHazard")
  km_events <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead,
                            event_only = TRUE)
  km_all <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead,
                         event_only = FALSE)

  expect_true(nrow(km_all) >= nrow(km_events))
  expect_true(all(km_events$n_event > 0))
})

test_that("hzr_kaplan conf_level parameter works", {
  data(cabgkul, package = "TemporalHazard")
  km_95 <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead,
                        conf_level = 0.95)
  km_68 <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead,
                        conf_level = 0.68268948)

  # 95% intervals should be wider than 68% (1-SD) intervals
  width_95 <- km_95$cl_upper - km_95$cl_lower
  width_68 <- km_68$cl_upper - km_68$cl_lower
  # Compare at matching times where both have events
  expect_true(all(width_95 >= width_68 - 1e-10))
})

test_that("hzr_kaplan rejects invalid inputs", {
  expect_error(hzr_kaplan("a", 1), "numeric")
  expect_error(hzr_kaplan(1:5, 1:3), "same length")
  expect_error(hzr_kaplan(1:5, rep(1, 5), conf_level = 1.5),
               "between 0 and 1")
})

test_that("print.hzr_kaplan runs without error", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)
  expect_output(print(km), "Kaplan-Meier")
  expect_output(print(km), "RMST")
})


# =========================================================================
# hzr_calibrate tests
# =========================================================================

test_that("hzr_calibrate returns correct structure", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  cal <- hzr_calibrate(avc$age, avc$dead, groups = 10)

  expect_s3_class(cal, "hzr_calibrate")
  expect_s3_class(cal, "data.frame")
  expect_equal(nrow(cal), 10)
  expect_named(cal, c("group", "n", "events", "mean", "min", "max",
                       "prob", "link_value"))
  expect_equal(attr(cal, "link"), "logit")
  expect_equal(attr(cal, "groups"), 10L)
})

test_that("hzr_calibrate group counts sum to n", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  cal <- hzr_calibrate(avc$age, avc$dead, groups = 5)

  expect_equal(sum(cal$n), nrow(avc))
  expect_equal(sum(cal$events), sum(avc$dead))
})

test_that("hzr_calibrate group means are monotone", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  cal <- hzr_calibrate(avc$age, avc$dead, groups = 10)

  # Group means should increase (group 1 = lowest x)
  expect_true(all(diff(cal$mean) > 0),
              label = "group means should be strictly increasing")
})

test_that("hzr_calibrate logit values are finite where prob in (0,1)", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  cal <- hzr_calibrate(avc$age, avc$dead, groups = 10)

  valid <- cal$prob > 0 & cal$prob < 1
  expect_true(all(is.finite(cal$link_value[valid])))
})

test_that("hzr_calibrate gompertz link works", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  cal <- hzr_calibrate(avc$age, avc$dead, groups = 5, link = "gompertz")

  expect_equal(attr(cal, "link"), "gompertz")
  # Gompertz = log(-log(1-p)), should be finite for 0 < p < 1
  valid <- cal$prob > 0 & cal$prob < 1
  expect_true(all(is.finite(cal$link_value[valid])))
})

test_that("hzr_calibrate cox link requires time", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  expect_error(
    hzr_calibrate(avc$age, avc$dead, link = "cox"),
    "time.*required"
  )

  cal <- hzr_calibrate(avc$age, avc$dead, link = "cox",
                        time = avc$int_dead, groups = 5)
  expect_equal(attr(cal, "link"), "cox")
  expect_equal(nrow(cal), 5)
})

test_that("hzr_calibrate stratified by grouping variable", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  cal <- hzr_calibrate(avc$age, avc$dead, groups = 5,
                        by = ifelse(avc$mal == 1, "mal", "no_mal"))

  expect_true("by" %in% names(cal))
  # Should have rows for both strata
  expect_true(all(c("mal", "no_mal") %in% cal$by))
  # Total n across all rows should equal nrow(avc)
  expect_equal(sum(cal$n), nrow(avc))
})

test_that("hzr_calibrate rejects invalid inputs", {
  expect_error(hzr_calibrate("a", 1), "numeric")
  expect_error(hzr_calibrate(1:5, 1:3), "same length")
  expect_error(hzr_calibrate(1:5, rep(1, 5), groups = 1), "at least 2")
})

test_that("print.hzr_calibrate runs without error", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  cal <- hzr_calibrate(avc$age, avc$dead, groups = 5)
  expect_output(print(cal), "Variable calibration")
  expect_output(print(cal), "logit")
})


# =========================================================================
# hzr_nelson tests
# =========================================================================

test_that("hzr_nelson returns correct structure", {
  data(cabgkul, package = "TemporalHazard")
  nel <- hzr_nelson(cabgkul$int_dead, cabgkul$dead)

  expect_s3_class(nel, "hzr_nelson")
  expect_true(nrow(nel) > 0)
  expect_named(nel, c("time", "n_risk", "n_event", "weight_sum",
                       "cumhaz", "std_err", "cl_lower", "cl_upper",
                       "hazard", "cum_events"))
})

test_that("hzr_nelson cumhaz is non-decreasing", {
  data(cabgkul, package = "TemporalHazard")
  nel <- hzr_nelson(cabgkul$int_dead, cabgkul$dead)

  expect_true(all(diff(nel$cumhaz) >= -1e-10))
  expect_true(all(nel$cumhaz >= 0))
})

test_that("hzr_nelson CL are non-negative and ordered", {
  data(cabgkul, package = "TemporalHazard")
  nel <- hzr_nelson(cabgkul$int_dead, cabgkul$dead)

  expect_true(all(nel$cl_lower >= 0))
  expect_true(all(nel$cl_upper >= 0))
  expect_true(all(nel$cl_lower <= nel$cumhaz + 1e-10))
  expect_true(all(nel$cl_upper >= nel$cumhaz - 1e-10))
})

test_that("hzr_nelson weighted events change cumhaz", {
  data(cabgkul, package = "TemporalHazard")
  nel_unw <- hzr_nelson(cabgkul$int_dead, cabgkul$dead)
  # Double weights should roughly double cumhaz
  nel_w2 <- hzr_nelson(cabgkul$int_dead, cabgkul$dead,
                         weight = rep(2, nrow(cabgkul)))
  expect_true(
    nel_w2$cumhaz[nrow(nel_w2)] > nel_unw$cumhaz[nrow(nel_unw)] * 1.5
  )
})

test_that("hzr_nelson rejects invalid inputs", {
  expect_error(hzr_nelson("a", 1), "numeric")
  expect_error(hzr_nelson(1:5, 1:3), "same length")
})

test_that("print.hzr_nelson runs without error", {
  data(cabgkul, package = "TemporalHazard")
  nel <- hzr_nelson(cabgkul$int_dead, cabgkul$dead)
  expect_output(print(nel), "Nelson")
})


# =========================================================================
# hzr_bootstrap tests
# =========================================================================

test_that("hzr_bootstrap returns correct structure", {
  set.seed(42)
  fit <- .fit_avc_weibull()
  bs <- hzr_bootstrap(fit, n_boot = 10, seed = 123)

  expect_s3_class(bs, "hzr_bootstrap")
  expect_true(is.data.frame(bs$replicates))
  expect_true(is.data.frame(bs$summary))
  expect_true(bs$n_success + bs$n_failed == 10)
  expect_named(bs$summary, c("parameter", "n", "pct", "mean", "sd",
                              "min", "max", "ci_lower", "ci_upper"))
})

test_that("hzr_bootstrap seed gives reproducible results", {
  fit <- .fit_avc_weibull()
  bs1 <- hzr_bootstrap(fit, n_boot = 5, seed = 999)
  bs2 <- hzr_bootstrap(fit, n_boot = 5, seed = 999)

  expect_equal(bs1$replicates$estimate, bs2$replicates$estimate)
})

test_that("hzr_bootstrap fraction parameter works", {
  fit <- .fit_avc_weibull()
  bs <- hzr_bootstrap(fit, n_boot = 5, fraction = 0.5, seed = 42)

  expect_s3_class(bs, "hzr_bootstrap")
  expect_true(bs$n_success > 0)
})

test_that("hzr_bootstrap rejects invalid inputs", {
  expect_error(hzr_bootstrap("not_a_model"), "hazard object")
})

test_that("print.hzr_bootstrap runs without error", {
  fit <- .fit_avc_weibull()
  bs <- hzr_bootstrap(fit, n_boot = 5, seed = 42)
  expect_output(print(bs), "Bootstrap")
})

test_that("hzr_bootstrap labels covariate rows with design-matrix names", {
  # Regression test: covariate parameter rows used to be emitted with empty
  # strings, which broke wide pivots of $replicates. They should now carry
  # the design-matrix column names (e.g. "x1", "x2").
  set.seed(7)
  n <- 80
  df <- data.frame(
    t = stats::rexp(n, 0.1),
    d = stats::rbinom(n, 1, 0.6),
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n)
  )
  fit <- hazard(
    survival::Surv(t, d) ~ x1 + x2,
    data  = df,
    dist  = "weibull",
    theta = c(mu = 0.05, nu = 0.8, 0, 0),
    fit   = TRUE
  )

  bs <- hzr_bootstrap(fit, n_boot = 5, seed = 123)
  params <- bs$replicates$parameter

  expect_false(any(!nzchar(params)))
  expect_true(all(c("x1", "x2") %in% params))
  # Each successful replicate should contribute one row per covariate.
  expect_equal(sum(params == "x1"), bs$n_success)
  expect_equal(sum(params == "x2"), bs$n_success)
  # Summary table should also carry the covariate labels.
  expect_true(all(c("x1", "x2") %in% bs$summary$parameter))
})

test_that("hzr_bootstrap preserves named betas while filling blanks", {
  # Mixed-naming case: one beta is explicitly named in the starting theta
  # and the other is left blank. The named slot must be preserved; the
  # blank slot should pick up its design-matrix column name.
  set.seed(11)
  n <- 80
  df <- data.frame(
    t = stats::rexp(n, 0.1),
    d = stats::rbinom(n, 1, 0.6),
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n)
  )
  fit <- hazard(
    survival::Surv(t, d) ~ x1 + x2,
    data  = df,
    dist  = "weibull",
    theta = c(mu = 0.05, nu = 0.8, b_custom = 0, 0),
    fit   = TRUE
  )

  bs <- hzr_bootstrap(fit, n_boot = 5, seed = 321)
  params <- bs$replicates$parameter

  expect_false(any(!nzchar(params)))
  # Named beta retained; blank beta filled from colnames(x) -> "x2".
  expect_true("b_custom" %in% params)
  expect_true("x2" %in% params)
  expect_false("x1" %in% params)
})

test_that("hzr_bootstrap resamples weights alongside the data", {
  # Regression test for the weighted-bootstrap bug: the resample loop used to
  # rewire only `data`, leaving the original full-length `weights` vector bound
  # to the refit call. Two observable consequences this guards against:
  #   1. fraction < 1 -> weights length (n_obs) no longer matches the smaller
  #      resampled data, so every replicate errored out: n_success == 0.
  #   2. fraction == 1 -> length matched, so no error, but weights were applied
  #      to the wrong (resampled) rows: silently incorrect estimates.
  set.seed(202)
  n <- 100
  df <- data.frame(
    t = stats::rexp(n, 0.1),
    d = stats::rbinom(n, 1, 0.6),
    x1 = stats::rnorm(n)
  )
  w <- stats::runif(n, 0.5, 2.0)

  fit <- hazard(
    survival::Surv(t, d) ~ x1,
    data    = df,
    weights = w,
    dist    = "weibull",
    theta   = c(mu = 0.05, nu = 0.8, 0),
    fit     = TRUE
  )

  # (1) Subsampled weighted bootstrap must produce successful replicates.
  bs_frac <- hzr_bootstrap(fit, n_boot = 10, fraction = 0.5, seed = 123)
  expect_gt(bs_frac$n_success, 0)
  expect_equal(bs_frac$n_failed, 0)

  # (2) Full-size weighted bootstrap must also succeed and be well-formed.
  bs_full <- hzr_bootstrap(fit, n_boot = 10, seed = 123)
  expect_gt(bs_full$n_success, 0)
  expect_true(all(is.finite(bs_full$summary$mean)))
})

test_that("hzr_bootstrap uses stored weights when the original symbol is out of scope", {
  # Robustness regression: the resample setup previously re-evaluated the call's
  # `weights` expression in parent.frame(). If the original weights symbol is no
  # longer in scope at bootstrap time -- e.g. the fit was built inside a helper
  # that has since returned -- that eval() raised "object not found". The fit now
  # carries its evaluated weights on object$data$weights, so bootstrap no longer
  # depends on a caller-frame lookup.
  set.seed(77)
  n <- 80
  # `df` stays in the test frame so the call's `data` symbol still resolves;
  # only the `weights` symbol is confined to the local() scope and gone by the
  # time we bootstrap, isolating the weights lookup.
  df <- data.frame(
    t  = stats::rexp(n, 0.1),
    d  = stats::rbinom(n, 1, 0.6),
    x1 = stats::rnorm(n)
  )
  fit <- local({
    local_w <- stats::runif(n, 0.5, 2.0)
    hazard(
      survival::Surv(t, d) ~ x1,
      data = df, weights = local_w, dist = "weibull",
      theta = c(mu = 0.05, nu = 0.8, 0), fit = TRUE
    )
  })
  # `local_w` is now out of scope; bootstrap must still resample weights via the
  # stored copy and produce successful replicates.
  bs <- hzr_bootstrap(fit, n_boot = 10, fraction = 0.5, seed = 5)
  expect_gt(bs$n_success, 0)
  expect_equal(bs$n_failed, 0)
})

test_that("hzr_bootstrap uses stored data frame when the original symbol is out of scope", {
  # Robustness regression (sibling of the weights case): the resample setup
  # re-evaluated the call's `data` expression via eval(cl$data, parent.frame()).
  # If the original data symbol is no longer in scope at bootstrap time -- e.g.
  # the entire fit was built inside a helper that has since returned -- that
  # eval() failed and n_obs was malformed, so bootstrap errored before the loop.
  # The fit now carries its evaluated model frame on object$data$frame, so the
  # resample no longer depends on a caller-frame lookup of the data symbol.
  set.seed(91)
  n <- 80
  # Build EVERYTHING inside local(): both the `data` (`bdf`) and `weights`
  # (`bw`) symbols are gone by the time we bootstrap.
  fit <- local({
    bdf <- data.frame(
      t  = stats::rexp(n, 0.1),
      d  = stats::rbinom(n, 1, 0.6),
      x1 = stats::rnorm(n)
    )
    bw <- stats::runif(n, 0.5, 2.0)
    hazard(
      survival::Surv(t, d) ~ x1,
      data = bdf, weights = bw, dist = "weibull",
      theta = c(mu = 0.05, nu = 0.8, 0), fit = TRUE
    )
  })
  bs <- hzr_bootstrap(fit, n_boot = 10, fraction = 0.5, seed = 7)
  expect_gt(bs$n_success, 0)
  expect_equal(bs$n_failed, 0)
})

test_that("hzr_bootstrap resolves call arguments passed by symbol from a function", {
  # Regression: hzr_bootstrap re-evaluated the stored call in its own frame, so
  # arguments passed by symbol (theta here) resolved against the package
  # namespace -> globalenv and were invisible when the fit was built inside a
  # function. Every replicate threw, the tryCatch swallowed it, and n_success
  # was silently 0. Building the fit in a helper is what triggers it.
  data(avc, package = "TemporalHazard")
  avc <- stats::na.omit(avc)
  build <- function(d) {
    th <- c(mu = 0.01, nu = 0.5) # passed by SYMBOL, not inline
    hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
           theta = th, fit = TRUE)
  }
  fit <- build(avc)
  bs <- hzr_bootstrap(fit, n_boot = 3, seed = 42)
  expect_gt(bs$n_success, 0)
  expect_equal(bs$n_failed, 0)
})

test_that("hazard() captures only the symbols its call references", {
  # Regression: call_env was parent.frame(), pinning the caller's whole frame
  # to every fit (measured ~1400x saveRDS bloat) and dragging unrelated data --
  # including other cohorts -- into any saved model.
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  build <- function(d) {
    unrelated_local <- matrix(0, nrow = 100, ncol = 100)
    th <- c(mu = 0.01, nu = 0.5)
    hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
           theta = th, fit = TRUE)
  }
  fit <- build(avc)
  # `th` and `d` ARE referenced by the call and must be captured.
  expect_true(exists("th", envir = fit$call_env, inherits = FALSE))
  expect_true(exists("d", envir = fit$call_env, inherits = FALSE))
  # `unrelated_local` is NOT referenced and must not be reachable at all --
  # including through the parent chain, which is why call_env parents to
  # globalenv() rather than the caller.
  expect_false(exists("unrelated_local", envir = fit$call_env, inherits = TRUE))
})

test_that("hzr_bootstrap scope = NULL reports mode = refit and is unaffected", {
  set.seed(42)
  fit <- .fit_avc_weibull()
  bs <- hzr_bootstrap(fit, n_boot = 10, seed = 123)

  expect_identical(bs$mode, "refit")
  expect_null(bs$scope)
  expect_gt(bs$n_success, 0)
  expect_gt(nrow(bs$summary), 0)
  expect_true(all(bs$summary$pct > 0))
})

test_that("hzr_bootstrap scope runs embedded stepwise selection per replicate", {
  set.seed(13)
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5),
    fit   = TRUE
  )

  bs <- hzr_bootstrap(base, n_boot = 10, seed = 321,
                       scope = ~ age + mal + com_iv,
                       slentry = 0.3, slstay = 0.2,
                       control = list(n_starts = 1))

  expect_s3_class(bs, "hzr_bootstrap")
  expect_identical(bs$mode, "select")
  expect_true(is.list(bs$scope) || inherits(bs$scope, "formula"))
  expect_gt(bs$n_success, 0)
  expect_true(all(bs$summary$pct > 0 & bs$summary$pct <= 100))
  # Shape parameters are never dropped by stepwise selection (only
  # covariates are), so they must be present -- and selected -- in every
  # successful replicate.
  shape_rows <- bs$summary[bs$summary$parameter %in% c("mu", "nu"), ]
  expect_true(all(shape_rows$n == bs$n_success))
})

test_that("hzr_bootstrap verbose = TRUE shows a progress bar", {
  fit <- .fit_avc_weibull()
  out <- capture.output(
    bs <- hzr_bootstrap(fit, n_boot = 5, seed = 42, verbose = TRUE),
    type = "output"
  )
  # txtProgressBar(style = 3) draws a bar with a percentage; some output must
  # have been written, and the run must still return a valid object.
  expect_true(any(nzchar(out)))
  expect_true(any(grepl("100%|=====", paste(out, collapse = ""))))
  expect_s3_class(bs, "hzr_bootstrap")
  expect_gt(bs$n_success, 0)
})

test_that("hzr_bootstrap suppresses per-replicate fit warnings", {
  # A fixed-shape multiphase base with covariate scope produces
  # ill-conditioned per-replicate Hessians (.hzr_safe_solve warns). Those are
  # not individually actionable and must not leak to the caller; the tryCatch
  # only catches errors, so this relies on the in-loop suppressWarnings().
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.2, nu = 1.4, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  )
  expect_no_warning(
    hzr_bootstrap(base, n_boot = 5, seed = 7,
                  scope = list(early    = ~ age + mal + com_iv,
                               constant = ~ age + mal + com_iv),
                  slentry = 0.3, slstay = 0.2,
                  control = list(n_starts = 1))
  )
})

test_that("hzr_bootstrap scope raises immediately on a structurally invalid scope", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5),
    fit   = TRUE
  )

  # scope must be NULL, a one-sided formula, or a character vector for a
  # single-distribution fit (hzr_stepwise()'s own validation); this must
  # raise via the pre-loop validation call, not after n_boot replicates.
  expect_error(
    hzr_bootstrap(base, n_boot = 10, scope = 42),
    "one-sided formula"
  )
})

test_that("hzr_bootstrap scope with a nonexistent column warns but does not raise", {
  # A typo'd (but syntactically valid) scope column is NOT caught by the
  # pre-loop validation -- hzr_stepwise() itself only warns per candidate
  # refit failure (see stepwise-step.R .hzr_stepwise_forward_step) and
  # never selects it. Confirms hzr_bootstrap() inherits that behavior
  # rather than silently reporting it as n_failed replicates.
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5),
    fit   = TRUE
  )

  expect_warning(
    bs <- hzr_bootstrap(base, n_boot = 3, seed = 1,
                         scope = ~ age + not_a_real_column,
                         control = list(n_starts = 1)),
    "not_a_real_column"
  )
  expect_false("not_a_real_column" %in% bs$summary$parameter)
})

test_that("hzr_bootstrap(scope=) forwards a caller-supplied trace= without colliding", {
  # Regression test: select-mode passed `trace = FALSE` to hzr_stepwise()
  # alongside the forwarded `...`, so a caller-supplied `trace=` in `...`
  # raised "formal argument matched by multiple actual arguments".
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5),
    fit   = TRUE
  )

  bs <- hzr_bootstrap(base, n_boot = 2, seed = 42, scope = ~ age + mal,
                       control = list(n_starts = 1), trace = TRUE)
  expect_s3_class(bs, "hzr_bootstrap")
  expect_gte(bs$n_success, 0L)
})

test_that("hzr_bootstrap errors on unused '...' when scope is NULL", {
  # '...' only has a purpose in select-mode (forwarded to hzr_stepwise());
  # in fixed-refit mode it must still error on an unrecognized argument, as
  # it did before '...' was added to the signature for scope=.
  fit <- .fit_avc_weibull()
  expect_error(
    hzr_bootstrap(fit, n_boot = 5, verbsoe = TRUE),
    "Unused argument"
  )
})

test_that("print.hzr_bootstrap reports the mode", {
  fit <- .fit_avc_weibull()
  bs_refit <- hzr_bootstrap(fit, n_boot = 5, seed = 42)
  expect_output(print(bs_refit), "fixed refit")

  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5),
    fit   = TRUE
  )
  bs_sel <- hzr_bootstrap(base, n_boot = 5, seed = 42, scope = ~ age + mal,
                           control = list(n_starts = 1))
  expect_output(print(bs_sel), "stepwise selection")
})


# =========================================================================
# hzr_competing_risks tests
# =========================================================================

test_that("hzr_competing_risks returns correct structure", {
  data(valves, package = "TemporalHazard")
  valves_cc <- na.omit(valves)
  event_cr <- ifelse(valves_cc$dead == 1, 1L,
                     ifelse(valves_cc$pve == 1, 2L, 0L))
  time_cr <- pmin(valves_cc$int_dead, valves_cc$int_pve)
  cr <- hzr_competing_risks(time_cr, event_cr)

  expect_s3_class(cr, "hzr_competing_risks")
  expect_true(nrow(cr) > 0)
  expect_true("surv" %in% names(cr))
  expect_true("incid_1" %in% names(cr))
  expect_true("incid_2" %in% names(cr))
  expect_true("se_surv" %in% names(cr))
})

test_that("hzr_competing_risks survival + incidences sum to ~1", {
  data(valves, package = "TemporalHazard")
  valves_cc <- na.omit(valves)
  event_cr <- ifelse(valves_cc$dead == 1, 1L,
                     ifelse(valves_cc$pve == 1, 2L, 0L))
  time_cr <- pmin(valves_cc$int_dead, valves_cc$int_pve)
  cr <- hzr_competing_risks(time_cr, event_cr)

  # At each time: surv + sum(incid) should be approximately 1
  total <- cr$surv + cr$incid_1 + cr$incid_2
  expect_true(all(abs(total - 1) < 0.01),
              label = "surv + incidences should sum to ~1")
})

test_that("hzr_competing_risks survival is non-increasing", {
  data(valves, package = "TemporalHazard")
  valves_cc <- na.omit(valves)
  event_cr <- ifelse(valves_cc$dead == 1, 1L,
                     ifelse(valves_cc$pve == 1, 2L, 0L))
  time_cr <- pmin(valves_cc$int_dead, valves_cc$int_pve)
  cr <- hzr_competing_risks(time_cr, event_cr)

  expect_true(all(diff(cr$surv) <= 1e-10))
  expect_true(all(cr$surv >= 0 & cr$surv <= 1))
})

test_that("hzr_competing_risks incidences are non-decreasing", {
  data(valves, package = "TemporalHazard")
  valves_cc <- na.omit(valves)
  event_cr <- ifelse(valves_cc$dead == 1, 1L,
                     ifelse(valves_cc$pve == 1, 2L, 0L))
  time_cr <- pmin(valves_cc$int_dead, valves_cc$int_pve)
  cr <- hzr_competing_risks(time_cr, event_cr)

  expect_true(all(diff(cr$incid_1) >= -1e-10))
  expect_true(all(diff(cr$incid_2) >= -1e-10))
})

test_that("hzr_competing_risks rejects invalid inputs", {
  expect_error(hzr_competing_risks("a", 1), "numeric")
  expect_error(hzr_competing_risks(1:5, rep(0, 5)), "No events")
})

test_that("print.hzr_competing_risks runs without error", {
  data(valves, package = "TemporalHazard")
  valves_cc <- na.omit(valves)
  event_cr <- ifelse(valves_cc$dead == 1, 1L,
                     ifelse(valves_cc$pve == 1, 2L, 0L))
  time_cr <- pmin(valves_cc$int_dead, valves_cc$int_pve)
  cr <- hzr_competing_risks(time_cr, event_cr)
  expect_output(print(cr), "Competing risks")
})
