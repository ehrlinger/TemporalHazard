# hzr_bootstrap() warns when a FREE parameter is identical in every
# successful replicate (#373).
#
# Resampled data moves every estimate a replicate actually makes, so an sd of
# exactly 0 on a free parameter means the replicates did not estimate it. The
# route #373 measured: a start at theta = 1e10 leaves a single-distribution
# fit at the optimizer's -1e10 sentinel, and every replicate reproduces it.
# That run reported n_success = 5, n_failed = 0 and no warning. A parameter
# held by `fixed =` is identical by design and must not be named.

zv_data_373 <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc[, c("int_dead", "dead", "age")])
}

# A replicate refit that returns `value` unchanged: the identical-fits shape
# of the 500-replicate incident in AGENTS.md, reached without the optimizer.
same_refit_373 <- function(fit, value) {
  env <- new.env(parent = fit$call_env %||% globalenv())
  assign("same_refit", function(...) value, envir = env)
  fit$call[[1L]] <- as.name("same_refit")
  fit$call_env <- env
  fit
}

mp_fit_373 <- function(d) {
  suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = list(n_starts = 1L, maxit = 300L)
  ))
}

zero_var_warnings <- function(w) grep("(sd = 0)", w, value = TRUE, fixed = TRUE)

run_373 <- function(expr) {
  w <- character()
  res <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(res = res, w = w)
}

test_that("a sentinel fit's replicates are named, not counted clean (#373)", {
  d <- zv_data_373()
  bad <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
    theta = c(1e10, 1e10), fit = TRUE
  ))
  # The mechanism: an objective at the sentinel, not a log-likelihood.
  expect_identical(bad$fit$objective, -1e10)
  out <- run_373(hzr_bootstrap(bad, n_boot = 5L, seed = 1L))
  expect_identical(out$res$n_success, 5L)
  expect_identical(out$res$summary$sd, c(0, 0))
  zw <- zero_var_warnings(out$w)
  expect_length(zw, 1L)
  expect_match(zw, "`param_1`, `param_2`", fixed = TRUE)
  expect_match(zw, "all 5 successful replicates", fixed = TRUE)
})

test_that("a working single-distribution bootstrap does not warn (#373)", {
  d <- zv_data_373()
  ok <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
               dist = "weibull", theta = c(0.1, 1), fit = TRUE)
  out <- run_373(hzr_bootstrap(ok, n_boot = 5L, seed = 1L))
  expect_identical(out$res$n_success, 5L)
  expect_true(all(out$res$summary$sd > 0))
  expect_length(zero_var_warnings(out$w), 0L)
})

test_that("a single-distribution fit's empty fixed_mask means nothing is fixed (#373)", {
  # fixed_mask is logical(0) on a single-distribution fit. Indexing theta
  # by an empty mask would check no parameter at all: the hollow shape this
  # guard exists to catch, reproduced inside it.
  d <- zv_data_373()
  ok <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
               dist = "weibull", theta = c(0.1, 1), fit = TRUE)
  expect_length(ok$fit$fixed_mask, 0L)
  out <- run_373(hzr_bootstrap(same_refit_373(ok, ok), n_boot = 3L,
                               seed = 1L))
  zw <- zero_var_warnings(out$w)
  expect_length(zw, 1L)
  expect_match(zw, "`param_1`, `param_2`", fixed = TRUE)
})

test_that("multiphase: free parameters are named, fixed ones are not (#373)", {
  skip_on_cran() # a multiphase fit
  d <- zv_data_373()
  mp <- mp_fit_373(d)
  # Known positive for the exemption: three shapes are fixed, and they sit in
  # theta beside the two free scales.
  expect_identical(mp$fit$fixed_mask, c(FALSE, TRUE, TRUE, TRUE, FALSE))
  out <- run_373(hzr_bootstrap(same_refit_373(mp, mp), n_boot = 3L,
                               seed = 1L))
  zw <- zero_var_warnings(out$w)
  expect_length(zw, 1L)
  expect_match(zw, "`early.log_mu`, `constant.log_mu`", fixed = TRUE)
  expect_no_match(zw, "log_t_half|early[.]nu|early[.]m`")
})

test_that("multiphase: a working bootstrap does not warn on its fixed shapes (#373)", {
  skip_on_cran() # a multiphase bootstrap
  d <- zv_data_373()
  out <- run_373(hzr_bootstrap(mp_fit_373(d), n_boot = 4L, seed = 1L))
  s <- out$res$summary
  # The fixed shapes are identical across replicates, as they must be ...
  expect_identical(s$sd[s$parameter == "early.nu"], 0)
  # ... and are not reported.
  expect_length(zero_var_warnings(out$w), 0L)
})

test_that("a parameter in only one replicate is not reported (#373)", {
  # With one replicate there is no spread to test; sd is NA, not 0.
  d <- zv_data_373()
  ok <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
               dist = "weibull", theta = c(0.1, 1), fit = TRUE)
  out <- run_373(hzr_bootstrap(same_refit_373(ok, ok), n_boot = 1L,
                               seed = 1L))
  expect_identical(out$res$n_success, 1L)
  expect_length(zero_var_warnings(out$w), 0L)
})
