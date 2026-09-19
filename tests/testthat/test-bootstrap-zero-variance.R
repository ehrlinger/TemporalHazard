# hzr_bootstrap() and replicates that did not estimate anything (#373).
#
# The route #373 measured: a start at theta = 1e10 leaves a single-
# distribution fit at the optimizer's -1e10 sentinel, which stands in for a
# likelihood that could not be evaluated, and every replicate reproduced it:
# n_success = 5, n_failed = 0, sd exactly 0, no warning. Two guards:
# a replicate at the sentinel is a failed replicate, and a free parameter
# that does not move across replicates is named in a warning. A parameter
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

zero_var_warnings <- function(w) grep("(sd = 0", w, value = TRUE, fixed = TRUE)

run_373 <- function(expr) {
  w <- character()
  res <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(res = res, w = w)
}

sentinel_reason_373 <-
  "objective at the optimizer's -1e10 sentinel (no log-likelihood)"

test_that("a replicate at the -1e10 sentinel is a failed replicate (#373)", {
  d <- zv_data_373()
  bad <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
    theta = c(1e10, 1e10), fit = TRUE
  ))
  # The mechanism: an objective at the sentinel, not a log-likelihood.
  expect_identical(bad$fit$objective, -1e10)
  out <- run_373(hzr_bootstrap(bad, n_boot = 5L, seed = 1L))
  expect_identical(out$res$n_success, 0L)
  expect_identical(out$res$failure_reasons,
                   stats::setNames(5L, sentinel_reason_373))
  expect_match(out$w, "no replicate succeeded", fixed = TRUE, all = FALSE)
})

test_that("replicates that drift before reaching the sentinel fail too (#373)", {
  # From theta = 20 the optimizer moves before the clamp stops it, so the
  # stuck replicates differ and an sd test alone would pass them. Measured:
  # four of five replicates end at the sentinel.
  d <- zv_data_373()
  drift <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
    theta = c(20, 20), fit = TRUE
  ))
  out <- run_373(hzr_bootstrap(drift, n_boot = 5L, seed = 1L))
  expect_identical(out$res$n_success, 1L)
  expect_identical(out$res$failure_reasons,
                   stats::setNames(4L, sentinel_reason_373))
})

test_that("identical replicates are named, to within rounding (#373)", {
  # From theta = 50 the objective is finite (about -3.6e196) and not the
  # sentinel, and every replicate stays put up to the last bits: sd is about
  # 1e-14 around a mean of 50. An exact-zero test would pass it.
  d <- zv_data_373()
  stuck <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
    theta = c(50, 50), fit = TRUE
  ))
  out <- run_373(hzr_bootstrap(stuck, n_boot = 5L, seed = 1L))
  expect_identical(out$res$n_success, 5L)
  expect_true(all(out$res$summary$sd > 0))
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

test_that("identical single-distribution replicates are named (#373)", {
  # fixed_mask is logical(0) on a single-distribution fit, so nothing is
  # exempt and both parameters are named.
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

test_that("a free parameter at exactly 0 in every replicate is named (#373)", {
  # The tolerance is relative to the mean, so at a mean of 0 it is sd == 0:
  # intended, since a free estimate that is exactly 0 in every replicate
  # did not move either.
  d <- zv_data_373()
  ok <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
               dist = "weibull", theta = c(0.1, 1), fit = TRUE)
  zero <- ok
  zero$fit$theta[2L] <- 0
  out <- run_373(hzr_bootstrap(same_refit_373(ok, zero), n_boot = 3L,
                               seed = 1L))
  s <- out$res$summary
  expect_identical(s$mean[s$parameter == "param_2"], 0)
  zw <- zero_var_warnings(out$w)
  expect_length(zw, 1L)
  expect_match(zw, "`param_1`, `param_2`", fixed = TRUE)
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
