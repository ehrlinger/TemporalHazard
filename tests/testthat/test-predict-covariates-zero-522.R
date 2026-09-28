# tests/testthat/test-predict-covariates-zero-522.R
# predict(newdata = <time only>) on a model WITH covariates evaluates every
# covariate at 0 -- the baseline, age 0 for avc -- and said nothing (#522).
# It now warns, once per call, with class "hzr_predict_covariates_zero",
# naming the covariates. A model without covariates does not warn.

.cz_avc <- local({
  data(avc, package = "TemporalHazard")
  stats::na.omit(avc)
})
.cz_nd <- data.frame(time = c(12, 60))

# Every warning a call gives, so a test can count them and read their class.
.cz_warnings <- function(expr) {
  ws <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = ws)
}

.cz_only_zero_warning <- function(res, covariates) {
  expect_length(res$warnings, 1L)
  w <- res$warnings[[1L]]
  expect_s3_class(w, "hzr_predict_covariates_zero")
  msg <- conditionMessage(w)
  expect_match(msg, "'newdata' has no covariate column", fixed = TRUE)
  expect_match(msg, "at 0", fixed = TRUE)
  for (cv in covariates) expect_match(msg, paste0("'", cv, "'"), fixed = TRUE)
}

test_that("single-distribution formula fit: time-only newdata warns (#522)", {
  d <- .cz_avc
  f <- hazard(survival::Surv(int_dead, dead) ~ age + mal, data = d,
              dist = "weibull", theta = c(0.1, 1, 0, 0), fit = TRUE)
  expect_true(isTRUE(f$fit$converged))
  for (ty in c("survival", "cumulative_hazard")) {
    res <- .cz_warnings(predict(f, newdata = .cz_nd, type = ty))
    .cz_only_zero_warning(res, c("age", "mal"))
    # The values are unchanged: exactly the age = 0, mal = 0 prediction.
    at_zero <- predict(f, newdata = data.frame(time = .cz_nd$time,
                                               age = 0, mal = 0), type = ty)
    expect_identical(res$value, at_zero)
  }
  # Supplying the covariates is silent (the warning is not unconditional).
  full <- .cz_warnings(predict(
    f, newdata = data.frame(time = .cz_nd$time, age = mean(d$age),
                            mal = mean(d$mal)),
    type = "survival"
  ))
  expect_length(full$warnings, 0L)
})

test_that("single-distribution x-matrix fit: time-only newdata warns (#522)", {
  d <- .cz_avc
  fx <- hazard(time = d$int_dead, status = d$dead,
               x = as.matrix(d[, c("age", "mal")]), dist = "weibull",
               theta = c(0.1, 1, 0, 0), fit = TRUE)
  res <- .cz_warnings(predict(fx, newdata = .cz_nd,
                              type = "cumulative_hazard"))
  .cz_only_zero_warning(res, c("age", "mal"))
})

test_that("known positive: a model without covariates does not warn (#522)", {
  f0 <- hazard(survival::Surv(int_dead, dead) ~ 1, data = .cz_avc,
               dist = "weibull", theta = c(0.1, 1), fit = TRUE)
  res <- .cz_warnings(predict(f0, newdata = .cz_nd, type = "survival"))
  expect_length(res$warnings, 0L)
  expect_length(res$value, 2L)
  expect_true(all(res$value > 0 & res$value < 1))
})

test_that("multiphase with phase covariates: time-only newdata warns (#522)", {
  skip_on_cran()
  d <- .cz_avc
  fm <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ age + mal),
      constant = hzr_phase("constant", formula = ~ age)
    ),
    fit = TRUE
  ))
  expect_true(isTRUE(fm$fit$converged))
  for (ty in c("survival", "cumulative_hazard", "hazard")) {
    res <- .cz_warnings(predict(fm, newdata = .cz_nd, type = ty))
    # Once per call, though two phases carry covariates.
    .cz_only_zero_warning(res, c("age", "mal"))
    at_zero <- predict(fm, newdata = data.frame(time = .cz_nd$time,
                                                age = 0, mal = 0), type = ty)
    expect_identical(res$value, at_zero)
  }
  full <- .cz_warnings(predict(
    fm, newdata = data.frame(time = .cz_nd$time, age = mean(d$age),
                             mal = mean(d$mal)),
    type = "survival"
  ))
  expect_length(full$warnings, 0L)
})

test_that("multiphase with global covariates: time-only newdata warns (#522)", {
  skip_on_cran()
  fg <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = .cz_avc,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  res <- .cz_warnings(predict(fg, newdata = .cz_nd, type = "survival"))
  .cz_only_zero_warning(res, "age")
})

test_that("known positive: a multiphase model without covariates does not warn", {
  skip_on_cran()
  fm0 <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = .cz_avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  )
  res <- .cz_warnings(predict(fm0, newdata = .cz_nd, type = "survival"))
  expect_length(res$warnings, 0L)
  expect_length(res$value, 2L)
  expect_true(all(res$value > 0 & res$value < 1))
})
