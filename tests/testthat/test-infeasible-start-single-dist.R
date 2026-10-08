# A single-distribution fit started where the likelihood is not defined
# (#486). The optimizer clamps a non-finite objective to 1e10, so from such a
# start it stops at once with convergence 0, and the fit read as converged
# with log-likelihood -1e10 and the start for estimates. Each family below
# reproduced that on main; each sane start is the known positive.

infeasible_486 <- list(
  exponential = 800,
  weibull     = c(1e10, 1e10),
  lognormal   = c(0, -800),
  loglogistic = c(0, 8)
)
sane_486 <- list(
  exponential = -3,
  weibull     = c(0.1, 1),
  lognormal   = c(1, 0.5),
  loglogistic = c(-1, 1)
)

fit_486 <- function(dist, theta) {
  data(avc, package = "TemporalHazard", envir = environment())
  classes <- character()
  fit <- withCallingHandlers(
    hazard(survival::Surv(int_dead, dead) ~ 1, data = avc, dist = dist,
           theta = theta, fit = TRUE),
    warning = function(w) {
      classes <<- c(classes, class(w)[1L])
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, classes = classes,
       printed = utils::capture.output(print(fit)),
       summarised = utils::capture.output(print(summary(fit))))
}

expect_infeasible_486 <- function(dist) {
  res <- fit_486(dist, infeasible_486[[dist]])
  f <- res$fit
  expect_identical(f$fit$converged, FALSE)
  expect_identical(f$fit$objective, NA_real_)
  expect_identical(sum(res$classes == "hzr_infeasible_start"), 1L)
  # print() shows no log-likelihood line without an objective; the record
  # names the cause. summary() shows the convergence flag.
  expect_match(res$printed,
               "standard_errors: the optimizer ended where the likelihood is not defined",
               fixed = TRUE, all = FALSE)
  expect_match(res$summarised, "converged:\\s+FALSE", all = FALSE)
  shown <- c(res$printed, res$summarised)
  expect_false(any(grepl("log-lik", shown, fixed = TRUE)))
  expect_false(any(grepl("-1e+10", shown, fixed = TRUE)))
}

expect_sane_486 <- function(dist) {
  res <- fit_486(dist, sane_486[[dist]])
  expect_identical(res$fit$fit$converged, TRUE)
  expect_true(is.finite(res$fit$fit$objective))
  expect_lt(abs(res$fit$fit$objective), 1e10)
  expect_identical(sum(res$classes == "hzr_infeasible_start"), 0L)
  # The probe reads real output: a converged fit prints both lines.
  expect_match(res$summarised, "converged:\\s+TRUE", all = FALSE)
  expect_match(res$printed, "log-lik:", fixed = TRUE, all = FALSE)
}

test_that("an exponential fit started at the clamp is not converged (#486)", {
  expect_infeasible_486("exponential")
  expect_sane_486("exponential")
})

test_that("a weibull fit started at the clamp is not converged (#486)", {
  expect_infeasible_486("weibull")
  expect_sane_486("weibull")
})

test_that("a lognormal fit started at the clamp is not converged (#486)", {
  expect_infeasible_486("lognormal")
  expect_sane_486("lognormal")
})

test_that("a loglogistic fit started at the clamp is not converged (#486)", {
  expect_infeasible_486("loglogistic")
  expect_sane_486("loglogistic")
})

test_that("a multiphase start at the clamp is still recorded as infeasible (#486)", {
  # The multiphase path opts out of the single-distribution marking: it asks
  # the likelihood itself and records the start in `starts`, and with one
  # start and none usable it stops, naming why.
  data(avc, package = "TemporalHazard", envir = environment())
  expect_error(
    suppressWarnings(hazard(
      survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
      phases = list(constant = hzr_phase("constant")),
      theta = 800, fit = TRUE, control = list(n_starts = 1L)
    )),
    "1 ended where the likelihood is not defined, 0 returned a non-finite",
    fixed = TRUE
  )
})
