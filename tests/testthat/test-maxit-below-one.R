# control$maxit below 1 was accepted without a word (#541). On some models
# the fit returned its starting values as `converged = TRUE`; on others the
# limit was ignored and the model optimised anyway. Neither said anything.
# A maxit that is not a whole number of at least 1 is now refused.

maxit_avc <- function() {
  data("avc", package = "TemporalHazard", envir = environment())
  stats::na.omit(avc)
}

fit_weibull <- function(maxit, fit = TRUE) {
  hazard(survival::Surv(int_dead, dead) ~ age, data = maxit_avc(),
         dist = "weibull", theta = c(0.1, 1, 0), fit = fit,
         control = list(maxit = maxit))
}

test_that("a maxit below 1 is refused, not fitted as a converged start (#541)", {
  # Known positive: a positive maxit fits, and maxit = 1 stops short.
  full <- fit_weibull(100)
  expect_true(full$fit$converged)
  one <- fit_weibull(1)
  expect_false(one$fit$converged)
  expect_lt(one$fit$objective, full$fit$objective)

  for (m in list(0, -1, 0L, -5L)) {
    expect_error(fit_weibull(m), "control\\$maxit", info = format(m))
  }
  msg <- tryCatch(fit_weibull(0), error = conditionMessage)
  # It says what to use instead of a zero-iteration "fit".
  expect_match(msg, "hzr_evaluate()", fixed = TRUE)
})

test_that("a maxit that is not one whole number is refused (#541)", {
  for (m in list(0.5, 2.5, NA_real_, Inf, "a", c(10, 20), numeric(0))) {
    expect_error(fit_weibull(m), "control\\$maxit", info = format(m))
  }
  # A whole number stored as a double is accepted.
  expect_true(fit_weibull(100)$fit$converged)
})

test_that("the refusal holds for a model that used to ignore maxit (#541)", {
  skip_on_cran()
  phases <- list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1),
                 constant = hzr_phase("constant"))
  mp <- function(maxit) {
    hazard(survival::Surv(int_dead, dead) ~ 1, data = maxit_avc(),
           dist = "multiphase", phases = phases, fit = TRUE,
           control = list(maxit = maxit, n_starts = 1))
  }
  expect_error(mp(0), "control\\$maxit")
  expect_error(mp(-1), "control\\$maxit")
})

test_that("hzr_stepwise() refuses a bad maxit once, before any refit (#541)", {
  skip_on_cran()
  base <- fit_weibull(100)
  expect_error(
    withCallingHandlers(
      hzr_stepwise(base, scope = c("mal", "nyha"), control = list(maxit = 0)),
      warning = function(w) invokeRestart("muffleWarning")),
    "control\\$maxit")
})

test_that("an unfitted specification with a bad maxit is refused too (#541)", {
  # A spec is refitted later by hzr_bootstrap() and hzr_stepwise(); refusing
  # it here keeps the refusal out of their candidate refits.
  expect_error(fit_weibull(0, fit = FALSE), "control\\$maxit")
})
