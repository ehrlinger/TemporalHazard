# Evaluating a model at parameters supplied from elsewhere -- SAS's converged
# estimates, say -- is the core operation of a parity check, and there was no
# way to ask for it: predict() on a fit = FALSE object died with "argument of
# length 0" from .hzr_split_theta (#144). predict() still refuses, by name,
# and hzr_evaluate() is the supported route.

eval_spec <- function(fit = FALSE, dist = "weibull") {
  data("avc", package = "TemporalHazard", envir = environment())
  hazard(survival::Surv(int_dead, dead) ~ 1, data = avc, dist = dist,
         theta = c(1, 1), fit = fit)
}

test_that("predict() names what an unfitted multiphase model lacks (#144)", {
  # Only multiphase fails: its per-phase designs are resolved at fit time.
  # The other families predict from supplied parameters perfectly well, and
  # the suite has long relied on that, so they are left alone.
  data("avc", package = "TemporalHazard", envir = environment())
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  theta <- c(log(0.05), log(0.15), 1.4, 1, log(0.03))
  spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = "multiphase", phases = phases, theta = theta,
                 fit = FALSE)
  expect_error(predict(spec, newdata = data.frame(time = 1), type = "hazard"),
               "no per-phase design matrices")
  expect_error(predict(spec, type = "hazard"), "hzr_evaluate\\(\\)")
  # Not the old failure, which named an internal helper's argument length.
  expect_false(grepl("argument of length 0",
                     tryCatch(predict(spec, type = "hazard"),
                              error = conditionMessage)))

  # A single-distribution model built with fit = FALSE still predicts.
  single <- eval_spec()
  expect_true(all(is.finite(predict(single, type = "hazard"))))
})

test_that("a fitted model still predicts (#144)", {
  fit <- eval_spec(fit = TRUE)
  expect_true(isTRUE(fit$fit$converged))
  expect_true(all(is.finite(predict(fit, type = "hazard"))))
  expect_true(all(is.finite(
    predict(fit, newdata = data.frame(time = c(1, 2)), type = "survival"))))
})

test_that("hzr_evaluate() scores the model's own likelihood (#144)", {
  # The oracle: at a fit's own estimates it must reproduce the objective the
  # optimizer reported, since it is the same likelihood.
  fit <- eval_spec(fit = TRUE)
  ev <- hzr_evaluate(fit, theta = fit$fit$theta)
  expect_equal(ev$logLik, fit$fit$objective, tolerance = 1e-12)
  expect_s3_class(ev, "hzr_evaluation")
  expect_identical(ev$dist, "weibull")
  expect_identical(ev$n_obs, length(fit$data$time))
  expect_identical(ev$n_events, sum(fit$data$status == 1))

  # An unfitted object evaluates too, and elsewhere in the parameter space
  # the likelihood is lower than at the optimum.
  spec <- eval_spec()
  away <- hzr_evaluate(spec, theta = fit$fit$theta + c(0.3, 0.3))
  expect_lt(away$logLik, ev$logLik)
  expect_identical(hzr_evaluate(spec, theta = fit$fit$theta)$logLik,
                   ev$logLik)
})

test_that("hzr_evaluate() scores a multiphase model at supplied parameters (#144)", {
  skip_on_cran()
  data("avc", package = "TemporalHazard", envir = environment())
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  fit <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1,
                                 data = avc, dist = "multiphase",
                                 phases = phases, fit = TRUE,
                                 control = list(n_starts = 1L,
                                                conserve = FALSE)))
  ev <- hzr_evaluate(fit, theta = fit$fit$theta)
  expect_equal(ev$logLik, fit$fit$objective, tolerance = 1e-10)

  # The unfitted specification gives the same value at the same parameters:
  # the designs come from the function the optimizer itself uses.
  spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = "multiphase", phases = phases,
                 theta = unname(fit$fit$theta), fit = FALSE)
  expect_equal(hzr_evaluate(spec, theta = fit$fit$theta)$logLik, ev$logLik,
               tolerance = 1e-12)
})

test_that("hzr_evaluate() returns the multiphase shape at supplied times (#144)", {
  data("avc", package = "TemporalHazard", envir = environment())
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  theta <- c(log(0.05), log(0.15), 1.4, 1, log(0.03))
  spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = "multiphase", phases = phases, theta = theta,
                 fit = FALSE)
  ev <- hzr_evaluate(spec, theta = theta, times = c(0.5, 2, 5))
  expect_identical(ev$curve$time, c(0.5, 2, 5))
  expect_true(all(is.finite(ev$curve$hazard)))
  expect_true(all(diff(ev$curve$cumulative_hazard) > 0))
  # The early phase is still declining towards the constant floor, so the
  # curve is the shape these parameters describe, not any fitted shape.
  expect_true(all(diff(ev$curve$hazard) < 0))
  expect_gt(ev$curve$hazard[1], exp(theta[5]))

  # The families without an internal shape function refuse rather than
  # duplicating predict()'s formulas.
  expect_error(hzr_evaluate(eval_spec(), theta = c(0.05, 0.9),
                            times = c(1, 2)),
               "supported for dist = .multiphase. only")
})

test_that("hzr_evaluate() refuses what it cannot evaluate (#144)", {
  spec <- eval_spec()
  expect_error(hzr_evaluate(spec, theta = c(1, 1, 1)), "has 3 parameters")
  expect_error(hzr_evaluate(spec, theta = c(NA_real_, 1)), "finite numeric")
  expect_error(hzr_evaluate(spec, theta = c(1, 1), times = -1),
               "non-negative times")
  expect_error(hzr_evaluate(structure(list(), class = "hazard"),
                            theta = c(1, 1)), "carries no data")
  expect_error(hzr_evaluate(list(a = 1), theta = c(1, 1)),
               "must be a 'hazard' object")
  # Only a multiphase model names its parameters, so that is where names can
  # be checked at all; the other families carry an unnamed theta.
  data("avc", package = "TemporalHazard", envir = environment())
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  theta <- c(log(0.05), log(0.15), 1.4, 1, log(0.03))
  mp <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
               dist = "multiphase", phases = phases, theta = theta,
               fit = FALSE)
  wrong <- stats::setNames(theta, paste0("wrong", seq_along(theta)))
  expect_error(hzr_evaluate(mp, theta = wrong), "not the model's parameter")
  # The model's own names are accepted.
  right <- stats::setNames(theta, hzr_theta_names(phases))
  expect_identical(hzr_evaluate(mp, theta = right)$logLik,
                   hzr_evaluate(mp, theta = unname(theta))$logLik)
})

test_that("printing an evaluation says it is not a fit (#144)", {
  spec <- eval_spec()
  out <- capture.output(print(hzr_evaluate(spec, theta = c(0.05, 0.9))))
  expect_match(out[1], "SUPPLIED parameters -- not a fit")
  expect_true(any(grepl("no standard errors, no convergence", out)))
  expect_true(any(grepl("logLik at the supplied parameters", out)))
})

test_that("hzr_evaluate() scores the covariates, not just the shape (#144)", {
  # Dropping the design would still return a number, and a plausible one:
  # the same model without its covariates. The oracle catches it because the
  # fit's own objective includes them.
  data("avc", package = "TemporalHazard", envir = environment())
  d <- na.omit(avc[, c("int_dead", "dead", "age")])
  fit <- hazard(survival::Surv(int_dead, dead) ~ age, data = d,
                dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  expect_false(isTRUE(all.equal(unname(fit$fit$theta)[3], 0)))
  expect_equal(hzr_evaluate(fit, theta = fit$fit$theta)$logLik,
               fit$fit$objective, tolerance = 1e-12)

  # And with the covariate coefficient zeroed the likelihood is different,
  # so the design is genuinely being used.
  no_beta <- fit$fit$theta
  no_beta[3] <- 0
  expect_false(isTRUE(all.equal(hzr_evaluate(fit, theta = no_beta)$logLik,
                                fit$fit$objective)))
})
