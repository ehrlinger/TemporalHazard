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
               "finite positive times")
  # 0 is refused too: every phase's cumulative hazard is 0 there, and a
  # shape's hazard at 0 need not be finite.
  expect_error(hzr_evaluate(spec, theta = c(1, 1), times = 0),
               "finite positive times")
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

test_that("hzr_evaluate() scores the window-expanded design (#144)", {
  skip_on_cran()
  # hazard() expands the design for time_windows before fitting, and the
  # object stores the UNEXPANDED x. Scoring that gave -26390.82 where the
  # fit reported -191.83, with no warning, because a 310x1 design against a
  # length-2 beta recycles cleanly (r-reviewer).
  data("avc", package = "TemporalHazard", envir = environment())
  d <- na.omit(avc[, c("int_dead", "dead", "age")])
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = d, dist = "weibull",
    theta = c(0.05, 1, 0, 0), time_windows = 1, fit = TRUE,
    control = list(n_starts = 1L)
  ))
  expect_length(fit$fit$theta, 4L)
  expect_equal(hzr_evaluate(fit, theta = fit$fit$theta)$logLik,
               fit$fit$objective, tolerance = 1e-10)

  # And multiphase, where the wrong design also splits theta at the wrong
  # offsets while the length check still passes.
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  mp <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = d, dist = "multiphase",
    phases = phases, time_windows = 1, fit = TRUE,
    control = list(n_starts = 1L, conserve = FALSE, maxit = 200)
  ))
  expect_equal(hzr_evaluate(mp, theta = mp$fit$theta)$logLik,
               mp$fit$objective, tolerance = 1e-8)
})

test_that("hzr_evaluate() counts the rows it scored (#144)", {
  skip_on_cran()
  # A phase design with an NA drops rows from the likelihood. Reporting the
  # object's row count beside that log-likelihood is a wrong denominator.
  data("avc", package = "TemporalHazard", envir = environment())
  d <- avc
  d$age[1:40] <- NA
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.2, nu = 1, m = -0.4, formula = ~ age),
    constant = hzr_phase("constant")
  )
  theta <- c(log(0.05), log(0.2), 1, -0.4, 0, log(0.03))
  spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                 dist = "multiphase", phases = phases, theta = theta,
                 fit = FALSE)
  ev <- hzr_evaluate(spec, theta = theta)
  kept <- !is.na(d$age)
  expect_identical(ev$n_obs, sum(kept))
  expect_identical(ev$n_events, sum(d$dead[kept] == 1))
  expect_lt(ev$n_obs, nrow(d))
})

test_that("hzr_evaluate() checks theta against the model, not against itself (#144)", {
  skip_on_cran()
  # An unfitted object built without `theta` has no fit$theta to compare
  # against, and that is exactly the object this feature exists for: the
  # count must come from the model's own parameter list.
  data("avc", package = "TemporalHazard", envir = environment())
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.2, nu = 1, m = -0.4),
    constant = hzr_phase("constant")
  )
  spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = "multiphase", phases = phases, fit = FALSE)
  expect_null(spec$fit$theta)
  theta <- c(log(0.05), log(0.2), 1, -0.4, log(0.03))
  expect_true(is.finite(hzr_evaluate(spec, theta = theta)$logLik))
  expect_error(hzr_evaluate(spec, theta = c(theta, 99, -99)),
               "has 7 parameters, but this multiphase model has 5")
  expect_error(hzr_evaluate(spec, theta = theta[1:4]), "has 4 parameters")

  # With covariates, the model's names include them, so a correctly named
  # theta is accepted rather than refused.
  d <- na.omit(avc[, c("int_dead", "dead", "age")])
  spec_cov <- hazard(survival::Surv(int_dead, dead) ~ age, data = d,
                     dist = "multiphase", phases = phases, fit = FALSE)
  nm <- hzr_theta_names(phases, covariates = list(early = "age",
                                                  constant = "age"))
  theta_cov <- c(log(0.05), log(0.2), 1, -0.4, 0, log(0.03), 0)
  expect_identical(
    hzr_evaluate(spec_cov, theta = stats::setNames(theta_cov, nm))$logLik,
    hzr_evaluate(spec_cov, theta = theta_cov)$logLik
  )
})

test_that("predict() still gives multiphase's own reason for linear_predictor (#144)", {
  data("avc", package = "TemporalHazard", envir = environment())
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  theta <- c(log(0.05), log(0.15), 1.4, 1, log(0.03))
  spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = "multiphase", phases = phases, theta = theta,
                 fit = FALSE)
  # The unfitted guard must not mask a refusal that is true of fitted models
  # too, and whose remedy is not "fit it".
  expect_error(predict(spec, type = "linear_predictor"),
               "linear_predictor")
  expect_false(grepl("fit = TRUE",
                     tryCatch(predict(spec, type = "linear_predictor"),
                              error = conditionMessage)))
})

test_that("hzr_evaluate() counts parameters without column names (#144)", {
  skip_on_cran()
  # The vector interface -- the parity interface -- has an unnamed design.
  # Deriving names from colnames() then fell back to the stored theta's
  # length, which is the vacuous check again: an 8-vector was accepted for a
  # 7-parameter model and its last entry silently dropped (r-reviewer).
  data("avc", package = "TemporalHazard", envir = environment())
  d <- na.omit(avc[, c("int_dead", "dead", "age")])
  x <- matrix(as.numeric(d$age), ncol = 1)
  expect_null(colnames(x))
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  theta7 <- c(log(0.05), log(0.15), 1.4, 1, 0, log(0.03), 0)
  spec <- hazard(time = d$int_dead, status = d$dead, x = x,
                 dist = "multiphase", phases = phases, theta = theta7,
                 fit = FALSE)
  expect_true(is.finite(hzr_evaluate(spec, theta = theta7)$logLik))
  expect_error(hzr_evaluate(spec, theta = c(theta7, 99)),
               "has 8 parameters, but this multiphase model has 7")
  # And the synthesised names are the fit's own, so a named theta works.
  nm <- names(hzr_evaluate(spec, theta = theta7)$theta)
  expect_length(nm, 7L)
  expect_identical(
    hzr_evaluate(spec, theta = stats::setNames(theta7, nm))$logLik,
    hzr_evaluate(spec, theta = theta7)$logLik
  )
  expect_error(hzr_evaluate(spec, theta = stats::setNames(theta7,
                                                          paste0("junk", 1:7))),
               "not the model's parameter names")
})

test_that("hzr_evaluate() applies a phase constraint to the supplied theta (#144)", {
  skip_on_cran()
  # A derived shape is a function of its sources. The fit re-applies the rule
  # on every evaluation, so honouring a contradictory value here would score
  # a model the package cannot fit (r-reviewer).
  data("avc", package = "TemporalHazard", envir = environment())
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    late = hzr_phase("g3", tau = 5, gamma = 1, alpha = 1, eta = 1,
                     constraint = "alpha_gamma_eta")
  )
  theta <- c(log(0.05), log(0.15), 1.4, 1, log(0.03), log(5), 1, 99, 1)
  spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = "multiphase", phases = phases, theta = theta,
                 fit = FALSE)
  # hazard() derived alpha = gamma * eta / 2 = 0.5, discarding the 99.
  expect_equal(spec$fit$theta[[8]], 0.5, tolerance = 1e-12)
  # The evaluation must agree with that, not with the 99.
  expect_identical(hzr_evaluate(spec, theta = theta)$logLik,
                   hzr_evaluate(spec, theta = unname(spec$fit$theta))$logLik)
})
