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

  # A single-distribution model built with fit = FALSE still predicts -- an
  # intended, tested capability -- but it says where the numbers come from.
  # The suite silences this warning wholesale (helper-unfitted-predictions.R);
  # switch it back on here, or the assertions below would have nothing to
  # catch.
  withr::local_options(TemporalHazard.warn_unfitted_prediction = TRUE)
  single <- eval_spec()
  expect_warning(p <- predict(single, type = "hazard"),
                 class = "hzr_unfitted_prediction")
  expect_true(all(is.finite(p)))
  expect_match(
    capture_warnings(predict(single, type = "hazard")),
    "come from the starting values"
  )
  # One warning per call, not one per row or per type.
  expect_length(capture_warnings(predict(single, type = "survival")), 1L)
  # A fitted model says nothing.
  expect_no_warning(predict(eval_spec(fit = TRUE), type = "hazard"))
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
  # hzr_phase() derives alpha at construction and says so.
  expect_warning(
    late <- hzr_phase("g3", tau = 5, gamma = 1, alpha = 1, eta = 1,
                      constraint = "alpha_gamma_eta"),
    "was replaced by"
  )
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    late = late
  )
  theta <- c(log(0.05), log(0.15), 1.4, 1, log(0.03), log(5), 1, 99, 1)
  # hazard() warns about the same replacement when it builds the object.
  expect_warning(
    spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                   dist = "multiphase", phases = phases, theta = theta,
                   fit = FALSE),
    "was replaced by"
  )
  # hazard() derived alpha = gamma * eta / 2 = 0.5, discarding the 99.
  expect_equal(spec$fit$theta[[8]], 0.5, tolerance = 1e-12)
  # The evaluation must agree with that, not with the 99, and it must warn
  # and report the vector it scored rather than the one handed in: the
  # parity user's SAS alpha will not be exactly gamma*eta/2 after rounding.
  expect_warning(ev <- hzr_evaluate(spec, theta = theta),
                 "the value its phase's constraint derives")
  expect_identical(ev$logLik,
                   suppressWarnings(
                     hzr_evaluate(spec,
                                  theta = unname(spec$fit$theta))$logLik))
  expect_equal(ev$theta[[8]], 0.5, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(ev$theta[[8]], 99)))
  # An unconstrained model is not warned about and is reported unchanged.
  plain <- list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                  fixed = "m"),
                constant = hzr_phase("constant"))
  th_plain <- c(log(0.05), log(0.15), 1.4, 1, log(0.03))
  spec2 <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                  dist = "multiphase", phases = plain, theta = th_plain,
                  fit = FALSE)
  expect_no_warning(ev2 <- hzr_evaluate(spec2, theta = th_plain))
  expect_identical(unname(ev2$theta), th_plain)
})

test_that("the curve matches predict() when the model has covariates (#144)", {
  skip_on_cran()
  # theta carries one slot per covariate per phase. Splitting it as though
  # the model had none read the constant phase's scale out of the early
  # phase's coefficient: 18 to 23 times wrong, monotone and finite, no
  # warning (r-reviewer). predict() at covariate 0 is the oracle.
  data("avc", package = "TemporalHazard", envir = environment())
  d <- na.omit(avc[, c("int_dead", "dead", "age")])
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = d, dist = "multiphase",
    phases = phases, fit = TRUE,
    control = list(n_starts = 1L, conserve = FALSE, maxit = 100)
  ))
  times <- c(1, 10)
  curve <- hzr_evaluate(fit, theta = fit$fit$theta, times = times)$curve
  at_zero <- data.frame(time = times, age = 0)
  expect_equal(curve$hazard,
               as.numeric(predict(fit, newdata = at_zero, type = "hazard")),
               tolerance = 1e-10)
  expect_equal(curve$cumulative_hazard,
               as.numeric(predict(fit, newdata = at_zero,
                                  type = "cumulative_hazard")),
               tolerance = 1e-10)
})

test_that("hzr_evaluate() refuses data whose rows were all dropped (#144)", {
  skip_on_cran()
  # A log-likelihood of 0 over no rows is the best value there is, and it
  # reads as a parity success.
  data("avc", package = "TemporalHazard", envir = environment())
  d <- avc[, c("int_dead", "dead", "age")]
  d$age <- NA_real_
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.2, nu = 1, m = -0.4, formula = ~ age),
    constant = hzr_phase("constant")
  )
  spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                 dist = "multiphase", phases = phases, fit = FALSE)
  theta <- c(log(0.05), log(0.2), 1, -0.4, 0, log(0.03))
  expect_error(hzr_evaluate(spec, theta = theta), "No rows are left")
})

test_that("the count error says which count is meant (#144)", {
  skip_on_cran()
  # hazard(fit = FALSE) does not resolve phase designs, so a specification
  # with covariates stores fewer parameters than a fit would use. Refusing
  # the object's own theta with a bare count is confusing; say why.
  data("avc", package = "TemporalHazard", envir = environment())
  d <- na.omit(avc[, c("int_dead", "dead", "age")])
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m",
                      formula = ~ age),
    constant = hzr_phase("constant")
  )
  spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                 dist = "multiphase", phases = phases,
                 theta = c(log(0.05), log(0.15), 1.4, 1, log(0.03)),
                 fit = FALSE)
  expect_length(spec$fit$theta, 5L)
  expect_error(hzr_evaluate(spec, theta = spec$fit$theta),
               "hazard\\(fit = FALSE\\) does not resolve phase designs")
  # The full vector is accepted and scores.
  full <- c(log(0.05), log(0.15), 1.4, 1, 0, log(0.03))
  expect_true(is.finite(hzr_evaluate(spec, theta = full)$logLik))
})

test_that("the refusals describe the model in front of them (#144)", {
  skip_on_cran()
  data("avc", package = "TemporalHazard", envir = environment())
  d <- na.omit(avc[, c("int_dead", "dead", "age")])
  x <- matrix(d$age, ncol = 1, dimnames = list(NULL, "age"))
  # A single-distribution model has no phases, and here the stored vector is
  # LONGER than the model's count: the multiphase explanation must not fire.
  w <- hazard(time = d$int_dead, status = d$dead, x = x, dist = "weibull",
              theta = c(0.05, 0.9, 0.01, 99), fit = FALSE)
  msg <- tryCatch(hzr_evaluate(w, theta = w$fit$theta),
                  error = conditionMessage)
  expect_match(msg, "has 4 parameters, but this weibull model has 3")
  expect_false(grepl("phase designs", msg))
  # A stored theta of the wrong length must not be pasted onto another
  # vector's names: that failed with base R's "'names' attribute [4] must be
  # the same length as the vector [3]". The correct theta evaluates, and its
  # names are simply not taken from the mismatched stored vector.
  wn <- hazard(time = d$int_dead, status = d$dead, x = x, dist = "weibull",
               theta = c(a = 0.05, b = 0.9, c = 0.01, d = 99), fit = FALSE)
  ev <- hzr_evaluate(wn, theta = c(0.05, 0.9, 0.01))
  expect_true(is.finite(ev$logLik))
  expect_null(names(ev$theta))

  # A model with no observations says that, rather than blaming a phase
  # design it does not have.
  z <- hazard(time = numeric(0), status = numeric(0), dist = "weibull",
              theta = c(0.05, 0.9), fit = FALSE)
  expect_error(hzr_evaluate(z, theta = c(0.05, 0.9)),
               "carries no observations")
})

test_that("the unfitted-prediction warning is on by DEFAULT (#144)", {
  # The suite switches this warning off wholesale
  # (helper-unfitted-predictions.R), and the tests that assert it switch it
  # back on for their own scope. Between those two, nothing would notice if
  # the shipped default flipped to FALSE: every test would stay green while
  # users stopped being told that a number came from a starting value. This
  # is the assertion that notices.
  withr::local_options(TemporalHazard.warn_unfitted_prediction = NULL)
  expect_null(getOption("TemporalHazard.warn_unfitted_prediction"))
  expect_warning(predict(eval_spec(), type = "hazard"),
                 class = "hzr_unfitted_prediction")

  # And the code's own fallback is TRUE, so "unset" means "warn" rather than
  # depending on something else having set it.
  src <- paste(deparse(args(predict.hazard)), collapse = " ")
  body_src <- paste(deparse(body(predict.hazard)), collapse = " ")
  body_src <- gsub("[[:space:]]+", " ", body_src)
  expect_match(
    body_src,
    'getOption("TemporalHazard.warn_unfitted_prediction", TRUE)',
    fixed = TRUE
  )
  expect_true(nzchar(src))
})

test_that("every single-distribution family predicts unfitted and warns once (#144)", {
  # The ruling for #144 is that only multiphase refuses: the other four
  # families keep predicting from supplied parameters and say where the
  # numbers came from. Until now that was asserted for weibull alone, via
  # eval_spec()'s default, even though eval_spec() takes a `dist` for
  # exactly this. The warning is raised on `converged` being unset rather
  # than on the family, so a family-specific regression would not be caught
  # by the weibull case.
  data("avc", package = "TemporalHazard", envir = environment())
  withr::local_options(TemporalHazard.warn_unfitted_prediction = TRUE)
  # exponential with an intercept-only formula takes ONE parameter; passing
  # two makes predict() look for a covariate and fail on its absence.
  thetas <- list(weibull = c(1, 1), exponential = 1,
                 lognormal = c(1, 1), loglogistic = c(1, 1))
  for (d in names(thetas)) {
    spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                   dist = d, theta = thetas[[d]], fit = FALSE)
    expect_warning(p <- predict(spec, type = "hazard"),
                   class = "hzr_unfitted_prediction",
                   info = d)
    expect_true(all(is.finite(p)), info = d)
    # One per call, for each family, not one per row.
    expect_length(capture_warnings(predict(spec, type = "hazard")), 1L)
  }
})

test_that("hzr_evaluate() does not nag that the model is unfitted (#144)", {
  # hzr_evaluate() IS the sanctioned way to score a model at supplied
  # parameters, so telling its caller that the model is not a fit is noise
  # about something they chose. It is quiet today because it computes the
  # curve itself rather than routing through predict(); this assertion is
  # what notices if it ever starts routing through predict() and inherits
  # the warning. The option is forced ON, or the suite-wide silencing in
  # helper-unfitted-predictions.R would make this pass over nothing.
  data("avc", package = "TemporalHazard", envir = environment())
  withr::local_options(TemporalHazard.warn_unfitted_prediction = TRUE)
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  theta <- c(log(0.05), log(0.15), 1.4, 1, log(0.03))
  spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = "multiphase", phases = phases, theta = theta,
                 fit = FALSE)
  ev <- expect_no_warning(hzr_evaluate(spec, theta = theta))
  expect_s3_class(ev, "hzr_evaluation")
  # With `times`, which is a DIFFERENT path: the curve is built only when
  # times is supplied, so an assertion that omits it leaves the branch a
  # refactor is most likely to route through predict() uncovered. A mutant
  # that warned inside .hzr_evaluate_curve() survived the times-free call.
  ev_t <- expect_no_warning(hzr_evaluate(spec, theta = theta, times = c(1, 5)))
  expect_false(is.null(ev_t$curve))
  # The known positive: predict() on the SAME unfitted object, with the
  # same option set, does warn. Without this the test above could pass
  # because nothing warns anywhere.
  single <- eval_spec()
  expect_warning(predict(single, type = "hazard"),
                 class = "hzr_unfitted_prediction")
})
