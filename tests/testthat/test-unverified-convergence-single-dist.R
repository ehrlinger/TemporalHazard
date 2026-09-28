# A single-distribution fit from a far start could report converged = TRUE
# over a point that is not a maximum (#518). optim()'s BFGS stops on the
# relative change in the objective, and the gradient it follows is zeroed
# wherever the score is not finite. On main 50a5ddd1:
#
# - loglogistic from c(-1e5, 1): converged TRUE, log-likelihood -4727490, the
#   score non-finite at the returned point; only Hessian warnings.
# - lognormal from c(1e5, 1): converged TRUE, log-likelihood -819.29 (the
#   maximum is -228.66), relative gradient 0.085 against SAS/C's 6.06e-06,
#   nlm code 2, which does not warn.
#
# The sane starts are the known positives: each family converges to the same
# log-likelihood it did before. Exponential from 50 must still reach -434.29.

fit_518 <- function(dist, theta) {
  data(avc, package = "TemporalHazard", envir = environment())
  classes <- list()
  fit <- withCallingHandlers(
    hazard(survival::Surv(int_dead, dead) ~ 1, data = avc, dist = dist,
           theta = theta, fit = TRUE),
    warning = function(w) {
      classes[[length(classes) + 1L]] <<- class(w)
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, classes = classes,
       summarised = utils::capture.output(print(summary(fit))))
}

n_class_518 <- function(classes, cls) {
  sum(vapply(classes, function(k) cls %in% k, logical(1)))
}

test_that("a loglogistic stop on a zeroed score is not converged (#518)", {
  res <- fit_518("loglogistic", c(-1e5, 1))
  f <- res$fit$fit
  # Premise: the probe reaches the case, a finite log-likelihood far below
  # the maximum with a score that is not finite there.
  expect_lt(f$objective, -1e6)
  expect_identical(f$rel_gradient_reason,
                   "the score has a non-finite component at the estimates")
  expect_identical(f$converged, FALSE)
  expect_identical(n_class_518(res$classes, "hzr_unverified_convergence"), 1L)
  expect_match(res$summarised, "converged:\\s+FALSE", all = FALSE)
})

test_that("a lognormal stop that fails SAS/C's gradient test is not converged (#518)", {
  res <- fit_518("lognormal", c(1e5, 1))
  f <- res$fit$fit
  # Premise: nlm() code 2 and a test failed by orders of magnitude.
  expect_identical(f$polish_code, 2L)
  expect_gt(f$rel_gradient, 1e3 * .Machine$double.eps^(1 / 3))
  expect_identical(f$converged, FALSE)
  expect_identical(n_class_518(res$classes, "hzr_unverified_convergence"), 1L)
  expect_match(res$summarised, "converged:\\s+FALSE", all = FALSE)
})

test_that("sane starts in every family still converge, unflagged (#518)", {
  sane <- list(exponential = 50, weibull = c(0.1, 1),
               lognormal = c(1, 1), loglogistic = c(1, 1))
  maxima <- c(exponential = -434.29, weibull = -234.327,
              lognormal = -228.663, loglogistic = -232.832)
  for (d in names(sane)) {
    res <- fit_518(d, sane[[d]])
    f <- res$fit$fit
    expect_identical(f$converged, TRUE, label = d)
    expect_equal(f$objective, maxima[[d]], tolerance = 1e-5, label = d)
    expect_lte(f$rel_gradient, .Machine$double.eps^(1 / 3))
    expect_identical(n_class_518(res$classes, "hzr_unverified_convergence"),
                     0L, label = d)
  }
})

test_that("a sound but poorly scaled fit that narrowly fails the test stays converged (#518)", {
  # The Weibull fit of avc on age and age^2 ends with mu near 0.002 and
  # fails SAS/C's test by about 17 times, because the test takes every
  # parameter's typical size to be 1. A scaled BFGS continuation from there
  # gains 1e-7 in log-likelihood: it is at its maximum. Only a failure by
  # more than 1000 times reports a fit as not converged.
  data(avc, package = "TemporalHazard", envir = environment())
  classes <- list()
  fit <- withCallingHandlers(
    hazard(survival::Surv(int_dead, dead) ~ age + I(age^2),
           data = stats::na.omit(avc), dist = "weibull",
           theta = c(mu = 0.01, nu = 0.5, 0, 0), fit = TRUE),
    warning = function(w) {
      classes[[length(classes) + 1L]] <<- class(w)
      invokeRestart("muffleWarning")
    }
  )
  gradtl <- .Machine$double.eps^(1 / 3)
  # Premise: the test fails here, by less than the margin.
  expect_gt(fit$fit$rel_gradient, 10 * gradtl)
  expect_lt(fit$fit$rel_gradient, 1000 * gradtl)
  expect_identical(fit$fit$converged, TRUE)
  expect_equal(fit$fit$objective, -219.2328881, tolerance = 1e-8)
  expect_identical(n_class_518(classes, "hzr_unverified_convergence"), 0L)
  # The failure is still shown to the reader.
  expect_output(print(fit), "not met, nlm code 3")
})

test_that("the multiphase path keeps its own reading of a zeroed score (#518)", {
  # Rosenbrock shifted by 1e4 with a NaN in the score: the wrapped gradient
  # zeroes it and BFGS stops far from the optimum at (1, 1). The
  # single-distribution path (mark_infeasible = TRUE) now says so; the
  # multiphase path opts out, as it does for #486 and #512, and records the
  # test as not evaluated.
  logl <- function(theta, ...) {
    -(1e4 + 100 * (theta[2] - theta[1]^2)^2 + (1 - theta[1])^2)
  }
  score <- function(theta, ...) c(NaN, -200 * (theta[2] - theta[1]^2))
  run <- function(mark) {
    classes <- character()
    fit <- withCallingHandlers(
      .hzr_optim_generic(
        logl_fn = logl, gradient_fn = score, time = 1, status = 1,
        theta_start = c(-1.2, 1), hessian_fn = function(theta) diag(2),
        mark_infeasible = mark
      ),
      warning = function(w) {
        classes <<- c(classes, class(w))
        invokeRestart("muffleWarning")
      }
    )
    list(fit = fit, classes = classes)
  }
  single <- run(TRUE)
  expect_identical(single$fit$convergence, 99L)
  expect_true("hzr_unverified_convergence" %in% single$classes)
  expect_match(single$fit$message, "gradient it had set to zero", fixed = TRUE)
  multi <- run(FALSE)
  expect_identical(multi$fit$convergence, 0L)
  expect_false("hzr_unverified_convergence" %in% multi$classes)
  expect_identical(multi$fit$rel_gradient_reason,
                   "the score has a non-finite component at the estimates")
})
