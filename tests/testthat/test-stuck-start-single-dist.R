# A single-distribution fit started where the log-likelihood is finite but
# below the optimizer's -1e10 penalty (#512). The optimizer clamps a
# non-finite objective to 1e10, so from such a start every trial point that
# leaves the finite region scores BETTER than the start, and the line search
# built on a gradient of order 1e177 or more cannot take a step. optim()
# stops where it began with convergence 0, and the fit read as converged,
# log-likelihood -3.55e+196, estimates equal to the start, and "Not done in
# this run: none". #486 caught only a non-finite likelihood at the end.
#
# Each stuck start reproduced that on main 97ff68c8. Loglogistic has no
# such start in the sweep: its likelihood turns non-finite first, which is
# #486's case. The sane starts are the known positives, and exponential from
# theta = 50 is the finite-but-huge START that must still fit (#373's
# -434.29): the criterion is where the optimizer ENDS, not where it began.

stuck_512 <- list(
  exponential = 400,
  weibull     = c(50, 50),
  lognormal   = c(1, -300)
)
sane_512 <- list(
  exponential = 50,
  weibull     = c(0.1, 1),
  lognormal   = c(1, 0.5)
)

fit_512 <- function(dist, theta) {
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
       printed = utils::capture.output(print(fit)),
       summarised = utils::capture.output(print(summary(fit))))
}

n_class_512 <- function(classes, cls) {
  sum(vapply(classes, function(k) cls %in% k, logical(1)))
}

expect_stuck_512 <- function(dist) {
  res <- fit_512(dist, stuck_512[[dist]])
  f <- res$fit
  expect_identical(f$fit$converged, FALSE)
  expect_identical(f$fit$objective, NA_real_)
  # One warning, of its own class, caught by a handler for #486's.
  expect_identical(n_class_512(res$classes, "hzr_start_past_penalty"), 1L)
  expect_identical(n_class_512(res$classes, "hzr_infeasible_start"), 1L)
  # What the user reads: the record names the cause, summary() the flag,
  # and no log-likelihood is shown.
  expect_match(res$printed,
               "standard_errors: the optimizer stopped where the log-likelihood is below its own penalty",
               fixed = TRUE, all = FALSE)
  expect_false(any(grepl("Not done in this run: none", res$printed, fixed = TRUE)))
  expect_match(res$summarised, "converged:\\s+FALSE", all = FALSE)
  shown <- c(res$printed, res$summarised)
  expect_false(any(grepl("log-lik:", shown, fixed = TRUE)))
  # The stuck log-likelihood (e+177 to e+263 here) is not shown anywhere.
  expect_false(any(grepl("e\\+[0-9]{3}", shown)))
}

expect_sane_512 <- function(dist) {
  res <- fit_512(dist, sane_512[[dist]])
  expect_identical(res$fit$fit$converged, TRUE)
  expect_true(is.finite(res$fit$fit$objective))
  expect_lt(abs(res$fit$fit$objective), 1e10)
  expect_identical(n_class_512(res$classes, "hzr_start_past_penalty"), 0L)
  # The probe reads real output: a converged fit prints both lines.
  expect_match(res$summarised, "converged:\\s+TRUE", all = FALSE)
  expect_match(res$printed, "log-lik:", fixed = TRUE, all = FALSE)
  res$fit
}

test_that("an exponential fit stuck below the penalty is not converged (#512)", {
  expect_stuck_512("exponential")
  # A finite-but-huge start that the optimizer leaves is a real fit.
  f <- expect_sane_512("exponential")
  expect_equal(f$fit$objective, -434.29, tolerance = 1e-5)
})

test_that("a weibull fit stuck below the penalty is not converged (#512)", {
  expect_stuck_512("weibull")
  expect_sane_512("weibull")
})

test_that("a lognormal fit stuck below the penalty is not converged (#512)", {
  expect_stuck_512("lognormal")
  expect_sane_512("lognormal")
})

test_that("a weighted optimum below -1e10 is still a converged fit (#512)", {
  # Size alone is not evidence of being stuck: with every weight 5e7 the
  # exponential MLE of avc is -2.17e10, and it passes the relative-gradient
  # test there. Started at the MLE, as a user refitting would.
  data(avc, package = "TemporalHazard", envir = environment())
  classes <- character()
  f <- withCallingHandlers(
    hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
           dist = "exponential", theta = -5.204146,
           weights = rep(5e7, nrow(avc)), fit = TRUE),
    warning = function(w) {
      classes <<- c(classes, class(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(f$fit$converged, TRUE)
  expect_equal(f$fit$objective, -2.17145e10, tolerance = 1e-5)
  expect_equal(unname(f$fit$theta), -5.204146, tolerance = 1e-6)
  expect_false("hzr_infeasible_start" %in% classes)
  expect_match(utils::capture.output(print(f)), "log-lik:", fixed = TRUE,
               all = FALSE)
})

test_that("a bootstrap of a stuck start counts no successes (#512)", {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc[, c("int_dead", "dead")])
  stuck <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
    theta = c(50, 50), fit = TRUE
  ))
  out <- suppressWarnings(hzr_bootstrap(stuck, n_boot = 3L, seed = 1L))
  expect_identical(out$n_success, 0L)
  expect_identical(out$n_failed, 3L)
})
