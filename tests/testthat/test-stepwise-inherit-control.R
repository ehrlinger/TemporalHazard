# hzr_stepwise() refits under the base fit's `control` (#688).
#
# Every step of a SAS SELECTION runs inside one PROC HAZARD, under its
# MAXITER=, CONSERVE and the rest. The refits used hazard()'s defaults
# instead, so a base fit run with conserve = FALSE or a raised maxit was
# stepped, and its selected model fitted, under settings it never had, and
# nothing said so. A `control` given to hzr_stepwise() still overrides the
# base, entry by entry.

ic_data <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc[, c("int_dead", "dead", "age", "mal", "com_iv")])
}

ic_screen <- function(base, d, ...) {
  suppressWarnings(hzr_stepwise(base, scope = ~ age + mal + com_iv, data = d,
                                direction = "forward", criterion = "wald",
                                trace = FALSE, ...))
}

test_that("the selected model is refitted under the base fit's control", {
  skip_on_cran() # a stepwise screen
  d <- ic_data()
  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
    theta = c(0.1, 1), fit = TRUE,
    control = list(maxit = 777L, reltol = 1e-9)
  ))
  sw <- ic_screen(base, d)
  # Known positive: the screen took a step, so the result is a refit and
  # not the base fit handed back.
  expect_gt(length(sw$fit$theta), length(base$fit$theta))
  expect_identical(sw$spec$control$maxit, 777L)
  expect_identical(sw$spec$control$reltol, 1e-9)
})

test_that("a control given to hzr_stepwise() overrides the base entry by entry", {
  skip_on_cran() # a stepwise screen
  d <- ic_data()
  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
    theta = c(0.1, 1), fit = TRUE,
    control = list(maxit = 777L, reltol = 1e-9)
  ))
  sw <- ic_screen(base, d, control = list(maxit = 300L))
  expect_gt(length(sw$fit$theta), length(base$fit$theta))
  expect_identical(sw$spec$control$maxit, 300L)
  expect_identical(sw$spec$control$reltol, 1e-9)
})

test_that("a multiphase base passes on conserve, not the state it applied", {
  skip_on_cran() # a multiphase screen
  d <- ic_data()
  ph <- list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                               fixed = "shapes"),
             constant = hzr_phase("constant"))
  mp <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = ph, fit = TRUE,
    control = list(n_starts = 1L, maxit = 300L, conserve = FALSE)
  ))
  # spec$control records what the fit applied beside what it was asked. A
  # refit given those would warn that it does not read them.
  expect_true("conserve_applied" %in% names(mp$spec$control))
  w <- character()
  sw <- withCallingHandlers(
    hzr_stepwise(mp, scope = list(constant = ~ mal + com_iv), data = d,
                 direction = "forward", criterion = "wald", slentry = 0.99,
                 trace = FALSE),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(grep("'control' element", w, fixed = TRUE), 0L)
  expect_gt(length(sw$fit$theta), length(mp$fit$theta))
  expect_false(sw$spec$control$conserve)
  expect_false(isTRUE(sw$spec$control$conserve_applied))
  expect_identical(sw$spec$control$n_starts, 1L)
})

test_that("an unnamed control element still warns when the base has a control", {
  skip_on_cran() # a stepwise screen
  d <- ic_data()
  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
    theta = c(0.1, 1), fit = TRUE, control = list(maxit = 777L)
  ))
  w <- character()
  withCallingHandlers(
    hzr_stepwise(base, scope = ~ age, data = d, direction = "forward",
                 criterion = "wald", trace = FALSE, control = list(500)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(grep("unnamed element", w, fixed = TRUE), 1L)
})

test_that("control = list(maxit = NULL) is refused whatever the base set", {
  d <- ic_data()
  fit_with <- function(control) {
    suppressWarnings(hazard(
      survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
      theta = c(0.1, 1), fit = TRUE, control = control
    ))
  }
  for (base in list(fit_with(list()), fit_with(list(maxit = 777L)),
                    fit_with(list(reltol = 1e-9)))) {
    expect_error(
      hzr_stepwise(base, scope = ~ age, data = d, trace = FALSE,
                   control = list(maxit = NULL)),
      "maxit"
    )
  }
})
