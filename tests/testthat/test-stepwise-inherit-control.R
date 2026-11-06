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

# hzr_bootstrap() rebuilds each replicate from the stored call, where
# `control` is an expression. hazard() snapshots a variable the call names
# from a function's frame, but not from the global environment, so there the
# replicates re-read whatever the variable holds by then. They must use the
# control the fit ran under.
boot_ctrl_runs <- function(change, ...) {
  d <- ic_data()
  assign("ctrl_711", list(maxit = 777L, reltol = 1e-9), envir = globalenv())
  on.exit(suppressWarnings(rm("ctrl_711", envir = globalenv())), add = TRUE)
  base <- suppressWarnings(eval(quote(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
    theta = c(0.1, 1), fit = TRUE, control = ctrl_711
  )), list(d = d), globalenv()))
  stopifnot(identical(base$spec$control$maxit, 777L))
  change()
  suppressWarnings(hzr_bootstrap(base, n_boot = 2L, seed = 1L, ...))
}

select_args <- list(scope = ~ age, direction = "forward", criterion = "wald")

test_that("bootstrap replicates use the fit's control, not the variable's new value", {
  skip_on_cran() # bootstrap refits
  # maxit = 1 cannot converge, so a replicate that read the changed
  # variable would fail.
  shrink <- function() assign("ctrl_711", list(maxit = 1L), envir = globalenv())
  for (args in list(list(), select_args)) {
    bs <- do.call(boot_ctrl_runs, c(list(shrink), args))
    expect_identical(bs$n_failed, 0L)
    expect_identical(bs$n_success, 2L)
  }
})

test_that("bootstrap replicates run when the control variable is gone", {
  skip_on_cran() # bootstrap refits
  drop <- function() rm("ctrl_711", envir = globalenv())
  for (args in list(list(), select_args)) {
    bs <- do.call(boot_ctrl_runs, c(list(drop), args))
    expect_identical(bs$n_failed, 0L)
    expect_identical(bs$n_success, 2L)
  }
})
