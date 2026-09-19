# hazard() checks `control` against the elements the fitter reads (#376).
# An element it would ignore used to be accepted silently: a typo such as
# `n_startz` left the default in force, and `control$fix` returned a fit
# identical to the unconstrained one, the "fixed" parameter having moved.
# Every such element now warns and the fit proceeds, as stats::optim() does
# for unknown control names. Nothing errors: an error inside a stepwise or
# bootstrap candidate refit is recorded as a failed candidate.

cv_data <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc)[1:150, ]
}

cv_weibull <- function(control, fit = FALSE) {
  d <- cv_data()
  hazard(survival::Surv(int_dead, dead) ~ age, data = d, dist = "weibull",
         theta = c(mu = 0.01, nu = 0.5, 0.004), fit = fit,
         control = control)
}

cv_multiphase <- function(control, fit = FALSE) {
  d <- cv_data()
  suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = fit, control = control
  ))
}

# cv_multiphase() suppresses warnings, so it cannot show a stray control
# warning; these calls do not fit (fit = FALSE), so any warning is the
# validator's.
cv_multiphase_raw <- function(control, dist = "multiphase") {
  hazard(
    survival::Surv(int_dead, dead) ~ 1, data = cv_data(), dist = dist,
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    control = control
  )
}

test_that("an unknown control element warns, is named, and the fit proceeds", {
  expect_warning(
    obj <- cv_weibull(list(n_startz = 99)),
    "control\\$n_startz \\(not an element any fit reads.*maxit, reltol, shape_param_count"
  )
  expect_s3_class(obj, "hazard")
  w <- character()
  withCallingHandlers(
    cv_multiphase_raw(list(n_startz = 99, nonsense = "x")),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(w, 1L)
  expect_match(w, "control\\$n_startz .*control\\$nonsense ")
  # The multiphase list names the multiphase elements too.
  expect_match(w, "n_starts, conserve")
  # The fit is the fit without it.
  want <- cv_weibull(list(maxit = 500), fit = TRUE)
  got <- suppressWarnings(cv_weibull(list(maxit = 500, n_startz = 9), fit = TRUE))
  expect_identical(got$fit$theta, want$fit$theta)
})
test_that("control$fix warns that it was never read, pointing at hzr_phase(fixed = )", {
  expect_warning(cv_weibull(list(fix = 2L)), "hzr_phase\\(fixed = \\)")
  expect_warning(cv_multiphase_raw(list(fix = 2L)), "hzr_phase\\(fixed = \\)")
  # It says what happened before, so a user of an old fit knows it was
  # unconstrained and that results obtained with it may be affected.
  expect_warning(cv_weibull(list(fix = 2L)), "never read.*unconstrained.*may be affected")
  # And the fit proceeds, as the unconstrained fit it always was.
  want <- cv_weibull(list(maxit = 500), fit = TRUE)
  got <- suppressWarnings(cv_weibull(list(maxit = 500, fix = 2L), fit = TRUE))
  expect_identical(got$fit$theta, want$fit$theta)
})
test_that("control$quasi warns that it was never read", {
  expect_warning(cv_multiphase_raw(list(quasi = TRUE)),
                 "control\\$quasi \\(hazard\\(\\) has never read it")
})
test_that("a documented name that no fit reads warns, and the fit is unchanged", {
  # cv_weibull(), not cv_multiphase(): the latter suppresses warnings.
  reasons <- c(abstol = "bounded optimizer", method = "BFGS",
               condition = "CONDITION= has no equivalent",
               nocov = "prints nothing", nocor = "prints nothing")
  for (nm in names(reasons)) {
    ctl <- stats::setNames(list(1), nm)
    expect_warning(
      obj <- cv_weibull(ctl),
      paste0("no effect.*control\\$", nm, " \\(.*", reasons[[nm]]),
      label = nm
    )
    expect_s3_class(obj, "hazard")
  }
  # Several at once share one warning, in the order given.
  expect_warning(cv_weibull(list(condition = 14, method = "bfgs")),
                 "control\\$condition.*control\\$method")
  # The fit is the fit without them: they are ignored, not applied.
  want <- cv_weibull(list(maxit = 500), fit = TRUE)
  expect_warning(
    got <- cv_weibull(list(maxit = 500, abstol = 1e-3, condition = 14),
                      fit = TRUE),
    "no effect"
  )
  expect_identical(got$fit$theta, want$fit$theta)
  expect_identical(got$fit$objective, want$fit$objective)
})

test_that("an unknown name and a documented one share one warning", {
  expect_warning(cv_weibull(list(condition = 14, n_startz = 1)),
                 "control\\$condition .*control\\$n_startz \\(not an element")
})
test_that("a multiphase element on a single distribution warns and fits", {
  expect_warning(
    obj <- cv_weibull(list(n_starts = 3)),
    "control\\$n_starts \\(it applies only to dist = \"multiphase\"\\)"
  )
  expect_s3_class(obj, "hazard")
})

test_that("a warned name forwarded by stepwise still lets the screen select", {
  # The case the warn-not-error choice exists for: hzr_stepwise() forwards
  # control to every candidate refit, and an error there fails each
  # candidate, emptying the screen as if nothing qualified.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  base <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                 dist = "weibull", fit = TRUE, theta = c(mu = 0.01, nu = 0.5))
  sw <- suppressWarnings(hzr_stepwise(
    base, scope = "age", data = d, direction = "forward",
    criterion = "wald", trace = FALSE, control = list(n_starts = 1L)
  ))
  expect_identical(sw$steps$variable, "age")
  expect_identical(sw$criteria$n_refit_failures, 0L)
})

test_that("an unknown name forwarded by stepwise still lets the screen select", {
  # The rejected policy made an unknown name an error inside every candidate
  # refit, which stepwise records as a failed candidate (#386): 0 steps.
  # A warning keeps the screen working; this fails if the error returns.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  base <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                 dist = "weibull", fit = TRUE, theta = c(mu = 0.01, nu = 0.5))
  sw <- suppressWarnings(hzr_stepwise(
    base, scope = "age", data = d, direction = "forward",
    criterion = "wald", trace = FALSE, control = list(n_startz = 1L)
  ))
  expect_identical(sw$steps$variable, "age")
  expect_identical(sw$criteria$n_refit_failures, 0L)
})
test_that("an unnamed control element warns and the fit proceeds", {
  # stats::optim() ignores an unnamed control element silently; hazard()
  # says so, because an element it cannot read is one the user meant.
  expect_warning(obj <- cv_weibull(list(200)), "1 unnamed element")
  expect_s3_class(obj, "hazard")
  expect_warning(cv_weibull(list(maxit = 500, 7)), "1 unnamed element")
})
test_that("the elements the fitter reads are accepted", {
  # Known positive: legitimate lists still construct and fit.
  expect_s3_class(cv_weibull(list(maxit = 200, reltol = 1e-8)), "hazard")
  expect_s3_class(cv_weibull(list(shape_param_count = 2L)), "hazard")
  expect_s3_class(
    cv_multiphase(list(n_starts = 1L, conserve = FALSE, phase_share_tol = 0,
                       start_seed = 3L, maxit = 200)),
    "hazard"
  )
  expect_s3_class(cv_weibull(list()), "hazard")
  fit <- cv_weibull(list(maxit = 500), fit = TRUE)
  expect_true(is.finite(fit$fit$objective))
})


test_that("the multiphase elements draw no warning on a multiphase fit", {
  expect_no_warning(cv_multiphase_raw(
    list(n_starts = 1L, conserve = FALSE, phase_share_tol = 0,
         start_seed = 3L, maxit = 200, reltol = 1e-8)
  ))
})

test_that("shape_param_count on a multiphase fit warns, and the fit is unchanged", {
  # Every multiphase path (the stepwise shape count, the refit, the score
  # test) derives its own theta layout, so nothing reads it there (#405).
  expect_warning(
    cv_multiphase_raw(list(shape_param_count = 2L)),
    "control\\$shape_param_count \\(it applies only to single-distribution"
  )
  plain <- cv_multiphase(list(n_starts = 1L), fit = TRUE)
  given <- cv_multiphase(list(n_starts = 1L, shape_param_count = 9L), fit = TRUE)
  expect_identical(given$fit$theta, plain$fit$theta)
})

test_that("a named dist scalar is classified by its value", {
  # hazard() accepts dist = c(model = "multiphase") and fits multiphase, so
  # the multiphase elements are read and must not be called off-path (#405).
  expect_no_warning(
    cv_multiphase_raw(list(n_starts = 1L), dist = c(model = "multiphase"))
  )
})

test_that("a named dist scalar still refuses a phase-scoped global term", {
  # With a visible function named like a phase, a skipped refusal is silent:
  # constant(age) became a global covariate and entered BOTH phases (#405).
  constant <- function(x) x
  expect_error(
    hazard(survival::Surv(int_dead, dead) ~ constant(age), data = cv_data(),
           dist = c(model = "multiphase"),
           phases = list(
             early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                               fixed = "shapes"),
             constant = hzr_phase("constant")
           )),
    "names phase 'constant' as a function"
  )
})

test_that("hazard() stores dist without names", {
  obj <- cv_multiphase_raw(list(), dist = c(model = "multiphase"))
  expect_identical(obj$spec$dist, "multiphase")
})

test_that("a warned name forwarded by bootstrap select mode still lets replicates select", {
  # hzr_bootstrap() catches each replicate's errors, so if a warned name
  # were an error, every candidate refit would fail and each replicate would
  # still count as a success, having selected nothing. Assert a selection.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  base <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                 dist = "weibull", theta = c(mu = 0.01, nu = 0.5), fit = TRUE)
  # An off-path name and an unknown one: both were errors under a rejected
  # policy, and either would empty every replicate's screen.
  for (ctl in list(list(n_starts = 1), list(n_startz = 1))) {
    bs <- suppressWarnings(hzr_bootstrap(
      base, n_boot = 10, seed = 321, scope = ~ age + mal + com_iv,
      slentry = 0.3, slstay = 0.2, control = ctl
    ))
    expect_gt(bs$n_success, 0)
    covariates <- bs$summary[!bs$summary$parameter %in% c("mu", "nu"), ]
    expect_gt(sum(covariates$n), 0, label = names(ctl))
  }
})

test_that("a warned name that extends a read one is ignored, not partial-matched", {
  # `$` partial-matches on lists, so control$n_starts would read
  # n_starts_extra. The warning says ignored, so the fit must be the default
  # one, not the one n_starts = 8 would give (#405 review).
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc[, c("int_dead", "dead")])
  mp_raw <- function(ctl) {
    hazard(
      survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
      phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1),
                    constant = hzr_phase("constant")),
      fit = TRUE, control = modifyList(list(start_seed = 1L), ctl)
    )
  }
  mp <- function(ctl) suppressWarnings(mp_raw(ctl))
  default <- mp(list())
  eight <- mp(list(n_starts = 8L))
  # Known positive: n_starts changes this fit, so the check can fail.
  expect_false(identical(default$fit$theta, eight$fit$theta))
  w <- character()
  mp_extra <- withCallingHandlers(
    mp_raw(list(n_starts_extra = 8L)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("n_starts_extra \\(not an element", w)))
  expect_identical(mp_extra$fit$theta, default$fit$theta)
})

test_that("a misspelled maxit is ignored, not partial-matched", {
  want <- cv_weibull(list(maxit = 1L), fit = TRUE)
  plain <- cv_weibull(list(), fit = TRUE)
  # Known positive: maxit = 1 changes this fit.
  expect_false(identical(want$fit$theta, plain$fit$theta))
  got <- suppressWarnings(cv_weibull(list(maxitt = 1L), fit = TRUE))
  expect_identical(got$fit$theta, plain$fit$theta)
})

test_that("the stored control keeps none of the ignored elements", {
  obj <- suppressWarnings(cv_weibull(
    list(maxit = 500, maxitt = 1, abstol = 1, n_starts = 2, 7)
  ))
  expect_identical(names(obj$spec$control), "maxit")
  # A fitted multiphase object also records what Conservation of Events
  # did, after the filter; the ignored name must be gone and that kept.
  mp <- suppressWarnings(cv_multiphase_raw(
    list(n_starts = 1L, n_starts_extra = 8L, abstol = 1)
  ))
  fitted <- suppressWarnings(cv_multiphase(
    list(n_starts = 1L, n_starts_extra = 8L, abstol = 1), fit = TRUE
  ))
  for (o in list(mp, fitted)) {
    expect_false(any(c("n_starts_extra", "abstol") %in% names(o$spec$control)))
    expect_true("n_starts" %in% names(o$spec$control))
  }
  expect_true("conserve_applied" %in% names(fitted$spec$control))
})

test_that("the warnings describe the fit as the documentation does", {
  # The unknown-name warning lists the names accepted without a warning:
  # shape_param_count is accepted but not read by the fit (#405 review).
  expect_warning(
    cv_weibull(list(n_startz = 1)),
    "accepted without a warning are maxit, reltol, shape_param_count\\)"
  )
  expect_warning(
    cv_multiphase_raw(list(n_startz = 1)),
    "accepted without a warning are maxit, reltol, n_starts, conserve"
  )
  # method and quasi name the optimizer as ?hazard does, warm-up included.
  for (nm in c("method", "quasi")) {
    expect_warning(
      cv_weibull(stats::setNames(list(TRUE), nm)),
      "BFGS, a quasi-Newton method \\(a multiphase fit may run a Nelder-Mead warm-up first",
      label = nm
    )
  }
})
