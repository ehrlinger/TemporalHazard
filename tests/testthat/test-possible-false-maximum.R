# A single-distribution fit can stop at a false maximum far below the best
# one and report converged = TRUE (#518, #531). The relative gradient does not
# separate those stops from good ones cleanly, so the maintainer's rule
# (2026-10-01) is a warning, not a refusal: a code 2 or 3 stop of the nlm
# polish with a relative gradient above 1e-3 warns, with class
# "hzr_possible_false_maximum", and `converged` is left as it is.
#
# Every case asserts its PREMISE first -- the stop code and the gradient the
# rule reads -- so a change that moves a fit off the rule's ground fails here
# rather than passing for the wrong reason.

pfm_avc <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc)
}

# Fit, collecting every warning instead of letting the first one escape.
pfm_fit <- function(expr) {
  ws <- list()
  f <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(fit = f, warnings = ws)
}

pfm_classed <- function(r) {
  Filter(function(w) inherits(w, "hzr_possible_false_maximum"), r$warnings)
}

test_that("the #518 false maxima that stop on code 2 warn, and stay converged", {
  skip_on_cran()
  d <- pfm_avc()
  one <- survival::Surv(int_dead, dead) ~ 1
  cov <- survival::Surv(int_dead, dead) ~ age + mal
  cases <- list(
    lognormal_far = list(one, "lognormal", c(1e5, 1)),
    weibull_far = list(one, "weibull", c(1e5, 0.5)),
    exponential_covariates = list(cov, "exponential", c(100, -1, -1)),
    lognormal_covariates = list(cov, "lognormal", c(-100, 1, 1, 1)),
    # An ordinary start, not a far one (#531, first comment).
    exponential_ordinary_start = list(cov, "exponential",
                                      c(15.693, -0.539, -9.526))
  )
  for (nm in names(cases)) {
    a <- cases[[nm]]
    r <- pfm_fit(hazard(a[[1L]], data = d, dist = a[[2L]], theta = a[[3L]],
                        fit = TRUE))
    f <- r$fit$fit
    # The premise: the stop the rule is about.
    expect_true(f$polish_code %in% c(2L, 3L), label = paste(nm, "code"))
    expect_gt(f$rel_gradient, 1e-3)
    w <- pfm_classed(r)
    expect_length(w, 1L)
    expect_match(conditionMessage(w[[1L]]), "may not be a maximum",
                 fixed = TRUE)
    # The figure is the fit's own, anchored to its phrase.
    expect_match(conditionMessage(w[[1L]]),
                 paste0("relative gradient of ", signif(f$rel_gradient, 3),
                        " (nlm code ", f$polish_code, ")"), fixed = TRUE)
    expect_match(conditionMessage(w[[1L]]), "theta", fixed = TRUE)
    expect_true(f$converged, label = paste(nm, "converged unchanged"))
  }
})

test_that("the gof-timefix fit stuck at a relative gradient of 15 warns", {
  skip_on_cran()
  d <- pfm_avc()
  r <- pfm_fit(hazard(time = d$int_dead + 1e4, status = d$dead,
                      dist = "weibull", theta = c(mu = 1e-4, nu = 1),
                      fit = TRUE))
  f <- r$fit$fit
  expect_true(f$polish_code %in% c(2L, 3L))
  expect_gt(f$rel_gradient, 1)
  expect_length(pfm_classed(r), 1L)
  expect_true(f$converged)
})

test_that("the #518 stops on code 4 keep their own warning, not this one", {
  skip_on_cran()
  # Code 4 (iteration limit) already warns, and says the estimates may not be
  # at the maximum. The new rule is for codes 2 and 3 only, so these get the
  # old warning and not a second one.
  d <- pfm_avc()
  for (th in list(c(1, -30), c(1000, 1))) {
    r <- pfm_fit(hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                        dist = "lognormal", theta = th, fit = TRUE))
    expect_identical(r$fit$fit$polish_code, 4L)
    msgs <- vapply(r$warnings, conditionMessage, "")
    expect_true(any(grepl("iteration limit", msgs, fixed = TRUE)))
    expect_length(pfm_classed(r), 0L)
  }
})
