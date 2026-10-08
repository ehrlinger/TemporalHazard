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
    # From c(-100, 1, 1, 1) this fit's false maximum depends on the platform:
    # code 3 at LL -9477.9 on Linux, code 4 at LL -1072.3 on Windows (#615).
    # The rule leaves codes 4 and 5 to their own warning, so there that one
    # fires and the classed one does not; the fit is warned about either way.
    if (nm == "lognormal_covariates" && isTRUE(f$polish_code %in% c(4L, 5L))) {
      msgs <- vapply(r$warnings, conditionMessage, character(1))
      expect_true(any(grepl("fail the relative-gradient test SAS/C HAZARD",
                            msgs, fixed = TRUE)), label = paste(nm, "code 4/5"))
      expect_length(pfm_classed(r), 0L)
      expect_gt(f$rel_gradient, 1e-3)
      expect_true(f$converged, label = paste(nm, "converged unchanged"))
      next
    }
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
    # Badly scaled covariates can trip it on a good fit, where rescaling,
    # not a restart, is the remedy.
    expect_match(conditionMessage(w[[1L]]), "centre or rescale",
                 fixed = TRUE)
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

# The #566 fixture, copied so this file stands alone.
pfm_w566_data <- function(alpha, b) {
  set.seed(1)
  n <- 400
  x <- 1000 + stats::rnorm(n, sd = 5)
  t <- (stats::rexp(n) / exp(alpha + b * x))^(1 / 0.2)
  cens <- stats::rexp(n)^5 * 3
  data.frame(time = pmin(t, cens), dead = as.integer(t <= cens),
             x = x, xc = x - 1000)
}

test_that("a stop with no polish code (no lower point found) warns too", {
  skip_on_cran()
  # The nlm continuation records a code only when it improves the point, so
  # a continuation that found no lower point leaves it NA -- the same stop as
  # code 2 or 3. This raw-x Weibull fit stops 0.04 below its centred twin,
  # the same model, with no code recorded.
  d <- pfm_w566_data(60, -0.06)
  r <- pfm_fit(hazard(survival::Surv(time, dead) ~ x, data = d,
                      dist = "weibull", theta = c(exp(300), 0.2, -0.06),
                      fit = TRUE))
  ctr <- suppressWarnings(hazard(survival::Surv(time, dead) ~ xc, data = d,
                                 dist = "weibull", theta = c(1, 0.2, -0.06),
                                 fit = TRUE))
  f <- r$fit$fit
  # The premises: no code, a large gradient, and genuinely short.
  expect_true(is.na(f$polish_code))
  expect_gt(f$rel_gradient, 1e-3)
  expect_gt(ctr$fit$objective - f$objective, 0.01)
  w <- pfm_classed(r)
  expect_length(w, 1L)
  # No code to name, so none is named.
  expect_false(grepl("nlm code", conditionMessage(w[[1L]]), fixed = TRUE))
  expect_true(f$converged)
})

test_that("good fits do not warn, and a stuck fit below 1e-3 is a known miss", {
  skip_on_cran()
  d <- pfm_avc()
  # A good fit: the #531 model from a start that reaches its maximum.
  r <- pfm_fit(hazard(survival::Surv(int_dead, dead) ~ age + I(age^2),
                      data = d, dist = "lognormal", theta = c(1, 1, 0, 0),
                      fit = TRUE))
  expect_lt(r$fit$fit$rel_gradient, 1e-3)
  expect_length(pfm_classed(r), 0L)

  # THE DOCUMENTED LIMITATION. An intercept-only Weibull on times near
  # 1e-170 stops at its starting shape, 1.5, against 1.72 from the same
  # data rescaled to O(1), with a relative gradient of about 5e-4: below the
  # threshold, and below the worst good fit in the suite (6.8e-4), so no
  # threshold separates it. The warning says "may not be"; its absence is
  # not a guarantee.
  set.seed(3)
  t <- stats::rweibull(200, 1.4, 1)
  cens <- stats::rexp(200, 0.3)
  tt <- pmin(t, cens)
  st <- as.integer(t <= cens)
  r <- pfm_fit(hazard(time = tt * 1e-170, status = st, dist = "weibull",
                      theta = c(mu = 1e170, nu = 1.5), fit = TRUE))
  ref <- suppressWarnings(hazard(time = tt, status = st, dist = "weibull",
                                 theta = c(mu = 1, nu = 1), fit = TRUE))
  f <- r$fit$fit
  # The premises: converged, stuck near its start, and below the line.
  expect_true(f$converged)
  expect_lt(abs(f$theta[[2L]] - 1.5), 1e-3)
  expect_gt(ref$fit$theta[[2L]] - f$theta[[2L]], 0.1)
  expect_gt(f$rel_gradient, 1e-4)
  expect_lt(f$rel_gradient, 1e-3)
  expect_length(pfm_classed(r), 0L)
})

test_that("multiphase fits are not flagged, whatever the stop", {
  # Left for 1.3.0: the relative gradient does not separate good from stuck
  # multiphase fits. The optimizer's real return, with the stop overridden.
  real <- .hzr_optim_multiphase
  testthat::local_mocked_bindings(.hzr_optim_multiphase = function(...) {
    r <- real(...)
    r$polish_code <- 2L
    r$rel_gradient <- 1e-2
    r$convergence <- 0L
    r
  })
  set.seed(5)
  tt <- stats::rexp(120, 0.3) + 0.05
  r <- pfm_fit(hazard(time = tt, status = rep(1L, 120), dist = "multiphase",
                      phases = list(a = hzr_phase("cdf"),
                                    b = hzr_phase("constant")),
                      fit = TRUE, control = list(n_starts = 1L)))
  expect_identical(r$fit$fit$polish_code, 2L)  # the premise: the mock took
  expect_equal(r$fit$fit$rel_gradient, 1e-2)
  expect_true(r$fit$fit$converged)
  expect_length(pfm_classed(r), 0L)
})

pfm_boot_base <- function() {
  set.seed(7)
  df <- data.frame(time = stats::rexp(80, 0.4),
                   status = rep(c(1, 1, 0), length.out = 80),
                   z = stats::rnorm(80))
  hazard(survival::Surv(time, status) ~ z, data = df, dist = "exponential",
         theta = c(log_rate = 0, z = 0), fit = TRUE)
}

pfm_boot <- function(base) {
  ws <- list()
  bs <- withCallingHandlers(
    hzr_bootstrap(base, n_boot = 4L, seed = 3L),
    warning = function(w) {
      ws[[length(ws) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  list(value = bs, classed = Filter(
    function(w) inherits(w, "hzr_possible_false_maximum"), ws))
}

test_that("hzr_bootstrap() counts replicates that meet the rule, and warns once", {
  # Replicates run with their warnings suppressed, so the per-fit warning
  # never reaches the user; the bootstrap reads each replicate's own fit and
  # reports the count once. Stuck BY CONSTRUCTION: the optimizer's real
  # return with the stop overridden, so every replicate meets the rule and
  # the count must equal n_success exactly.
  base <- pfm_boot_base()
  real <- .hzr_optim_exponential
  testthat::local_mocked_bindings(.hzr_optim_exponential = function(...) {
    r <- real(...)
    r$polish_code <- 2L
    r$rel_gradient <- 1e-2
    r$convergence <- 0L
    r
  })
  b <- pfm_boot(base)
  n_ok <- b$value$n_success
  expect_gt(n_ok, 0L)
  expect_length(b$classed, 1L)
  expect_match(conditionMessage(b$classed[[1L]]),
               paste0("^", n_ok, " of ", n_ok, " successful replicates"))
  # Still pooled: no replicate was turned into a failure by the rule.
  expect_identical(b$value$n_failed, 0L)
})

pfm_select_data <- function() {
  set.seed(11)
  n <- 150
  z <- stats::rnorm(n)
  t <- stats::rexp(n, 0.3 * exp(0.9 * z))
  data.frame(time = pmin(t, 6), status = as.integer(t <= 6), z = z)
}

test_that("only the replicate that met the rule is counted", {
  # One stuck call among clean ones: the count is per replicate, so it must
  # be exactly 1, not every replicate after it.
  base <- pfm_boot_base()
  calls <- new.env()
  calls$n <- 0L
  real <- .hzr_optim_exponential
  testthat::local_mocked_bindings(.hzr_optim_exponential = function(...) {
    r <- real(...)
    calls$n <- calls$n + 1L
    if (calls$n == 2L) {
      r$polish_code <- 2L
      r$rel_gradient <- 1e-2
      r$convergence <- 0L
    }
    r
  })
  b <- pfm_boot(base)
  n_ok <- b$value$n_success
  expect_gt(n_ok, 2L)                  # the premise: clean ones followed
  expect_gt(calls$n, 2L)
  expect_length(b$classed, 1L)
  expect_match(conditionMessage(b$classed[[1L]]),
               paste0("^1 of ", n_ok, " successful replicates"))
})

test_that("a stepwise bootstrap counts a stuck BASE refit, not only the final fit", {
  skip_on_cran()
  # The reviewer's shape: in select mode each replicate refits the base model
  # and then runs a stepwise screen. The base refit is stuck BY CONSTRUCTION
  # (the intercept-only exponential's return overridden to code 2 and a
  # relative gradient of 0.01); every refit with z in it is left clean, so
  # the final fit is clean. The replicate's base drove its candidate scores,
  # so the replicate counts.
  df <- pfm_select_data()
  base <- hazard(survival::Surv(time, status) ~ 1, data = df,
                 dist = "exponential", theta = c(log_rate = 0), fit = TRUE)
  calls <- new.env()
  calls$stuck <- 0L
  calls$clean_with_z <- 0L
  real <- .hzr_optim_exponential
  testthat::local_mocked_bindings(.hzr_optim_exponential = function(...) {
    r <- real(...)
    if (length(r$par) == 1L) {
      r$polish_code <- 2L
      r$rel_gradient <- 1e-2
      r$convergence <- 0L
      calls$stuck <- calls$stuck + 1L
    } else {
      calls$clean_with_z <- calls$clean_with_z + 1L
    }
    r
  })
  ws <- list()
  bs <- withCallingHandlers(
    hzr_bootstrap(base, n_boot = 3L, seed = 2L, scope = ~ z,
                  direction = "forward", criterion = "wald"),
    warning = function(w) {
      ws[[length(ws) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  # The premises: base refits were stuck, and refits with z ran clean.
  expect_gt(calls$stuck, 0L)
  expect_gt(calls$clean_with_z, 0L)
  hit <- Filter(function(w) inherits(w, "hzr_possible_false_maximum"), ws)
  n_ok <- bs$n_success
  expect_gt(n_ok, 0L)
  expect_length(hit, 1L)
  expect_match(conditionMessage(hit[[1L]]),
               paste0("^", n_ok, " of ", n_ok, " successful replicates"))
  expect_match(conditionMessage(hit[[1L]]), "in the base fit or a refit",
               fixed = TRUE)
})

test_that("a stepwise bootstrap with every fit clean raises no such count", {
  skip_on_cran()
  df <- pfm_select_data()
  base <- hazard(survival::Surv(time, status) ~ 1, data = df,
                 dist = "exponential", theta = c(log_rate = 0), fit = TRUE)
  ws <- list()
  bs <- withCallingHandlers(
    hzr_bootstrap(base, n_boot = 3L, seed = 2L, scope = ~ z,
                  direction = "forward", criterion = "wald"),
    warning = function(w) {
      ws[[length(ws) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_gt(bs$n_success, 0L)
  expect_false(any(vapply(ws, inherits, TRUE, "hzr_possible_false_maximum")))
})

test_that("a good fit's bootstrap raises no such count", {
  b <- pfm_boot(pfm_boot_base())
  expect_gt(b$value$n_success, 0L)
  expect_length(b$classed, 0L)
})
