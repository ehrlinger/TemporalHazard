# A wrapper that forwards its formula argument by variable makes a vector
# fit whose call still names a formula (#406).
#
#   wrap <- function(dat) { fml <- NULL; hazard(formula = fml, time = ...) }
#
# stores `call$formula = quote(fml)`: a symbol, not NULL. Every test of
# `is.null(call$formula)` therefore took the fit for a formula-interface one
# and tried to resolve `fml`, so the refit failed with a message about
# resolving a name rather than about the interface.
#
# This recognises ONE shape, the one #406 reported: `time =` in the call, the
# name STILL bound to NULL where the fit was made, and no data frame beside
# any stored design -- which a formula fit must have, since hazard() requires
# `data` with a formula. Every other shape keeps the behaviour it has on
# main, and the tests below pin those too, by vintage, so what was left
# alone is visible.
# Deciding the rest needs hazard() to record its interface at fit time, which
# is the 1.3.1 issue.

wv_data_406 <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc[, c("int_dead", "dead", "age", "mal", "com_iv")])
}

# The wrapper of the issue: a vector fit, no `x`, `fml` NULL in its frame.
wv_wrap_406 <- function(dat) {
  fml <- NULL
  hazard(formula = fml, time = dat$int_dead, status = dat$dead,
         dist = "weibull", theta = c(0.1, 1), fit = TRUE)
}

wv_wrap_x_406 <- function(dat) {
  fml <- NULL
  hazard(formula = fml, time = dat$int_dead, status = dat$dead,
         x = cbind(age = dat$age), dist = "weibull", theta = c(0.1, 1, 0),
         fit = TRUE)
}

wv_msg_406 <- function(expr) tryCatch(expr, error = conditionMessage)

test_that("the wrapper's vector fit resamples as the plain fit does (#406)", {
  d <- wv_data_406()
  plain <- hazard(time = d$int_dead, status = d$dead, dist = "weibull",
                  theta = c(0.1, 1), fit = TRUE)
  want <- hzr_bootstrap(plain, n_boot = 3L, seed = 1L)
  # Known positive: the plain fit's replicates are real ones, so an identical
  # result is evidence of resampling and not of two empty tables.
  expect_identical(nrow(want$replicates), 6L)
  expect_true(all(want$summary$sd > 0))

  wf <- wv_wrap_406(d)
  expect_identical(wf$call$formula, quote(fml))
  expect_null(wf$data$x)
  expect_null(wf$data$x_design)
  got <- hzr_bootstrap(wf, n_boot = 3L, seed = 1L)
  expect_identical(got$replicates, want$replicates)
  expect_identical(got$n_success, want$n_success)
})

test_that("the wrapper's vector fit is refused as the vector fit it is (#406)", {
  # A single-distribution vector fit cannot be screened either way; what
  # changes is that the refusal now names the interface instead of blaming a
  # name it could not resolve.
  d <- wv_data_406()
  plain <- hazard(time = d$int_dead, status = d$dead, dist = "weibull",
                  theta = c(0.1, 1), fit = TRUE)
  screen <- function(f) {
    wv_msg_406(hzr_stepwise(f, scope = ~ mal, data = d, direction = "forward",
                            criterion = "wald", trace = FALSE))
  }
  expect_match(screen(plain), "built via the vector interface", fixed = TRUE)
  expect_match(screen(wv_wrap_406(d)), "built via the vector interface",
               fixed = TRUE)
})

test_that("a multiphase wrapper fit screens as the plain fit does (#406)", {
  skip_on_cran() # multiphase screens
  # A multiphase vector fit CAN be refit: its scope lives in the phase
  # formulas and the refit rebuilds the response from the stored vectors.
  d <- wv_data_406()
  ph <- function() {
    list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ age),
         constant = hzr_phase("constant"))
  }
  ctl <- list(n_starts = 1L, maxit = 300L)
  plain <- suppressWarnings(hazard(
    time = d$int_dead, status = d$dead, data = d, dist = "multiphase",
    phases = ph(), fit = TRUE, control = ctl
  ))
  wrap_m <- function(dat) {
    fml <- NULL
    suppressWarnings(hazard(
      formula = fml, time = dat$int_dead, status = dat$dead, data = dat,
      dist = "multiphase", phases = ph(), fit = TRUE, control = ctl
    ))
  }
  screen <- function(f) {
    sw <- suppressWarnings(hzr_stepwise(
      f, scope = list(constant = ~ mal), data = d, direction = "forward",
      criterion = "wald", slentry = 0.99, trace = FALSE
    ))
    list(steps = sw$steps[, c("action", "variable", "phase")],
         theta = sw$fit$theta)
  }
  want <- screen(plain)
  # Known positive: the screen takes a step, so a refit ran.
  expect_identical(nrow(want$steps), 1L)
  expect_identical(screen(wrap_m(d)), want)
})

test_that("the wrapper's vector fit with `x` matches the plain one (#406)", {
  # #406's own reproduction passes `x =`. Both paths must answer exactly as
  # they do for the same fit made without the wrapper -- including the
  # bootstrap refusal, which exists because a design passed as `x` cannot be
  # resampled with the rows. Skipping that refusal would pair resampled
  # outcomes with the original design.
  d <- wv_data_406()
  plain <- hazard(time = d$int_dead, status = d$dead,
                  x = cbind(age = d$age), dist = "weibull",
                  theta = c(0.1, 1, 0), fit = TRUE)
  wf <- wv_wrap_x_406(d)
  expect_false(is.null(wf$data$x))
  expect_null(wf$data$frame)
  screen <- function(f) {
    wv_msg_406(hzr_stepwise(f, scope = ~ mal, data = d, direction = "forward",
                            criterion = "wald", trace = FALSE))
  }
  boot <- function(f) wv_msg_406(hzr_bootstrap(f, n_boot = 2L, seed = 1L))
  expect_match(screen(plain), "built via the vector interface", fixed = TRUE)
  expect_identical(screen(wf), screen(plain))
  expect_match(boot(plain), "passed directly as `x`", fixed = TRUE)
  expect_identical(boot(wf), boot(plain))
})

test_that("a fit storing a design beside a frame keeps its behaviour (#406)", {
  # hazard() requires `data` with a formula, so a fit holding both a design
  # and a frame could be either interface. The rule declines it and its
  # message stays main's -- about resolving the name. Changing that needs
  # the interface recorded at fit time (#432).
  d <- wv_data_406()
  wf <- local({
    fml <- NULL
    hazard(formula = fml, time = d$int_dead, status = d$dead,
           x = cbind(age = d$age), data = d, dist = "weibull",
           theta = c(0.1, 1, 0), fit = TRUE)
  })
  expect_false(is.null(wf$data$x))
  expect_false(is.null(wf$data$frame))
  expect_match(
    wv_msg_406(hzr_stepwise(wf, scope = ~ mal, data = d,
                            direction = "forward", criterion = "wald",
                            trace = FALSE)),
    "did not resolve to a formula", fixed = TRUE
  )
  expect_match(wv_msg_406(hzr_bootstrap(wf, n_boot = 2L, seed = 1L)),
               "did not resolve to a formula", fixed = TRUE)
})

test_that("a fit saved before `call_env` keeps the behaviour it had (#406)", {
  # Vintage: saved before v1.2.2, which is when `call_env` was first stored.
  # Without it the binding cannot be read, so the fit is left alone.
  d <- wv_data_406()
  old <- wv_wrap_406(d)
  old$call_env <- NULL
  expect_match(
    wv_msg_406(hzr_stepwise(old, scope = ~ mal, data = d,
                            direction = "forward", criterion = "wald",
                            trace = FALSE)),
    "could not be resolved", fixed = TRUE
  )
})

test_that("a name bound to something else keeps the behaviour it had (#406)", {
  # The rule reads the binding only to confirm it is still NULL. A name
  # bound to a formula is not classified, so nothing about such a fit
  # changes -- including, on main and here, that its screen runs.
  d <- wv_data_406()
  wf <- wv_wrap_406(d)
  assign("fml", survival::Surv(int_dead, dead) ~ mal, envir = wf$call_env)
  sw <- hzr_stepwise(wf, scope = ~ mal, data = d, direction = "forward",
                     criterion = "wald", trace = FALSE)
  expect_identical(nrow(sw$steps), 0L)
  expect_match(wv_msg_406(hzr_bootstrap(wf, n_boot = 2L, seed = 1L)),
               "cannot count the rows to resample", fixed = TRUE)
})

test_that("a fit whose call names no `time =` is left alone (#406)", {
  # The rule reads the binding only for a call that names `time =`, which is
  # what a vector fit records. A covariate-free formula fit passed by
  # variable does not, so a NULL binding must not turn it into a vector fit:
  # its message stays main's.
  #
  # The NULL is written into `call_env` itself. hazard() captures the
  # caller's bindings by copy (.hzr_capture_call_env), so rebinding the
  # name in this frame would not reach the fit -- except for a fit made at
  # top level, whose captured environment is empty and parented to the
  # global environment, where the binding stays live.
  d <- wv_data_406()
  fml <- survival::Surv(int_dead, dead) ~ 1
  f <- hazard(formula = fml, data = d, dist = "weibull", theta = c(0.1, 1),
              fit = TRUE)
  expect_false("time" %in% names(f$call))
  expect_null(f$data$x_design)
  assign("fml", NULL, envir = f$call_env)
  expect_match(
    wv_msg_406(hzr_stepwise(f, scope = ~ mal, data = d, direction = "forward",
                            criterion = "wald", trace = FALSE)),
    "did not resolve to a formula", fixed = TRUE
  )
})

test_that("a formula fit saved before `x_design` keeps its guards (#406)", {
  # Vintage: saved by 1.2.10 or earlier, which stored no `x_design`. It has
  # `x`, so the rule declines it and #278's guard still runs; classifying it
  # as a vector fit would resample a model whose variable is not a column of
  # the data, reporting full success.
  d <- wv_data_406()
  extra <- d$age * 2
  legacy <- hazard(formula = survival::Surv(int_dead, dead) ~ extra, data = d,
                   time = d$int_dead, dist = "weibull",
                   theta = c(0.1, 1, 0), fit = TRUE)
  legacy$data$x_design <- NULL
  expect_identical(colnames(legacy$data$x), "extra")
  expect_match(wv_msg_406(hzr_bootstrap(legacy, n_boot = 3L, seed = 1L)),
               "which is not a column of it", fixed = TRUE)
})
