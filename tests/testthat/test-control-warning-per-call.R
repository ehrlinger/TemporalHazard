# A `control` name the fit does not read is warned about ONCE per selection
# call, not once per candidate refit (#410).
#
# Since #376 an unread `control` element draws a warning and hazard()
# proceeds. hzr_stepwise() and hzr_bootstrap() pass `control` to every
# candidate refit, each of which calls hazard(), so one warning became six in
# a three-step screen and three in a select-mode bootstrap, one per candidate
# refit of the up-front screen it runs before resampling. (The replicate
# screens run muffled, so the bootstrap's count never grew with n_boot.) The
# result is unaffected; the harm is to OTHER warnings. Past 50, R prints only
# "There were 50 or more warnings", so repeated control warnings can bury the
# ill-conditioned-Hessian and gradient-test warnings that say a fit is not
# trustworthy.

cw_data_410 <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc[, c("int_dead", "dead", "age", "mal", "com_iv")])
}

cw_base_410 <- function(d) {
  suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                          dist = "weibull", theta = c(0.1, 1), fit = TRUE))
}

# Every warning, kept apart by whether it is the control one.
cw_run_410 <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(value = v, control = grep("'control' element", w, fixed = TRUE,
                                 value = TRUE),
       other = grep("'control' element", w, fixed = TRUE, value = TRUE,
                    invert = TRUE))
}

test_that("one hazard() call warns once: the figure to match (#410)", {
  d <- cw_data_410()
  one <- cw_run_410(hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                           dist = "weibull", theta = c(0.1, 1), fit = TRUE,
                           control = list(n_starts = 1)))
  expect_length(one$control, 1L)
  expect_match(one$control, "control$n_starts", fixed = TRUE)
})

test_that("a screen warns once, not once per candidate refit (#410)", {
  skip_on_cran() # two stepwise screens
  d <- cw_data_410()
  base <- cw_base_410(d)
  screen <- function(...) {
    hzr_stepwise(base, scope = ~ age + mal + com_iv, data = d,
                 direction = "forward", criterion = "wald", trace = FALSE,
                 ...)
  }
  got <- cw_run_410(screen(control = list(n_starts = 1)))
  expect_length(got$control, 1L)
  # The screen itself is unaffected: same steps as with no `control` at all.
  plain <- cw_run_410(screen())
  expect_identical(got$value$steps$variable, plain$value$steps$variable)
  expect_identical(nrow(got$value$steps), 3L)
  # This screen raises no other warnings either way, so the muffling check
  # that can fail is the multiphase one below; this only pins that the
  # control warning was the ONLY warning.
  expect_identical(got$other, character())
  expect_identical(plain$other, character())
})

test_that("a select-mode bootstrap warns once for the whole call (#410)", {
  skip_on_cran() # bootstrap replicates with a screen each
  d <- cw_data_410()
  base <- cw_base_410(d)
  boot <- function(...) {
    hzr_bootstrap(base, n_boot = 3L, seed = 1L, scope = ~ age + mal,
                  criterion = "wald", ...)
  }
  got <- cw_run_410(boot(control = list(n_starts = 1)))
  expect_length(got$control, 1L)
  plain <- cw_run_410(boot())
  # Known positive: the replicates really selected, so "unchanged" is not
  # two empty tables agreeing.
  expect_true(any(c("age", "mal") %in% plain$value$summary$parameter))
  expect_identical(got$value$replicates, plain$value$replicates)
})

test_that("a genuine diagnostic still surfaces beside it (#410)", {
  skip_on_cran() # a multiphase screen
  # The point of the fix is the other warnings. A multiphase screen that
  # raises its own diagnostics must raise exactly as many of them whether or
  # not an unread `control` name is passed.
  d <- cw_data_410()
  ph <- list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                               fixed = "shapes", formula = ~ age),
             constant = hzr_phase("constant"))
  mp <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = ph, fit = TRUE, control = list(n_starts = 1L, maxit = 300L)
  ))
  screen <- function(...) {
    hzr_stepwise(mp, scope = list(constant = ~ mal + com_iv), data = d,
                 direction = "forward", criterion = "wald", slentry = 0.99,
                 trace = FALSE, ...)
  }
  plain <- cw_run_410(screen())
  # Known positive: this screen does raise diagnostics of its own.
  expect_gt(length(plain$other), 0L)
  got <- cw_run_410(screen(control = list(zzz_unread = 1)))
  expect_length(got$control, 1L)
  expect_identical(sort(got$other), sort(plain$other))
})

test_that("a `control` element the fit reads is still forwarded (#410)", {
  skip_on_cran() # a stepwise screen
  # Validating once must not stop the cleaned list reaching the refits: the
  # final model is the last accepted refit, and a fit records its control.
  d <- cw_data_410()
  base <- cw_base_410(d)
  sw <- cw_run_410(
    hzr_stepwise(base, scope = ~ age + mal + com_iv, data = d,
                 direction = "forward", criterion = "wald", trace = FALSE,
                 control = list(maxit = 77L, n_starts = 1))
  )
  expect_length(sw$control, 1L)
  expect_identical(sw$value$spec$control$maxit, 77L)
  expect_identical(nrow(sw$value$steps), 3L)
})

test_that("a screen that refits nothing still warns once (#410)", {
  skip_on_cran() # a score screen
  # Before, the warning came from the candidate refits, so a screen that
  # never refit anything -- the score criterion with a threshold nothing can
  # clear -- reported the ignored name NOT AT ALL. Validating at the call
  # means the caller is told once whether or not a refit happens.
  d <- cw_data_410()
  base <- cw_base_410(d)
  got <- cw_run_410(
    hzr_stepwise(base, scope = ~ age, data = d, direction = "forward",
                 criterion = "score", slentry = 1e-12, trace = FALSE,
                 control = list(n_starts = 1L))
  )
  expect_identical(nrow(got$value$steps), 0L)
  expect_identical(got$value$criteria$n_refit_failures, 0L)
  expect_length(got$control, 1L)
})
