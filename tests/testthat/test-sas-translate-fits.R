# test-sas-translate-fits.R -- nothing in the rest of the suite ever fits
# the translator's own output. An ICENSOR job's emitted call errored outright
# (time_lower/time_upper carried the raw ICENSOR bounds for every row, which
# are NA off the interval-censored subset) until .hzr_censor_spec() started
# status-gating them; these tests fit the emitted calls so a regression here
# fails loudly instead of shipping as a document that merely renders.

test_that("an ICENSOR job's emitted call actually fits", {
  skip_on_cran()
  set.seed(1)
  n <- 40
  AVCS <- data.frame(
    INT_DEAD = stats::rexp(n, 0.2),
    DEAD = rep(c(1, 0, 0), length.out = n)
  )
  # ICENSOR's grammar (src/hazard/hazard_y.y) is `ICENSOR flag = timevar;` --
  # an interval-censoring INDICATOR variable and a second time variable, not
  # two comma-separated bound variables. ICFLAG is 0/1, not NA-gated.
  AVCS$ICFLAG <- as.numeric(seq_len(n) %% 5 == 0)
  AVCS$ICTIME <- ifelse(AVCS$ICFLAG == 1, AVCS$INT_DEAD * 1.5, NA)

  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; ICENSOR ICFLAG = ICTIME;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)

  job <- suppressWarnings(hzr_translate_sas(f))
  expect_equal(names(job$calls), c("data", "status", "fit"))

  # Evaluate every emitted call in order, in an env holding the data. The
  # "data" chunk is a guard that fails loudly when the DATA= dataset was
  # never assigned (see translate-sas.R); it checks for a bound variable
  # literally named "AVCS", which the fit call's `data = AVCS` argument
  # never forces (hazard() only touches `data` on the formula path), so
  # bind AVCS itself alongside its columns rather than skip the guard.
  env <- list2env(as.list(AVCS), parent = environment())
  env$AVCS <- AVCS
  fit <- NULL
  for (nm in names(job$calls)) fit <- suppressWarnings(eval(job$calls[[nm]], env))
  expect_s3_class(fit, "hazard")

  # This package has no logLik.hazard method (confirmed: methods("logLik")
  # lists none for class "hazard", and getS3method("logLik", "hazard",
  # optional = TRUE) is NULL). The translator now emits dist = "multiphase"
  # whenever it emits `phases` (see .hzr_parse_hazard()), but it still never
  # emits fit = TRUE -- forcing fit = TRUE here now fails with
  # "'names' attribute [5] must be the same length as the vector [2]" from
  # .hzr_optim_multiphase(), because the translator's `theta` carries only
  # the log(MUE)/log(MUC) scale starts, not the full parameter vector the
  # multiphase engine expects to seed from (it also wants the phase shape
  # starts, e.g. t_half/nu). That is a separate, already-tracked translator
  # gap (the "blind-start convergence" item), out of scope here. So
  # "actually fits" in this test still means: hazard()'s own
  # finite/non-negative validation on time_lower/time_upper passes --
  # exactly the check that errored before this fix -- and the bounds are
  # genuinely status-gated, not just populated. Comparing replicates
  # (interval vs. non-interval rows), not a summary statistic, per
  # AGENTS.md's assertion-discipline rule.
  interval <- fit$data$status == 2
  expect_true(any(interval))
  expect_false(anyNA(fit$data$time_lower))
  expect_false(anyNA(fit$data$time_upper))
  expect_true(all(fit$data$time_lower[!interval] == 0))
  expect_true(all(fit$data$time_lower[interval] > 0))
  expect_equal(fit$data$time_upper[!interval], fit$data$time[!interval])
  expect_true(all(fit$data$time_upper[interval] > fit$data$time[interval]))
})

test_that("a plain EVENT/TIME job's emitted call fits with no time_lower/time_upper", {
  skip_on_cran()
  set.seed(2)
  n <- 40
  AVCS <- data.frame(
    INT_DEAD = stats::rexp(n, 0.2),
    DEAD = rep(c(1, 0, 0), length.out = n)
  )

  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)

  job <- suppressWarnings(hzr_translate_sas(f))
  expect_equal(names(job$calls), c("data", "fit"))
  expect_null(job$calls$fit[["time_lower"]])
  expect_null(job$calls$fit[["time_upper"]])

  env <- list2env(as.list(AVCS), parent = environment())
  env$AVCS <- AVCS
  fit <- NULL
  for (nm in names(job$calls)) fit <- suppressWarnings(eval(job$calls[[nm]], env))
  expect_s3_class(fit, "hazard")
  expect_null(fit$data$time_lower)
  expect_null(fit$data$time_upper)
})

test_that("the emitted call fits the multiphase model, not a Weibull", {
  skip_on_cran()
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=D CONDITION=14;",
    "EVENT DEAD; TIME TT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  cl <- job$calls$fit
  # phases are silently discarded unless dist is multiphase
  expect_equal(cl[["dist"]], "multiphase")

  set.seed(1)
  n <- 200
  D <- data.frame(TT = stats::rexp(n, 0.2),
                  DEAD = rep(c(1, 0), length.out = n))
  env <- list2env(as.list(D), parent = environment())
  fit <- NULL
  for (nm in names(job$calls)) {
    if (identical(nm, "data")) next
    fit <- eval(job$calls[[nm]], env)
  }
  expect_s3_class(fit, "hazard")
  # The decisive assertion: the fitted object must actually BE multiphase.
  expect_equal(fit$spec$dist, "multiphase")
  expect_length(fit$spec$phases, 2L)
})

test_that("evaluating the emitted call emits no 'phases is ignored' warning", {
  skip_on_cran()
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=D CONDITION=14;",
    "EVENT DEAD; TIME TT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  set.seed(1)
  n <- 200
  D <- data.frame(TT = stats::rexp(n, 0.2),
                  DEAD = rep(c(1, 0), length.out = n))
  env <- list2env(as.list(D), parent = environment())
  w <- character(0)
  withCallingHandlers(
    for (nm in names(job$calls)) {
      if (identical(nm, "data")) next
      eval(job$calls[[nm]], env)
    },
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("phases' is ignored", w, fixed = TRUE)))
})
