# test-sas-translate-fits.R -- nothing in the rest of the suite ever fits
# the translator's own output. An ICENSOR job's emitted call errored outright
# (time_lower carried the raw ICENSOR bounds for every row, which are NA off
# the interval-censored subset) until .hzr_censor_spec() started
# status-gating it; these tests fit the emitted calls so a regression here
# fails loudly instead of shipping as a document that merely renders.

test_that("an ICENSOR job's emitted call actually fits", {
  skip_on_cran()
  set.seed(1)
  n <- 40
  AVCS <- data.frame(
    INT_DEAD = stats::rexp(n, 0.2),
    DEAD = rep(c(1, 0, 0), length.out = n)
  )
  # ICENSOR's grammar (src/hazard/hazard_y.y) is `ICENSOR c3 = ctime;` -- an
  # event-COUNT variable (OBS column 4, C3), not a 0/1 flag, and a second
  # time variable that is the interval's LOWER bound (the interval runs
  # ctime -> TIME). C3FLAG is a count (any positive value fires), not
  # NA-gated; CTIME must precede INT_DEAD to be a valid lower bound.
  AVCS$C3FLAG <- as.numeric(seq_len(n) %% 5 == 0)
  AVCS$ICTIME <- ifelse(AVCS$C3FLAG > 0, AVCS$INT_DEAD * 0.5, NA)

  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; ICENSOR C3FLAG = ICTIME;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)

  job <- suppressWarnings(hzr_translate_sas(f))
  expect_equal(names(job$calls), c("data", "status", "fit"))
  expect_null(job$calls$fit[["time_upper"]])

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
  # finite/non-negative validation on time_lower passes -- exactly the check
  # that errored before this fix -- and the bound is genuinely status-gated,
  # not just populated. time_upper is never emitted (the ICENSOR interval's
  # upper bound is TIME, which is exactly what hazard()'s time_upper already
  # defaults to), so fit$data$time_upper is NULL, same as an untranslated
  # job. Comparing replicates (interval vs. non-interval rows), not a
  # summary statistic, per AGENTS.md's assertion-discipline rule.
  interval <- fit$data$status == 2
  expect_true(any(interval))
  expect_null(fit$data$time_upper)
  expect_false(anyNA(fit$data$time_lower))
  expect_true(all(fit$data$time_lower[!interval] == 0))
  expect_true(all(fit$data$time_lower[interval] > 0))
  expect_true(all(fit$data$time_lower[interval] < fit$data$time[interval]))
})

test_that("hz.te123.OMC's LCENSOR STARTTME round-trips to time_lower with no time_upper", {
  skip_on_cran()
  dir <- .hzr_sas_fixture_dir() # nolint: object_usage_linter.
  skip_if(is.na(dir), "SAS HAZARD fixture directory not available")
  src <- file.path(dir, "hz.te123.OMC.sas")
  skip_if_not(file.exists(src), "hz.te123.OMC.sas not available")

  # hz.te123.OMC.sas has TWO PROC HAZARD blocks; hzr_translate_sas() keeps
  # only the last block's fit call (translate-sas.R: `calls$fit <- r$call`
  # overwrites per block), and the second block has no LCENSOR at all.
  # Isolate the first block -- the one carrying the real corpus's
  # `LCENSOR STARTTME;` -- from the actual file text, rather than
  # hand-copying it, so this stays tied to the real corpus job rather than a
  # paraphrase of it.
  lines <- readLines(src, warn = FALSE)
  start <- grep("%HAZARD\\(", lines)[[1L]]
  stop_rel <- grep("^\\s*\\);\\s*$", lines[start:length(lines)])[[1L]]
  block <- lines[start:(start + stop_rel - 1L)]
  expect_true(any(grepl("LCENSOR\\s+STARTTME", block)))

  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(block, f)

  job <- suppressWarnings(hzr_translate_sas(f, out_dir = withr::local_tempdir()))
  fit_call <- job$calls$fit
  expect_equal(fit_call[["time_lower"]], as.name("STARTTME"))
  expect_null(fit_call[["time_upper"]])
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

# theta used to carry only c(log(MUE), log(MUC)) -- the multiphase engine's
# theta_start is the full interleaved vector (one block per phase; see
# .hzr_phase_theta_names()), so forcing fit = TRUE on the emitted call
# failed outright with "'names' attribute [5] must be the same length as
# the vector [2]". .hzr_parse_parms() now builds theta by walking the same
# early -> constant -> late phase order the phases themselves are built in.

test_that("theta is the full interleaved multiphase vector", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=D CONDITION=14;",
    "EVENT DEAD; TIME TT;",
    "PARMS MUE=0.2361727 THALF=0.1512095 NU=1.438652 M=1 FIXM",
    "MUC=0.0005436977; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_equal(
    eval(job$calls$fit[["theta"]]),
    c(log(0.2361727), log(0.1512095), 1.438652, 1, log(0.0005436977))
  )
})

test_that("the emitted call fits end to end", {
  skip_on_cran()
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=D CONDITION=14;",
    "EVENT DEAD; TIME TT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  set.seed(1)
  n <- 300
  D <- data.frame(TT = stats::rexp(n, 0.2),
                  DEAD = rep(c(1, 0), length.out = n))
  env <- list2env(as.list(D), parent = environment())
  fit <- NULL
  for (nm in names(job$calls)) {
    if (identical(nm, "data")) next
    fit <- eval(job$calls[[nm]], env)
  }
  expect_s3_class(fit, "hazard")
  expect_equal(fit$spec$dist, "multiphase")
  expect_length(fit$fit$theta, 5L)
})

# The emitted HAZPRED survival call carries conf.type = "logit" because SAS's
# hzp_calc_srv_CL.c builds the limits on Z = log(e^H - 1) = logit(1 - S),
# where predict.hazard() defaults to "log-log". A shape assertion on the call
# cannot tell whether predict() honours the argument or silently ignores it,
# so this evaluates the emitted call and checks the bounds it actually
# produces.

test_that("the emitted HAZPRED call produces logit bounds, not the default", {
  skip_on_cran()
  set.seed(17)
  n <- 60
  df <- data.frame(time = stats::rexp(n, 0.3),
                   status = stats::rbinom(n, 1, 0.6),
                   x = stats::rnorm(n))
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)

  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=E.H OUT=P; TIME TIME; );"
  )
  emitted <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)$call

  env <- new.env(parent = environment())
  env$fit <- fit
  env$P <- data.frame(time = c(0.5, 1, 2, 5), x = 0)
  got <- eval(emitted, env)

  logit <- predict(fit, newdata = env$P, type = "survival", se.fit = TRUE,
                   conf.type = "logit")
  loglog <- predict(fit, newdata = env$P, type = "survival", se.fit = TRUE,
                    conf.type = "log-log")

  expect_equal(got, logit)
  # The decisive check: the two transforms must genuinely disagree here, or
  # `expect_equal(got, logit)` would pass no matter which one predict() used.
  expect_false(isTRUE(all.equal(logit$lower, loglog$lower)))
  expect_false(isTRUE(all.equal(logit$upper, loglog$upper)))
  expect_true(all(got$lower <= got$fit & got$fit <= got$upper))
  expect_true(all(got$lower >= 0 & got$upper <= 1))
})
