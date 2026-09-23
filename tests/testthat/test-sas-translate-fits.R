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
  # Gated on DEAD == 0: a row where EVENT and C3 both fire is refused at
  # fit time (setlik.c sums the two contributions and hazard() carries one
  # status per row), so a fixture that mixes them tests the refusal, not
  # the fit.
  AVCS$C3FLAG <- as.numeric(seq_len(n) %% 5 == 0 & AVCS$DEAD == 0)
  AVCS$ICTIME <- ifelse(AVCS$C3FLAG > 0, AVCS$INT_DEAD * 0.5, NA)

  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; ICENSOR C3FLAG = ICTIME;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)

  job <- suppressWarnings(hzr_translate_sas(f))
  expect_equal(names(job$calls), c("data", "status", "fit"))
  # job$calls$fit is now `fit <- hazard(...)`; index into the hazard() call.
  expect_null(job$calls$fit[[3L]][["time_upper"]])

  # Evaluate every emitted call in order, in an env holding the data. The
  # "data" chunk is a guard that fails loudly when the DATA= dataset was
  # never assigned (see translate-sas.R); it checks for a bound variable
  # literally named "AVCS", which the fit call's `data = AVCS` argument
  # never forces (hazard() only touches `data` on the formula path), so
  # bind AVCS itself. The "status" chunk this job emits
  # (.hzr_status <- with(AVCS, ifelse(...))) now wraps its SAS-column
  # expression in with(AVCS, ...), so it resolves those columns from AVCS
  # itself and needs no separate per-column binding.
  env <- new.env(parent = environment())
  env$AVCS <- AVCS
  for (nm in names(job$calls)) suppressWarnings(eval(job$calls[[nm]], env))
  fit <- env$fit
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
  # job$calls$fit is now `fit <- hazard(...)`; the hazard() call itself is
  # the assignment's right-hand side.
  fit_call <- job$calls$fit[[3L]]
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
  # The status chunk is emitted for every job now, not only ICENSOR or
  # multi-count ones: it carries the missing/negative guard on EVENT, which
  # has to run and be seen to run ahead of the fit.
  expect_equal(names(job$calls), c("data", "status", "fit"))
  # job$calls$fit is now `fit <- hazard(...)`; index into the hazard() call.
  expect_null(job$calls$fit[[3L]][["time_lower"]])
  expect_null(job$calls$fit[[3L]][["time_upper"]])

  env <- new.env(parent = environment())
  env$AVCS <- AVCS
  for (nm in names(job$calls)) suppressWarnings(eval(job$calls[[nm]], env))
  fit <- env$fit
  expect_s3_class(fit, "hazard")
  expect_null(fit$data$time_lower)
  expect_null(fit$data$time_upper)
})

test_that("a PARMS that builds no usable phase emits a stop(), not a fit", {
  # A template's `?` builds no phase and is not refused: the parser cannot
  # read it. (A MUE without a shape operand used to be here too; it now
  # builds on SAS's defaults, #345, and is executed in the test below.) The fit chunk
  # used to be hazard(fit = TRUE, theta = c()) under the default Weibull,
  # which rendered an unfitted object. The test above is the paired case: a
  # usable PARMS still emits hazard().
  # (A template's `?` used to be the example; PROC HAZARD's lexer rejects `?`,
  # so it is now a syntax stop (U1). An operand written with spaces around `=`
  # is joined and read now (#421). A macro-only PARMS is the parser's own
  # limit: SAS expands it and runs the job, this parser cannot read it.)
  for (parms in c("PARMS &ALLPARMS;")) {
    f <- withr::local_tempfile(fileext = ".sas")
    writeLines(paste(
      "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
      "EVENT DEAD; TIME INT_DEAD;", parms, ");"
    ), f)
    job <- suppressWarnings(hzr_translate_sas(f))
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"), info = parms)
    msg <- tryCatch(eval(job$calls$fit, new.env()), error = conditionMessage)
    expect_match(msg, "builds no phase this translator could use", info = parms)
    # A MUE or MUL with no shape operand now builds its phase (#345), so the
    # message must not offer it as a cause.
    expect_no_match(msg, "with no shape operand", fixed = TRUE, info = parms)
    # Not the reference's own refusal, which needs a PARMS the parser read.
    expect_false(any(grepl("modterm.c", job$untranslated$reason, fixed = TRUE)),
                 info = parms)
  }
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
  # job$calls$fit is now `fit <- hazard(...)`; index into the hazard() call.
  cl <- job$calls$fit[[3L]]
  # phases are silently discarded unless dist is multiphase
  expect_equal(cl[["dist"]], "multiphase")

  set.seed(1)
  n <- 200
  D <- data.frame(TT = stats::rexp(n, 0.2),
                  DEAD = rep(c(1, 0), length.out = n))
  env <- new.env(parent = environment())
  for (nm in names(job$calls)) {
    if (identical(nm, "data")) next
    eval(job$calls[[nm]], env)
  }
  fit <- env$fit
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
  env <- new.env(parent = environment())
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
  # job$calls$fit is now `fit <- hazard(...)`; index into the hazard() call.
  expect_equal(
    eval(job$calls$fit[[3L]][["theta"]]),
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
  env <- new.env(parent = environment())
  for (nm in names(job$calls)) {
    if (identical(nm, "data")) next
    eval(job$calls[[nm]], env)
  }
  fit <- env$fit
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

  txt <- .hzr_sas_normalise(paste(
    "DATA P; DO MONTHS=1 TO 12 BY 1; OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=E.H OUT=P; TIME TIME; );"
  ))
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

test_that("an orphan MUE and MUL translate to the fit of their written defaults (#345)", {
  skip_on_cran()
  # Execute the emitted chunks, not their text: PARMS names both phases by
  # scale only, and PROC HAZARD runs them on its shape defaults. The fit must
  # be the one the same defaults written out produce -- same model, same
  # start -- and a real likelihood, not the optimizer's 1e10 penalty.
  translate_and_fit <- function(parms, data) {
    f <- withr::local_tempfile(fileext = ".sas")
    writeLines(paste(
      "%HAZARD( PROC HAZARD DATA=D NOCONSERVE;",
      "EVENT DEAD; TIME TT;", parms, ");"
    ), f)
    job <- suppressWarnings(hzr_translate_sas(f))
    env <- new.env(parent = asNamespace("TemporalHazard"))
    env$D <- data
    for (nm in names(job$calls)) suppressWarnings(eval(job$calls[[nm]], env))
    list(job = job, fit = env$fit)
  }
  withr::local_seed(345)
  n <- 400
  dat <- data.frame(TT = stats::rweibull(n, 0.8, 5),
                    DEAD = stats::rbinom(n, 1, 0.7))
  orphan <- translate_and_fit("PARMS MUE=0.2 MUL=0.05;", dat)
  written <- translate_and_fit(
    "PARMS MUE=0.2 THALF=1 NU=2 M=1 MUL=0.05 GAMMA=1 ALPHA=1 ETA=2;", dat
  )
  expect_s3_class(orphan$fit, "hazard")
  expect_equal(length(orphan$fit$spec$phases), 2L)
  expect_true(orphan$fit$fit$converged)
  expect_gt(orphan$fit$fit$objective, -1e9)
  expect_identical(orphan$fit$fit$theta, written$fit$fit$theta)
  expect_identical(orphan$fit$fit$objective, written$fit$fit$objective)
  expect_equal(orphan$job$untranslated$reason, written$job$untranslated$reason)
})

test_that("every covariate after a '/' option reaches the fitted model (#342)", {
  skip_on_cran()
  # Executed, not shape-asserted: the fit must carry OPMOS and Y, and must
  # not carry the excluded X. Before #342 the chunk fitted ~AGE + MAL only.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=D CONDITION=14;",
    "EVENT DEAD; TIME TT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005;",
    "EARLY AGE, MAL/I, OPMOS;",
    "CONSTANT X/E, Y; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  set.seed(3)
  n <- 150
  D <- data.frame(TT = stats::rexp(n, 0.2),
                  DEAD = rep(c(1, 0, 1), length.out = n),
                  AGE = stats::rnorm(n), MAL = rep(0:1, length.out = n),
                  OPMOS = stats::runif(n), X = stats::rnorm(n),
                  Y = stats::rnorm(n))
  res <- suppressWarnings(render_sim(job, list(D = D)))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  cf <- names(stats::coef(res$env$fit))
  expect_true(all(c("phase_1.AGE", "phase_1.MAL", "phase_1.OPMOS",
                    "phase_2.Y") %in% cf), info = paste(cf, collapse = " "))
  expect_false(any(grepl("\\.X$", cf)), info = paste(cf, collapse = " "))
})

test_that("a phase variable outside the model still deletes its missing rows (#342)", {
  skip_on_cran()
  # getrisk.c collects every phase-statement variable and readobs.c drops a
  # row where any is missing, whether or not the variable is estimated. An
  # /E variable is not in hazard()'s formula, so hazard() cannot drop those
  # rows: the status chunk has to stop and say so rather than fit more rows
  # than SAS did.
  job_for <- function(stmts) {
    f <- withr::local_tempfile(fileext = ".sas", .local_envir = parent.frame())
    writeLines(paste(
      "%HAZARD( PROC HAZARD DATA=D CONDITION=14; EVENT DEAD; TIME TT;",
      "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005;", stmts, ");"), f)
    suppressWarnings(hzr_translate_sas(f))
  }
  set.seed(4)
  n <- 150
  D <- data.frame(TT = stats::rexp(n, 0.2),
                  DEAD = rep(c(1, 0, 1), length.out = n),
                  AGE = stats::rnorm(n), X = stats::rnorm(n),
                  Y = stats::rnorm(n))
  D_na <- D
  D_na$X[1:40] <- NA

  job <- job_for("EARLY AGE; CONSTANT X/E, Y;")
  res <- suppressWarnings(render_sim(job, list(D = D_na)))
  expect_false(res$ok)
  expect_match(res$results[["status"]], "X")
  expect_match(res$results[["status"]], "missing")
  expect_false(exists("fit", envir = res$env, inherits = FALSE))
  # Complete X: nothing to delete, and the model still leaves X out.
  res <- suppressWarnings(render_sim(job, list(D = D)))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))

  # X missing only where the modelled AGE is missing too: hazard() drops
  # those rows itself, so they are exactly the rows SAS deletes, and there is
  # nothing to stop for (#340 item 8). The guard used to stop here.
  D_same <- D
  D_same$AGE[1:40] <- NA
  D_same$X[1:40] <- NA
  res <- suppressWarnings(render_sim(job, list(D = D_same)))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  # The fit used exactly the rows SAS keeps: the same chunks rendered on the
  # 110 complete rows give the same estimates. (fit$data keeps all rows; the
  # multiphase path drops incomplete design rows while fitting.)
  sub <- suppressWarnings(render_sim(job, list(D = D_same[-(1:40), ])))
  expect_true(sub$ok, info = paste(sub$results, collapse = "; "))
  expect_equal(stats::coef(res$env$fit), stats::coef(sub$env$fit))

  # A second statement for a phase adds to the first (hazard_y.y appends),
  # it does not replace it.
  job <- job_for("EARLY AGE; EARLY Y;")
  res <- suppressWarnings(render_sim(job, list(D = D)))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  expect_true(all(c("phase_1.AGE", "phase_1.Y") %in%
                    names(stats::coef(res$env$fit))))
})

test_that("a job with no DATA= and phase covariates emits a stop(), not a fit (#311)", {
  # With no DATA= the fit chunk has no `data`, and hazard() refuses a phase
  # formula with nothing to evaluate it in (#299): the chunk stopped with
  # advice to pass `data =`, an argument the SAS job never had.
  translate <- function(src) {
    f <- withr::local_tempfile(fileext = ".sas", .local_envir = parent.frame())
    writeLines(src, f)
    suppressWarnings(hzr_translate_sas(f))
  }
  parms <- "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005;"
  set.seed(1)
  n <- 60
  D <- data.frame(TT = stats::rexp(n, 0.2),
                  DEAD = rep(c(1, 0), length.out = n),
                  MAL = rep(c(0, 1, 1), length.out = n))

  for (cov in c("EARLY MAL=0;", "CONSTANT MAL=0;")) {
    job <- translate(paste("%HAZARD( PROC HAZARD CONDITION=14;",
                           "EVENT DEAD; TIME TT;", parms, cov, ");"))
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"), info = cov)
    expect_true(any(job$untranslated$construct == "DATA="), info = cov)
    # Executed where the variables exist as vectors, the way a reader who
    # binds them would render it: the chunk must fail with the translator's
    # reason, not hazard()'s advice to pass `data =`.
    res <- render_sim(job, as.list(D))
    expect_false(res$ok, info = cov)
    expect_match(res$results[["fit"]], "no DATA=", info = cov)
  }

  # Paired controls. No covariates: no DATA= is fine, the chunk fits from
  # the bound vectors.
  job <- translate(paste("%HAZARD( PROC HAZARD CONDITION=14;",
                         "EVENT DEAD; TIME TT;", parms, ");"))
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
  expect_false(any(job$untranslated$construct == "DATA="))
  res <- suppressWarnings(render_sim(job, as.list(D)))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  expect_s3_class(res$env$fit, "hazard")

  # DATA= with the same covariate: still a fit, and the phase carries MAL.
  skip_on_cran()
  job <- translate(paste("%HAZARD( PROC HAZARD DATA=D CONDITION=14;",
                         "EVENT DEAD; TIME TT;", parms, "EARLY MAL=0; );"))
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
  expect_false(any(job$untranslated$construct == "DATA="))
  res <- suppressWarnings(render_sim(job, list(D = D)))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  expect_true("phase_1.MAL" %in% names(stats::coef(res$env$fit)))
})

test_that("a job PROC HAZARD rejects at parse emits a stop(), not a fit (#340)", {
  job_for <- function(stmt) {
    f <- withr::local_tempfile(fileext = ".sas", .local_envir = parent.frame())
    writeLines(paste("%HAZARD( PROC HAZARD DATA=D CONDITION=14; EVENT DEAD; TIME TT;",
                     "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005;", stmt, ");"), f)
    suppressWarnings(hzr_translate_sas(f))
  }
  cases <- c("EARLY AGE/EI, Y;" = "hazard_l.l:176",
             "EARLY AGE/E ORDER=2, Y;" = "przconc.c:45-53",
             "EARLY AGE/E/I, Y;" = "syntax error",
             "SELECTION; EARLY AGE/EI, Y;" = "hazard_l.l:176")
  for (stmt in names(cases)) {
    job <- job_for(stmt)
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"), info = stmt)
    expect_error(eval(job$calls$fit), cases[[stmt]], fixed = TRUE, info = stmt)
    expect_error(eval(job$calls$fit), "PROC HAZARD does not run this job",
                 fixed = TRUE, info = stmt)
    expect_false("fit_base" %in% names(job$calls), info = stmt)
  }
  # SAS stops at parse before any semantic check, so the parse refusal must
  # win even over a job that also lacks EVENT, which otherwise stops the
  # whole translation first (Copilot, #396).
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines("%HAZARD( PROC HAZARD DATA=D; TIME TT; PARMS MUE=0.2 THALF=0.15 NU=1; EARLY AGE/EI; );", f)
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_error(eval(job$calls$fit), "hazard_l.l:176", fixed = TRUE)
  # Control: the same options written the way SAS accepts them still fit.
  job <- job_for("EARLY AGE/E I, Y;")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
})

test_that("a phase variable missing from the data is named, not 'object not found' (#340)", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste("%HAZARD( PROC HAZARD DATA=D CONDITION=14; EVENT DEAD; TIME TT;",
                   "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005;",
                   "EARLY AGE, ZZ; CONSTANT QQ/E; );"), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  set.seed(9)
  D <- data.frame(TT = stats::rexp(60), DEAD = rep(c(1, 0), 30), AGE = stats::rnorm(60))
  res <- suppressWarnings(render_sim(job, list(D = D)))
  expect_false(res$ok)
  expect_match(res$results[["status"]], "phase statements name ZZ, QQ", fixed = TRUE)
  expect_match(res$results[["status"]], "not columns of D", fixed = TRUE)
  expect_false(exists("fit", envir = res$env, inherits = FALSE))
  # Control: with the columns present, the chunk runs.
  D$ZZ <- stats::rnorm(60)
  D$QQ <- stats::rnorm(60)
  res <- suppressWarnings(render_sim(job, list(D = D)))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  # A column named like the dataset must not mask it (Copilot, #396). Inside
  # transform(), `D` resolved to that column, whose names() are NULL, so the
  # check refused a job whose variables were all present.
  D$D <- seq_len(60)
  res <- suppressWarnings(render_sim(job, list(D = D)))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  # A phase variable that is not numeric: PROC HAZARD refuses the job
  # (vfynvar.c:22-26 "VARIABLE NOT NUMERIC" sets semerr; hazard.c:249-251
  # exits), while hazard() would dummy-code it and fit.
  D$ZZ <- rep(c("a", "b"), 30)
  res <- suppressWarnings(render_sim(job, list(D = D)))
  expect_false(res$ok)
  expect_match(res$results[["status"]], "These phase variables are not numeric: ZZ.",
               fixed = TRUE)
  expect_false(exists("fit", envir = res$env, inherits = FALSE))
})

test_that("every reachable SETG3 refusal reaches the reader (#359)", {
  # The whole reachable set, not a sample. Twelve codes can fire through
  # .hzr_parse_parms(): the nine the exhaustive search in
  # test-sas-parse-parms.R pins for the .hzr_setg3_notes() trace, plus the
  # three the FIXGE2/FIXGAE2 constraint block raises itself (SETG3940,
  # SETG3990, SETG31000). Each is RENDERED: the fit chunk must fail and bind
  # no `fit`, which is the consequence a reader meets, rather than the
  # message text we happen to emit.
  refusals <- c(
    # site: .hzr_setg3_notes() -- entry checks, setg3.c:269-284
    SETG3900 = "PARMS MUL=0.2 TAU=0 FIXTAU GAMMA=2 ETA=1;",
    SETG3910 = "PARMS MUL=0.2 TAU=1 GAMMA=0 FIXGAMMA ETA=1;",
    SETG3920 = "PARMS MUL=0.2 TAU=1 GAMMA=2 ALPHA=-1 FIXALPHA ETA=1;",
    SETG3930 = "PARMS MUL=0.2 TAU=1 GAMMA=2 ETA=0 FIXETA;",
    # site: .hzr_setg3_notes() -- SETG3_weibull and the alpha fixup
    SETG3960 = "PARMS MUL=0.2 TAU=1 GAMMA=0 ETA=0.25 FIXGE2 WEIBULL;",
    SETG3970 = "PARMS MUL=0.2 TAU=1 GAMMA=2 ETA=0 WEIBULL;",
    SETG3980 = "PARMS MUL=0.2 TAU=1 GAMMA=4 ETA=0.25 ALPHA=0 WEIBULL;",
    SETG31020 = "PARMS MUL=0.2 TAU=1 GAMMA=1 ETA=1 FIXGAMMA FIXETA;",
    SETG31040 = "PARMS MUL=0.2 TAU=1 GAMMA=1 ETA=1 ALPHA=3 FIXALPHA;",
    # site: the constraint block's own alpha = 0 case, which the trace cannot
    # see because it takes ga_two as FALSE
    SETG3980_gae2 = paste("PARMS MUL=0.2 TAU=1 GAMMA=4 ETA=0.25 ALPHA=0",
                          "FIXALPHA FIXGAE2 WEIBULL;"),
    # site: the constraint block, FIXGE2 with both shapes fixed off 2
    SETG3990 = paste("PARMS MUL=0.2 TAU=1 GAMMA=4 ETA=0.25 FIXGAMMA FIXETA",
                     "FIXGE2 WEIBULL;"),
    # site: the constraint block, FIXGAE2 with ALPHA fixed off the constraint
    SETG31000 = paste("PARMS MUL=0.2 TAU=1 GAMMA=4 ETA=0.25 ALPHA=3 FIXALPHA",
                      "FIXGAE2 WEIBULL;"),
    # site: SETG3_ignore_tau under both flags, ALPHA fixed away from 1
    SETG3940 = paste("PARMS MUL=0.2 TAU=1 GAMMA=4 ETA=0.5 ALPHA=3 FIXALPHA",
                     "FIXGAE2 FIXGE2 WEIBULL;")
  )
  translate <- function(parms) {
    f <- withr::local_tempfile(fileext = ".sas", .local_envir = parent.frame())
    writeLines(paste("%HAZARD( PROC HAZARD DATA=D CONDITION=14;",
                     "EVENT DEAD; TIME TT;", parms, ");"), f)
    suppressWarnings(hzr_translate_sas(f))
  }
  set.seed(359)
  n <- 80
  D <- data.frame(TT = stats::rexp(n, 0.2), DEAD = rep(c(1, 0), length.out = n))

  # Since 2026-09-22 EVERY such job emits the fit and warns: a rendered
  # document completes and carries the reason instead of halting on it.
  #
  # Three of these (SETG3910/3920/3930) still halt further down, at
  # hzr_phase(), because SAS refuses them for a shape value that is out of
  # range and hzr_phase() will not build a phase from that same value. The
  # warning is emitted in its own chunk ABOVE the fit so that the SETG3 code
  # and the operand are stated before that happens. Measured caveat, recorded
  # here because it is easy to assume otherwise: when the fit chunk errors,
  # Quarto writes no document, and knitr captures warnings INTO the document,
  # so that warning does not reach the render console either -- it is in the
  # emitted .qmd source above the failing chunk, and in $untranslated.
  for (label in names(refusals)) {
    code <- sub("_gae2$", "", label)
    job <- translate(refusals[[label]])
    warn_nm <- grep("^refusal", names(job$calls), value = TRUE)
    expect_length(warn_nm, 1L)
    msgs <- character(0)
    withCallingHandlers(eval(job$calls[[warn_nm[[1L]]]], new.env()),
                        warning = function(x) {
                          msgs <<- c(msgs, conditionMessage(x))
                          invokeRestart("muffleWarning")
                        })
    # The code reaches the reader, the construct is listed, and the fit is
    # emitted rather than replaced. All three, for every class.
    expect_match(paste(msgs, collapse = " "), code, fixed = TRUE, info = label)
    expect_gt(NROW(job$untranslated), 0L)
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"), info = label)
  }

  # The paired control: a job PROC HAZARD runs still renders a fit. Without
  # it, a fix that refused everything would pass every assertion above.
  ok <- translate("PARMS MUL=0.2 TAU=1 GAMMA=4 ETA=0.25 WEIBULL;")
  ok_res <- suppressWarnings(render_sim(ok, list(D = D)))
  expect_true(ok_res$ok)
  expect_s3_class(get("fit", envir = ok_res$env, inherits = FALSE), "hazard")

  # The second control, and the one this change nearly got wrong: a constraint
  # flag without WEIBULL reaches SETG3 down a path the trace does not model,
  # so the trace's verdict is not PROC HAZARD's. PROC HAZARD fits both of
  # these, and the old narrowing let them render. Under U1 (and after two
  # review passes found the hand-derived non-WEIBULL path wrong both ways)
  # the translation stops on them, saying it cannot tell -- never calling
  # them refused.
  for (ops in c("PARMS MUL=0.2 TAU=1 GAMMA=4 ETA=0.5 FIXGAMMA FIXETA FIXGE2;",
                "PARMS MUL=0.2 TAU=1 GAMMA=4 ETA=1 ALPHA=2 FIXALPHA FIXGAE2;")) {
    job <- translate(ops)
    # These now warn and still fit, so the reason reaches the reader through
    # the refusal chunk rather than through an error.
    warn_nm <- grep("^refusal", names(job$calls), value = TRUE)
    expect_length(warn_nm, 1L)
    msgs <- character(0)
    withCallingHandlers(eval(job$calls[[warn_nm[[1L]]]], new.env()),
                        warning = function(x) {
                          msgs <<- c(msgs, conditionMessage(x))
                          invokeRestart("muffleWarning")
                        })
    msg <- paste(msgs, collapse = " ")
    expect_match(msg, "cannot tell whether PROC HAZARD refuses", fixed = TRUE, info = ops)
    expect_no_match(msg, "refused before any fit", fixed = TRUE, info = ops)
  }
})

# --- U1: a job PROC HAZARD refuses, or fits differently, warns and fits -----
# John's decision (2026-09-19), as amended on 2026-09-22: when the translator
# knows PROC HAZARD refuses a job, or fits a model other than the one it would
# emit, the document EMITS the fit with a loud warning above it and a row in
# $untranslated, rather than either stopping or fitting with a row the reader
# may never see. The first form of the decision stopped; that halted Quarto
# before it wrote any output, including for the jobs preceding the refused
# one, so the warning carries the message instead.

.u1_job <- function(proc = "", parms, env = parent.frame()) {
  f <- withr::local_tempfile(fileext = ".sas", .local_envir = env)
  writeLines(paste0("%HAZARD( PROC HAZARD DATA=D", proc,
                    "; EVENT DEAD; TIME TT; PARMS ", parms, "; );"), f)
  suppressWarnings(hzr_translate_sas(f))
}
# Since 2026-09-22 a job PROC HAZARD would refuse, or fit differently, is
# EMITTED with a loud warning and an $untranslated row, so a rendered
# document completes and carries the reason instead of halting on it. These
# helpers assert that whole contract, not just one half of it.
.u1_refusal_chunk <- function(job) {
  nm <- grep("^refusal", names(job$calls), value = TRUE)
  if (!length(nm)) NULL else job$calls[[nm[[1L]]]]
}
.u1_stops <- function(job) {
  identical(job$calls$fit[[3L]][[1L]], as.name("stop"))
}
# The warn route: the fit IS emitted, a refusal chunk warns, and the
# construct is listed. All three, because any one alone is a half-contract.
.u1_warns_and_fits <- function(job) {
  !is.null(.u1_refusal_chunk(job)) &&
    NROW(job$untranslated) > 0L &&
    identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
}
# A refusal must reach the reader by ONE of two routes and never neither:
# it warns and still fits, or -- where the value PROC HAZARD refuses is also
# one hzr_phase() will not build -- there is no fit to emit and it stops.
.u1_refuses <- function(job) .u1_warns_and_fits(job) || .u1_stops(job)
# The reason text, by whichever route carried it.
.u1_msg <- function(job) {
  ch <- .u1_refusal_chunk(job)
  if (is.null(ch)) {
    if (!.u1_stops(job)) {
      return("no refusal")
    }
    return(tryCatch({
      eval(job$calls$fit, new.env())
      "no error"
    }, error = conditionMessage))
  }
  # Evaluate the emitted chunk rather than reading its text: what a reader of
  # the rendered document receives is the warning, not the source.
  msgs <- character(0)
  withCallingHandlers(eval(ch, new.env()), warning = function(x) {
    msgs <<- c(msgs, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  if (length(msgs)) paste(msgs, collapse = " ") else "no warning"
}

test_that("a PARMS syntax error warns and still fits (U1, #421)", {
  for (p in c("MUE=0.2 THALF=0.5 NU=1E-3", "MUE=0.2 THALF=0.5 NU",
              "MUE=0.2 THALF=0.5 NU = ABC", "MUE=0.2 THALF=0.5 FIXG1")) {
    job <- .u1_job(parms = p)
    expect_true(.u1_refuses(job), info = p)
    expect_match(.u1_msg(job), "PROC HAZARD does not run this job", fixed = TRUE,
                 info = p)
  }
  # Controls. An ordinary job is not refused at all. A macro piece is
  # INDETERMINATE, not a syntax error: SAS expands it before PROC HAZARD
  # reads the statement, so it must never be reported as a job PROC HAZARD
  # does not run. (It is still carried as an unread operand, which is a
  # different claim and is asserted elsewhere.)
  job <- .u1_job(parms = "MUE=0.2 THALF=0.5 NU=1")
  expect_false(.u1_refuses(job))
  job <- .u1_job(parms = "MUE=0.2 THALF=0.5 NU=1 &X")
  expect_false(grepl("PROC HAZARD does not run this job", .u1_msg(job),
                     fixed = TRUE))
})

test_that("a PROC-line value the lexer rejects warns and still fits (U1, #403)", {
  for (proc in c(" MAXITER=1E5", " CONDITION=5.")) {
    job <- .u1_job(proc = proc, parms = "MUE=0.2 THALF=1 NU=1")
    expect_true(.u1_refuses(job), info = proc)
    msg <- .u1_msg(job)
    expect_match(msg, "PROC HAZARD does not run this job", fixed = TRUE, info = proc)
    expect_match(msg, "hazard_l.l:33-38", fixed = TRUE, info = proc)
  }
  job <- .u1_job(proc = " MAXITER=200 CONDITION=14", parms = "MUE=0.2 THALF=1 NU=1")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
})

test_that("MAXITER or CONDITION with no value warns and still fits (U1, #433)", {
  # hazard_y.y:63-64 are `MAXITER '=' NUMBER` and `CONDITION '=' NUMBER`.
  # No NUMBER, no rule: the option falls to `hazardopt : error` (:77),
  # yyerror latches yysynerr (yyerror.c:19) and initprz.c:75-77 terminates
  # the procedure. So an empty value is a syntax error exactly as a
  # non-numeric one is, and a bare keyword with no `=` is too.
  for (proc in c(" MAXITER=", " MAXITER =", " CONDITION=", " CONDITION =",
                 " MAXITER", " CONDITION")) {
    job <- .u1_job(proc = proc, parms = "MUE=0.2 THALF=1 NU=1")
    expect_true(.u1_refuses(job), info = proc)
    msg <- .u1_msg(job)
    expect_match(msg, "PROC HAZARD does not run this job", fixed = TRUE, info = proc)
    expect_match(msg, "hazard_y.y:63-64", fixed = TRUE, info = proc)
  }
  # A macro is still exempt: SAS expands it before PROC HAZARD reads the
  # statement, so whether a NUMBER arrives is not knowable here.
  job <- .u1_job(proc = " MAXITER=&N", parms = "MUE=0.2 THALF=1 NU=1")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
  # And a value that is present and numeric is untouched. MAXITER alone
  # leaves no row; CONDITION always leaves one, because hazard() reads no
  # `condition` and that is recorded rather than emitted (#384). So the
  # assertion is that neither is REJECTED, not that nothing is recorded.
  job <- .u1_job(proc = " MAXITER=200", parms = "MUE=0.2 THALF=1 NU=1")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
  expect_identical(NROW(job$untranslated), 0L)
  job <- .u1_job(proc = " MAXITER=200 CONDITION=14", parms = "MUE=0.2 THALF=1 NU=1")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
  expect_false(any(grepl("syntax error", job$untranslated$reason, fixed = TRUE)))
})

test_that("a spaced PROC-line option is read, not dropped (U1 review 4, #421)", {
  # SAS's lexer skips whitespace (hazard_l.l:32), so MAXITER = 250 is one
  # option. Read as three tokens it is dropped, and the fit then runs on
  # hazard()'s own iteration limit rather than the job's 250 -- a different
  # model with no refusal. Asserted on the emitted call, not on the parse.
  for (p in c(" MAXITER=250", " MAXITER = 250", " MAXITER= 250",
              " MAXITER =250")) {
    job <- .u1_job(proc = p, parms = "MUE=0.2 THALF=1 NU=1")
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"), info = p)
    expect_true(any(grepl("maxit = 250", deparse(job$calls$fit), fixed = TRUE)),
                info = p)
    expect_identical(NROW(job$untranslated), 0L, info = p)
  }
  # The PROC HAZPRED line joins the same way: its grid name survives.
  hp <- .hzr_parse_hazpred(list(text = "PROC HAZPRED DATA = G INHAZ = HZ"), "")
  expect_false(any(hp$untranslated$reason == "unknown PROC HAZPRED option"))
})

test_that("FIXMNU1 on an active early phase warns and still fits (U1, #358)", {
  job <- .u1_job(parms = "MUE=0.2 THALF=1 NU=2 M=0.5 FIXMNU1")
  expect_true(.u1_refuses(job))
  msg <- .u1_msg(job)
  expect_match(msg, "FIXMNU1", fixed = TRUE)
  expect_match(msg, "|M*NU| = 1", fixed = TRUE)
  expect_match(msg, "#358", fixed = TRUE)
  # With no early phase there is nothing for FIXMNU1 to constrain: SAS fits.
  job <- .u1_job(parms = "MUL=0.2 TAU=1 GAMMA=2 ETA=1 FIXMNU1")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
})

# --- U1, r-reviewer pass 1 on the branch: false refusals and missed ones ---

test_that("a macro value or call is never refused as a syntax error (U1 review)", {
  # SAS expands `&X` and `%CALL` before PROC HAZARD reads the statement.
  for (cs in list(list(proc = " MAXITER=&MX", parms = "MUE=0.2 THALF=1 NU=1"),
                  list(proc = " CONDITION=&C", parms = "MUE=0.2 THALF=1 NU=1"),
                  list(proc = "", parms = "MUE=0.2 THALF=1 NU=1 %FIXOPTS"),
                  list(proc = "", parms = "MUE=0.2 THALF=1 NU=&V"))) {
    job <- .u1_job(proc = cs$proc, parms = cs$parms)
    info <- paste(cs$proc, cs$parms)
    expect_false(grepl("does not run", .u1_msg(job), fixed = TRUE), info = info)
  }
})

test_that("an unspaced non-numeric PARMS value is a syntax error (U1 review)", {
  job <- .u1_job(parms = "MUE=0.2 THALF=0.5 NU=ABC")
  expect_true(.u1_refuses(job))
  expect_match(.u1_msg(job), "PROC HAZARD does not run this job", fixed = TRUE)
})

test_that("FIXMNU1 warns honestly, even where it might constrain nothing (U1 review 2)", {
  # The vacuous-FIXMNU1 exemption was dropped after its sign cases went wrong
  # (M*NU = -1 is vacuous too; M = NU = -1 is SETG1920). The stop claims
  # neither that PROC HAZARD runs the job nor that it refuses it.
  for (p in c("MUE=0.2 THALF=1 NU=1 M=1 FIXM FIXNU FIXMNU1",
              "MUE=0.2 THALF=0.5 M=-1 NU=1 FIXM FIXNU FIXMNU1",
              "MUE=0.2 THALF=0.5 M=-1 NU=-1 FIXM FIXNU FIXMNU1")) {
    job <- .u1_job(parms = p)
    expect_true(.u1_refuses(job), info = p)
    msg <- .u1_msg(job)
    expect_match(msg, "cannot emit PROC HAZARD's model", fixed = TRUE, info = p)
    expect_no_match(msg, "runs this job", fixed = TRUE, info = p)
  }
})

test_that("SETG3's entry refusals warn on the constraint path too (U1 review)", {
  # setg3.c:269-284 run before any constraint or WEIBULL logic.
  for (p in c("MUL=0.2 TAU=0 FIXTAU GAMMA=2 ETA=1 FIXGE2",
              "MUL=0.2 TAU=1 GAMMA=0 FIXGAMMA ETA=1 FIXGAE2",
              "MUL=0.2 TAU=1 GAMMA=2 ALPHA=-1 FIXALPHA ETA=1 FIXGE2")) {
    job <- .u1_job(parms = p)
    expect_true(.u1_refuses(job), info = p)
    expect_match(.u1_msg(job), "SETG39[0-3]0", info = p)
  }
})

test_that("FIXGE2/FIXGAE2 without WEIBULL warns, saying it cannot tell (U1 review 2)", {
  # The non-WEIBULL constraint path (SETG3_all_gt_0(), SETG3_alpha_le_0(), ...)
  # is not modelled here, and hand-deriving it went wrong twice. Every such
  # job stops and claims neither a refusal nor a run -- including the jobs
  # where the emitted phase happens to be SAS's model (loud, not wrong).
  for (p in c("MUL=0.2 TAU=1 GAMMA=4 ETA=0.25 FIXGAMMA FIXETA FIXGE2",
              "MUL=0.2 TAU=1 GAMMA=3 ETA=1 FIXGE2",
              "MUL=0.2 TAU=1 GAMMA=4 ETA=0.5 FIXGAMMA FIXETA FIXGE2",
              "MUL=0.2 TAU=1 GAMMA=1 ETA=1 ALPHA=0 FIXALPHA FIXGAMMA FIXETA FIXGE2",
              "MUL=0.2 TAU=1 GAMMA=4 ETA=1 ALPHA=0 FIXALPHA FIXGAMMA FIXETA FIXGAE2",
              "MUL=0.2 TAU=1 GAMMA=4 ETA=1 ALPHA=3 FIXGAMMA FIXETA FIXGAE2")) {
    job <- .u1_job(parms = p)
    expect_true(.u1_refuses(job), info = p)
    msg <- .u1_msg(job)
    expect_match(msg, "cannot tell whether PROC HAZARD refuses", fixed = TRUE, info = p)
    expect_no_match(msg, "runs this job", fixed = TRUE, info = p)
  }
  # The same shapes with WEIBULL are the mirrored path, and fit.
  expect_false(.u1_refuses(.u1_job(parms = "MUL=0.2 TAU=1 GAMMA=3 ETA=1 FIXGE2 WEIBULL")))
  expect_false(.u1_refuses(.u1_job(parms = "MUL=0.2 TAU=1 GAMMA=4 ETA=0.5 FIXGAMMA FIXETA FIXGE2 WEIBULL")))
})

test_that("DELTA != 0 and a FIXTAU with no TAU written warn (U1 review)", {
  job <- .u1_job(parms = "MUE=0.2 THALF=1 NU=1 DELTA=0.5")
  expect_true(.u1_refuses(job))
  expect_match(.u1_msg(job), "DELTA", fixed = TRUE)
  job <- .u1_job(parms = "MUL=0.2 GAMMA=2 ETA=1 FIXTAU")
  expect_true(.u1_refuses(job))
  expect_match(.u1_msg(job), "0.75*Tmax", fixed = TRUE)
  # Controls: DELTA = 0 is R's model; a written positive TAU is fixed at it.
  expect_false(.u1_refuses(.u1_job(parms = "MUE=0.2 THALF=1 NU=1 DELTA=0")))
  # DELTA is read only by SETG1 (setg1.c:306), which runs only for an active
  # early phase (shape.c:19-21): a late-only job ignores it.
  expect_false(.u1_refuses(.u1_job(parms = "MUL=0.2 TAU=1 GAMMA=2 ETA=1 DELTA=0.5")))
  expect_false(.u1_refuses(.u1_job(parms = "MUL=0.2 TAU=2 GAMMA=2 ETA=1 FIXTAU")))
})

test_that("a template placeholder or a bare % is a syntax error, not filled in (U1 review 2)", {
  # `NU=?` used to fit with NU filled from SAS's default; PROC HAZARD's lexer
  # rejects `?` (hazard_l.l:178). `50%` is no macro: % with no name after it.
  for (p in c("MUE=0.2 THALF=1 NU=?", "MUE=0.2 THALF=1 NU=50%")) {
    job <- .u1_job(parms = p)
    expect_true(.u1_refuses(job), info = p)
    expect_match(.u1_msg(job), "PROC HAZARD does not run this job", fixed = TRUE, info = p)
  }
  expect_match(.u1_msg(.u1_job(parms = "MUE=0.2 THALF=1 NU=?")), "fill it in", fixed = TRUE)
})

test_that("U1 review 3: the last DELTA wins, and more known-unfittable jobs warn", {
  # hazard_y.y:138 is last-wins, so DELTA=0.5 DELTA=0 runs at delta = 0 --
  # exactly what is emitted.
  expect_false(.u1_refuses(.u1_job(parms = "MUE=0.2 DELTA=0.5 DELTA=0 NU=1 M=1 THALF=1")))
  # An unknown PROC option IS a lexer catch-all in PROC HAZARD, but this
  # parser's block text can carry another step's keywords (a %repeat call
  # brings a DATA step through), so it is recorded rather than refused.
  job <- .u1_job(proc = " FOO", parms = "MUE=0.2 THALF=1 NU=1")
  expect_false(.u1_refuses(job))
  expect_true("FOO" %in% job$untranslated$construct)
  # FIXTAU whose TAU this parser could not read: PROC HAZARD fixes TAU at the
  # written value or at 0.75*Tmax, never at the 1 the emitted phase pins.
  job <- .u1_job(parms = "MUL=0.2 GAMMA=1 TAU=&T FIXTAU")
  expect_true(.u1_refuses(job))
  expect_match(.u1_msg(job), "FIXTAU", fixed = TRUE)
  # An active MU whose phase this parser could not build: PROC HAZARD fits
  # that phase, so the emitted model is short of one.
  job <- .u1_job(parms = "MUE=0.2 THALF=0.5 NU=1 M=1 MUL=0.3 &SHAPE FIXGE2")
  expect_true(.u1_refuses(job))
  expect_match(.u1_msg(job), "MUL", fixed = TRUE)
  # Control: the same job with the late shape written builds both phases.
  expect_false(.u1_refuses(.u1_job(parms = "MUE=0.2 THALF=0.5 NU=1 M=1 MUL=0.3 GAMMA=2 ETA=1 WEIBULL")))
})

test_that("U1 review 4: spaces around `=` are SAS's job, not a refusal (#421)", {
  # hazard_l.l:32 skips whitespace, so `MAXITER = 50` and `THALF = 0.3` are
  # the same jobs as their unspaced forms. They used to be split apart here:
  # the PROC line called them a syntax error (a false refusal), and PARMS
  # filled the operand from SAS's default and fitted a model PROC HAZARD does
  # not fit. Operands are joined before parsing.
  job <- .u1_job(proc = " MAXITER = 50", parms = "MUE=0.2 THALF=1 NU=1")
  expect_false(.u1_refuses(job))
  # `fit_base` exists only for a SELECTION job, and this one has none, so an
  # earlier `expect_equal(job$calls$fit_base %||% job$calls$fit,
  # job$calls$fit)` reduced to expect_equal(x, x) and could not fail (#433
  # review). Assert the property that was meant instead: no base-fit chunk is
  # emitted, which a SELECTION regression here WOULD break.
  expect_false("fit_base" %in% names(job$calls))
  # The value is read, not defaulted.
  job <- .u1_job(parms = "MUE=0.2 THALF = 0.3 NU = 1 M=1")
  expect_false(.u1_refuses(job))
  src <- paste(deparse(job$calls$fit), collapse = " ")
  expect_match(src, "t_half = 0.3", fixed = TRUE)
  expect_match(src, "nu = 1", fixed = TRUE)
  # Each spelling joins.
  for (p in c("MUE=0.2 THALF= 0.3 NU=1 M=1", "MUE=0.2 THALF =0.3 NU=1 M=1")) {
    src <- paste(deparse(.u1_job(parms = p)$calls$fit), collapse = " ")
    expect_match(src, "t_half = 0.3", fixed = TRUE, info = p)
  }
  # A joined operand SAS still rejects is still a syntax stop.
  expect_true(.u1_refuses(.u1_job(parms = "MUE=0.2 THALF = ABC NU=1")))
  # An unknown HAZARD statement keyword is recorded, for the same reason.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste("%HAZARD( PROC HAZARD DATA=D; EVENT DEAD; TIME TT; FOO BAR;",
                   "PARMS MUL=0.2 TAU=1 GAMMA=2 ETA=1 WEIBULL; );"), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_false(.u1_refuses(job))
  expect_true("FOO" %in% job$untranslated$construct)
})

test_that("every refusal class warns AND lists a row, and a clean job does neither", {
  # The contract John set on 2026-09-22, asserted as a whole and over the
  # whole class list rather than one example: a job PROC HAZARD would refuse,
  # or would fit differently, EMITS the fit, warns loudly, and names the
  # construct in $untranslated. Any one of the three alone is a half-contract:
  # a warning nobody can grep afterwards, a row nobody reads at render, or a
  # fit that quietly disappeared.
  classes <- list(
    "PARMS value the lexer rejects"  = list("", "MUE=0.2 THALF=1 NU=1E-3"),
    "PARMS value keyword, no number" = list("", "MUE=0.2 THALF=1 NU"),
    "PARMS keyword outside grammar"  = list("", "MUE=0.2 THALF=1 NU=1 FIXG1"),
    "PARMS flag given a value"       = list("", "MUE=0.2 THALF=1 NU=1 FIXNU=1"),
    "PROC value the lexer rejects"   = list(" MAXITER=1E5", "MUE=0.2 THALF=1 NU=1"),
    "PROC option with no value"      = list(" MAXITER=", "MUE=0.2 THALF=1 NU=1"),
    "template ? placeholder"         = list("", "MUE=0.2 THALF=? NU=1"),
    "FIXMNU1 on an active early"     = list("", "MUE=0.2 THALF=1 NU=2 M=0.5 FIXMNU1"),
    "DELTA != 0 on an active early"  = list("", "MUE=0.2 THALF=0.3 NU=1 DELTA=0.5"),
    "FIXTAU with no TAU written"     = list("", "MUL=0.2 GAMMA=2 ETA=1 ALPHA=1 FIXTAU"),
    "FIXGAE2 without WEIBULL"        = list("", "MUL=0.1 TAU=8 ALPHA=2 GAMMA=5 ETA=1 FIXGAE2"),
    "an operand read as a macro"     = list("", "MUE=0.2 THALF=0.3 NU=1 M=0 &FLAGS")
  )
  for (nm in names(classes)) {
    job <- .u1_job(proc = classes[[nm]][[1L]], parms = classes[[nm]][[2L]])
    expect_true(.u1_warns_and_fits(job), info = nm)
    # Spelled out, so a failure says which half broke.
    expect_false(is.null(.u1_refusal_chunk(job)), info = nm)
    expect_gt(NROW(job$untranslated), 0L)
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"), info = nm)
    expect_no_match(.u1_msg(job), "no warning", fixed = TRUE, info = nm)
  }

  # KNOWN NEGATIVE: a job PROC HAZARD runs, translated faithfully, must raise
  # neither half. Without it a change that warned on everything would satisfy
  # every assertion above.
  clean <- .u1_job(parms = "MUE=0.2 THALF=1 NU=1")
  expect_null(.u1_refusal_chunk(clean))
  expect_identical(NROW(clean$untranslated), 0L)
  expect_identical(clean$calls$fit[[3L]][[1L]], as.name("hazard"))
})

test_that("SETG3's entry refusals are raised in the C's own order (#433 review)", {
  # setg3.c:269-284 (src/model/, pin dad7978) checks TAU, then GAMMA, then
  # ALPHA, then ETA, and EACH RETURNS before the next -- and all four return
  # before SETG3_ignore_tau() at :309-323 and before the WEIBULL branch. So a
  # job tripping two of them is refused by the FIRST, and a message naming the
  # later code names a refusal PROC HAZARD never reaches.
  #
  # One case per adjacent pair, plus the entry-vs-WEIBULL case that prompted
  # this: it named SETG3980 while SAS returns SETG3910.
  pairs <- list(
    # TAU before GAMMA
    list(parms = "MUL=0.2 TAU=0 FIXTAU GAMMA=0 FIXGAMMA ETA=1",
         want = "SETG3900", notwant = "SETG3910"),
    # GAMMA before ALPHA
    list(parms = "MUL=0.2 TAU=1 GAMMA=0 FIXGAMMA ALPHA=-1 FIXALPHA ETA=1",
         want = "SETG3910", notwant = "SETG3920"),
    # ALPHA before ETA
    list(parms = "MUL=0.2 TAU=1 GAMMA=2 ALPHA=-1 FIXALPHA ETA=0 FIXETA",
         want = "SETG3920", notwant = "SETG3930"),
    # ETA before the WEIBULL branch
    list(parms = "MUL=0.2 TAU=1 GAMMA=2 ETA=0 FIXETA ALPHA=0 FIXALPHA WEIBULL",
         want = "SETG3930", notwant = "SETG3980"),
    # the reported case: an entry refusal before the FIXGAE2 alpha rule
    list(parms = paste("MUL=0.2 TAU=1 GAMMA=0 FIXGAMMA ETA=1 ALPHA=0",
                       "FIXALPHA FIXGAE2 WEIBULL"),
         want = "SETG3910", notwant = "SETG3980")
  )
  for (p in pairs) {
    job <- .u1_job(parms = p$parms)
    msg <- .u1_msg(job)
    expect_match(msg, p$want, fixed = TRUE, info = p$parms)
    # The point is not only that the right code appears, but that the LATER
    # one does not: flag_refusal() keeps the first reason, so a wrong order
    # shows up as the later code appearing instead.
    expect_no_match(msg, p$notwant, fixed = TRUE, info = p$parms)
  }
})

test_that("the refusal message's claim about hzr_phase() matches what happens (#433 review)", {
  skip_on_cran()
  # The sentence used to assert that hzr_phase() accepts the shape, from a
  # hand-maintained idea of which codes were "shape" refusals. It was false
  # for FIVE of the seven classes, not the three the text claimed: SAS
  # refuses several of these precisely BECAUSE a shape is out of range, and
  # the same value is out of range for hzr_phase().
  #
  # So the message is now derived by CONSTRUCTING the phase, and this test
  # executes each class's emitted chunks and requires the message and the
  # outcome to agree. No list of codes appears in either, so neither can
  # drift from the other.
  set.seed(3)
  n <- 80
  D <- data.frame(TT = stats::rexp(n, 0.2),
                  DEAD = rep(c(1, 1, 0), length.out = n))
  classes <- c(
    SETG3900 = "MUL=0.1 TAU=0 GAMMA=2 ETA=1 ALPHA=1 FIXTAU",
    SETG3910 = "MUL=0.1 TAU=8 GAMMA=0 ETA=1 ALPHA=1 FIXGAMMA",
    SETG3920 = "MUL=0.1 TAU=8 GAMMA=2 ETA=1 ALPHA=-1 FIXALPHA",
    SETG3930 = "MUL=0.1 TAU=8 GAMMA=2 ETA=0 ALPHA=1 FIXETA",
    SETG3960 = "MUL=0.2 TAU=1 GAMMA=0 ETA=0.25 FIXGE2 WEIBULL",
    SETG3970 = "MUL=0.2 TAU=1 GAMMA=2 ETA=0 WEIBULL",
    SETG3980 = "MUL=0.2 TAU=1 GAMMA=4 ETA=0.25 ALPHA=0 WEIBULL"
  )
  agreed <- 0L
  says_yes <- 0L
  says_no <- 0L
  for (nm in names(classes)) {
    job <- .u1_job(parms = classes[[nm]])
    msg <- .u1_msg(job)
    claims_accepts <- grepl("accepts this shape", msg, fixed = TRUE)
    env <- new.env(parent = environment())
    env$D <- D
    completes <- tryCatch({
      for (k in names(job$calls)) suppressWarnings(eval(job$calls[[k]], env))
      TRUE
    }, error = function(e) FALSE)
    expect_identical(claims_accepts, completes, info = nm)
    agreed <- agreed + 1L
    if (claims_accepts) says_yes <- says_yes + 1L else says_no <- says_no + 1L
  }
  expect_identical(agreed, length(classes))
  # BOTH outcomes must occur, or an implementation that always said one thing
  # would satisfy every assertion above.
  expect_gt(says_yes, 0L)
  expect_gt(says_no, 0L)
})

test_that("a macro call whose arguments contain spaces stays one operand (#433 review)", {
  # Operands are split on whitespace, so `%FLAGS(A, B)` became `%FLAGS(A,`
  # and `B)`. Only the first looked like a macro; the remainder was judged on
  # its own, giving a job SAS runs a false $untranslated row and a false
  # warning. SAS expands the whole call before PROC HAZARD reads any operand,
  # so it must travel as one token and stay indeterminate.
  one_row <- function(parms) {
    job <- .u1_job(parms = parms)
    job$untranslated
  }
  # The reported two-argument case, and its no-space control.
  u <- one_row("MUE=0.2 THALF=1 NU=1 %FLAGS(A, B)")
  expect_identical(NROW(u), 1L)
  expect_identical(u$construct, "%FLAGS(A, B)")
  u <- one_row("MUE=0.2 THALF=1 NU=1 %FLAGS(A,B)")
  expect_identical(NROW(u), 1L)
  expect_identical(u$construct, "%FLAGS(A,B)")
  # Three arguments, so the joiner is not special-cased to one space.
  u <- one_row("MUE=0.2 THALF=1 NU=1 %F(A, B, C)")
  expect_identical(NROW(u), 1L)
  expect_identical(u$construct, "%F(A, B, C)")

  # OVER-REACH CONTROL, the direction a joiner fails in: an operand AFTER the
  # macro must still be read. Without this the test could not tell a correct
  # join from one that swallowed the rest of the statement.
  job <- .u1_job(parms = "MUE=0.2 THALF=1 %F(A, B) NU=1")
  expect_identical(NROW(job$untranslated), 1L)
  expect_true(any(grepl("nu = 1", deparse(job$calls$fit), fixed = TRUE)))

  # A non-macro token carrying parentheses is NOT joined: it is not a macro,
  # and reading unreadable phase text as a variable is tracked by #440.
  u <- one_row("MUE=0.2 THALF=1 NU=1 LOG(A, B)")
  expect_gt(NROW(u), 1L)
})

test_that("a rejected PROC option gets exactly one applicable row (#433 review)", {
  # check_number() recorded the rejection, then the option's own switch arm
  # continued and added a second row. `CONDITION=5.` said both that PROC
  # HAZARD's lexer rejects the number AND what its optimizer does with the
  # value, although a rejected job never runs. `MAXITER =` also left an
  # operand whose key was the empty string, reported as an unknown option
  # with a blank name.
  rows <- function(proc) {
    .u1_job(proc = proc, parms = "MUE=0.2 THALF=1 NU=1")$untranslated
  }
  for (proc in c(" MAXITER=", " MAXITER =", " CONDITION=5.", " MAXITER=1E5",
                 " CONDITION=")) {
    u <- rows(proc)
    expect_identical(NROW(u), 1L, info = proc)
    # No construct may be blank: that is the dangling half of a spaced
    # assignment, not an option anybody wrote.
    expect_true(all(nzchar(u$construct)), info = proc)
    # And no row may describe what the optimizer does with a value in a job
    # PROC HAZARD does not run.
    expect_false(any(grepl("stops PROC HAZARD's optimizer", u$reason,
                           fixed = TRUE)), info = proc)
  }

  # CONTROLS, so the test cannot pass by suppressing rows generally.
  # An ACCEPTED CONDITION still records its one explanatory row, because
  # hazard() reads no `condition` (#384).
  u <- rows(" CONDITION=14")
  expect_identical(NROW(u), 1L)
  expect_match(u$reason, "stops PROC HAZARD's optimizer", fixed = TRUE)
  # An accepted MAXITER records nothing and reaches the emitted call.
  job <- .u1_job(proc = " MAXITER=250", parms = "MUE=0.2 THALF=1 NU=1")
  expect_identical(NROW(job$untranslated), 0L)
  expect_true(any(grepl("maxit = 250", deparse(job$calls$fit), fixed = TRUE)))
})

test_that("a valueless option does not swallow the option after it (#433 review)", {
  # The spaced-operand joiner reads `MUE= 0.2` as one operand. It must not
  # read `DATA= MAXITER=50` the same way: the second token is another OPTION,
  # not this one's value. Joining them emitted
  # hazard(data = `MAXITER=50`) with NO row and NO warning, and silently lost
  # MAXITER -- a fit for a job PROC HAZARD rejects (`dsfield : NAME`,
  # hazard_y.y:82-84), which is the defect this whole branch exists to stop.
  ops <- .hzr_sas_join_spaced(c("DATA=", "MAXITER=50"))
  expect_equal(ops, c("DATA=", "MAXITER=50"))
  # The same guard on the bare-`=` spelling.
  expect_equal(.hzr_sas_join_spaced(c("DATA", "=", "MAXITER=50")),
               c("DATA=", "MAXITER=50"))
  # KNOWN NEGATIVE: a genuine spaced value must still join, or the guard has
  # simply disabled the feature it is protecting.
  expect_equal(.hzr_sas_join_spaced(c("MUE=", "0.2")), "MUE=0.2")
  expect_equal(.hzr_sas_join_spaced(c("MUE", "=", "0.2")), "MUE=0.2")

  # End to end, the property that actually matters: the job must NOT silently
  # emit a fit. On this branch before the guard it produced a clean
  # hazard(data = `MAXITER=50`) with no row and no warning. It now errors,
  # which is what main does -- loud, and therefore acceptable. Note the
  # error is RESTORED by this fix, not unchanged across it: the silent fit
  # was this branch's own regression, so "as it does on main" was true of
  # main and false of this branch's base (#433 review 2). Turning it into a
  # proper U1 refusal is a separate, LOUD leftover (see the leftovers issue):
  # `DATA=` with no NAME is a syntax error at `dsfield : NAME`
  # (hazard_y.y:80-82), so the job is one PROC HAZARD rejects.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste0("%HAZARD( PROC HAZARD DATA= MAXITER=50; EVENT DEAD;",
                    " TIME TT; PARMS MUE=0.2 THALF=1 NU=1; );"), f)
  out <- tryCatch(suppressWarnings(hzr_translate_sas(f)), error = function(e) e)
  silent_fit <- !inherits(out, "error") &&
    identical(out$calls$fit[[3L]][[1L]], as.name("hazard")) &&
    NROW(out$untranslated) == 0L &&
    !length(grep("^refusal", names(out$calls)))
  expect_false(silent_fit)
})

test_that("an entry refusal suppresses the 'cannot tell' verdict (#433 review)", {
  # SETG3's entry checks (setg3.c:269-284) each `return` immediately, so the
  # untraced non-WEIBULL dispatch (setg3.c:359-374) is NEVER reached once one
  # of them fires. A job carrying both therefore has exactly one true verdict:
  # PROC HAZARD refuses it at entry. Emitting "cannot tell whether PROC HAZARD
  # refuses this job" alongside "refused before any fit is computed" told the
  # reader both that SAS produces nothing and that we cannot say.
  job <- .u1_job(parms = "MUL=0.2 TAU=0 FIXTAU GAMMA=2 ETA=1 FIXGE2")
  m <- .u1_msg(job)
  # The entry refusal must SURVIVE -- suppressing the contradiction must not
  # suppress the verdict with it.
  expect_match(m, "(SETG3900)", fixed = TRUE)
  expect_match(m, "refused before any fit is computed", fixed = TRUE)
  expect_no_match(m, "cannot tell whether PROC HAZARD refuses", fixed = TRUE)

  # KNOWN NEGATIVE: with no entry refusal, the untraced path must still say
  # it cannot tell. Otherwise the fix has deleted the honest verdict too.
  plain <- .u1_job(parms = "MUL=0.1 TAU=8 ALPHA=2 GAMMA=5 ETA=1 FIXGAE2")
  pm <- .u1_msg(plain)
  expect_match(pm, "cannot tell whether PROC HAZARD refuses", fixed = TRUE)
  expect_no_match(pm, "refused before any fit is computed", fixed = TRUE)
})

test_that("two refusal reasons are separated in the emitted warning (#433 review)", {
  # warning(a, b) pastes its arguments with NO separator, so a job carrying
  # two refusal classes rendered as "...fit the model by hand.This translation
  # cannot emit...". The test helper hid it: it collapses the captured
  # messages itself, and a single warning() call yields ONE message however
  # many pieces were pasted into it.
  job <- .u1_job(parms = "MUE=0.2 THALF=1 NU=1E-3 M=1 FIXMNU1")
  m <- .u1_msg(job)
  expect_no_match(m, "[a-z]\\.[A-Z]")
})

test_that("a valueless option does not swallow a SPACED following option (#433 review 2)", {
  # Round 1's guard asked whether the next token CONTAINS `=`. That is the
  # wrong thing to index on: written fully spaced, the next token is a bare
  # keyword with no `=` in it, so the guard passed and the joiner swallowed
  # the option AND its value. Vary the spacing, which is what the mechanism
  # actually depends on, not just the content.
  expect_equal(.hzr_sas_join_spaced(c("DATA", "=", "MAXITER", "=", "50")),
               c("DATA=", "MAXITER=50"))
  expect_equal(.hzr_sas_join_spaced(c("DATA=", "MAXITER", "=", "50")),
               c("DATA=", "MAXITER=50"))
  expect_equal(.hzr_sas_join_spaced(c("DATA", "=", "MAXITER=", "50")),
               c("DATA=", "MAXITER=50"))
  # KNOWN NEGATIVES: every genuine spaced value must still join, at the
  # start, middle and END of the list (the end is where `i < n` stops
  # applying, so it is its own case).
  expect_equal(.hzr_sas_join_spaced(c("MUE", "=", "0.2", "THALF", "=", "0.3")),
               c("MUE=0.2", "THALF=0.3"))
  expect_equal(.hzr_sas_join_spaced(c("MUE=", "0.2", "NU", "=", "1")),
               c("MUE=0.2", "NU=1"))
  expect_equal(.hzr_sas_join_spaced(c("FIXNU", "MUE", "=", "0.2")),
               c("FIXNU", "MUE=0.2"))

  # End to end: the fully-spaced job must not emit a clean fit.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste0("%HAZARD( PROC HAZARD DATA = MAXITER = 50; EVENT DEAD;",
                    " TIME TT; PARMS MUE=0.2 THALF=1 NU=1; );"), f)
  out <- tryCatch(suppressWarnings(hzr_translate_sas(f)), error = function(e) e)
  silent_fit <- !inherits(out, "error") &&
    identical(out$calls$fit[[3L]][[1L]], as.name("hazard")) &&
    NROW(out$untranslated) == 0L &&
    !length(grep("^refusal", names(out$calls)))
  expect_false(silent_fit)
})
