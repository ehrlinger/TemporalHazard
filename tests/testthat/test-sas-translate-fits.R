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
              "MUE=0.2 THALF=0.5 M=-1 NU=1 FIXM FIXNU FIXMNU1")) {
    job <- .u1_job(parms = p)
    expect_true(.u1_refuses(job), info = p)
    msg <- .u1_msg(job)
    expect_match(msg, "cannot emit PROC HAZARD's model", fixed = TRUE, info = p)
    expect_no_match(msg, "runs this job", fixed = TRUE, info = p)
  }
  # M = NU = -1, both fixed, IS SETG1920, and the HAZARD binary refuses it
  # (fixtures/setg1-oracle.csv), so it now says so rather than that it
  # cannot tell (#424).
  msg <- .u1_msg(.u1_job(
    parms = "MUE=0.2 THALF=0.5 M=-1 NU=-1 FIXM FIXNU FIXMNU1"))
  expect_match(msg, "(SETG1920)", fixed = TRUE)
  expect_no_match(msg, "cannot emit PROC HAZARD's model", fixed = TRUE)
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

test_that("operand joining is INVARIANT under spacing (#433 review 2)", {
  # PROC HAZARD's lexer is whitespace-insensitive (hazard_l.l:32, :55), so
  # every spacing of one statement is the same token stream to SAS. Two
  # earlier rounds of this review fixed one spelling each and each time
  # another slipped through, because the set of spellings cannot be
  # enumerated by a fix. Assert the PROPERTY instead: generate every spacing
  # of each statement and require one answer.
  spacings <- function(pairs) {
    # every `=` with or without a space on each side
    grid <- expand.grid(rep(list(c("", " ")), 2L * length(pairs)),
                        stringsAsFactors = FALSE)
    out <- character(nrow(grid))
    for (r in seq_len(nrow(grid))) {
      s <- ""
      for (k in seq_along(pairs)) {
        l <- grid[[2L * k - 1L]][r]
        rgt <- grid[[2L * k]][r]
        s <- paste0(s, if (nzchar(s)) " " else "",
                    pairs[[k]][[1L]], l, "=", rgt, pairs[[k]][[2L]])
      }
      out[[r]] <- s
    }
    unique(out)
  }
  statements <- list(
    list(c("DATA", "MAXITER"), c("", "50")),       # the defect: key as value
    list(c("MUE", "0.2"), c("THALF", "0.3")),      # two genuine operands
    list(c("MAXITER", "250")),                     # one genuine operand
    # A MACRO operand. SAS expands `&N` before the lexer runs, so every
    # spelling is one job that RUNS -- but a macro is opaque to the
    # tokeniser, and when the `=` is glued to it (`MAXITER =&N`) the whole
    # `=&N` tested as a macro and was never split. That was the last branch
    # keying on spacing, and it produced a FALSE refusal (#433 review 3).
    list(c("MAXITER", "&N")),
    list(c("DATA", "&LIB")),
    list(c("MUE", "%N(1)"))
  )
  for (st in statements) {
    variants <- spacings(st)
    expect_gt(length(variants), 1L)                # the generator must vary
    joined <- lapply(variants, function(v) {
      .hzr_sas_join_spaced(strsplit(trimws(v), "[ \t]+")[[1L]])
    })
    # ONE answer for all spellings of this statement.
    expect_length(unique(joined), 1L)
  }

  # And the consequence that matters: no spelling may emit a `data` argument
  # that is itself an option. This is the assertion the two spelling-specific
  # tests were each half of.
  for (v in spacings(list(c("DATA", "MAXITER"), c("", "50")))) {
    f <- withr::local_tempfile(fileext = ".sas")
    writeLines(paste0("%HAZARD( PROC HAZARD ", v, "; EVENT DEAD; TIME TT;",
                      " PARMS MUE=0.2 THALF=1 NU=1; );"), f)
    out <- tryCatch(suppressWarnings(hzr_translate_sas(f)), error = function(e) e)
    if (inherits(out, "error")) next          # loud is acceptable
    d <- out$calls$fit[[3L]]$data
    expect_false(is.name(d) && grepl("=", as.character(d), fixed = TRUE),
                 info = v)
    silent <- identical(out$calls$fit[[3L]][[1L]], as.name("hazard")) &&
      NROW(out$untranslated) == 0L &&
      !length(grep("^refusal", names(out$calls)))
    expect_false(silent, info = v)
  }
})

test_that("a name-valued PROC option with no value is refused (#433 review 2)", {
  # `DATA '=' dsfield` and `OUTHAZ '=' dsfield`, dsfield : NAME | LIBMEM
  # (hazard_y.y:61-62, :80-81). Neither has a form without a name, so an
  # empty value is a syntax error and the job does not run. Only MAXITER and
  # CONDITION had a presence check; OUTHAZ= was dropped silently and the job
  # fitted, and DATA= surfaced as an internal R error naming neither.
  #
  # The option must be LAST to be genuinely valueless. `OUTHAZ= MAXITER=50`
  # is NOT this case: <HZRP>OUTHAZ switches the lexer to DSNM, where MAXITER
  # lexes as a NAME (hazard_l.l:60,80), so SAS reads OUTHAZ=MAXITER and then
  # a stray `=`. An earlier draft of this test asserted that spelling and was
  # wrong about SAS, not about the code.
  for (opt in c("DATA", "OUTHAZ")) {
    f <- withr::local_tempfile(fileext = ".sas")
    writeLines(paste0("%HAZARD( PROC HAZARD DATA=D MAXITER=50 ", opt, "=;",
                      " EVENT DEAD; TIME TT; PARMS MUE=0.2 THALF=1 NU=1; );"), f)
    job <- suppressWarnings(hzr_translate_sas(f))
    expect_false(is.null(.u1_refusal_chunk(job)), info = opt)
    expect_true(any(grepl(opt, job$untranslated$construct, fixed = TRUE)),
                info = opt)
    # An earlier option on the same line is still read: a refusal must not
    # eat the whole statement.
    expect_identical(job$calls$fit[[3L]]$control$maxit, 50, info = opt)
  }
  # KNOWN NEGATIVE: real values refuse nothing.
  f2 <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste0("%HAZARD( PROC HAZARD DATA=D OUTHAZ=H MAXITER=50;",
                    " EVENT DEAD; TIME TT; PARMS MUE=0.2 THALF=1 NU=1; );"), f2)
  clean <- suppressWarnings(hzr_translate_sas(f2))
  expect_null(.u1_refusal_chunk(clean))
})

test_that("a libref with no member is refused, not left empty (#433 review 3)", {
  # `WORK.` is neither NAME nor LIBMEM (hazard_l.l:39-40), so SAS rejects it.
  # The presence check ran BEFORE the WORK. strip, so "WORK." passed it and
  # the strip then left an empty name, surfacing as an internal
  # "attempt to use zero-length variable name" that named neither the option
  # nor the reason.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste0("%HAZARD( PROC HAZARD DATA=WORK. MAXITER=50; EVENT DEAD;",
                    " TIME TT; PARMS MUE=0.2 THALF=1 NU=1; );"), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_false(is.null(.u1_refusal_chunk(job)))
  expect_true(any(grepl("DATA", job$untranslated$construct, fixed = TRUE)))
  # KNOWN NEGATIVE: a real WORK-qualified name still translates, stripped.
  f2 <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste0("%HAZARD( PROC HAZARD DATA=WORK.D MAXITER=50; EVENT DEAD;",
                    " TIME TT; PARMS MUE=0.2 THALF=1 NU=1; );"), f2)
  ok <- suppressWarnings(hzr_translate_sas(f2))
  expect_null(.u1_refusal_chunk(ok))
  expect_identical(ok$calls$fit[[3L]]$data, as.name("D"))
})

test_that("PROC HAZPRED refuses a stray `=` and an empty name (#433 review 3)", {
  # The HAZPRED caller shares the joiner but had neither the stray-`=` block
  # nor the presence check, so `DATA=G = INHAZ=H` recorded a BLANK-keyword
  # "unknown option" row and emitted the prediction calls anyway, for a job
  # SAS rejects.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste0("%HAZARD( PROC HAZARD DATA=D; EVENT DEAD; TIME TT;",
                    " PARMS MUE=0.2 THALF=1 NU=1; );\n",
                    "%HAZPRED( PROC HAZPRED DATA=PGRID = INHAZ=H OUT=P;",
                    " TIME TT; );"), f)
  # A %HAZPRED block is folded into the SAME hzr_sas_job as the %HAZARD one
  # it predicts from, so this is one job, not two.
  hp <- suppressWarnings(hzr_translate_sas(f))
  expect_s3_class(hp, "hzr_sas_job")
  expect_false(any(!nzchar(hp$untranslated$construct)))
  expect_true(any(grepl("stray", hp$untranslated$reason, fixed = TRUE)))
})

# --- #411: a SAS covariate name that begins with an underscore ---

test_that("an underscore-named covariate fits end to end (#411)", {
  skip_on_cran()
  set.seed(7)
  n <- 120
  # check.names = FALSE: the column really is named `_X1`, as the SAS dataset
  # named it. data.frame()'s default would rename it (asserted below).
  D <- data.frame(T = stats::rexp(n, 0.2), E = rep(c(1, 1, 0), length.out = n),
                  AGE = stats::rnorm(n), check.names = FALSE)
  D[["_X1"]] <- stats::rnorm(n)
  expect_true("_X1" %in% names(D))

  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(c("%HAZARD( PROC HAZARD DATA=D; TIME T; EVENT E;",
               " PARMS MUE=0.2 THALF=1 NU=1 M=1 MUC=0.01;",
               " EARLY AGE=0.1, _X1=0.1; );"), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_true(any(grepl("`_X1`", deparse(job$calls$fit), fixed = TRUE)))

  # Execute the emitted chunks, not their text.
  env <- new.env(parent = environment())
  env$D <- D
  for (nm in names(job$calls)) suppressWarnings(eval(job$calls[[nm]], env))
  expect_s3_class(env$fit, "hazard")
  # The covariate is actually IN the model, not merely in the formula text.
  # model.matrix() backquotes a non-syntactic name, so the coefficient is
  # named with the backquotes: assert the name the fit really carries.
  expect_true("phase_1.`_X1`" %in% names(stats::coef(env$fit)))
})

test_that("a renamed underscore column fails loudly, it does not fit a smaller model (#411)", {
  skip_on_cran()
  set.seed(7)
  n <- 120
  # data.frame()'s default check.names = TRUE turns `_X1` into `X_X1`. The
  # danger is a fit that quietly drops the covariate; it must refuse instead.
  D <- data.frame(T = stats::rexp(n, 0.2), E = rep(c(1, 1, 0), length.out = n),
                  AGE = stats::rnorm(n), "_X1" = stats::rnorm(n))
  expect_false("_X1" %in% names(D))
  expect_true("X_X1" %in% names(D))

  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(c("%HAZARD( PROC HAZARD DATA=D; TIME T; EVENT E;",
               " PARMS MUE=0.2 THALF=1 NU=1 M=1 MUC=0.01;",
               " EARLY AGE=0.1, _X1=0.1; );"), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  env <- new.env(parent = environment())
  env$D <- D
  err <- tryCatch({
    for (nm in names(job$calls)) suppressWarnings(eval(job$calls[[nm]], env))
    "no error"
  }, error = conditionMessage)
  expect_match(err, "_X1", fixed = TRUE)
  expect_match(err, "not columns of D", fixed = TRUE)
})

test_that("a SELECTION job carrying a non-syntactic name is refused, not screened (#411)", {
  skip_on_cran()
  # The phase formulas now carry `_X1`, but hzr_stepwise() spells such a name
  # two ways at once: backquoted in its terms() candidate labels, bare in
  # force_in. Measured consequences, which is why this is a refusal and not a
  # screen: a /I pin never matches, so BACKWARD DROPS a variable SAS holds in
  # with no warning naming it; and the score criterion, the only one this
  # translator emits, indexes `data` by the backquoted label and skips the
  # candidate. Both would be wrong models from a populated result.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(c("%HAZARD( PROC HAZARD DATA=D; TIME T; EVENT E;",
               " PARMS MUE=0.2 THALF=1 NU=1 M=1;",
               " EARLY AGE=0.1 /I, _X1=0.1 /I;",
               " SELECTION BACKWARD SLS=0.05; );"), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"))
  msg <- tryCatch({
    eval(job$calls$fit, new.env())
    "no error"
  }, error = conditionMessage)
  expect_match(msg, "_X1", fixed = TRUE)
  expect_match(msg, "not a syntactic R name", fixed = TRUE)
  expect_match(msg, "hazard_l.l:39", fixed = TRUE)

  # Control, so the refusal is not blanket: ordinary names still screen.
  f2 <- withr::local_tempfile(fileext = ".sas")
  writeLines(c("%HAZARD( PROC HAZARD DATA=D; TIME T; EVENT E;",
               " PARMS MUE=0.2 THALF=1 NU=1 M=1;",
               " EARLY AGE, SEX;",
               " SELECTION FORWARD SLE=0.3; );"), f2)
  job2 <- suppressWarnings(hzr_translate_sas(f2))
  expect_identical(job2$calls$fit[[3L]][[1L]], as.name("hzr_stepwise"))
})

test_that(".hzr_sas_covar_formula() refuses an empty name vector (#411)", {
  # Reduce() over an empty list is NULL and `~NULL` is a hollow formula: a
  # phase fitted with no covariates, and nothing would error.
  expect_error(.hzr_sas_covar_formula(character(0)), "at least one name")
  expect_identical(.hzr_sas_covar_formula("AGE"), quote(~AGE))
})

test_that("a SELECTION name PROC HAZARD accepts is refused; one it rejects is not (#411)", {
  # Only a name PROC HAZARD ACCEPTS is refused here. `_X1` is in its grammar
  # (hazard_l.l:39) and only R objects to it, and that job already died on
  # main, so refusing it is not a new stop.
  #
  # Text PROC HAZARD REJECTS at parse (`AGE*SEX`, `LOG(AGE)`) is NOT refused
  # here: stopping it would be a NEW stop for a job that translated on main.
  # It warns, with a row, and is left out of the screen (#440; tested below).
  job <- function(early) {
    f <- withr::local_tempfile(fileext = ".sas", .local_envir = parent.frame())
    writeLines(c("%HAZARD( PROC HAZARD DATA=D; TIME T; EVENT E;",
                 " PARMS MUE=0.2 THALF=1 NU=1 M=1;",
                 paste0(" EARLY ", early, ";"),
                 " SELECTION FORWARD SLE=0.3; );"), f)
    suppressWarnings(hzr_translate_sas(f))
  }
  msg <- function(j) {
    tryCatch({
      eval(j$calls$fit, new.env())
      "no error"
    }, error = conditionMessage)
  }

  # In the grammar: refused, and told the truth about why.
  for (nm in c("_X1", "NA", "TRUE")) {
    j <- job(paste0("AGE /I, ", nm))
    expect_identical(j$calls$fit[[3L]][[1L]], as.name("stop"), info = nm)
    expect_match(msg(j), "lexer accepts this name", fixed = TRUE, info = nm)
    expect_match(msg(j), "hazard_l.l:39", fixed = TRUE, info = nm)
  }

  # NOT in the grammar: NOT refused. The emitted call is still a screen --
  # no stop is added here.
  for (nm in c("AGE*SEX", "LOG(AGE)")) {
    j <- job(paste0("AGE /I, ", nm))
    expect_identical(j$calls$fit[[3L]][[1L]], as.name("hzr_stepwise"), info = nm)
    expect_false(grepl("not a syntactic R name", msg(j), fixed = TRUE), info = nm)
  }

  # Control: an ordinary pair still screens.
  expect_identical(job("AGE /I, SEX")$calls$fit[[3L]][[1L]], as.name("hzr_stepwise"))
})


# --- #440: a phase variable PROC HAZARD cannot lex as a NAME ---------------
# hazard_l.l:39 is `name ([_A-Z][_A-Z0-9]*)` and hazard_y.y:213 is
# `phasevar : NAME`, so `AGE*SEX`, `LOG(AGE)`, `B SEX` and `1AGE` are not
# phase variables: PROC HAZARD rejects the job at parse. On main the parser
# passed each through as a column name, with no row and no warning. Under U1
# such a job now warns, records one row, and emits a fit WITHOUT the operand,
# so the emitted formula and theta agree.
.p440_classes <- c(interaction = "AGE*SEX", call = "LOG(AGE)",
                   spaced = "B SEX", digit = "1AGE")
.p440_job <- function(early, selection = FALSE, env = parent.frame()) {
  f <- withr::local_tempfile(fileext = ".sas", .local_envir = env)
  writeLines(paste0("%HAZARD( PROC HAZARD DATA=D; TIME T; EVENT E;",
                    " PARMS MUE=0.2 THALF=1 NU=1 M=1 MUC=0.01;",
                    " EARLY ", early, ";",
                    if (selection) " SELECTION SLE=0.3 SLS=0.2;",
                    " );"), f)
  w <- character(0)
  withCallingHandlers(job <- hzr_translate_sas(f), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  job$translate_warnings <- w
  job
}
# The warning the DOCUMENT raises, by evaluating the emitted refusal chunk.
.p440_chunk_warning <- function(job) {
  nm <- grep("^refusal", names(job$calls), value = TRUE)
  if (!length(nm)) return(character(0))
  msgs <- character(0)
  withCallingHandlers(eval(job$calls[[nm[[1L]]]], new.env()),
                      warning = function(x) {
                        msgs <<- c(msgs, conditionMessage(x))
                        invokeRestart("muffleWarning")
                      })
  msgs
}
.p440_rows <- function(job) {
  u <- job$untranslated
  u[grepl("^(not a PROC HAZARD variable name|follows a `[(]`)", u$reason), ,
    drop = FALSE]
}
# Covariate TERMS in the emitted early-phase formula against covariate STARTS
# in the emitted theta (the early block is log_mu, log_t_half, nu, m; the
# constant phase adds log_mu). A formula that expands to more terms than it
# has starts is the defect #440 describes.
.p440_terms_vs_starts <- function(fit_call) {
  ph <- fit_call$phases[[2L]]
  n_terms <- if (is.null(ph$formula)) 0L else
    length(attr(stats::terms(eval(ph$formula)), "term.labels"))
  n_starts <- length(fit_call$theta) - 1L - 4L - 1L
  c(terms = n_terms, starts = n_starts)
}

test_that("a phase variable that is not a PROC HAZARD NAME warns, rows and drops (#440)", {
  base <- .p440_job("AGE=0.1")
  base_sel <- .p440_job("AGE=0.1", selection = TRUE)
  for (cls in names(.p440_classes)) {
    x <- .p440_classes[[cls]]
    for (sel in c(FALSE, TRUE)) {
      info <- paste(cls, if (sel) "with SELECTION" else "plain")
      b <- if (sel) base_sel else base
      job <- .p440_job(paste0("AGE=0.1, ", x, "=0.2"), selection = sel)
      rows <- .p440_rows(job)
      # Exactly one row, naming the construct and the grammar.
      expect_identical(NROW(rows), 1L, info = info)
      expect_identical(rows$construct, x, info = info)
      expect_match(rows$reason, "hazard_l.l:39", fixed = TRUE, info = info)
      expect_match(rows$reason, "hazard_y.y:213", fixed = TRUE, info = info)
      # And nothing else changed: every other row is the baseline's.
      expect_identical(NROW(job$untranslated), NROW(b$untranslated) + 1L,
                       info = info)
      # The translation warns, naming the construct.
      expect_true(any(grepl(x, job$translate_warnings, fixed = TRUE)),
                  info = info)
      # The document warns, naming the construct and the grammar reason,
      # and the warning carries the ROW's own construct and reason, so the
      # two cannot disagree about what was refused.
      w <- .p440_chunk_warning(job)
      expect_length(w, 1L)
      expect_match(w, "PROC HAZARD does not run this job", fixed = TRUE,
                   info = info)
      expect_match(w, paste0(rows$construct, ": ", rows$reason),
                   fixed = TRUE, info = info)
      # The fit is emitted, not a stop, and it is the baseline's model: the
      # operand is gone from the formula AND from theta.
      if (sel) {
        expect_identical(job$calls$fit[[3L]][[1L]], as.name("hzr_stepwise"),
                         info = info)
        expect_identical(job$calls$fit_base, b$calls$fit_base, info = info)
        expect_identical(job$calls$fit, b$calls$fit, info = info)
      } else {
        expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"),
                         info = info)
        expect_identical(job$calls$fit, b$calls$fit, info = info)
        tv <- .p440_terms_vs_starts(job$calls$fit[[3L]])
        expect_identical(tv[["terms"]], tv[["starts"]], info = info)
      }
    }
  }
})

test_that("names PROC HAZARD accepts do not warn (#440 known negatives)", {
  # From hazard_l.l:39, clause by clause: a leading underscore, a bare
  # underscore, digits after the first character, and words that are
  # keywords elsewhere but lex as NAME in the phase-variable state (PHVR has
  # no keyword rules, hazard_l.l:151-152 and :174).
  ok <- .hzr_parse_phase_covars("_X1, _, A1_2, E, I, S, EARLY, PARMS, age")
  expect_identical(ok$names, c("_X1", "_", "A1_2", "E", "I", "S", "EARLY",
                               "PARMS", "age"))
  expect_length(ok$untranslated_construct, 0L)
  expect_length(ok$not_a_name, 0L)
  # `(` returns no token (hazard_l.l:56) and `)` is whitespace (:32), so a
  # last item `LOG()` is the variable LOG, with or without a start value.
  for (s in c("LOG()", "LOG ( ) = 0.2", "LOG)", "LOG((")) {
    p <- .hzr_parse_phase_covars(paste0("AGE, ", s))
    expect_identical(p$names, c("AGE", "LOG"), info = s)
    expect_length(p$not_a_name, 0L)
  }
  expect_identical(.hzr_parse_phase_covars("AGE, LOG()=0.2")$values,
                   c(NA, 0.2))
  # Each clause of the regex, violated.
  bad <- .hzr_parse_phase_covars("AGE, 1AGE, A.B, A-B, A$B, A*B, F(X), B SEX")
  expect_identical(bad$names, "AGE")
  expect_identical(bad$untranslated_construct,
                   c("1AGE", "A.B", "A-B", "A$B", "A*B", "F(X)", "B SEX"))
  expect_length(bad$not_a_name, 7L)
  # Not phase-statement syntax errors of the #340 kind, so they do not stop.
  expect_length(bad$rejected, 0L)

  # End to end: no row, no refusal chunk, no warning.
  for (x in c("AGE=0.1", "AGE=0.1, _X1=0.2")) {
    job <- .p440_job(x)
    expect_identical(NROW(job$untranslated), 0L, info = x)
    expect_length(grep("^refusal", names(job$calls)), 0L)
    expect_length(job$translate_warnings, 0L)
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
  }
  # `_X1` under SELECTION keeps its own R-side refusal (#411); this check
  # adds nothing to it.
  job <- .p440_job("AGE=0.1, _X1=0.2", selection = TRUE)
  expect_identical(NROW(.p440_rows(job)), 0L)
  expect_match(job$untranslated$reason[job$untranslated$construct ==
                                          "SELECTION"],
               "not a syntactic R name", fixed = TRUE)
})

test_that("an item with no variable stops, as it did on main (#440)", {
  # `EARLY AGE, =0.2;` errored inside R on main ("attempt to use
  # zero-length variable name"). It is a parse error to PROC HAZARD too,
  # and it stays a stop -- now the #340 one, which names the source.
  p <- .hzr_parse_phase_covars("AGE, =0.2")
  expect_identical(p$names, "AGE")
  expect_length(p$not_a_name, 0L)
  expect_length(p$rejected, 1L)
  expect_match(p$rejected, "hazard_y.y:210-213", fixed = TRUE)
  job <- .p440_job("AGE, =0.2")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"))
})

test_that("a macro reference in a phase statement is not judged (#440)", {
  # SAS expands `&V` before PROC HAZARD's lexer sees it, so whether it is a
  # NAME cannot be known here.
  p <- .hzr_parse_phase_covars("AGE, &V")
  expect_length(p$not_a_name, 0L)
  expect_identical(p$names, c("AGE", "&V"))
})

test_that("the #440 verdict is invariant under spacing", {
  # PROC HAZARD's lexer skips whitespace (hazard_l.l:32), so a space changes
  # the token stream only where it separates two word-rule characters
  # (`B SEX` is two NAMEs, `BSEX` one). Generate every spelling of each
  # token stream: each free gap with and without a space, each required gap
  # with one space or two. This generalises #433's generator, which varies
  # only the gaps around `=`.
  word <- function(t) grepl("^[-._A-Z0-9]+$", t)
  spellings <- function(toks) {
    n <- length(toks) - 1L
    gaps <- lapply(seq_len(n), function(k) {
      if (word(toks[[k]]) && word(toks[[k + 1L]])) c(" ", "  ") else c("", " ")
    })
    grid <- expand.grid(gaps, stringsAsFactors = FALSE)
    apply(grid, 1L, function(g) {
      paste0(c(rbind(toks[-length(toks)], g), toks[[length(toks)]]),
             collapse = "")
    })
  }
  streams <- list(
    c("AGE", "*", "SEX", "=", "0.2"),
    c("LOG", "(", "AGE", ")", "=", "0.2"),
    c("B", "SEX", "=", "0.2"),
    c("1AGE", "=", "0.2")
  )
  for (toks in streams) {
    v <- unique(spellings(toks))
    expect_gt(length(v), 1L)          # the generator must vary
    # Parser level: every spelling, one answer.
    parsed <- lapply(v, function(s) {
      .hzr_parse_phase_covars(paste0("AGE=0.1, ", s))
    })
    expect_length(unique(parsed), 1L)
    expect_length(parsed[[1L]]$not_a_name, 1L)
    # Translation level: rows, calls and the document's warning.
    out <- lapply(v, function(s) {
      job <- .p440_job(paste0("AGE=0.1, ", s))
      list(rows = job$untranslated[, c("construct", "reason")],
           calls = job$calls[names(job$calls) != "data"],
           warning = .p440_chunk_warning(job))
    })
    expect_length(unique(out), 1L)
    expect_length(out[[1L]]$warning, 1L)
  }
})

test_that("the emitted #440 fit runs, with one start per covariate term", {
  skip_on_cran()
  set.seed(440)
  n <- 150
  D <- data.frame(AGE = stats::rnorm(n), SEX = stats::rbinom(n, 1, 0.5))
  D$T <- stats::rexp(n, 0.2)
  D$E <- stats::rbinom(n, 1, 0.7)
  run <- function(job) {
    env <- new.env(parent = environment())
    env$D <- D
    for (k in names(job$calls)) suppressWarnings(eval(job$calls[[k]], env))
    env
  }
  for (cls in names(.p440_classes)) {
    job <- .p440_job(paste0("AGE=0.1, ", .p440_classes[[cls]], "=0.2"))
    env <- run(job)
    fit <- env$fit
    expect_s3_class(fit, "hazard")
    # The fitted parameter vector is exactly as long as the emitted start,
    # and carries one covariate: a formula with more terms than starts
    # cannot reach here.
    expect_identical(length(fit$fit$theta),
                     length(eval(job$calls$fit[[3L]]$theta)), info = cls)
    expect_identical(sum(fit$fit$covariate_counts), 1L, info = cls)
  }
  # The SELECTION jobs emit the baseline's calls (asserted above for every
  # class), so one run of the screen covers them; this one is the #411 case.
  job <- .p440_job("AGE=0.1, AGE*SEX=0.2", selection = TRUE)
  utils::capture.output(env <- suppressMessages(run(job)))
  expect_s3_class(env$fit, "hzr_stepwise")
  expect_false(any(grepl("SEX", unlist(lapply(env$fit$scope, deparse)),
                         fixed = TRUE)))
})

test_that("the phase-name verdict matches the HAZARD binary on a grid (#440)", {
  # The oracle is the C binary itself, not a reading of its lexer: each row
  # records whether PROC HAZARD ran one EARLY statement and which phase
  # variables its listing reports. Provenance (binary version, file date,
  # source checkout) is in the fixture's header; regenerate it with
  # data-raw/phase-name-oracle.R.
  path <- test_path("fixtures", "phase-name-oracle.csv")
  oracle <- utils::read.csv(path, comment.char = "#", stringsAsFactors = FALSE)
  # Coverage before comparison: every verdict class is present, so the test
  # cannot pass over a grid that exercises only one side.
  expect_setequal(unique(oracle$verdict),
                  c("runs", "rejected", "runs_after_reset"))
  expect_gte(nrow(oracle), 20L)
  for (k in seq_len(nrow(oracle))) {
    o <- oracle[k, ]
    job <- .p440_job(o$early)
    rows <- .p440_rows(job)
    w <- .p440_chunk_warning(job)
    fc <- job$calls$fit[[3L]]
    ph <- fc$phases[[2L]]
    emitted <- if (is.null(ph$formula)) character(0) else
      all.vars(eval(ph$formula))
    if (o$verdict == "runs") {
      # SAS runs it: no warning, and the SAME variables in the SAME order.
      expect_identical(NROW(rows), 0L, info = o$early)
      expect_length(w, 0L)
      expect_identical(emitted, strsplit(o$early_vars, " ")[[1L]],
                       info = o$early)
    } else {
      # SAS rejects it (or runs it only because a later `(` cleared the
      # error): warned, recorded, and never an operand SAS did not read.
      expect_gte(NROW(rows), 1L)
      expect_length(w, 1L)
      expect_identical(as.character(fc[[1L]]), "hazard", info = o$early)
      expect_true(all(emitted %in% c("AGE", "SEX", "LOG")), info = o$early)
      if (o$verdict == "runs_after_reset") {
        expect_match(w, "clears its syntax-error flag", fixed = TRUE,
                     info = o$early)
      }
    }
  }
})

# --- SETG1: the early phase's refusals, rewrites and no-result case (#424) --
# The oracle is the HAZARD binary itself (data-raw/setg1-oracle.R): each row
# records whether PROC HAZARD refused one early-phase PARMS statement, ran it,
# or produced no result, with the data staged as PROC HAZARD stages it.

.p424_oracle <- function() {
  path <- test_path("fixtures", "setg1-oracle.csv")
  utils::read.csv(path, comment.char = "#", stringsAsFactors = FALSE)
}
# The jobs PROC HAZARD runs only after SETG1 has moved a starting value
# (setg1.c:342-349, :614-617, :686-691, :763-767): one row, no warning.
.p424_moved <- c("thalf_neg_free", "thalf_zero_free", "thalf_half_free",
                 "mnu_zero_free", "nu_zero_fixm")
.p424_class <- function(job) {
  msg <- .u1_msg(job)
  if (grepl("refused before any fit is computed", msg, fixed = TRUE)) {
    "refused"
  } else if (grepl("cannot emit PROC HAZARD's model", msg, fixed = TRUE)) {
    # A different model, which SAS RUNS. Checked before "no result",
    # because this message quotes the row's reason, and a no-result reason
    # routed here would otherwise read as a no-result verdict.
    "runs"
  } else if (grepl("produces no result", msg, fixed = TRUE)) {
    "no_result"
  } else {
    "runs"
  }
}
.p424_warnings <- function(job) {
  ch <- .u1_refusal_chunk(job)
  if (is.null(ch)) return(character(0))
  msgs <- character(0)
  withCallingHandlers(eval(ch, new.env()), warning = function(x) {
    msgs <<- c(msgs, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  msgs
}

test_that("the SETG1 verdict matches the HAZARD binary on a grid (#424)", {
  oracle <- .p424_oracle()
  # Coverage before comparison: every class, and the known positive, present.
  expect_setequal(unique(oracle$binary), c("runs", "refused", "no_result"))
  expect_gte(nrow(oracle), 20L)
  expect_identical(oracle$binary[oracle$id == "known_positive"], "runs")
  for (k in seq_len(nrow(oracle))) {
    o <- oracle[k, ]
    job <- .u1_job(parms = o$parms)
    expect_identical(.p424_class(job), o$binary, info = o$id)
    w <- .p424_warnings(job)
    n_rows <- NROW(job$untranslated)
    expect_identical(as.character(job$calls$fit[[3L]][[1L]]), "hazard",
                     info = o$id)
    if (o$binary == "runs") {
      expect_length(w, 0L)
      expect_identical(n_rows, if (o$id %in% .p424_moved) 1L else 0L,
                       info = o$id)
    } else {
      # Warned once, recorded once.
      expect_length(w, 1L)
      expect_identical(n_rows, 1L, info = o$id)
    }
    if (startsWith(o$id, "SETG")) {
      expect_match(paste(w, collapse = " "),
                   paste0("(", sub("_.*$", "", o$id), ")"), fixed = TRUE,
                   info = o$id)
    }
  }
})

test_that("a free non-positive THALF starts at 1, as SETG1 does (#421 item 3)", {
  for (p in c("MUE=0.2 THALF=-1 NU=1 M=1", "MUE=0.2 THALF=0 NU=1 M=1",
              "MUE=0.2 THALF=-.5")) {
    job <- .u1_job(parms = p)
    ph <- job$calls$fit[[3L]]$phases[[2L]]
    expect_identical(ph$t_half, 1, info = p)
    # The theta block is built from the same value: log(1), finite.
    th <- eval(job$calls$fit[[3L]]$theta)
    expect_true(all(is.finite(th)), info = p)
    expect_identical(NROW(job$untranslated), 1L)
    expect_match(job$untranslated$reason, "setg1.c:343-349", fixed = TRUE)
    expect_null(.u1_refusal_chunk(job))
  }
  # Fixed, the same value is SETG1910, and the warning says so.
  job <- .u1_job(parms = "MUE=0.2 THALF=-1 FIXTHALF NU=1 M=1")
  expect_match(.u1_msg(job), "(SETG1910)", fixed = TRUE)
})

test_that("SETG1's g1flag 4 with a free M is NO RESULT, not another model (#424)", {
  # PROC HAZARD does not fit a different model here: SETG1 selects the
  # limiting positive generic case (setg1.c:763-776) and the fit then stops
  # on a domain error (DLG1980), with no estimates. Calling it not_mirrored
  # would be a false statement about SAS.
  for (p in c("MUE=0.2 THALF=1 M=1 NU=0", "MUE=0.2 THALF=1 M=-1 NU=0",
              "MUE=0.2 THALF=1 M=0 NU=0 FIXNU")) {
    job <- .u1_job(parms = p)
    msg <- .u1_msg(job)
    expect_match(msg, "produces no result", fixed = TRUE, info = p)
    expect_no_match(msg, "cannot emit PROC HAZARD's model", fixed = TRUE,
                    info = p)
    expect_no_match(msg, "refused before any fit", fixed = TRUE, info = p)
  }
  # M fixed: SAS runs it, so nothing is said.
  job <- .u1_job(parms = "MUE=0.2 THALF=1 M=-1 NU=0 FIXM FIXNU")
  expect_null(.u1_refusal_chunk(job))
})

test_that("the SETG1 verdict is invariant under spacing (#424)", {
  oracle <- .p424_oracle()
  oracle <- oracle[oracle$binary != "runs" | oracle$id %in% .p424_moved, ]
  expect_gte(nrow(oracle), 10L)
  for (k in seq_len(nrow(oracle))) {
    o <- oracle[k, ]
    base <- .u1_job(parms = o$parms)
    for (sp in c(" = ", " =", "= ")) {
      v <- gsub("=", sp, o$parms, fixed = TRUE)
      expect_false(identical(v, o$parms))          # the spelling varies
      job <- .u1_job(parms = v)
      expect_identical(.p424_class(job), o$binary, info = v)
      expect_identical(job$untranslated, base$untranslated, info = v)
      expect_identical(.p424_warnings(job), .p424_warnings(base), info = v)
    }
  }
})

test_that("the SETG1 documents render past the warning (#424)", {
  skip_on_cran()
  e <- new.env()
  utils::data("avc", package = "TemporalHazard", envir = e)
  a <- e$avc[stats::complete.cases(e$avc), ]
  D <- data.frame(TT = a$int_dead, DEAD = a$dead)
  oracle <- .p424_oracle()
  # Rows whose operands are ALSO outside what hzr_phase() or the likelihood
  # accepts: the fit halts after the warning, as for SETG3910-3930 above.
  halts <- c("SETG1910_neg", "SETG1910_zero", "SETG1940", "SETG1950",
             "SETG1960", "SETG1920", "SETG1930", "g1flag4_mpos_fixnu")
  for (id in c(halts, "SETG1900", "SETG1901", .p424_moved, "g1flag4_mpos",
               "g1flag4_mneg", "g1flag4_mzero_fixnu")) {
    job <- .u1_job(parms = oracle$parms[oracle$id == id])
    res <- suppressWarnings(render_sim(job, list(D = D)))
    expect_identical(res$ok, !(id %in% halts), info = id)
    # The warning chunk itself always ran.
    rn <- grep("^refusal", names(res$results), value = TRUE)
    if (length(rn)) expect_identical(unname(res$results[[rn]]), "ok", info = id)
  }
})
