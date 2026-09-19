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
  # so it is now a syntax stop (U1). An operand written with spaces around
  # `=` is the parser's own limit: SAS accepts it, this parser cannot read it.)
  for (parms in c("PARMS MUE = 0.2 THALF = 1 NU = 1;")) {
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

test_that("every reachable SETG3 refusal renders a stop(), not a fit (#359)", {
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

  for (label in names(refusals)) {
    code <- sub("_gae2$", "", label)
    job <- translate(refusals[[label]])
    res <- suppressWarnings(render_sim(job, list(D = D)))
    expect_false(res$ok, info = label)
    expect_match(res$results[["fit"]], code, fixed = TRUE, info = label)
    expect_false(exists("fit", envir = res$env, inherits = FALSE), info = label)
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
    msg <- tryCatch({
      eval(job$calls$fit, new.env())
      "no error"
    }, error = conditionMessage)
    expect_match(msg, "cannot tell whether PROC HAZARD refuses", fixed = TRUE, info = ops)
    expect_no_match(msg, "refused before any fit", fixed = TRUE, info = ops)
  }
})

# --- U1: a job PROC HAZARD refuses, or fits differently, emits stop() -------
# John's decision (2026-09-19): when the translator knows PROC HAZARD refuses
# a job, or fits a model other than the one it would emit, the document stops
# rather than fitting with a row the reader may never see.

.u1_job <- function(proc = "", parms, env = parent.frame()) {
  f <- withr::local_tempfile(fileext = ".sas", .local_envir = env)
  writeLines(paste0("%HAZARD( PROC HAZARD DATA=D", proc,
                    "; EVENT DEAD; TIME TT; PARMS ", parms, "; );"), f)
  suppressWarnings(hzr_translate_sas(f))
}
.u1_msg <- function(job) {
  tryCatch({
    eval(job$calls$fit, new.env())
    "no error"
  }, error = conditionMessage)
}

test_that("a PARMS syntax error stops the document (U1, #421)", {
  for (p in c("MUE=0.2 THALF=0.5 NU=1E-3", "MUE=0.2 THALF=0.5 NU",
              "MUE=0.2 THALF=0.5 NU = ABC", "MUE=0.2 THALF=0.5 FIXG1")) {
    job <- .u1_job(parms = p)
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"), info = p)
    expect_match(.u1_msg(job), "PROC HAZARD does not run this job", fixed = TRUE,
                 info = p)
  }
  # Controls: a job SAS runs still fits, including a macro piece this parser
  # cannot judge, which is indeterminate and so not a known refusal.
  for (p in c("MUE=0.2 THALF=0.5 NU=1", "MUE=0.2 THALF=0.5 NU=1 &X")) {
    job <- .u1_job(parms = p)
    expect_false(identical(job$calls$fit[[3L]][[1L]], as.name("stop")) &&
                   grepl("does not run", .u1_msg(job), fixed = TRUE), info = p)
  }
})

test_that("a PROC-line value the lexer rejects stops the document (U1, #403)", {
  for (proc in c(" MAXITER=1E5", " CONDITION=5.")) {
    job <- .u1_job(proc = proc, parms = "MUE=0.2 THALF=1 NU=1")
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"), info = proc)
    msg <- .u1_msg(job)
    expect_match(msg, "PROC HAZARD does not run this job", fixed = TRUE, info = proc)
    expect_match(msg, "hazard_l.l:34-38", fixed = TRUE, info = proc)
  }
  job <- .u1_job(proc = " MAXITER=200 CONDITION=14", parms = "MUE=0.2 THALF=1 NU=1")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
})

test_that("FIXMNU1 on an active early phase stops the document (U1, #358)", {
  job <- .u1_job(parms = "MUE=0.2 THALF=1 NU=2 M=0.5 FIXMNU1")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"))
  msg <- .u1_msg(job)
  expect_match(msg, "FIXMNU1", fixed = TRUE)
  expect_match(msg, "|M*NU| = 1", fixed = TRUE)
  expect_match(msg, "#358", fixed = TRUE)
  # With no early phase there is nothing for FIXMNU1 to constrain: SAS fits.
  job <- .u1_job(parms = "MUL=0.2 TAU=1 GAMMA=2 ETA=1 FIXMNU1")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
})

# --- U1, r-reviewer pass 1 on the branch: false refusals and missed ones ---
.u1_stops <- function(job) identical(job$calls$fit[[3L]][[1L]], as.name("stop"))

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
  expect_true(.u1_stops(job))
  expect_match(.u1_msg(job), "PROC HAZARD does not run this job", fixed = TRUE)
})

test_that("FIXMNU1 stops honestly, even where it might constrain nothing (U1 review 2)", {
  # The vacuous-FIXMNU1 exemption was dropped after its sign cases went wrong
  # (M*NU = -1 is vacuous too; M = NU = -1 is SETG1920). The stop claims
  # neither that PROC HAZARD runs the job nor that it refuses it.
  for (p in c("MUE=0.2 THALF=1 NU=1 M=1 FIXM FIXNU FIXMNU1",
              "MUE=0.2 THALF=0.5 M=-1 NU=1 FIXM FIXNU FIXMNU1",
              "MUE=0.2 THALF=0.5 M=-1 NU=-1 FIXM FIXNU FIXMNU1")) {
    job <- .u1_job(parms = p)
    expect_true(.u1_stops(job), info = p)
    msg <- .u1_msg(job)
    expect_match(msg, "cannot emit PROC HAZARD's model", fixed = TRUE, info = p)
    expect_no_match(msg, "runs this job", fixed = TRUE, info = p)
  }
})

test_that("SETG3's entry refusals stop on the constraint path too (U1 review)", {
  # setg3.c:269-284 run before any constraint or WEIBULL logic.
  for (p in c("MUL=0.2 TAU=0 FIXTAU GAMMA=2 ETA=1 FIXGE2",
              "MUL=0.2 TAU=1 GAMMA=0 FIXGAMMA ETA=1 FIXGAE2",
              "MUL=0.2 TAU=1 GAMMA=2 ALPHA=-1 FIXALPHA ETA=1 FIXGE2")) {
    job <- .u1_job(parms = p)
    expect_true(.u1_stops(job), info = p)
    expect_match(.u1_msg(job), "SETG39[0-3]0", info = p)
  }
})

test_that("FIXGE2/FIXGAE2 without WEIBULL stops, saying it cannot tell (U1 review 2)", {
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
    expect_true(.u1_stops(job), info = p)
    msg <- .u1_msg(job)
    expect_match(msg, "cannot tell whether PROC HAZARD refuses", fixed = TRUE, info = p)
    expect_no_match(msg, "runs this job", fixed = TRUE, info = p)
  }
  # The same shapes with WEIBULL are the mirrored path, and fit.
  expect_false(.u1_stops(.u1_job(parms = "MUL=0.2 TAU=1 GAMMA=3 ETA=1 FIXGE2 WEIBULL")))
  expect_false(.u1_stops(.u1_job(parms = "MUL=0.2 TAU=1 GAMMA=4 ETA=0.5 FIXGAMMA FIXETA FIXGE2 WEIBULL")))
})

test_that("DELTA != 0 and a FIXTAU with no TAU written stop (U1 review)", {
  job <- .u1_job(parms = "MUE=0.2 THALF=1 NU=1 DELTA=0.5")
  expect_true(.u1_stops(job))
  expect_match(.u1_msg(job), "DELTA", fixed = TRUE)
  job <- .u1_job(parms = "MUL=0.2 GAMMA=2 ETA=1 FIXTAU")
  expect_true(.u1_stops(job))
  expect_match(.u1_msg(job), "0.75*Tmax", fixed = TRUE)
  # Controls: DELTA = 0 is R's model; a written positive TAU is fixed at it.
  expect_false(.u1_stops(.u1_job(parms = "MUE=0.2 THALF=1 NU=1 DELTA=0")))
  # DELTA is read only by SETG1 (setg1.c:306), which runs only for an active
  # early phase (shape.c:19-21): a late-only job ignores it.
  expect_false(.u1_stops(.u1_job(parms = "MUL=0.2 TAU=1 GAMMA=2 ETA=1 DELTA=0.5")))
  expect_false(.u1_stops(.u1_job(parms = "MUL=0.2 TAU=2 GAMMA=2 ETA=1 FIXTAU")))
})

test_that("a template placeholder or a bare % is a syntax error, not filled in (U1 review 2)", {
  # `NU=?` used to fit with NU filled from SAS's default; PROC HAZARD's lexer
  # rejects `?` (hazard_l.l:178). `50%` is no macro: % with no name after it.
  for (p in c("MUE=0.2 THALF=1 NU=?", "MUE=0.2 THALF=1 NU=50%")) {
    job <- .u1_job(parms = p)
    expect_true(.u1_stops(job), info = p)
    expect_match(.u1_msg(job), "PROC HAZARD does not run this job", fixed = TRUE, info = p)
  }
  expect_match(.u1_msg(.u1_job(parms = "MUE=0.2 THALF=1 NU=?")), "fill it in", fixed = TRUE)
})
