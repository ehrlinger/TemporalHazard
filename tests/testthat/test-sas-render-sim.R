test_that("render_sim reports per-chunk success and failure", {
  job <- list(calls = list(
    good = quote(x <- 1 + 1),
    bad  = quote(stop("boom")),
    uses = quote(y <- x * 2)
  ))
  got <- render_sim(job)
  expect_false(got$ok)
  expect_equal(unname(got$results[["good"]]), "ok")
  expect_match(got$results[["bad"]], "^ERROR: boom")
  expect_equal(got$env$y, 4)
})

test_that("render_sim binds supplied data and nothing else", {
  job <- list(calls = list(a = quote(n <- nrow(D)), b = quote(z <- MISSING_VAR)))
  got <- render_sim(job, data = list(D = data.frame(x = 1:3)))
  expect_equal(got$env$n, 3L)
  expect_match(got$results[["b"]], "object 'MISSING_VAR' not found")
})

test_that("render_sim does not report ok when nothing was evaluated", {
  # all(character(0) == "ok") is TRUE -- a job with zero emitted calls must
  # not be indistinguishable from one where every call succeeded.
  got <- render_sim(list(calls = list()))
  expect_false(got$ok)
})

test_that("a minimal translated job renders end to end", {
  skip_on_cran()
  set.seed(1)
  AVCS <- data.frame(INT_DEAD = stats::rexp(200, 0.2),
                     DEAD = rep(c(1, 0), length.out = 200))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  expect_s3_class(got$env$fit, "hazard")
  expect_true(isTRUE(got$env$fit$fit$converged))
})

test_that("a censoring job's status chunk renders", {
  skip_on_cran()
  set.seed(2)
  n <- 120
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     DEAD = rep(c(1, 0, 0), length.out = n))
  AVCS$C3FLAG <- as.numeric(seq_len(n) %% 5 == 0 & AVCS$DEAD == 0)
  AVCS$ICTIME <- ifelse(AVCS$C3FLAG > 0, AVCS$INT_DEAD * 0.5, NA)
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; ICENSOR C3FLAG = ICTIME;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  # The status is derived INTO the dataset, not bound locally, so that
  # hazard()'s data mask cannot shadow it with a column of the same name.
  expect_true(any(got$env$AVCS$.hzr_status == 2))
})

test_that("a HAZARD + HAZPRED job renders, predictions included", {
  skip_on_cran()
  set.seed(3)
  AVCS <- data.frame(INT_DEAD = stats::rexp(200, 0.2),
                     DEAD = rep(c(1, 0), length.out = 200))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(c(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14 OUTHAZ=E.H;",
    "EVENT DEAD; TIME INT_DEAD; PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );",
    "DATA PGRID; DO MONTHS = 1 TO 12 BY 1; OUTPUT; END; RUN;",
    "%HAZPRED( PROC HAZPRED DATA=PGRID INHAZ=E.H OUT=P; TIME MONTHS; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  expect_true("time" %in% names(got$env$PGRID))
})

test_that("the ICENSOR count reaches the fitted object as weights", {
  skip_on_cran()
  set.seed(5)
  n <- 120
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     DEAD = rep(c(1, 0, 0), length.out = n))
  AVCS$C3FLAG <- ifelse(seq_len(n) %% 5 == 0 & AVCS$DEAD == 0, 3, 0)
  AVCS$ICTIME <- ifelse(AVCS$C3FLAG > 0, AVCS$INT_DEAD * 0.5, NA)
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; ICENSOR C3FLAG = ICTIME;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  # The decisive check: a count of 3 must not arrive as a weight of 1 -- and
  # a count of 0 on a non-interval row must not arrive as a weight of 0,
  # which would delete that row from the likelihood.
  status <- ifelse(AVCS$DEAD > 0, 1, ifelse(AVCS$C3FLAG > 0, 2, 0))
  expect_equal(got$env$fit$data$status, status)
  expect_equal(got$env$fit$data$weights,
               ifelse(AVCS$DEAD > 0, AVCS$DEAD,
                      ifelse(AVCS$C3FLAG > 0, AVCS$C3FLAG, 1)))
  expect_true(any(got$env$fit$data$weights == 3))
  expect_false(any(got$env$fit$data$weights == 0))
})

test_that("a job with no ICENSOR fits with unit weights, unchanged", {
  skip_on_cran()
  set.seed(6)
  AVCS <- data.frame(INT_DEAD = stats::rexp(200, 0.2),
                     DEAD = rep(c(1, 0), length.out = 200))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  # EVENT is a count, so a weights argument is always emitted now (#157).
  # For the 0/1 DEAD this job carries it must evaluate to unit weights on
  # every row -- the same fit the job got before, asserted on the fitted
  # object rather than on the presence or absence of an argument.
  expect_equal(got$env$fit$data$weights, rep(1, 200))
  expect_equal(got$env$fit$data$status, AVCS$DEAD)
})

test_that("a repeat-event count reaches the fitted object as weight and status 1", {
  # #157: EVENT is a COUNT (readc1.c keeps any C1 >= 0; setlik.c weights the
  # row's whole contribution by c1w = C1 * WT). Mapped straight onto status,
  # DEAD = 2 becomes INTERVAL-CENSORED in this package's coding -- a
  # different likelihood branch -- and DEAD = 3 falls outside the coding
  # entirely, which also disables Conservation of Events.
  skip_on_cran()
  set.seed(11)
  n <- 120
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     DEAD = rep(c(0, 1, 2, 3), length.out = n))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  # Every count > 0 is an event, never an interval or an out-of-range code.
  expect_equal(got$env$fit$data$status, as.numeric(AVCS$DEAD > 0))
  expect_false(any(got$env$fit$data$status == 2))
  # The magnitude survives as the weight, and no row is deleted by a zero.
  expect_equal(got$env$fit$data$weights, ifelse(AVCS$DEAD > 0, AVCS$DEAD, 1))
  expect_equal(sort(unique(got$env$fit$data$weights)), c(1, 2, 3))
})

test_that("WEIGHT leaves right-censored rows at 1, as c2 enters setlik.c unweighted", {
  # #158: setlik.c accumulates c1c2c3 = c1w + c2 + c3w, and c2 -- which
  # readc2.c sets to ONE on exactly the rows where neither other count fires
  # -- is the term NOT multiplied by the weight. A bare WT on censored rows
  # over-weights them, and a WT of 0 there deletes them from the likelihood
  # silently.
  skip_on_cran()
  set.seed(12)
  n <- 120
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     DEAD = rep(c(1, 0, 0), length.out = n))
  # MORBID is 0 on censored rows -- the shape of hz.tm123.OMC.sas's WEIGHT
  # MORBID, and the case that silently deletes rows.
  AVCS$MORBID <- ifelse(AVCS$DEAD > 0, rep(c(2, 4), length.out = n), 0)
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; WEIGHT MORBID;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  w <- got$env$fit$data$weights
  expect_equal(w, ifelse(AVCS$DEAD > 0, AVCS$MORBID, 1))
  # The decisive one: not a single censored row was deleted by a zero
  # weight, and there are as many weight-1 rows as there are censored rows.
  expect_equal(sum(w == 0), 0L)
  expect_equal(sum(w == 1), sum(AVCS$DEAD == 0))
  expect_equal(got$env$fit$data$status, as.numeric(AVCS$DEAD > 0))
})

test_that("a row that is both an event and interval-censored is refused, not guessed", {
  # setlik.c SUMS the contributions (c1c2c3 = c1w + c2 + c3w), so such a row
  # is an event of weight C1*WT AND an interval observation of weight C3*WT
  # at once. hazard() has one status and one weight per row. Picking the
  # event branch converges and reports plausibly, which is the shape this
  # package exists to refuse.
  skip_on_cran()
  set.seed(13)
  n <- 60
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     DEAD = rep(c(1, 0, 0), length.out = n))
  AVCS$C3FLAG <- ifelse(seq_len(n) %% 5 == 0, 2, 0)
  AVCS$ICTIME <- ifelse(AVCS$C3FLAG > 0, AVCS$INT_DEAD * 0.5, NA)
  expect_true(any(AVCS$DEAD > 0 & AVCS$C3FLAG > 0))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; ICENSOR C3FLAG = ICTIME;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_false(got$ok)
  expect_match(got$results[["status"]], "both non-zero")
  # And no fit was produced: a refusal that still leaves a fitted object
  # behind is not a refusal.
  expect_null(got$env$fit)
})

test_that("the RCENSOR count reaches the fitted object as weights", {
  # #162: RCENSOR names C2 (hazard_y.y; rcnsprc.c sets c2name from statement
  # field 13), and C2 is a COUNT of censored individuals at T (setlik.c
  # header). With c2name non-blank readc2.c reads it straight from the data
  # -- the `C2 = ONE` derivation is skipped entirely -- so four censored
  # individuals must arrive as weight 4, not as one observation.
  skip_on_cran()
  set.seed(21)
  n <- 120
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     DEAD = rep(c(1, 0, 0), length.out = n))
  # A genuine count above 1, and only on rows where DEAD does not fire --
  # the both-fire case is refused, and is tested separately below.
  AVCS$NCENS <- ifelse(AVCS$DEAD > 0, 0, rep(c(1, 4), length.out = n))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; RCENSOR NCENS;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  # The decisive check, on the fitted object: a count of 4 must not arrive
  # as a weight of 1, and the event rows keep their own count.
  expect_equal(got$env$fit$data$weights,
               ifelse(AVCS$DEAD > 0, AVCS$DEAD, AVCS$NCENS))
  expect_true(any(got$env$fit$data$weights == 4))
  # C2 changes the weight, never the status: it is right-censoring, code 0.
  expect_equal(got$env$fit$data$status, as.numeric(AVCS$DEAD > 0))
  expect_false(any(got$env$fit$data$status == 2))
  # And the counts really are being read -- a fit whose weights are all 1
  # would satisfy a weaker assertion while ignoring RCENSOR entirely.
  expect_gt(sum(got$env$fit$data$weights), n)
})

test_that("a WEIGHT variable never multiplies the RCENSOR count", {
  # c2 is the one term setlik.c leaves out of the WT product
  # (c1c2c3 = c1w + c2 + c3w). Naming C2 explicitly does not change that.
  skip_on_cran()
  set.seed(22)
  n <- 90
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     DEAD = rep(c(1, 0, 0), length.out = n))
  AVCS$NCENS <- ifelse(AVCS$DEAD > 0, 0, 3)
  # MORBID is 0 on censored rows, the hz.tm123.OMC.sas shape: multiplied in,
  # it would delete every censored row from the likelihood.
  AVCS$MORBID <- ifelse(AVCS$DEAD > 0, 2, 0)
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; RCENSOR NCENS; WEIGHT MORBID;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  expect_equal(got$env$fit$data$weights,
               ifelse(AVCS$DEAD > 0, AVCS$DEAD * AVCS$MORBID, AVCS$NCENS))
  expect_equal(sum(got$env$fit$data$weights == 0), 0L)
})

test_that("a row that is both an event and right-censored is refused, not guessed", {
  # readobs.c deletes a row only when BOTH ic10 and ic20 fire, so a row with
  # C1 > 0 AND C2 > 0 reaches setlik.c and contributes c1w + c2 -- an event
  # of weight C1*WT and a right-censored observation of weight C2 at once.
  # Fitting it as an event alone converges and reports plausibly (#162).
  skip_on_cran()
  set.seed(23)
  n <- 60
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     DEAD = rep(c(1, 0, 0), length.out = n))
  AVCS$NCENS <- ifelse(seq_len(n) %% 5 == 0, 2, 0)
  expect_true(any(AVCS$DEAD > 0 & AVCS$NCENS > 0))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; RCENSOR NCENS;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_false(got$ok)
  expect_match(got$results[["status"]], "both non-zero")
  # A refusal that still leaves a fitted object behind is not a refusal.
  expect_null(got$env$fit)
})

test_that("a 0/1 RCENSOR flag fits identically with and without the statement", {
  # The two corpus jobs that carry RCENSOR (hz.te123.OMC.sas,
  # hz.tm123.OMC.sas) set a mutually-exclusive 0/1 CENSORED. Their fits must
  # be bit-for-bit what they were: this fix is for counts above 1. Compared
  # as two real fits, not as two call shapes.
  skip_on_cran()
  set.seed(24)
  n <- 150
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     TE = rep(c(1, 0, 0), length.out = n))
  AVCS$CENSORED <- 1 - AVCS$TE
  head <- "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14; EVENT TE; TIME INT_DEAD;"
  tail <- "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  fits <- lapply(c("", " RCENSOR CENSORED;"), function(stmt) {
    f <- withr::local_tempfile(fileext = ".sas")
    writeLines(paste0(head, stmt, " ", tail), f)
    job <- suppressWarnings(hzr_translate_sas(f))
    # .hzr_optim_multiphase()'s random restart draws from the global RNG
    # stream (R/likelihood-multiphase.R), so the SAME call twice returns
    # different theta. Seed each fit identically or this comparison tests
    # the RNG, not the translation.
    set.seed(99)
    got <- render_sim(job, data = list(AVCS = AVCS))
    expect_true(got$ok, info = paste(names(got$results), got$results,
                                     collapse = " | "))
    got$env$fit
  })
  expect_equal(fits[[2L]]$data$weights, fits[[1L]]$data$weights)
  expect_equal(fits[[2L]]$data$status, fits[[1L]]$data$status)
  expect_equal(fits[[2L]]$fit$theta, fits[[1L]]$fit$theta)
  expect_equal(fits[[2L]]$fit$objective, fits[[1L]]$fit$objective)
  # Guard against the comparison being vacuous. The fitted parameters live
  # under fit$fit, not fit -- `fits[[1L]]$theta` is NULL, and every
  # all()/is.finite() assertion over NULL passes while comparing nothing.
  expect_equal(length(fits[[1L]]$data$weights), n)
  expect_true(is.finite(fits[[1L]]$fit$objective))
  expect_true(all(is.finite(fits[[1L]]$fit$theta)))
  expect_gt(length(fits[[1L]]$fit$theta), 0L)
  expect_equal(fits[[1L]]$data$weights, rep(1, n))
})

test_that("a data column named .hzr_status cannot change the fitted status", {
  skip_on_cran()
  set.seed(7)
  n <- 120
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     DEAD = rep(c(1, 0, 0), length.out = n))
  AVCS$C3FLAG <- as.numeric(seq_len(n) %% 5 == 0)
  AVCS$ICTIME <- ifelse(AVCS$C3FLAG > 0, AVCS$INT_DEAD * 0.5, NA)
  # The hostile column: every row an exact event, which is a different
  # censoring structure from the one the job specifies. `.hzr_status` is
  # protected against SAS-name collisions, but nothing stops a column of the
  # user's R data frame carrying it -- and hazard()'s vector path gives
  # columns precedence over the calling frame, so a locally bound temporary
  # would lose to this silently, for `status =` and the `time_lower` gate
  # alike.
  AVCS$.hzr_status <- rep(1, n)

  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; ICENSOR C3FLAG = ICTIME;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))

  want <- ifelse(AVCS$DEAD == 0 & AVCS$C3FLAG > 0, 2, AVCS$DEAD)
  # The assertion that can fail: the hostile column claims 1 everywhere, so
  # `want` must not be all 1s or this proves nothing.
  expect_true(any(want == 2))
  expect_false(all(want == 1))
  expect_equal(got$env$fit$data$status, want)
  # The gate reads the same classification, so the interval rows enter at
  # ICTIME rather than 0.
  expect_equal(got$env$fit$data$time_lower,
               ifelse(want == 2, AVCS$ICTIME, 0))
})

test_that("a refused prediction grid stops rather than predicting", {
  skip_on_cran()
  set.seed(8)
  AVCS <- data.frame(INT_DEAD = stats::rexp(200, 0.2),
                     DEAD = rep(c(1, 0), length.out = 200))
  f <- withr::local_tempfile(fileext = ".sas")
  # A SET-derived grid: .hzr_parse_grid() refuses it, so no chunk builds
  # PGRID. Emitting predict(fit, newdata = PGRID) would fail on an unbound
  # name -- or, worse, predict over whatever else is named PGRID in the
  # rendering session and report it.
  writeLines(c(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14 OUTHAZ=E.H;",
    "EVENT DEAD; TIME INT_DEAD; PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );",
    "DATA PGRID; SET AVCS; KEEP INT_DEAD; RUN;",
    "%HAZPRED( PROC HAZPRED DATA=PGRID INHAZ=E.H OUT=P; TIME INT_DEAD; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))

  shape <- sas_job_shape(job)
  # Both prediction families are refused: no predict() chunk survives, and
  # each one became a stop().
  expect_equal(length(shape$preds), 0L)
  expect_equal(length(shape$refusals), 2L)
  # The refusal is recorded, not merely emitted.
  expect_true(any(grepl("PGRID", job$untranslated$construct, fixed = TRUE)))

  # A hostile object of the grid's name in the rendering session must not
  # turn the refusal into a prediction.
  got <- suppressWarnings(render_sim(
    job, data = c(sas_synth_data(job),
                  list(PGRID = data.frame(time = c(1, 6, 12))))
  ))
  for (nm in shape$refusals) {
    expect_match(got$results[[nm]], "^ERROR: ")
    expect_match(got$results[[nm]], "PGRID", fixed = TRUE)
  }
  expect_false(got$ok)
})

test_that("render_sim cannot see the caller's globals", {
  # The helper is a test oracle, not a Quarto simulator: a global left over
  # from another test that happens to supply a missing symbol turns a broken
  # translation into "ok". Parenting the render environment on globalenv()
  # made exactly that possible.
  assign("HZR_RENDER_PROBE", 99, envir = globalenv())
  withr::defer(rm("HZR_RENDER_PROBE", envir = globalenv()))
  # The probe really is visible from the test's own frame, so a failure here
  # is about the render environment and not about the assign() above.
  expect_equal(HZR_RENDER_PROBE, 99)

  got <- render_sim(list(calls = list(a = quote(z <- HZR_RENDER_PROBE))))
  expect_match(got$results[["a"]], "object 'HZR_RENDER_PROBE' not found")
  # The package's own exports still resolve, or nothing could render.
  ok <- render_sim(list(calls = list(a = quote(f <- hzr_decompos(1, 1, 1, 1)))))
  expect_true(ok$ok)
})
