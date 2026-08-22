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
  AVCS$C3FLAG <- as.numeric(seq_len(n) %% 5 == 0)
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
  AVCS$C3FLAG <- ifelse(seq_len(n) %% 5 == 0, 3, 0)
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
  status <- ifelse(AVCS$DEAD == 0 & AVCS$C3FLAG > 0, 2, AVCS$DEAD)
  expect_equal(got$env$fit$data$weights, ifelse(status == 2, AVCS$C3FLAG, 1))
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
  expect_null(got$env$fit$data$weights)
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
