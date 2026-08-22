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
  expect_true(any(got$env$.hzr_status == 2))
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
