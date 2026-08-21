test_that("a HAZPRED block becomes a predict() call with se.fit", {
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=EX.HZD OUT=PREDICT; TIME MONTHS; );"
  )
  b <- .hzr_sas_blocks(txt)[[1L]]
  got <- .hzr_parse_hazpred(b, txt)
  expect_equal(got$inhaz, "EX.HZD")
  expect_equal(got$call[["se.fit"]], TRUE)
  expect_equal(got$call[["newdata"]], as.name("PREDICT"))
})

test_that("a log-spaced DO grid becomes an exp(seq(...)) call", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; MAX=180; LN_MAX=LOG(MAX); INC=(5+LN_MAX)/99.9;",
    "DO LN_TIME=-5 TO LN_MAX BY INC, LN_MAX; MONTHS=EXP(LN_TIME); OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_equal(got$grid,
               quote(data.frame(MONTHS = exp(seq(-5, log(180), length.out = 100)))))
})

test_that("a grid built by SET is untranslated, not guessed at", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; SET COHORT; ",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_null(got$grid)
  expect_true(any(grepl("grid", got$untranslated$reason)))
})

test_that("the default HAZPRED block emits both call and call_haz with se.fit", {
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=EX.HZD OUT=PREDICT; TIME MONTHS; );"
  )
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_equal(got$call[["type"]], "survival")
  expect_equal(got$call[["se.fit"]], TRUE)
  expect_equal(got$call_haz[["type"]], "hazard")
  expect_equal(got$call_haz[["se.fit"]], TRUE)
})

test_that("NOHAZ makes call_haz NULL", {
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=EX.HZD OUT=PREDICT NOHAZ; TIME MONTHS; );"
  )
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_equal(got$call[["type"]], "survival")
  expect_null(got$call_haz)
})

test_that("NOSURV makes call the hazard prediction", {
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=EX.HZD OUT=PREDICT NOSURV; TIME MONTHS; );"
  )
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_equal(got$call[["type"]], "hazard")
  expect_null(got$call_haz)
})

test_that("NOSURV and NOHAZ together yield no predict() call at all", {
  # The ternary that picks `call` tests only want_surv, so without a guard
  # this degenerate input would silently produce a hazard predict() nobody
  # asked for. Both call and call_haz must be NULL, and the suppression must
  # be recorded, not dropped.
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=EX.HZD OUT=PREDICT NOSURV NOHAZ; TIME MONTHS; );"
  )
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_null(got$call)
  expect_null(got$call_haz)
  expect_true(any(grepl("NOSURV", got$untranslated$construct)))
})
