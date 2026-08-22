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
               quote(data.frame(
                 time = exp(-5 + seq(0, 99) * ((5 + log(180)) / 99.9))
               )))
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

test_that("an explicit DO grid resolves DATA-step constants, incl. 1*DTY", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; DIGITAL=0;",
    "DTY=12/365.2425;",
    "DO MONTHS=1*DTY,2*DTY,24 TO 180 BY 12;",
    "OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_false(is.null(got$grid))

  dty <- 12 / 365.2425
  expected <- c(1 * dty, 2 * dty, seq(24, 180, by = 12))
  expect_equal(eval(got$grid)$time, expected)
})

test_that("a DO list referencing an unknown name still refuses and records", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; DO MONTHS=1*FOO,2*FOO; OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_null(got$grid)
  expect_true(any(grepl("grid", got$untranslated$reason)))
})

test_that("a constant defined in terms of an earlier constant resolves", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; A=10; B=A*2; DO MONTHS=1*B,2*B; OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_false(is.null(got$grid))
  expect_equal(eval(got$grid)$time, c(20, 40))
})

test_that("a DO list constant defined via a function call refuses, not evaluates", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; A=SQRT(4); DO MONTHS=1*A,2*A; OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_null(got$grid)
  expect_true(any(grepl("grid", got$untranslated$reason)))
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

test_that("the survival predict call sets conf.type = logit for SAS parity", {
  # SAS HAZPRED's survival confidence limits are logit-scale
  # (hzp_calc_srv_CL.c); predict.hazard() defaults to "log-log". Emitting the
  # default reproduces the job with silently different bounds.
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=E.H OUT=P; TIME MONTHS; );"
  )
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_equal(got$call[["conf.type"]], "logit")
  # Hazard CLs are log scale in both engines, so that call must NOT be steered.
  expect_null(got$call_haz[["conf.type"]])
})

test_that("conf.type is omitted when NOCL suppresses confidence limits", {
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=E.H OUT=P NOCL; TIME MONTHS; );"
  )
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_equal(got$call[["se.fit"]], FALSE)
  expect_null(got$call[["conf.type"]])
})

test_that("the log grid matches SAS's INC = (5 + LN_MAX)/99.9 step", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PGRID; MAX = 180; LN_MAX = LOG(MAX);",
    "INC = (5 + LN_MAX)/99.9;",
    "DO LN_TIME = -5 TO LN_MAX BY INC; MONTHS = EXP(LN_TIME);",
    "OUTPUT; END; RUN;"
  ))
  cl <- .hzr_parse_grid(txt, "PGRID")
  grid <- eval(cl)
  expect_equal(nrow(grid), 100L)
  # SAS's last point is exp(-5 + 99*(log(180)+5)/99.9) = 164.2, not 180.
  expect_equal(max(grid$time), exp(-5 + 99 * (log(180) + 5) / 99.9),
               tolerance = 1e-8)
  expect_false(isTRUE(all.equal(max(grid$time), 180)))
})

test_that("a directly-assigned log bound is not mistaken for MAX", {
  # No separate MAX variable is ever assigned here -- LN_MAX is set directly,
  # as #153 describes. The unrelated LOG(2) is only there so the body still
  # trips the log-grid branch's own " DO "/"LOG("/EXP() detection heuristic.
  txt <- .hzr_sas_normalise(paste(
    "DATA PGRID; LN_MAX = 5.2; X = LOG(2); INC = (5 + LN_MAX)/99.9;",
    "DO LN_TIME = -5 TO LN_MAX BY INC; MONTHS = EXP(LN_TIME);",
    "OUTPUT; END; RUN;"
  ))
  cl <- .hzr_parse_grid(txt, "PGRID")
  grid <- if (is.null(cl)) NULL else eval(cl)
  # Must not silently produce a grid ending at t = 5.2 by reading LN_MAX=5.2
  # as MAX=5.2 and taking its log (#153). There is no MAX to log here, so
  # refusing to translate (NULL) is the correct outcome.
  expect_true(is.null(grid) || max(grid$time) > 100)
})
