test_that("a canonical AVC-style block becomes a hazard() call", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONSERVE P OUTHAZ=EX.HZD",
    "STEEPEST QUASI CONDITION=14 MI=200;",
    "EVENT DEAD; TIME INT_DEAD;",
    "PARMS MUE=0.2361727 THALF=0.1512095 NU=1.438652 M=1 FIXM MUC=0.0005436977; );"
  ))
  b <- .hzr_sas_blocks(txt)[[1L]]
  got <- .hzr_parse_hazard(b)

  expect_equal(got$outhaz, "EX.HZD")
  expect_equal(got$call[["data"]], as.name("AVCS"))
  expect_equal(got$call[["time"]], as.name("INT_DEAD"))
  expect_equal(got$call[["control"]],
               quote(list(maxit = 200, condition = 14, conserve = TRUE,
                          method = "bfgs")))
  # STEEPEST has no R equivalent and must be surfaced, not dropped.
  expect_true("STEEPEST" %in% got$untranslated$construct)
})

test_that("NOCONSERVE is emitted explicitly rather than defaulted", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A NOCONSERVE CONDITION=14;",
    "EVENT D; TIME T; PARMS MUE=1 THALF=1 NU=1; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(got$call[["control"]][["conserve"]], FALSE)
})

test_that("a LATE statement with comma-separated VAR=VALUE operands parses", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A;",
    "EVENT D; TIME T;",
    "PARMS MUL=0.01 TAU=2 GAMMA=1.5;",
    "LATE NOPREVTE=2.653528, NOTEE=-0.3680751; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(
    got$call[["phases"]],
    quote(list(hzr_phase("g3", tau = 2, gamma = 1.5,
                          formula = ~NOPREVTE + NOTEE)))
  )
  expect_true(any(got$untranslated$reason ==
    "phase covariate starting values are not yet mapped to theta"))
})
