test_that("SELECTION FORWARD maps to direction = both, not forward", {
  # The obvious mapping is wrong and fails silently on every stepwise job.
  got <- .hzr_selection_spec(c("FORWARD", "SLENTRY=0.05", "SLSTAY=0.10"))
  expect_equal(got$direction, "both")
  expect_equal(got$slentry, 0.05)
  expect_equal(got$slstay, 0.10)
})

test_that("every forward/stepwise spelling reaches the same direction", {
  for (kw in c("FORWARD", "FW", "SW", "SELECT", "STEPWISE")) {
    expect_equal(.hzr_selection_spec(kw)$direction, "both",
                 info = paste("spelling:", kw))
  }
})

test_that("BACKWARD and ONEWAY are distinct from stepwise", {
  expect_equal(.hzr_selection_spec("BACKWARD")$direction, "backward")
  expect_equal(.hzr_selection_spec("BW")$direction, "backward")
  # ONEWAY means no stepwise at all: the caller emits a plain hazard() fit.
  expect_null(.hzr_selection_spec("NOSTEPWISE")$direction)
})

test_that("a SELECTION block emits hzr_stepwise, not hazard", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT D; TIME T; EARLY X1=1, X2=1, X3=1;",
    "PARMS MUE=1 THALF=1 NU=1;",
    "SELECTION FORWARD SLENTRY=0.05 SLSTAY=0.1; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(got$call[[1L]], as.name("hzr_stepwise"))
  expect_equal(got$call[["direction"]], "both")
})

test_that("a non-numeric SLENTRY value is untranslated, not coerced to NA silently", {
  got <- .hzr_selection_spec(c("STEPWISE", "SLENTRY=oops"))
  expect_null(got$slentry)
  expect_true(any(grepl("SLENTRY", got$untranslated$construct)))
})

test_that("an unknown SELECTION option is recorded, not dropped", {
  got <- .hzr_selection_spec("BOGUS")
  expect_true(any(grepl("BOGUS", got$untranslated$construct)))
})

test_that("a SELECTION statement with no direction still enables stepwise", {
  # stepwisestmt : STEPWISE stepwiseopts { setopt(33); } -- stepwiseopts may be
  # empty, so the statement itself turns stepwise on.
  got <- .hzr_selection_spec(c("SLENTRY=0.05", "SLSTAY=0.1"))
  expect_true(got$stepwise)
  expect_equal(got$direction, "both")
  expect_equal(got$slentry, 0.05)
  expect_equal(got$slstay, 0.1)
})

test_that("a bare SELECTION with no options at all enables stepwise", {
  got <- .hzr_selection_spec(character(0))
  expect_true(got$stepwise)
  expect_equal(got$direction, "both")
})

test_that("NOSTEPWISE disables it and records any orphaned thresholds", {
  got <- .hzr_selection_spec(c("NOSTEPWISE", "SLENTRY=0.05"))
  expect_false(got$stepwise)
  expect_true(any(grepl("SLENTRY", got$untranslated$construct)))
})

test_that("a directionless SELECTION block emits hzr_stepwise, not hazard", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT D; TIME T; EARLY X1=0.5, X2=0.25;",
    "PARMS MUE=1 THALF=1 NU=1;",
    "SELECTION SLENTRY=0.05 SLSTAY=0.1; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(got$call[[1L]], as.name("hzr_stepwise"))
  expect_equal(got$call[["slentry"]], 0.05)
  expect_equal(got$call[["slstay"]], 0.1)
})
