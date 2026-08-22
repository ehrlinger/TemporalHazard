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

test_that("a SELECTION block is refused, not emitted as hzr_stepwise", {
  # This test used to assert got$call[[1L]] == as.name("hzr_stepwise") and
  # got$call[["direction"]] == "both". It does not any more: hzr_stepwise()'s
  # refit path (.hzr_refit_with_scope(), R/stepwise-refit.R) requires a
  # formula-interface base fit, and this translator only ever emits the
  # vector interface, so every candidate refit would error, the forward step
  # downgrades those errors to warnings, and the run would silently report
  # zero steps -- indistinguishable from "nothing met slentry" (#159, #160).
  # Refuse loudly instead, mirroring .hzr_censor_spec()'s LCENSOR + ICENSOR
  # refusal.
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT D; TIME T; EARLY X1=1, X2=1, X3=1;",
    "PARMS MUE=1 THALF=1 NU=1;",
    "SELECTION FORWARD SLENTRY=0.05 SLSTAY=0.1; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_true(any(grepl("SELECTION", got$untranslated$construct, fixed = TRUE)))
  # The callout alone is not enough -- a reader who renders past it would
  # still get a fit. The emitted call must itself refuse.
  expect_identical(got$call[[1L]], as.name("stop"))
  expect_error(eval(got$call), "SELECTION")
  expect_null(got$call[["direction"]])
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

test_that("a directionless SELECTION block is refused too", {
  # This test used to assert got$call[[1L]] == as.name("hzr_stepwise") and
  # got$call[["slentry"]]/["slstay"] == 0.05/0.1 -- i.e. that a bare
  # SELECTION (stepwisestmt with empty stepwiseopts, still enables stepwise)
  # translated cleanly. It does not any more, for the same reason as the
  # directed case above: refused, not translated (#159, #160).
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT D; TIME T; EARLY X1=0.5, X2=0.25;",
    "PARMS MUE=1 THALF=1 NU=1;",
    "SELECTION SLENTRY=0.05 SLSTAY=0.1; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_true(any(grepl("SELECTION", got$untranslated$construct, fixed = TRUE)))
  expect_identical(got$call[[1L]], as.name("stop"))
  expect_null(got$call[["slentry"]])
  expect_null(got$call[["slstay"]])
})

test_that("a SELECTION job is refused, not translated into a no-op screen", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; EARLY X1, X2, X3;",
    "SELECTION STEPWISE SLENTRY=0.05 SLSTAY=0.1;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_true(any(grepl("SELECTION", job$untranslated$construct, fixed = TRUE)))
  # No hzr_stepwise() call may be emitted at all: one that runs and screens
  # nothing is the defect this refusal exists to avoid (#159, #160).
  deparsed <- paste(vapply(job$calls, function(c0) {
    paste(deparse(c0), collapse = " ")
  }, character(1)), collapse = " ")
  expect_false(grepl("hzr_stepwise", deparsed, fixed = TRUE))
  expect_true(grepl("stop(", deparsed, fixed = TRUE))
})

test_that("a job with no SELECTION statement is unaffected", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; EARLY X1, X2, X3;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_false(any(grepl("SELECTION", job$untranslated$construct, fixed = TRUE)))
  expect_equal(job$calls$fit[[3L]][[1L]], as.name("hazard"))
})
