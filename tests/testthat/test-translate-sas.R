test_that("a job with no recognisable block errors and writes nothing", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(c("DATA X; SET Y;", "PROC PRINT;"), f)
  out <- withr::local_tempdir()
  expect_error(hzr_translate_sas(f, out_dir = out), "no SAS statements|no HAZARD")
  expect_length(list.files(out, pattern = "[.]qmd$"), 0L)
})

test_that("librefs resolves an INHAZ that no translated job wrote", {
  f <- withr::local_tempfile(fileext = ".sas")
  # The grid DATA step is load bearing: a DATA= grid the translator cannot
  # build is refused, and the refusal is itself a stop() chunk, so without it
  # the "no stop()" assertion below would be about the wrong stop().
  writeLines(paste(
    "DATA P; DO MONTHS=1 TO 12 BY 1; OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.HZD OUT=P; TIME MONTHS; );"
  ), f)
  out <- withr::local_tempdir()
  job <- hzr_translate_sas(f, out_dir = out, librefs = c(EX = "estimates"))
  qmd <- readLines(file.path(out, sub("[.]sas$", ".qmd", basename(f))))
  expect_true(any(grepl("hzr_read_outhaz", qmd, fixed = TRUE)))
  expect_false(any(grepl("stop(", qmd, fixed = TRUE)))
})

test_that("without librefs an unresolved INHAZ warns and emits stop()", {
  f <- withr::local_tempfile(fileext = ".sas")
  # Translatable grid, so the stop() this asserts on can only be the
  # unresolved INHAZ -- a refused grid emits a stop() of its own.
  writeLines(paste(
    "DATA P; DO MONTHS=1 TO 12 BY 1; OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.HZD OUT=P; TIME MONTHS; );"
  ), f)
  out <- withr::local_tempdir()
  expect_warning(hzr_translate_sas(f, out_dir = out), "unresolved")
  qmd <- readLines(file.path(out, sub("[.]sas$", ".qmd", basename(f))))
  expect_true(any(grepl("stop(", qmd, fixed = TRUE)))
  expect_true(any(grepl("EX.HZD", qmd, fixed = TRUE)))
})

test_that("a malformed librefs is rejected with a clear error", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines("%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.HZD OUT=P; TIME MONTHS; );", f)
  expect_error(hzr_translate_sas(f, librefs = "estimates"), "named character vector")
  expect_error(hzr_translate_sas(f, librefs = c(EX = 1)), "named character vector")
  # An NA name used to pass: nzchar() defaults to keepNA = FALSE, so
  # nzchar(NA_character_) is TRUE and any(!nzchar(...)) never sees it. The
  # validation now rejects it with the same message the empty-name case gets.
  na_named <- c(EX = "estimates", "other")
  names(na_named)[2] <- NA
  expect_error(hzr_translate_sas(f, librefs = na_named),
               "named character vector")
  expect_error(hzr_translate_sas(f, librefs = c(EX = "estimates", "")),
               "named character vector")
})

test_that("a block missing both EVENT and ICENSOR fails cleanly, naming the file", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(
    "%HAZARD( PROC HAZARD DATA=X; TIME T; PARMS MUE=1 THALF=1; );",
    f
  )
  err <- tryCatch(hzr_translate_sas(f), error = function(e) e)
  expect_s3_class(err, "error")
  expect_true(grepl(basename(f), conditionMessage(err), fixed = TRUE))
  expect_true(grepl("EVENT", conditionMessage(err), fixed = TRUE))
})

test_that("a hazard() call with data= gets a loud missing-data guard chunk, first", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(
    "%HAZARD( PROC HAZARD DATA=AVCS; EVENT DEAD; TIME T; PARMS MUE=1 THALF=1; );",
    f
  )
  out <- withr::local_tempdir()
  job <- hzr_translate_sas(f, out_dir = out)

  expect_identical(names(job$calls)[[1L]], "data")
  expect_true(any(grepl("exists(\"AVCS\")",
                         deparse(job$calls$data, width.cutoff = 500L),
                         fixed = TRUE)))

  qmd <- readLines(file.path(out, sub("[.]sas$", ".qmd", basename(f))))
  expect_true(any(grepl("AVCS", qmd, fixed = TRUE)))
  data_line <- grep("label: data", qmd, fixed = TRUE)
  fit_line <- grep("label: fit", qmd, fixed = TRUE)
  expect_length(data_line, 1L)
  expect_length(fit_line, 1L)
  expect_true(data_line < fit_line)
})

test_that("hzr_translate_sas parses the packaged example without writing when out_dir is NULL", {
  f <- system.file("extdata", "hz-example.sas", package = "TemporalHazard")
  job <- hzr_translate_sas(f)
  expect_s3_class(job, "hzr_sas_job")
  expect_true(job$coverage$tokens_seen > 0L)
})

test_that("a librefs-resolved INHAZ emits a path with a real extension", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.HZD OUT=P; TIME MONTHS; );", f)
  out <- withr::local_tempdir()
  suppressWarnings(hzr_translate_sas(f, out_dir = out, librefs = c(EX = "est")))
  qmd <- readLines(list.files(out, pattern = "[.]qmd$", full.names = TRUE))
  line <- grep("hzr_read_outhaz", qmd, value = TRUE)
  expect_length(line, 1L)
  # The emitted path must be one hzr_read_outhaz() can actually dispatch on.
  expect_match(line, "\\.(sas7bdat|rds)")
})

test_that("a librefs value that already carries an .rds extension is used verbatim", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.HZD OUT=P; TIME MONTHS; );", f)
  out <- withr::local_tempdir()
  suppressWarnings(hzr_translate_sas(
    f, out_dir = out, librefs = c(EX = "estimates/hzdeath.rds")))
  qmd <- readLines(list.files(out, pattern = "[.]qmd$", full.names = TRUE))
  line <- grep("hzr_read_outhaz", qmd, value = TRUE)
  expect_length(line, 1L)
  expect_true(any(grepl("estimates/hzdeath.rds", line, fixed = TRUE)))
  # No file.path() re-derivation from the member name when the librefs
  # value already names the file directly -- no separate ".sas7bdat" form.
  expect_false(any(grepl("sas7bdat", line, fixed = TRUE)))
})

test_that("a two-PROC-HAZARD job yields two distinct fit calls, not one dropped", {
  # This is the package's signature defect: a second block of the same kind
  # silently overwrote the first. Both fits must survive, under distinct,
  # stable names, and both must actually appear in the emitted .qmd.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=X OUTHAZ=A; EVENT DEAD; TIME T;",
    "PARMS MUE=1 THALF=1; );",
    "%HAZARD( PROC HAZARD DATA=X OUTHAZ=B; EVENT DEAD; TIME T;",
    "PARMS MUE=2 THALF=2; );"
  ), f)
  out <- withr::local_tempdir()
  job <- hzr_translate_sas(f, out_dir = out)

  expect_true(all(c("fit", "fit_2") %in% names(job$calls)))
  expect_equal(job$outhaz, c("A", "B"))
  # The two fits must not be identical calls -- otherwise the second block
  # would still be silently dropped in substance, even under a new name.
  expect_false(identical(job$calls$fit, job$calls$fit_2))

  qmd <- readLines(file.path(out, sub("[.]sas$", ".qmd", basename(f))))
  expect_length(grep("^#\\| label: fit$", qmd), 1L)
  expect_length(grep("^#\\| label: fit_2$", qmd), 1L)
  expect_true(any(grepl("t_half = 1", qmd, fixed = TRUE)))
  expect_true(any(grepl("t_half = 2", qmd, fixed = TRUE)))
})

test_that("a three-HAZPRED job yields three predict calls", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.A OUT=P NOHAZ; TIME MONTHS; );",
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.B OUT=P NOHAZ; TIME MONTHS; );",
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.C OUT=P NOHAZ; TIME MONTHS; );"
  ), f)
  out <- withr::local_tempdir()
  job <- suppressWarnings(hzr_translate_sas(f, out_dir = out))

  expect_equal(names(job$calls), c("pred", "pred_2", "pred_3"))

  qmd <- readLines(file.path(out, sub("[.]sas$", ".qmd", basename(f))))
  expect_length(grep("^#\\| label: pred$", qmd), 1L)
  expect_length(grep("^#\\| label: pred_2$", qmd), 1L)
  expect_length(grep("^#\\| label: pred_3$", qmd), 1L)
})

test_that("several fits with the same OUTHAZ resolve predict() positionally", {
  # Two PROC HAZARD blocks that both write OUTHAZ=OUTEST -- SAS semantics:
  # the second write overwrites the dataset the first one wrote. A HAZPRED
  # block that appears between them must reference the first fit; one after
  # both must reference the second. This is not "ambiguous": document order
  # of the blocks fully determines it.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "DATA P; DO MONTHS=1 TO 12 BY 1; OUTPUT; END;",
    "%HAZARD( PROC HAZARD DATA=X OUTHAZ=OUTEST; EVENT DEAD; TIME T;",
    "PARMS MUE=1 THALF=1; );",
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=OUTEST OUT=P NOHAZ; TIME MONTHS; );",
    "%HAZARD( PROC HAZARD DATA=X OUTHAZ=OUTEST; EVENT DEAD; TIME T;",
    "PARMS MUE=2 THALF=2; );",
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=OUTEST OUT=P NOHAZ; TIME MONTHS; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))

  expect_equal(job$calls$pred[[2L]], as.name("fit"))
  expect_equal(job$calls$pred_2[[2L]], as.name("fit_2"))
  # Neither is a genuine ambiguity: the same-job resolution above is
  # deterministic, so no "ambiguous" row is expected here.
  expect_false(any(grepl("ambiguous", job$untranslated$reason)))
})

test_that("an INHAZ that matches none of several local OUTHAZ is recorded ambiguous", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "DATA P; DO MONTHS=1 TO 12 BY 1; OUTPUT; END;",
    "%HAZARD( PROC HAZARD DATA=X OUTHAZ=A; EVENT DEAD; TIME T;",
    "PARMS MUE=1 THALF=1; );",
    "%HAZARD( PROC HAZARD DATA=X OUTHAZ=B; EVENT DEAD; TIME T;",
    "PARMS MUE=2 THALF=2; );",
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.C OUT=P NOHAZ; TIME MONTHS; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))

  # Falls back to `fit` rather than guessing which local fit was meant.
  expect_equal(job$calls$pred[[2L]], as.name("fit"))
  expect_true(any(grepl("ambiguous", job$untranslated$reason)))
})

test_that("a corpus round-trip on hz.te123.OMC.sas preserves both HAZARD and both HAZPRED blocks", {
  skip_on_cran()
  repo <- Sys.getenv("HAZARD_REPO", "~/Documents/GitHub/hazard")
  f <- file.path(path.expand(repo), "examples", "hz.te123.OMC.sas")
  skip_if_not(file.exists(f), "hz.te123.OMC.sas not available")

  job <- suppressWarnings(hzr_translate_sas(f))

  # HAZARD 2, HAZPRED 2 (measured 2026-08-20): both fits and both predict()
  # families -- four calls of the expected kinds -- must all survive.
  expect_true(all(c("fit", "fit_2", "pred", "pred_2") %in% names(job$calls)))
  expect_equal(job$outhaz, c("OUTEST", "OUTEST"))
  expect_false(identical(job$calls$fit, job$calls$fit_2))
  # Both PROC HAZPRED blocks in this file omit NOSURV/NOHAZ, so each also
  # emits a paired hazard predict() -- these must survive too, not just the
  # four minimum kinds.
  expect_true(all(c("pred_haz", "pred_haz_2") %in% names(job$calls)))
})
