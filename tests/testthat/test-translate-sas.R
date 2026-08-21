test_that("a job with no recognisable block errors and writes nothing", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(c("DATA X; SET Y;", "PROC PRINT;"), f)
  out <- withr::local_tempdir()
  expect_error(hzr_translate_sas(f, out_dir = out), "no SAS statements|no HAZARD")
  expect_length(list.files(out, pattern = "[.]qmd$"), 0L)
})

test_that("librefs resolves an INHAZ that no translated job wrote", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
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
  writeLines("%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.HZD OUT=P; TIME MONTHS; );", f)
  out <- withr::local_tempdir()
  expect_warning(hzr_translate_sas(f, out_dir = out), "unresolved")
  qmd <- readLines(file.path(out, sub("[.]sas$", ".qmd", basename(f))))
  expect_true(any(grepl("stop(", qmd, fixed = TRUE)))
})

test_that("a malformed librefs is rejected with a clear error", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines("%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.HZD OUT=P; TIME MONTHS; );", f)
  expect_error(hzr_translate_sas(f, librefs = "estimates"), "named character vector")
  expect_error(hzr_translate_sas(f, librefs = c(EX = 1)), "named character vector")
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
