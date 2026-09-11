# test-sas-translate-repeat.R -- hzr_translate_sas() on jobs that call the
# SAS macro %repeat. Every test evaluates the emitted chunks; expected values
# are derived by hand from ~/Documents/macro.library/repeat.sas (see
# inst/dev/REPEAT-TRANSLATOR-PLAN.md, "The synthetic fixture"), never from
# hzr_repeated_events() itself.

translate_lines <- function(lines) {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(lines, f)
  suppressWarnings(hzr_translate_sas(f))
}

test_that("a job with %repeat but no PROC HAZARD or HAZPRED is still refused", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines("%repeat(in=bd, out=events);", f)
  expect_error(hzr_translate_sas(f), "no HAZARD or HAZPRED block found")
})
