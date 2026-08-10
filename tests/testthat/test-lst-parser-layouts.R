# Layout variants in SAS HAZARD .lst output.
#
# The parsers were built against the AVC/KUL captures, which happen to print
# one particular shape of the observation-count block and one particular set
# of nomogram columns. Production listings from other studies print others.
# Both variants below were found in a 2006-vintage CCF listing
# (hz.dead_JR.lst, AVR/LV-function survival study); the fixture here
# reproduces the two layouts synthetically, with no patient data.

fx <- testthat::test_path("fixtures", "layout-variants.lst")

test_that("observation counts parse when the events line ends in a colon", {
  # AVC prints "545 events"; this study prints "1032 events:" with the
  # per-censoring-type breakdown nested underneath. The old anchor
  # \\bevents\\s*$ rejected the colon and returned NA.
  counts <- .hzr_extract_obs_counts(.hzr_read_lst(fx))
  expect_equal(counts$n_obs, 3049L)
  expect_equal(counts$n_events, 1032L)
  expect_equal(counts$n_censored, 2017L)
})

test_that("the nested Uncensored / Interval Censored lines are not mistaken for the total", {
  # "1029 Uncensored" and "3 Interval Censored" sit between the events line
  # and the right-censored line. A looser regex could pick one of them up.
  counts <- .hzr_extract_obs_counts(.hzr_read_lst(fx))
  expect_false(identical(counts$n_events, 1029L))
  expect_false(identical(counts$n_events, 3L))
})

test_that("the nomogram parses without a MONTHS column", {
  # hz.dead_JR prints Obs YEARS _SURVIV _CLLSURV _CLUSURV _HAZARD _CLLHAZ
  # _CLUHAZ -- no MONTHS. The old header regex required YEARS\\s+MONTHS and
  # returned NULL, and the row reader was hardcoded to toks[2:9].
  nom <- .hzr_parse_sas_nomogram(fx)
  expect_s3_class(nom, "data.frame")
  expect_equal(nrow(nom), 3L)
  expect_equal(names(nom),
               c("YEARS", "SURVIV", "CLLSURV", "CLUSURV",
                 "HAZARD", "CLLHAZ", "CLUHAZ"))
  expect_false("MONTHS" %in% names(nom))
})

test_that("nomogram values are read from the right columns", {
  nom <- .hzr_parse_sas_nomogram(fx)
  # Row 1 is 30 days (30 / 365.2425 = 0.0821).
  expect_equal(nom$YEARS[1], 0.0821, tolerance = 1e-6)
  expect_equal(nom$SURVIV[1], 0.97106, tolerance = 1e-6)
  expect_equal(nom$CLLSURV[1], 0.96812, tolerance = 1e-6)
  expect_equal(nom$CLUHAZ[1], 0.26891, tolerance = 1e-6)
  expect_equal(nom$SURVIV[3], 0.91050, tolerance = 1e-6)
})

test_that("a nomogram WITH a MONTHS column still parses", {
  # Regression guard: the AVC fixtures print MONTHS, and must keep working.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(
    "                    Obs     YEARS   MONTHS    _SURVIV   _CLLSURV   _CLUSURV   _HAZARD   _CLLHAZ   _CLUHAZ",
    "",
    "                      1    0.0821     1.00    0.97106    0.96812    0.97373   0.24456   0.22241   0.26891",
    ""
  ), tmp)
  nom <- .hzr_parse_sas_nomogram(tmp)
  expect_equal(names(nom),
               c("YEARS", "MONTHS", "SURVIV", "CLLSURV", "CLUSURV",
                 "HAZARD", "CLLHAZ", "CLUHAZ"))
  expect_equal(nom$MONTHS[1], 1.00, tolerance = 1e-6)
  expect_equal(nom$SURVIV[1], 0.97106, tolerance = 1e-6)
})

test_that("the nomogram parser still returns NULL when no table is present", {
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c("nothing here", "no nomogram at all"), tmp)
  expect_null(.hzr_parse_sas_nomogram(tmp))
})
