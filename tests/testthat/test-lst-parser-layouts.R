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

test_that("SAS missing markers become NA without a coercion warning", {
  # SAS prints a bare "." for missing. Coercing it with as.numeric() yields NA
  # anyway, but warns once per row -- noise that buries real warnings.
  # .hzr_parse_sas_lifetable() already normalises this; the nomogram must too.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(
    "          Obs     YEARS     _SURVIV    _CLLSURV    _CLUSURV    _HAZARD    _CLLHAZ    _CLUHAZ",
    "",
    "            1    0.0821     0.97106     0.96812     0.97373    0.24456    0.22241    0.26891",
    "",
    "            2    0.2500     0.94683     0.94276     0.95062          .          .          .",
    ""
  ), tmp)

  expect_no_warning(nom <- .hzr_parse_sas_nomogram(tmp))
  expect_equal(nrow(nom), 2L)
  expect_equal(nom$SURVIV[2], 0.94683, tolerance = 1e-6)
  expect_true(is.na(nom$HAZARD[2]))
  expect_true(is.na(nom$CLLHAZ[2]))
  expect_true(is.na(nom$CLUHAZ[2]))
  # A missing marker in one column must not disturb the others on that row.
  expect_equal(nom$CLUSURV[2], 0.95062, tolerance = 1e-6)
})

test_that("the nomogram parser still returns NULL when no table is present", {
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c("nothing here", "no nomogram at all"), tmp)
  expect_null(.hzr_parse_sas_nomogram(tmp))
})

test_that("the life table parses whatever the time variable is called", {
  # The second column is the study's time variable: INT_DEAD in the AVC
  # captures, iv_dead in the CCF AVR/LV listings. Hardcoding one study's name
  # made every other study's life table return NULL.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(
    "       Obs iv_dead NUMBER CENSORED dead CUM_SURV   SE_EXACT CL_LOWER CL_UPPER",
    "",
    "         1  0.0821   3049       12   35  0.98852    0.00193  0.98659  0.99045",
    "         2  0.2500   3002       31   48  0.97281    0.00297  0.96984  0.97578",
    ""
  ), tmp)

  lt <- .hzr_parse_sas_lifetable(tmp, which = "kaplan")
  expect_s3_class(lt, "data.frame")
  expect_equal(nrow(lt), 2L)
  expect_true("iv_dead" %in% names(lt))
  expect_equal(lt$CUM_SURV[1], 0.98852, tolerance = 1e-6)
  expect_equal(lt$SE_EXACT[2], 0.00297, tolerance = 1e-6)
})

test_that("the life table still parses the AVC INT_DEAD layout", {
  # Regression guard: the fixtures this parser was built against must keep
  # working with the loosened header match.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(
    "       Obs INT_DEAD NUMBER CENSORED DEAD CUM_SURV   SE_EXACT CL_LOWER CL_UPPER",
    "",
    "         1   0.0821    310        2    3  0.99032    0.00556  0.98476  0.99588",
    ""
  ), tmp)

  lt <- .hzr_parse_sas_lifetable(tmp, which = "kaplan")
  expect_equal(nrow(lt), 1L)
  expect_true("INT_DEAD" %in% names(lt))
  expect_equal(lt$CUM_SURV[1], 0.99032, tolerance = 1e-6)
})

test_that("a paginated life table is read past the first page break", {
  # SAS paginates long life tables: form feed, titles, a _CATG rule, a
  # "(continued)" marker, and a repeated header, with the Obs counter running
  # on. Stopping at the first page break returned page one only -- 97 rows of
  # an 844-row table on the AVR/LV listing, across 56 page headers.
  #
  # This is the worst failure shape for a parity harness: it returns REAL
  # DATA, so no absence guard fires, the comparison covers a fraction of the
  # table, and the badge reads PASS.
  fx <- testthat::test_path("fixtures", "paginated-lifetable.lst")
  lt <- .hzr_parse_sas_lifetable(fx, which = "kaplan")

  expect_s3_class(lt, "data.frame")
  expect_equal(nrow(lt), 15L)            # 3 pages x 5 rows, not 5
  expect_true("iv_dead" %in% names(lt))
  # Rows from the 2nd and 3rd pages are present and in order.
  expect_equal(lt$iv_dead[6],  0.06, tolerance = 1e-8)
  expect_equal(lt$iv_dead[15], 0.15, tolerance = 1e-8)
  expect_false(is.unsorted(lt$iv_dead))
})

test_that("a following BY group is not absorbed into the previous one", {
  # The stratified tables share the overall table's columns, so "same header"
  # cannot be the continuation test on its own. The _CATG rule is the boundary.
  fx <- testthat::test_path("fixtures", "paginated-lifetable.lst")
  lt <- .hzr_parse_sas_lifetable(fx, which = "kaplan")

  expect_equal(nrow(lt), 15L)
  # The _CATG=1 group contributes 4 more rows; they must be absent here.
  expect_equal(max(lt$iv_dead), 0.15, tolerance = 1e-8)
  expect_true(all(lt$NUMBER > 1000))     # group 1 rows carry NUMBER < 1000
})

test_that("a listing with no fit yields no fits, not one hollow fit", {
  # .hzr_split_fits() fell back to the whole listing as one block when no
  # "Initial Summary:" anchor was present, so a %KAPLAN job produced ONE fit
  # with every field NA. That survives is.null() and a length check and reads
  # downstream as a result.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c("Life Table Analyses", "no fit anywhere in this listing",
               "       Obs iv_dead NUMBER", "         1  0.0821   3049"), tmp)

  out <- .hzr_parse_sas_lst(tmp)
  expect_false(is.null(out))
  expect_equal(length(out$fits), 0L)
})

test_that("a listing with a fit but no anchor still parses as one fit", {
  # The fallback is right when there IS a fit; only the no-fit case changed.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c("Some listing without the usual anchor",
               "        Log likelihood =     -3659.01",
               "      There are  3049 observations available for analysis with:",
               "                           1032 events:",
               "                           2017 Right Censored Observations"), tmp)

  out <- .hzr_parse_sas_lst(tmp)
  expect_equal(length(out$fits), 1L)
  expect_equal(out$fits[[1]]$loglik, -3659.01, tolerance = 1e-8)
  expect_equal(out$fits[[1]]$n_events, 1032L)
})

test_that("a listing that only MENTIONS log likelihood yields no fits", {
  # The no-anchor fallback must ask the same question the extractor asks. A
  # bare mention with no parsable value would pass a phrase test, then produce
  # loglik = NA -- rebuilding the hollow fit the guard exists to prevent.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c("Notes on the Log likelihood reported by PROC HAZARD",
               "(no value is printed in this listing)"), tmp)

  out <- .hzr_parse_sas_lst(tmp)
  expect_equal(length(out$fits), 0L)
})

test_that("a header on the last line returns NULL rather than walking backwards", {
  # (pick + 1L):length(lines) is a DECREASING sequence when pick is the last
  # line, so the loop would count down from an NA index.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c("Life Table Analyses",
               "       Obs iv_dead NUMBER CENSORED dead CUM_SURV"), tmp)

  expect_null(.hzr_parse_sas_lifetable(tmp, which = "kaplan"))
})

test_that("a header with only blank lines after it returns NULL", {
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c("       Obs iv_dead NUMBER CENSORED dead CUM_SURV", "", "  "), tmp)

  expect_null(.hzr_parse_sas_lifetable(tmp, which = "kaplan"))
})
