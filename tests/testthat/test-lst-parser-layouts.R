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

test_that("the nomogram parses when the counter column is labelled OBS", {
  # #184. SAS's casing for PROC PRINT's counter is not stable across the
  # corpus. Dropping it with cols[cols != "Obs"] left "OBS" in the header
  # names, which made the row-width guard reject every data row -- so the
  # parser returned NULL as though the job had printed no nomogram at all.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(
    "        OBS    MONTHS     YEARS     _SURVIV   _CLLSURV   _CLUSURV   _HAZARD   _CLLHAZ   _CLUHAZ",
    "",
    "          1     0.986    0.0821     0.99983    0.99961    0.99993   0.00208   0.00091   0.00477",
    "",
    "          2     3.000    0.2500     0.99946    0.99879    0.99976   0.00238   0.00112   0.00504",
    ""
  ), tmp)
  nom <- .hzr_parse_sas_nomogram(tmp)
  expect_s3_class(nom, "data.frame")
  expect_equal(nrow(nom), 2L)
  expect_equal(names(nom),
               c("MONTHS", "YEARS", "SURVIV", "CLLSURV", "CLUSURV",
                 "HAZARD", "CLLHAZ", "CLUHAZ"))
  # The counter must be gone, not shifted into the first data column.
  expect_false("OBS" %in% names(nom))
  expect_equal(nom$MONTHS[1], 0.986, tolerance = 1e-6)
  expect_equal(nom$YEARS[2], 0.2500, tolerance = 1e-6)
  expect_equal(nom$CLUHAZ[2], 0.00504, tolerance = 1e-6)
})

test_that("the counter is dropped under any of its known labels", {
  # The fix for #184 must not trade one hardcoded label for two, so the set is
  # case-insensitive and holds every spelling the corpus has shown. The label
  # is what decides: a column the header names as a counter is dropped whatever
  # values it carries, which is why a paginated counter starting at 41 is
  # dropped too. The 1..n run is only the fallback for an unrecognised name,
  # and it is not sufficient alone -- a whole-year YEARS grid runs 1..n and was
  # being deleted (#212).
  # Every member of counter_labels, so the set cannot be silently trimmed.
  # "_N_" is written as SAS prints it; the parser strips the leading
  # underscore before matching, which is why the set holds "n_".
  for (label in c("Obs", "OBS", "Observation", "Row", "_N_", "#")) {
    tmp <- withr::local_tempfile(fileext = ".lst")
    writeLines(c(
      paste0("        ", label,
             "     YEARS     _SURVIV   _CLLSURV   _CLUSURV   _HAZARD   _CLLHAZ   _CLUHAZ"),
      "",
      "          1    0.0821     0.97106    0.96812    0.97373   0.24456   0.22241   0.26891",
      "          2    0.2500     0.94683    0.94276    0.95062   0.09424   0.08597   0.10332",
      ""
    ), tmp)
    nom <- .hzr_parse_sas_nomogram(tmp)
    expect_equal(names(nom),
                 c("YEARS", "SURVIV", "CLLSURV", "CLUSURV",
                   "HAZARD", "CLLHAZ", "CLUHAZ"),
                 info = label)
    expect_equal(nom$YEARS[1], 0.0821, tolerance = 1e-6, info = label)
  }
})

test_that("a counter emitted but absent from the header is dropped on width", {
  # The other arm of the same block, which nothing covered: PROC PRINT can
  # print the counter while the header names only the measurements, so the
  # data row is one token wider than the header. That arm is decided by width
  # alone and needs no label, which is why it survives an unknown spelling.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(
    "           YEARS     _SURVIV   _CLLSURV   _CLUSURV   _HAZARD   _CLLHAZ   _CLUHAZ",
    "",
    "      1    0.0821     0.97106    0.96812    0.97373   0.24456   0.22241   0.26891",
    "      2    0.2500     0.94683    0.94276    0.95062   0.09424   0.08597   0.10332",
    ""
  ), tmp)

  nom <- .hzr_parse_sas_nomogram(tmp)
  expect_equal(names(nom),
               c("YEARS", "SURVIV", "CLLSURV", "CLUSURV",
                 "HAZARD", "CLLHAZ", "CLUHAZ"))
  expect_equal(nom$YEARS, c(0.0821, 0.2500), tolerance = 1e-6)
})

test_that("a labelled counter is dropped even when it does not start at 1", {
  # Page two of a paginated PROC PRINT starts at 41. Requiring a 1..n run as
  # well as the label left that counter in the data as a measurement column,
  # silently -- the same paginated-listing family as the life table that was
  # read only to page one.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(
    "             Obs     YEARS     _SURVIV   _CLLSURV   _CLUSURV   _HAZARD   _CLLHAZ   _CLUHAZ",
    "",
    "              41    0.0821     0.97106    0.96812    0.97373   0.24456   0.22241   0.26891",
    "              42    0.2500     0.94683    0.94276    0.95062   0.09424   0.08597   0.10332",
    ""
  ), tmp)

  nom <- .hzr_parse_sas_nomogram(tmp)
  expect_equal(names(nom),
               c("YEARS", "SURVIV", "CLLSURV", "CLUSURV",
                 "HAZARD", "CLLHAZ", "CLUHAZ"))
  expect_equal(nom$YEARS, c(0.0821, 0.2500), tolerance = 1e-6)
})

test_that("an unlabelled column that is not a 1..n run is kept in silence", {
  # The third branch is only for the ambiguous case. A leading measurement
  # that looks nothing like a counter must not warn, or the warning becomes
  # noise a reader learns to skip.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(
    "           YEARS     _SURVIV   _CLLSURV   _CLUSURV   _HAZARD   _CLLHAZ   _CLUHAZ",
    "",
    "          0.0821     0.97106    0.96812    0.97373   0.24456   0.22241   0.26891",
    "          0.2500     0.94683    0.94276    0.95062   0.09424   0.08597   0.10332",
    ""
  ), tmp)

  expect_no_warning(nom <- .hzr_parse_sas_nomogram(tmp))
  expect_true("YEARS" %in% names(nom))
})

test_that("a whole-year time grid is kept, not deleted as a counter", {
  # #212. YEARS on a whole-year grid runs 1.0000, 2.0000, 3.0000 -- a gapless
  # 1-based integer run -- so the structural test alone deleted the time key.
  # `ncol == length(cols)` afterwards, so the shape guard could not catch it:
  # the listing parsed, the badge was green, and the comparison ran against a
  # table with no time column.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(
    "           YEARS     _SURVIV   _CLLSURV   _CLUSURV   _HAZARD   _CLLHAZ   _CLUHAZ",
    "",
    "          1.0000     0.97106    0.96812    0.97373   0.24456   0.22241   0.26891",
    "          2.0000     0.94683    0.94001    0.95102   0.21001   0.19004   0.23555",
    "          3.0000     0.91050    0.90112    0.92004   0.18775   0.16888   0.20901",
    ""
  ), tmp)

  expect_warning(nom <- .hzr_parse_sas_nomogram(tmp), "does not name it as a counter")
  expect_equal(names(nom),
               c("YEARS", "SURVIV", "CLLSURV", "CLUSURV",
                 "HAZARD", "CLLHAZ", "CLUHAZ"))
  expect_equal(nom$YEARS, c(1, 2, 3), tolerance = 1e-9)
  expect_equal(nom$SURVIV, c(0.97106, 0.94683, 0.91050), tolerance = 1e-6)
})

test_that("an unknown leading label on a 1..n run warns rather than guessing", {
  # If SAS starts spelling the counter something new, the ambiguity is real
  # and unresolvable from the listing alone. Keeping the column errs toward an
  # extra column, which a caller's shape check catches, rather than toward
  # deleting a real one, which nothing catches.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(
    "         SEQNO     YEARS     _SURVIV   _CLLSURV   _CLUSURV   _HAZARD   _CLLHAZ   _CLUHAZ",
    "",
    "             1    0.0821     0.97106    0.96812    0.97373   0.24456   0.22241   0.26891",
    "             2    0.2500     0.94683    0.94276    0.95062   0.09424   0.08597   0.10332",
    ""
  ), tmp)

  expect_warning(nom <- .hzr_parse_sas_nomogram(tmp), "SEQNO")
  expect_true("SEQNO" %in% names(nom))
})

test_that("a leading measurement column is not mistaken for the counter", {
  # The structural test must fire on a 1..n run and nothing else. Here the
  # counter is suppressed and YEARS leads; dropping it would silently shift
  # every column left and corrupt the whole table.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(
    "           YEARS     _SURVIV   _CLLSURV   _CLUSURV   _HAZARD   _CLLHAZ   _CLUHAZ",
    "",
    "          0.0821     0.97106    0.96812    0.97373   0.24456   0.22241   0.26891",
    "          0.2500     0.94683    0.94276    0.95062   0.09424   0.08597   0.10332",
    ""
  ), tmp)
  nom <- .hzr_parse_sas_nomogram(tmp)
  expect_equal(names(nom),
               c("YEARS", "SURVIV", "CLLSURV", "CLUSURV",
                 "HAZARD", "CLLHAZ", "CLUHAZ"))
  expect_equal(nom$YEARS, c(0.0821, 0.2500), tolerance = 1e-6)
  expect_equal(nom$SURVIV, c(0.97106, 0.94683), tolerance = 1e-6)
})

test_that("an unparseable table warns instead of returning NULL in silence", {
  # #184's real damage was the silence: a caller cannot tell "no table" from
  # "table I failed to read", so a sweep reports its own parse failures as
  # properties of the SAS corpus. Once a header has matched, NULL must talk.
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(
    "        Obs     YEARS     _SURVIV   _CLLSURV   _CLUSURV   _HAZARD   _CLLHAZ   _CLUHAZ",
    "",
    "          1     0.0821    not-a-number   0.96812   0.97373   0.24456   0.22241   0.26891",
    ""
  ), tmp)
  expect_warning(res <- .hzr_parse_sas_nomogram(tmp), "no data rows parsed")
  expect_null(res)
})

test_that("a header on the final line warns rather than erroring", {
  # (h + 1):length(lines) is a DECREASING sequence when the header is last,
  # which would walk off the end of the file. The scan range is guarded, so
  # this lands on the same "no data rows parsed" warning as any other
  # unreadable table -- not an error, and not a silent NULL.
  for (tail_line in list(character(0), "")) {
    tmp <- withr::local_tempfile(fileext = ".lst")
    writeLines(c(
      "preamble",
      "     Obs     YEARS     _SURVIV   _CLLSURV   _CLUSURV   _HAZARD   _CLLHAZ   _CLUHAZ",
      tail_line
    ), tmp)
    expect_warning(res <- .hzr_parse_sas_nomogram(tmp), "no data rows parsed")
    expect_null(res)
  }
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

# DELTA in a .lst natural-estimates table (#181).
#
# DELTA is unimplemented in R, not "absorbed by decompos()" as two comments
# used to claim. It enters rho, the time argument and the density Jacobian
# separately, and R computes the delta = 0 branch of all three. A listing with
# a non-zero DELTA therefore describes a fit against a different function, and
# comparing R to it measures that difference while reporting it as a parity
# discrepancy.

.hzr_delta_lst <- function(delta) {
  c(
    "                    Estimates for Model Parameters",
    "",
    "        Phase       Parameter     Fixed?     Estimate",
    "        ---------------------------------------------",
    paste0("        Early:      DELTA         Yes       ", format(delta)),
    "                    THALF         Yes                0.2",
    "                    NU            No                 1.4",
    "                    MUE           No          0.02283304",
    "",
    "        Asymptotic Variance-Covariance Matrix"
  )
}

test_that("a .lst with DELTA = 0 parses without a warning", {
  # The branch R actually implements. Warning here would be a false alarm, and
  # a warning that fires on the safe case gets ignored on the unsafe one.
  expect_silent(nat <- .hzr_extract_natural(.hzr_delta_lst(0)))
  # Guard the guard: assert the table really was parsed, or "silent" is
  # satisfied by a parser that returned NULL without reading anything.
  expect_true(is.data.frame(nat))
  expect_true("DELTA" %in% nat$name)
  expect_equal(nat$estimate[nat$name == "DELTA"], 0)
})

test_that("a .lst with DELTA != 0 warns that R fits a different function", {
  expect_warning(nat <- .hzr_extract_natural(.hzr_delta_lst(0.35)),
                 "different function")
  expect_warning(.hzr_extract_natural(.hzr_delta_lst(0.35)), "DELTA = 0")
  # The table is still returned -- this is a flag on the comparison, not a
  # refusal to read the listing.
  expect_true(is.data.frame(nat))
  expect_equal(nat$estimate[nat$name == "DELTA"], 0.35)
  # And a negative one is caught too, not just a positive.
  expect_warning(.hzr_extract_natural(.hzr_delta_lst(-0.2)),
                 "different function")
})

# A listing holding more than one nomogram (#183).
#
# `h <- grep(...)` may match several headers; the parser reads h[1] and never
# looks at the rest, with no warning and nothing in the return value to say a
# second existed. Same shape as the three layout defects fixed in #110 --
# silent, and found only by pointing the parser at a second study. Every file
# in the resilia corpus prints exactly one, so the code gets the right answer
# there by luck of the corpus, not by construction.

.hzr_nomo_block <- function(surv) {
  c(
    "         Obs     YEARS    _SURVIV   _CLLSURV   _CLUSURV",
    "",
    paste0("           1    0.0821    ", format(surv, nsmall = 5),
           "    0.96812    0.97373"),
    ""
  )
}

test_that("a single nomogram parses silently and reports n_found = 1", {
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(.hzr_nomo_block(0.97106), tmp)
  expect_no_warning(nom <- .hzr_parse_sas_nomogram(tmp))
  # n_found is present whatever the count, so a caller can test it without
  # knowing whether it is set only in the plural case.
  expect_identical(attr(nom, "n_found"), 1L)
  expect_equal(nom$SURVIV[1], 0.97106, tolerance = 1e-6)
})

test_that("a listing with two nomograms warns and says how many", {
  tmp <- withr::local_tempfile(fileext = ".lst")
  writeLines(c(.hzr_nomo_block(0.97106), "Initial Summary:",
               .hzr_nomo_block(0.88231)), tmp)
  expect_warning(nom <- .hzr_parse_sas_nomogram(tmp), "2 tables found")
  expect_identical(attr(nom, "n_found"), 2L)
  # Still the FIRST one -- the behaviour is unchanged, only no longer silent.
  # Asserting the value pins which table was returned; asserting only n_found
  # would pass if the parser had started returning the second.
  expect_equal(suppressWarnings(
    .hzr_parse_sas_nomogram(tmp))$SURVIV[1], 0.97106, tolerance = 1e-6)
})

test_that("a caller can parse one fit's block without a tempfile round-trip", {
  # The workaround this replaces: split the listing by fit, write each block
  # back out to a tempfile, and re-read it. A multi-fit listing needs the
  # nomogram attributed to the fit whose block contains it, not to the file.
  all_lines <- c(.hzr_nomo_block(0.97106), "Initial Summary:",
                 .hzr_nomo_block(0.88231))
  split_at <- grep("Initial Summary:", all_lines)
  second <- all_lines[split_at:length(all_lines)]

  expect_no_warning(nom2 <- .hzr_parse_sas_nomogram(lines = second))
  expect_identical(attr(nom2, "n_found"), 1L)
  # The SECOND fit's value, which the file-level call cannot reach at all.
  expect_equal(nom2$SURVIV[1], 0.88231, tolerance = 1e-6)
})

test_that("the lines argument is validated rather than silently misread", {
  expect_error(.hzr_parse_sas_nomogram(), "Supply either")
  expect_error(.hzr_parse_sas_nomogram(lines = 42), "character vector")
})
