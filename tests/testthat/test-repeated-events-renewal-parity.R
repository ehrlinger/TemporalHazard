# Parity with SAS for ac.reintervention (achalasia): the renewal column.
#
# The maze cardioversion job (test-repeated-events-parity.R) never reads
# renewal. This job does -- it prints ev_rein by renewal -- and it saved the
# macro's output, unchanged, as library.bdrein. So every column
# hzr_repeated_events() returns can be compared with SAS's own, row by row.
# Every expected value below comes from the job's .log, its .lst or bdrein --
# none from the R code.
#
# The data are PHI and live only on the study volume. This test skips when the
# volume is not mounted, which is every CI runner, so a green CI run is NOT
# evidence of parity. Failure output is counts only: no row and no ccfid is
# ever printed.

skip_if_no_achalasia <- function() {
  d <- .hzr_achalasia_dir() # nolint: object_usage_linter.
  testthat::skip_if(is.na(d), "achalasia study datasets not mounted")
  d
}

test_that("ac.reintervention: %repeat reproduces the SAS log, the .lst and library.bdrein", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("haven")
  root <- skip_if_no_achalasia()

  r <- .hzr_derive_bdmult1(root)
  # The rebuild, against the .log: WORK.PNDIL1 29 x 7, WORK.REINT1 15 x 7,
  # WORK.CMB 44 x 9, WORK.BDMULT 425 rows and 218 columns, WORK.BDMULT1 425 x 219.
  expect_equal(unname(vapply(r, dim, integer(2))),
               matrix(c(29L, 7L, 15L, 7L, 44L, 9L, 425L, 218L, 425L, 219L), nrow = 2))
  bdmult1 <- r$bdmult1

  id <- "ccfid"
  tm <- "iv_rein"
  fu <- "iv_fup"
  ind <- "ev_rein"
  s1 <- .hzr_re_stage1(bdmult1)
  s2 <- .hzr_re_stage2(s1, id, tm, ind)
  s3 <- .hzr_re_stage3(s2, id, tm, fu, ind)
  s4 <- .hzr_re_stage4(s3, id, tm, fu, ind)
  s5 <- .hzr_re_stage5(s4, ind)
  s6 <- .hzr_re_stage6(s5, id, tm, fu)
  s7 <- .hzr_re_stage7(s6, id, tm, ind)

  # Row coverage first: the macro's thirteen DATA and SORT steps (.log lines
  # 349-443) give 425 x 220 through stage 3, 464 x 222 at stages 4 and 5, and
  # 464 x 228 at stages 6 and 7. From stage 5 R is one column wider: the job
  # passes event=ev_rein, the same column as eventype=, so SAS overwrites its
  # indicator where R adds `event`. Stages 6 and 7 also drop lag_iv and number,
  # the macro's loop counters.
  shapes <- rbind(dim(s1), dim(s2), dim(s3), dim(s4), dim(s5), dim(s6), dim(s7))
  expect_equal(shapes[, 1], c(425L, 425L, 425L, 464L, 464L, 464L, 464L))
  expect_equal(shapes[, 2], c(220L, 220L, 220L, 222L, 222L + 1L, 228L - 2L + 1L, 228L - 2L + 1L))

  # This data trips none of the function's input warnings. Those warnings name
  # subjects, so they are counted here, never printed; and the comparison with
  # the stages reports a bare FALSE rather than a diff of patient rows.
  warned <- testthat::capture_warnings(events <- hzr_repeated_events(bdmult1, id, tm, fu, ind))
  expect_length(warned, 0L)
  expect_true(isTRUE(all.equal(events, `row.names<-`(s7, NULL))))

  # The .lst's "Table of ev_rein by renewal".
  tab <- table(events$event, events$renewal)
  expect_equal(unname(dimnames(tab)), list(c("0", "1"), c("1", "2", "3", "4", "5")))
  expect_equal(as.vector(tab), c(381L, 39L, 36L, 3L, 2L, 1L, 0L, 1L, 1L, 0L))

  # library.bdrein, row by row. Coverage before any value: the .log's shape
  # (LIBRARY.BDREIN has 464 observations and 228 variables), then the same
  # subjects in the same order.
  sas <- haven::read_sas(file.path(root, "datasets", "bdrein.sas7bdat"))
  sas <- as.data.frame(haven::zap_labels(haven::zap_formats(sas)))
  names(sas) <- tolower(names(sas))
  expect_equal(dim(sas), c(464L, 228L))
  expect_equal(nrow(events), nrow(sas))
  expect_equal(sum(events$ccfid != sas$ccfid), 0L)

  # R's `event` and `rcensor` are the job's ev_rein and cn_rein.
  pairs <- c(iv_rein = "iv_rein", iv_fup = "iv_fup", iv_start = "iv_start", iv_seg = "iv_seg",
             event = "ev_rein", rcensor = "cn_rein", event_no = "event_no",
             renewal = "renewal", first = "first", last = "last")
  for (v in names(pairs)) {
    got <- events[[v]]
    want <- sas[[pairs[[v]]]]
    # A missing column would make both sums below zero-length, and so zero.
    expect_length(got, nrow(sas))
    # A constant column would match whatever the code did.
    expect_gt(length(unique(want)), 1L, label = paste(v, "distinct values in SAS"))
    expect_equal(sum(is.na(got) != is.na(want)), 0L, label = paste(v, "missingness mismatches"))
    # Elementwise and relative, so a small value is not waved through by an
    # absolute tolerance; a zero must be exactly zero.
    bad <- abs(got - want) > 1e-9 * abs(want)
    expect_equal(sum(bad, na.rm = TRUE), 0L, label = paste(v, "value mismatches"))
  }

  # bdrein's first and last are the stage-4 flags, carried stale into the
  # renewal step. On 75 rows they are not the subject's first or last row of
  # the output, and R matches SAS on every one of them above: SAS evidence that
  # the flags are carried forward, not recomputed.
  #
  # What this data cannot show is renewal reading them. The 36 stale `first`
  # rows are the censored rows appended after a subject's only event, which
  # copy that event's first = 1; they have event_no = 1, and the first-row bump
  # needs event_no = 0. The 39 stale `last` rows are each subject's last event;
  # they have event = 1, and the last-row bump needs event = 0. So renewal is
  # the same from stale flags or fresh ones here, and that branch is pinned
  # only by the hand-traced fixtures in test-repeated-events.R.
  fresh <- .hzr_re_flags(events$ccfid)
  stale_first <- fresh$first != (events$first == 1)
  stale_last <- fresh$last != (events$last == 1)
  expect_equal(c(sum(stale_first), sum(stale_last)), c(36L, 39L))
  expect_equal(unique(events$event_no[stale_first]), 1)
  expect_equal(unique(events$event[stale_last]), 1)
})
