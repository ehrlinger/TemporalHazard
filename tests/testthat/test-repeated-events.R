test_that(".hzr_re_flags marks group boundaries", {
  f <- .hzr_re_flags(c("a", "a", "b", "c", "c", "c"))
  expect_equal(f$first, c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(f$last, c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
})

test_that(".hzr_re_flags handles a single row and zero rows", {
  expect_equal(.hzr_re_flags("a"), list(first = TRUE, last = TRUE))
  expect_equal(.hzr_re_flags(character(0)), list(first = logical(0), last = logical(0)))
})

test_that(".hzr_re_order sorts by id then time, with missing times first, as SAS does", {
  d <- data.frame(id = c(2, 1, 1, 1), t = c(5, 3, NA, 1))
  expect_equal(.hzr_re_order(d, "id", "t"), c(3L, 4L, 2L, 1L))
})

test_that(".hzr_re_order is stable at tied times", {
  d <- data.frame(id = c("a", "a", "a"), t = c(1, 1, 1), tag = 1:3)
  expect_equal(d$tag[.hzr_re_order(d, "id", "t")], 1:3)
})

test_that(".hzr_re_nonevent follows SAS (&eventype=0 or &eventype=.) and excludes other codes", {
  expect_equal(.hzr_re_nonevent(c(0, 1, NA, 2)), c(TRUE, FALSE, TRUE, FALSE))
})

test_that(".hzr_re_validate rejects a missing column", {
  d <- data.frame(id = 1, t = 1, fu = 1)
  expect_error(.hzr_re_validate(d, "id", "t", "fu", "ev"), "not found")
})

test_that(".hzr_re_validate rejects an input column that would be overwritten by an output column", {
  d <- data.frame(id = 1, t = 1, fu = 1, ev = 1, rcensor = 0)
  expect_error(.hzr_re_validate(d, "id", "t", "fu", "ev"), "rcensor")
})

test_that(".hzr_re_validate rejects a missing id", {
  d <- data.frame(id = c(1, NA), t = 1, fu = 1, ev = 1)
  expect_error(.hzr_re_validate(d, "id", "t", "fu", "ev"), "missing")
})

test_that("hzr_repeated_events rejects a non-numeric, non-logical indicator instead of silently dropping events", {
  # "Y"/"N" match neither `== 1` nor `== 0`, so before this fix every row
  # was dropped and censored rows backfilled: a populated, plausible frame
  # with zero events and no error. This test proves the guard fires.
  d <- data.frame(id = c("s1", "s1", "s2"), t = c(1, 3, 2), fu = 10, ev = c("Y", "Y", "N"),
    stringsAsFactors = FALSE)
  expect_error(hzr_repeated_events(d, "id", "t", "fu", "ev"), "character")
})

test_that("hzr_repeated_events warns, but does not error, when the indicator has no event at all", {
  d <- data.frame(id = c("s1", "s2"), t = c(1, 2), fu = c(10, 10), ev = c(0, 0), stringsAsFactors = FALSE)
  expect_warning(hzr_repeated_events(d, "id", "t", "fu", "ev"), "no events")
})

test_that("hzr_repeated_events rejects a missing followup value", {
  # A missing followup used to sort the appended censored row BEFORE the
  # event it terminates (NA maps to -Inf), scrambling first/last flags.
  d3 <- data.frame(id = c("s1", "s2"), t = c(1, NA), fu = c(NA, NA), ev = c(1, 0), stringsAsFactors = FALSE)
  expect_error(hzr_repeated_events(d3, "id", "t", "fu", "ev"), "followup")
})

# A four-subject fixture covering the branches stages 2-3 discriminate:
#   s1 -- two events, no leading non-event
#   s2 -- a leading non-event row then one event
#   s3 -- a single non-event row (no first event at all)
#   s4 -- one event followed by a trailing non-event row
re_fixture <- function() {
  data.frame(
    id = c("s1", "s1", "s2", "s2", "s3", "s4", "s4"),
    t  = c(1, 3, 0, 2, NA, 4, 6),
    fu = c(10, 10, 10, 10, 10, 10, 10),
    ev = c(1, 1, 0, 1, 0, 1, 0),
    stringsAsFactors = FALSE
  )
}

test_that("hzr_repeated_events rejects a factor id rather than sorting by level order", {
  # A factor sorts by level order, not value; with levels reversed relative
  # to the subject labels this would silently give a different subject
  # ordering, and therefore different `first` rows, than SAS's proc sort.
  d <- re_fixture()
  d$id <- factor(d$id, levels = c("s4", "s3", "s2", "s1"))
  expect_error(hzr_repeated_events(d, "id", "t", "fu", "ev"), "factor")
})

test_that("stage 1 adds rcensor as all zero and changes nothing else", {
  d <- re_fixture()
  out <- .hzr_re_stage1(d)
  expect_equal(out$rcensor, rep(0, 7))
  expect_equal(out[names(d)], d)
})

test_that("stage 2 keeps every first row and every event, dropping non-first non-events", {
  out <- .hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev")
  # s2's leading non-event is at t=0 so it sorts first and is KEPT as the
  # group's first row; s4's trailing non-event is dropped.
  expect_equal(out$id, c("s1", "s1", "s2", "s2", "s3", "s4"))
  expect_equal(out$t, c(1, 3, 0, 2, NA, 4))
})

test_that("stage 2 does not treat a non-1 non-0 indicator code as a non-event", {
  d <- data.frame(id = c("a", "a"), t = c(1, 2), fu = c(9, 9), ev = c(1, 2), stringsAsFactors = FALSE)
  out <- .hzr_re_stage2(.hzr_re_stage1(d), "id", "t", "ev")
  expect_equal(nrow(out), 2L)
})

test_that("stage 3 pads a subject that never had a first event", {
  out <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  s3 <- out[out$id == "s3", ]
  expect_equal(nrow(s3), 1L)
  expect_equal(s3$t, 10)
  expect_equal(s3$rcensor, 1)
})

test_that("stage 3 drops a leading non-event when the subject has later rows", {
  out <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  expect_equal(out$id, c("s1", "s1", "s2", "s3", "s4"))
  expect_equal(out$t, c(1, 3, 2, 10, 4))
})

test_that("stage 4 records first and last as numeric data columns", {
  d <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  out <- .hzr_re_stage4(d, "id", "t", "fu", "ev")
  expect_type(out$first, "double")
  expect_type(out$last, "double")
})

test_that("stage 4 appends a terminal censored row only where the last row is not already censored", {
  d <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  out <- .hzr_re_stage4(d, "id", "t", "fu", "ev")
  # s1, s2 and s4 gain a terminal row; s3 was already rcensor=1 at stage 3.
  expect_equal(as.vector(table(out$id)[c("s1", "s2", "s3", "s4")]), c(3L, 2L, 1L, 2L))
  expect_equal(nrow(out), 8L)
})

test_that("stage 4's appended row carries followup as the time and a zero indicator", {
  d <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  out <- .hzr_re_stage4(d, "id", "t", "fu", "ev")
  added <- out[out$id == "s1" & out$rcensor == 1, ]
  expect_equal(nrow(added), 1L)
  expect_equal(added$t, 10)
  expect_equal(added$ev, 0)
})

test_that("stage 4 puts the appended row after the row it was copied from", {
  d <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  out <- .hzr_re_stage4(d, "id", "t", "fu", "ev")
  s1 <- out[out$id == "s1", ]
  expect_equal(s1$t, c(1, 3, 10))
})

test_that("stage 4's appended row inherits first and last from its source row", {
  d <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  out <- .hzr_re_stage4(d, "id", "t", "fu", "ev")
  s1 <- out[out$id == "s1", ]
  # s1 has two original rows (t=1 and t=3) plus one appended row (t=10).
  # Row 1 (t=1) is first in its group: first=1, last=0.
  # Row 2 (t=3) is the last original row: first=0, last=1.
  # Row 3 (appended from row 2) must inherit those flags: first=0, last=1.
  # Stage 6 reads the carried values, so this inheritance is load-bearing.
  expect_equal(s1$first, c(1, 0, 0))
  expect_equal(s1$last, c(0, 1, 1))
})

re_through_stage4 <- function(d = re_fixture()) {
  .hzr_re_stage4(
    .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(d), "id", "t", "ev"), "id", "t", "fu", "ev"),
    "id", "t", "fu", "ev"
  )
}

test_that("stage 5 derives event from eventype=1 exactly, not from the non-event complement", {
  d <- data.frame(
    id = c("a", "a", "a"), t = c(1, 2, 3), fu = c(9, 9, 9), ev = c(1, 2, 0),
    rcensor = c(0, 1, 1), first = c(1, 0, 0), last = c(0, 0, 1), stringsAsFactors = FALSE
  )
  out <- .hzr_re_stage5(d, "ev")
  expect_equal(out$ev, c(1, 2, 0))
  expect_equal(out$event, c(1, 0, 0))
})

test_that("stage 5 keeps only censored rows and events", {
  # Row 3 (ev = 2, rcensor = 0) is neither an event nor censored, and SAS
  # drops it. The complement-of-nonevent form (!.hzr_re_nonevent(ev)) would
  # treat ev = 2 as event-like and keep it, since it is not 0 or missing --
  # this row is what makes the test able to fail against that mutant.
  d <- data.frame(
    id = c("a", "a", "a"), t = c(1, 2, 3), fu = c(9, 9, 9), ev = c(0, 1, 2),
    rcensor = c(0, 0, 0), first = c(1, 0, 0), last = c(0, 0, 1), stringsAsFactors = FALSE
  )
  out <- .hzr_re_stage5(d, "ev")
  expect_equal(nrow(out), 1L)
  expect_equal(out$t, 2)
})

test_that("stage 6 lags the event time within subject to give iv_start and iv_seg", {
  out <- .hzr_re_stage6(.hzr_re_stage5(re_through_stage4(), "ev"), "id", "t", "fu")
  s1 <- out[out$id == "s1", ]
  expect_equal(s1$t, c(1, 3, 10))
  expect_equal(s1$iv_start, c(0, 1, 3))
  expect_equal(s1$iv_seg, c(1, 2, 7))
})

test_that("stage 6 counts event repeats within subject and restarts at each subject", {
  out <- .hzr_re_stage6(.hzr_re_stage5(re_through_stage4(), "ev"), "id", "t", "fu")
  expect_equal(out$event_no[out$id == "s1"], c(1, 2, 2))
  expect_equal(out$event_no[out$id == "s2"], c(1, 1))
  expect_equal(out$event_no[out$id == "s3"], 0)
})

test_that("stage 6 sets rcensor where the event time equals end of follow-up", {
  d <- data.frame(
    id = c("a"), t = 9, fu = 9, ev = 1, rcensor = 0, first = 1, last = 1, stringsAsFactors = FALSE
  )
  out <- .hzr_re_stage6(.hzr_re_stage5(d, "ev"), "id", "t", "fu")
  expect_equal(out$rcensor, 1)
})

test_that("stage 6 sets rcensor at end of follow-up when both time and followup are missing", {
  # SAS numeric missing compares equal to itself, so `if . = . then` is TRUE:
  # a row with both &iv_event and &iv_end missing is at end of follow-up.
  # This is still SAS-faithful and correct at the stage level, but a missing
  # `followup` is now refused by hzr_repeated_events() itself (see the test
  # "hzr_repeated_events rejects a missing followup value" below), so this
  # input is unreachable through the exported function. The stages are
  # called directly here to keep the behaviour covered.
  d <- data.frame(
    id = c("a"), t = NA_real_, fu = NA_real_, ev = 1, rcensor = 0, first = 1, last = 1,
    stringsAsFactors = FALSE
  )
  out <- .hzr_re_stage6(.hzr_re_stage5(d, "ev"), "id", "t", "fu")
  expect_equal(out$rcensor, 1)
})

test_that("stage 6 bumps renewal for a subject with no events and for a trailing censored row", {
  out <- .hzr_re_stage6(.hzr_re_stage5(re_through_stage4(), "ev"), "id", "t", "fu")
  # s3 never had an event: first == 1 and event_no == 0, so renewal is bumped to 1.
  expect_equal(out$renewal[out$id == "s3"], 1)
  # s1's trailing censored row: last == 1, rcensor == 1, event == 0, so 2 -> 3.
  expect_equal(out$renewal[out$id == "s1"], c(1, 2, 3))
})

test_that("stage 6 reads the STALE first/last from stage 4, not flags recomputed after stage 5", {
  # Subject "a" has three rows at stage 4, with last == 1 on the third.
  # Stage 5 removes that third row (rcensor 0 and not an event), so the
  # carried flags leave NO row with last == 1, while flags recomputed
  # after the subset would put last == 1 on the second row.  That second
  # row satisfies the rest of the renewal bump (rcensor == 1, event == 0),
  # so the two readings give different renewal values -- which is what
  # makes this test able to fail.
  #
  #   carried    -> last c(0, 0) -> no bump -> renewal c(1, 1)
  #   recomputed -> last c(0, 1) ->    bump -> renewal c(1, 2)
  d <- data.frame(
    id      = c("a", "a", "a"),
    t       = c(1, 5, 7),
    fu      = c(9, 9, 9),
    ev      = c(1, 0, 0),
    rcensor = c(0, 1, 0),
    first   = c(1, 0, 0),
    last    = c(0, 0, 1),
    stringsAsFactors = FALSE
  )
  out <- .hzr_re_stage6(.hzr_re_stage5(d, "ev"), "id", "t", "fu")
  expect_equal(nrow(out), 2L)
  expect_equal(out$last, c(0, 0))
  expect_equal(out$renewal, c(1, 1))
})

test_that("stage 7 drops a zero-duration non-event row that is not the subject's first", {
  d <- data.frame(
    id = c("a", "a"), t = c(0, 5), fu = c(9, 9), ev = c(0, 0),
    rcensor = c(1, 1), first = c(1, 0), last = c(0, 1),
    event = c(0, 0), event_no = c(0, 0), iv_start = c(0, 5), iv_seg = c(0, 0),
    renewal = c(1, 1), stringsAsFactors = FALSE
  )
  out <- .hzr_re_stage7(d, "id", "t", "ev")
  expect_equal(nrow(out), 1L)
  expect_equal(out$iv_start, 0)
})

test_that("stage 7 tests eventype=0 exactly and keeps a zero-duration row with a missing indicator", {
  d <- data.frame(
    id = c("a", "a"), t = c(0, 5), fu = c(9, 9), ev = c(0, NA),
    rcensor = c(1, 1), first = c(1, 0), last = c(0, 1),
    event = c(0, 0), event_no = c(0, 0), iv_start = c(0, 5), iv_seg = c(0, 0),
    renewal = c(1, 1), stringsAsFactors = FALSE
  )
  expect_equal(nrow(.hzr_re_stage7(d, "id", "t", "ev")), 2L)
})

test_that("stage 7 keeps a zero-duration non-first row when iv_seg is NA rather than injecting a phantom row", {
  # iv_seg can legitimately be NA (time - iv_start with a missing event
  # time). SAS's `if &iv_seg=0 ...` is FALSE when iv_seg is missing (a
  # missing never equals 0), so the row is kept. Without a !is.na() guard
  # on the iv_seg comparison, `drop` evaluates to NA for this row and
  # `data[!drop, ]` injects an all-NA phantom row instead of keeping or
  # dropping it.
  # Row 2 is not the subject's first row, has indicator 0 (not missing), and
  # an iv_seg of NA: every conjunct but the iv_seg comparison is TRUE.
  d <- data.frame(
    id = c("a", "a"), t = c(0, 5), fu = c(9, 9), ev = c(0, 0),
    rcensor = c(1, 1), first = c(1, 0), last = c(0, 1),
    event = c(0, 0), event_no = c(0, 0), iv_start = c(0, 5), iv_seg = c(0, NA),
    renewal = c(1, 1), stringsAsFactors = FALSE
  )
  out <- .hzr_re_stage7(d, "id", "t", "ev")
  expect_equal(nrow(out), 2L)
  expect_false(anyNA(out$id))
})

test_that("hzr_repeated_events returns the documented columns and drops SAS loop state", {
  out <- hzr_repeated_events(re_fixture(), "id", "t", "fu", "ev")
  expect_true(all(c("event", "event_no", "rcensor", "iv_start", "iv_seg", "renewal", "first", "last")
                  %in% names(out)))
  expect_false(any(c("lag_iv", "number") %in% names(out)))
  expect_equal(names(out)[1:4], c("id", "t", "fu", "ev"))
})

test_that("hzr_repeated_events gives every subject at least one row", {
  out <- hzr_repeated_events(re_fixture(), "id", "t", "fu", "ev")
  expect_setequal(unique(out$id), c("s1", "s2", "s3", "s4"))
})

test_that("hzr_repeated_events produces segments that tile each subject's follow-up without gaps", {
  # The original assertions (iv_start == c(0, head(t, -1)), iv_seg == t -
  # iv_start) restated exactly how stage 6 computes those columns from `t`,
  # and on this fixture stage 7 drops no rows so they always held -- a
  # mutant that broke the tiling would still pass them if it also broke
  # them the same way. These instead constrain the OUTPUT: within each
  # subject, ordered by iv_start, each segment must start exactly where the
  # previous one ended (contiguous, no gap or overlap) and the segments
  # together must cover [0, followup] -- a property that would legitimately
  # break if stage 7 ever dropped a row, unlike the arithmetic restatement.
  out <- hzr_repeated_events(re_fixture(), "id", "t", "fu", "ev")
  for (subject in unique(out$id)) {
    rows <- out[out$id == subject, ]
    rows <- rows[order(rows$iv_start), ]
    seg_end <- rows$iv_start + rows$iv_seg
    expect_equal(rows$iv_start[1], 0, info = subject)
    if (nrow(rows) > 1L) {
      expect_equal(rows$iv_start[-1], utils::head(seg_end, -1), info = subject)
    }
    expect_equal(utils::tail(seg_end, 1), 10, info = subject)
  }
})

test_that("hzr_repeated_events rejects a column name clash rather than overwriting", {
  d <- re_fixture()
  d$renewal <- 0
  expect_error(hzr_repeated_events(d, "id", "t", "fu", "ev"), "renewal")
})

test_that("rcensor and event are not mutually exclusive: a genuine event at end of follow-up sets both", {
  # Confirmed directly: with the last event coinciding with end of follow-up,
  # the returned row carries event == 1 AND rcensor == 1 together, exactly the
  # overlap documented in \value and @note. A "fix" that made them exclusive
  # (e.g. clearing rcensor whenever event == 1) would fail this.
  d <- data.frame(id = c("a", "a"), t = c(2, 9), fu = c(9, 9), ev = c(1, 1))
  out <- hzr_repeated_events(d, "id", "t", "fu", "ev")
  last_row <- out[out$id == "a" & out$t == 9, ]
  expect_equal(nrow(last_row), 1L)
  expect_equal(last_row$event, 1)
  expect_equal(last_row$rcensor, 1)
})

test_that(".hzr_re_validate gives the intended error for a multi-element argument, not 'subscript out of bounds'", {
  # Before the fix, c(id = id, ...) FLATTENED a multi-element `id` into
  # entries named id1/id2, so length(value) != 1L could never fire here: the
  # bad call instead died later with "subscript out of bounds" from deep in
  # the pipeline. This pins the intended, earlier, clearer error.
  d <- data.frame(id = 1, t = 1, fu = 1, ev = 1)
  expect_error(.hzr_re_validate(d, c("id", "t"), "t", "fu", "ev"), "single column name")
})

test_that(".hzr_re_validate rejects a non-character argument instead of silently coercing it", {
  # c()'s type coercion, not just its flattening, made !is.character(value)
  # dead too: c(id = 5, time = "t", ...) coerces 5 to the string "5" before
  # the guard ever sees it. list() keeps the original type so the guard fires.
  d <- data.frame(id = 1, t = 1, fu = 1, ev = 1)
  expect_error(.hzr_re_validate(d, 5, "t", "fu", "ev"), "single column name")
})
