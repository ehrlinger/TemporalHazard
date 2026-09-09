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
