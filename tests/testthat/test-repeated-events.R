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
