test_that("ICENSOR produces interval code 2 and both bounds", {
  st <- list(EVENT = "DEAD", TIME = "T", ICENSOR = c("LO", "HI"))
  got <- .hzr_censor_spec(st)
  expect_equal(got$time_lower, as.name("LO"))
  expect_equal(got$time_upper, as.name("HI"))
  # Assert the exact codes. `all(status %in% c(0,1,2,3))` cannot fail and is
  # the assertion shape that hid this class of bug for years.
  expect_equal(eval(got$status_expr, list(DEAD = c(1, 0), LO = c(1, 1), HI = c(2, 2))),
               c(1, 2))
})

test_that("LCENSOR produces -1, not survival's 2", {
  st <- list(EVENT = "DEAD", TIME = "T", LCENSOR = "LFLAG")
  got <- .hzr_censor_spec(st)
  expect_equal(eval(got$status_expr, list(DEAD = c(0, 0), LFLAG = c(1, 0))),
               c(-1, 0))
})

test_that("RCENSOR produces 0", {
  st <- list(EVENT = "DEAD", TIME = "T", RCENSOR = "RFLAG")
  got <- .hzr_censor_spec(st)
  expect_equal(eval(got$status_expr, list(DEAD = c(0, 1), RFLAG = c(1, 0))),
               c(0, 1))
})

test_that("no censoring statement leaves status as the bare EVENT variable", {
  st <- list(EVENT = "DEAD", TIME = "T")
  got <- .hzr_censor_spec(st)
  expect_equal(got$status_expr, as.name("DEAD"))
  expect_null(got$time_lower)
  expect_null(got$time_upper)
})

test_that("time_lower/time_upper are NULL, not names, when ICENSOR is absent", {
  st <- list(EVENT = "DEAD", TIME = "T", RCENSOR = "RFLAG")
  got <- .hzr_censor_spec(st)
  expect_null(got$time_lower)
  expect_null(got$time_upper)
})

test_that("untranslated is a well-formed empty frame, not merely non-NULL", {
  st <- list(EVENT = "DEAD", TIME = "T")
  got <- .hzr_censor_spec(st)
  expect_equal(nrow(got$untranslated), 0L)
  expect_equal(names(got$untranslated), c("line", "construct", "reason"))
})
