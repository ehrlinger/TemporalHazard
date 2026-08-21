test_that("ICENSOR produces interval code 2 and status-gated bounds", {
  # ICENSOR's grammar is `ICENSOR flag = timevar;` (src/hazard/hazard_y.y):
  # an interval-censoring INDICATOR, not a bound variable -- ICFLAG is 0/1,
  # not NA-gated like a bound would be.
  st <- list(EVENT = "DEAD", TIME = "T", ICENSOR = c("ICFLAG", "ICTIME"))
  got <- .hzr_censor_spec(st)
  expect_equal(got$status_name, as.name(".hzr_status"))
  # time_lower/time_upper are gated on status, not unconditional -- passing
  # bounds for every row is wrong for every non-interval row, and the
  # interval's orientation (TIME vs. ICTIME) is unverified, so both are
  # emitted as pmin()/pmax() of the two, correct under either reading.
  env <- list(.hzr_status = c(1, 2), T = c(3, 8), ICTIME = c(99, 5))
  expect_equal(eval(got$time_lower, env), c(0, 5))
  expect_equal(eval(got$time_upper, env), c(3, 8))
  # Assert the exact codes. `all(status %in% c(0,1,2,3))` cannot fail and is
  # the assertion shape that hid this class of bug for years.
  expect_equal(
    eval(got$status_expr, list(DEAD = c(1, 0), ICFLAG = c(1, 1), ICTIME = c(2, 2))),
    c(1, 2)
  )
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

test_that("an ICENSOR-only job needs no EVENT statement", {
  # HAZARD terminates only when BOTH EVENT and ICENSOR are missing
  # (src/hazard/varterm.c), so ICENSOR alone is a legitimate job.
  st <- list(TIME = "T", ICENSOR = c("ICFLAG", "ICTIME"))
  got <- .hzr_censor_spec(st)
  expect_equal(eval(got$status_expr, list(ICFLAG = c(1, 0), ICTIME = c(2, NA))),
               c(2, 0))
})

test_that("a job with neither EVENT nor ICENSOR is rejected", {
  expect_error(.hzr_censor_spec(list(TIME = "T")), "EVENT|ICENSOR")
})

test_that("an event outranks a left-censoring flag", {
  st <- list(EVENT = "DEAD", TIME = "T", LCENSOR = "LFLAG")
  got <- .hzr_censor_spec(st)
  expect_equal(eval(got$status_expr, list(DEAD = c(1, 0), LFLAG = c(1, 1))),
               c(1, -1))
})

test_that("ICENSOR wins over LCENSOR where both fire", {
  st <- list(EVENT = "DEAD", TIME = "T", LCENSOR = "LFLAG",
             ICENSOR = c("ICFLAG", "ICTIME"))
  got <- .hzr_censor_spec(st)
  expect_equal(
    eval(got$status_expr,
         list(DEAD = c(0, 0), LFLAG = c(1, 1), ICFLAG = c(1, 0), ICTIME = c(5, NA))),
    c(2, -1)
  )
})
