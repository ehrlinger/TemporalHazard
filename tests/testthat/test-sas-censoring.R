test_that("ICENSOR produces interval code 2 and a status-gated time_lower", {
  # ICENSOR's grammar is `ICENSOR c3 = ctime;` (src/hazard/hazard_y.y): C3 is
  # an event COUNT (OBS column 4), not a 0/1 flag -- a row is
  # interval-censored where C3 > 0.
  st <- list(EVENT = "DEAD", TIME = "T", ICENSOR = c("C3", "CTIME"))
  got <- .hzr_censor_spec(st)
  expect_equal(got$status_name, as.name(".hzr_status"))
  # time_lower is gated on status, not unconditional -- passing a bound for
  # every row is wrong for every non-interval row. No time_upper is emitted:
  # the interval's upper bound is TIME, which is exactly what hazard()'s
  # time_upper already defaults to.
  expect_null(got$time_upper)
  env <- list(.hzr_status = c(1, 2), CTIME = c(99, 5))
  expect_equal(eval(got$time_lower, env), c(0, 5))
  # Assert the exact codes. `all(status %in% c(0,1,2,3))` cannot fail and is
  # the assertion shape that hid this class of bug for years.
  expect_equal(
    eval(got$status_expr, list(DEAD = c(1, 0), C3 = c(1, 3), CTIME = c(2, 2))),
    c(1, 2)
  )
})

test_that("LCENSOR is left-truncation: status is untouched, time_lower is the entry time", {
  # LCENSOR is the counting-process entry time, not left-censoring -- HAZARD
  # has no left-censoring in this package's sense (FIXTURE-GAP-LIST.md Q1).
  # status = -1 must never be produced by this translator.
  st <- list(EVENT = "DEAD", TIME = "T", LCENSOR = "STARTTME")
  got <- .hzr_censor_spec(st)
  expect_null(got$status_name)
  expect_equal(eval(got$status_expr, list(DEAD = c(1, 0))), c(1, 0))
  expect_equal(got$time_lower, as.name("STARTTME"))
  expect_null(got$time_upper)
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

test_that("time_lower/time_upper are NULL, not names, when neither ICENSOR nor LCENSOR is given", {
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
  st <- list(TIME = "T", ICENSOR = c("C3", "CTIME"))
  got <- .hzr_censor_spec(st)
  expect_equal(eval(got$status_expr, list(C3 = c(1, 0), CTIME = c(2, NA))),
               c(2, 0))
})

test_that("a job with neither EVENT nor ICENSOR is rejected", {
  expect_error(.hzr_censor_spec(list(TIME = "T")), "EVENT|ICENSOR")
})

test_that("an event outranks interval censoring", {
  st <- list(EVENT = "DEAD", TIME = "T", ICENSOR = c("C3", "CTIME"))
  got <- .hzr_censor_spec(st)
  expect_equal(eval(got$status_expr, list(DEAD = c(1, 0), C3 = c(1, 1), CTIME = c(2, 2))),
               c(1, 2))
})

test_that("LCENSOR and ICENSOR together: time_lower picks CTIME on interval rows, LCENSOR's var elsewhere", {
  st <- list(EVENT = "DEAD", TIME = "T", LCENSOR = "STARTTME",
             ICENSOR = c("C3", "CTIME"))
  got <- .hzr_censor_spec(st)
  # status is unaffected by LCENSOR -- it stays a pure EVENT/ICENSOR mapping.
  expect_equal(
    eval(got$status_expr, list(DEAD = c(0, 0), C3 = c(1, 0))),
    c(2, 0)
  )
  env <- list(.hzr_status = c(2, 0, 1), CTIME = c(5, 5, 5), STARTTME = c(1, 1, 1))
  expect_equal(eval(got$time_lower, env), c(5, 1, 1))
})

test_that("ICENSOR's event count is carried as a status-gated weight", {
  # setlik.c: C3 is a COUNT of interval-censored individuals at TIME, and it
  # multiplies the log-likelihood contribution in BOTH terms
  # (c3w = c3 * weight; llike = -(c1w + c2 + c3w) * (cumhaz - cumhst);
  # if (c3 > ZERO) llike += c3w * lct). Discarding the magnitude fits a
  # 3-event row as a single observation (#154).
  st <- list(EVENT = "DEAD", TIME = "T", ICENSOR = c("C3", "CTIME"))
  got <- .hzr_censor_spec(st)
  # Gated on status, never the bare count. C3 is 0 on every non-interval row,
  # and hazard() multiplies each row's log-likelihood contribution by its
  # weight, so a bare C3 would delete every event and right-censored row from
  # the fit outright. readc2.c sets C2 = 1 on exactly those rows when no
  # RCENSOR statement is given, so their weight is 1.
  env <- list(.hzr_status = c(1, 2, 0), C3 = c(0, 3, 0))
  expect_equal(eval(got$weights_expr, env), c(1, 3, 1))
  # An event outranks interval censoring, so a row that is both carries the
  # event's weight of 1, not C3.
  expect_equal(eval(got$weights_expr, list(.hzr_status = 1, C3 = 3)), 1)
})

test_that("a job with no ICENSOR emits no weights expression", {
  expect_null(.hzr_censor_spec(list(EVENT = "DEAD", TIME = "T"))$weights_expr)
  expect_null(
    .hzr_censor_spec(list(EVENT = "DEAD", TIME = "T",
                          LCENSOR = "STARTTME"))$weights_expr
  )
  expect_null(
    .hzr_censor_spec(list(EVENT = "DEAD", TIME = "T",
                          RCENSOR = "RFLAG"))$weights_expr
  )
})

test_that("the emitted call carries ICENSOR's count as weights", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; ICENSOR C3 = CT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(
    eval(got$call[["weights"]], list(.hzr_status = c(2, 0), C3 = c(3, 0))),
    c(3, 1)
  )
})

test_that("ICENSOR and WEIGHT together multiply, as c3w = c3 * weight does", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; ICENSOR C3 = CT; WEIGHT WT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(
    eval(got$call[["weights"]],
         list(.hzr_status = c(2, 0), C3 = c(3, 0), WT = c(2, 2))),
    c(6, 2)
  )
})

test_that("a job with no ICENSOR emits no weights argument at all", {
  # The non-ICENSOR path must be untouched by #154: no weights argument, so
  # hazard() uses unit weights exactly as it did before.
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_false("weights" %in% names(as.list(got$call)))
  # A WEIGHT statement without ICENSOR still emits the bare weight variable.
  txt2 <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; WEIGHT WT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got2 <- .hzr_parse_hazard(.hzr_sas_blocks(txt2)[[1L]])
  expect_equal(got2$call[["weights"]], as.name("WT"))
})
