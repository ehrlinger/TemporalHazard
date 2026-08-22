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
    eval(got$status_expr, list(DEAD = c(1, 0), C3 = c(0, 3), CTIME = c(2, 2))),
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

test_that("no censoring statement still derives status from EVENT > 0, not from EVENT", {
  # EVENT is a COUNT (readc1.c keeps any C1 >= 0). Mapped straight onto
  # status, DEAD = 2 reads as INTERVAL-CENSORED in this package's -1/0/1/2
  # coding and DEAD = 3 falls outside it altogether (#157).
  st <- list(EVENT = "DEAD", TIME = "T")
  got <- .hzr_censor_spec(st)
  expect_equal(eval(got$status_expr, list(DEAD = c(0, 1, 2, 3))),
               c(0, 1, 1, 1))
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

test_that("a row that is both an event and interval-censored stops, rather than picking one", {
  # setlik.c SUMS the contributions of one OBS record
  # (c1c2c3 = c1w + c2 + c3w; llike gains c1w*log h AND c3w*lct), so such a
  # row is an event of weight C1*WT and an interval observation of weight
  # C3*WT at once. hazard() has one status and one weight per row and cannot
  # say that. This translator used to pick the event branch and drop c3w --
  # a fit that converges and reports plausibly over a model it is not.
  st <- list(EVENT = "DEAD", TIME = "T", ICENSOR = c("C3", "CTIME"))
  got <- .hzr_censor_spec(st)
  expect_error(
    eval(got$status_expr,
         list(DEAD = c(1, 0), C3 = c(1, 1), CTIME = c(2, 2))),
    "both non-zero"
  )
  # Non-mixed rows are unaffected: the guard is not a blanket refusal of
  # EVENT + ICENSOR jobs.
  expect_equal(
    eval(got$status_expr, list(DEAD = c(1, 0), C3 = c(0, 1), CTIME = c(2, 2))),
    c(1, 2)
  )
})

test_that("the both-fire guard names the two weights the split-out rows need", {
  st <- list(EVENT = "DEAD", TIME = "T", ICENSOR = c("C3", "CTIME"),
             WEIGHT = "WT")
  got <- .hzr_censor_spec(st)
  err <- tryCatch(
    eval(got$status_expr, list(DEAD = 1, C3 = 1, WT = 2)),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "DEAD \\* WT", fixed = FALSE)
  expect_match(err, "C3 \\* WT", fixed = FALSE)
})

test_that("LCENSOR and ICENSOR together are refused, not translated", {
  # This test used to assert the opposite -- time_lower = ifelse(status == 2,
  # CTIME, STARTTME) -- which reads as if it preserved the entry time. It does
  # not: on the interval rows the entry time is gone, so a left-truncated
  # interval-censored subject is fitted as at risk from time 0. There is no
  # expression that serves, because time_lower is one column carrying two
  # meanings. SAS carries three times (TIME, CTIME, STIME) and subtracts
  # H(STIME) for every row, so a faithful translation needs a hazard()
  # argument that does not exist yet (#155). Refuse until it does.
  st <- list(EVENT = "DEAD", TIME = "T", LCENSOR = "STARTTME",
             ICENSOR = c("C3", "CTIME"))
  got <- .hzr_censor_spec(st)
  expect_true(isTRUE(got$refused))
  expect_equal(got$untranslated$construct, "LCENSOR + ICENSOR")
  expect_match(got$untranslated$reason, "entry-time")
  # Nothing is emitted: a partial spec is how a wrong fit gets built anyway.
  expect_null(got$status_expr)
  expect_null(got$status_name)
  expect_null(got$time_lower)
  expect_null(got$weights_expr)
  # Both returns carry the same element names, so a caller reading one never
  # meets an unexpected missing field.
  expect_equal(
    sort(names(got)),
    sort(names(.hzr_censor_spec(list(EVENT = "DEAD", TIME = "T"))))
  )
})

test_that("LCENSOR combined with ICENSOR is refused, not silently mis-fitted", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; LCENSOR ST; ICENSOR C3 = CT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_true(any(grepl("LCENSOR", got$untranslated$construct)))
  # time_lower cannot carry both an entry time and an interval lower bound,
  # so the job must not translate to a silently wrong fit (#155).
  expect_null(got$call[["time_lower"]])
  # The callout alone is not enough -- a reader who renders past it would
  # still get a fit. The emitted call must itself refuse.
  expect_identical(got$call[[1L]], as.name("stop"))
  expect_error(eval(got$call), "LCENSOR")
  expect_null(got$status_call)
})

test_that("the refusal does not leak into LCENSOR-only or ICENSOR-only jobs", {
  # ICENSOR appears in 64 of 93 production PROC HAZARD jobs and LCENSOR in 10;
  # a refusal that fired on either alone would be far worse than the bug it
  # replaces.
  lonly <- .hzr_censor_spec(list(EVENT = "DEAD", TIME = "T",
                                 LCENSOR = "STARTTME"))
  expect_false(isTRUE(lonly$refused))
  expect_equal(lonly$time_lower, as.name("STARTTME"))
  expect_equal(nrow(lonly$untranslated), 0L)

  ionly <- .hzr_censor_spec(list(EVENT = "DEAD", TIME = "T",
                                 ICENSOR = c("C3", "CTIME")))
  expect_false(isTRUE(ionly$refused))
  expect_equal(
    eval(ionly$time_lower, list(.hzr_status = c(2, 0), CTIME = c(5, 5))),
    c(5, 0)
  )
  expect_equal(nrow(ionly$untranslated), 0L)

  # And through the whole block parser, not just the spec helper: both still
  # emit a real hazard() call.
  for (stmt in c("LCENSOR ST;", "ICENSOR C3 = CT;")) {
    txt <- .hzr_sas_normalise(paste(
      "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
      "EVENT DEAD; TIME T;", stmt,
      "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
    ))
    got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
    expect_identical(got$call[[1L]], as.name("hazard"))
    expect_false(any(grepl("LCENSOR", got$untranslated$construct)))
    expect_false(is.null(got$call[["time_lower"]]))
  }
})

test_that("the weight is each row's own count, and 1 on right-censored rows", {
  # setlik.c: c1c2c3 = c1w + c2 + c3w, with c1w = C1 * WT and c3w = C3 * WT.
  # C2 -- which readc2.c sets to ONE on exactly the rows where neither other
  # count fires -- is the ONE term not multiplied by the weight (#158).
  # Discarding a magnitude fits a 3-event row as a single observation (#154,
  # #157); emitting a bare count gives every other row a weight of ZERO,
  # deleting it from the likelihood without a word.
  st <- list(EVENT = "DEAD", TIME = "T", ICENSOR = c("C3", "CTIME"))
  got <- .hzr_censor_spec(st)
  env <- list(DEAD = c(2, 0, 0), C3 = c(0, 3, 0))
  expect_equal(eval(got$weights_expr, env), c(2, 3, 1))
  expect_false(any(eval(got$weights_expr, env) == 0))
})

test_that("EVENT's count is carried as a weight even with no ICENSOR", {
  # #157: c1w = C1 * weight, so the event COUNT weights the row's whole
  # contribution whether or not the job has an ICENSOR statement.
  got <- .hzr_censor_spec(list(EVENT = "DEAD", TIME = "T"))
  expect_equal(eval(got$weights_expr, list(DEAD = c(0, 1, 2, 3))),
               c(1, 1, 2, 3))
})

test_that("WEIGHT multiplies the event count but leaves censored rows at 1", {
  # c1w = C1 * WT, but c2 enters the sum unweighted (#158). A WEIGHT that is
  # 0 on censored rows -- the shape of hz.tm123.OMC.sas's WEIGHT MORBID --
  # deletes them from the fit outright when emitted bare.
  got <- .hzr_censor_spec(list(EVENT = "DEAD", TIME = "T", WEIGHT = "WT"))
  w <- eval(got$weights_expr, list(DEAD = c(1, 0, 2), WT = c(4, 0, 3)))
  expect_equal(w, c(4, 1, 6))
  expect_false(any(w == 0))
})

test_that("an ICENSOR-only job weights interval rows by C3 * WT and the rest by 1", {
  got <- .hzr_censor_spec(list(TIME = "T", ICENSOR = c("C3", "CTIME"),
                               WEIGHT = "WT"))
  expect_equal(eval(got$weights_expr, list(C3 = c(3, 0), WT = c(2, 5))),
               c(6, 1))
})

test_that("the emitted call carries ICENSOR's count as weights", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; ICENSOR C3 = CT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(
    eval(got$call[["weights"]], list(DEAD = c(0, 0), C3 = c(3, 0))),
    c(3, 1)
  )
})

test_that("ICENSOR and WEIGHT multiply on interval rows, as c3w = c3 * weight does", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; ICENSOR C3 = CT; WEIGHT WT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(
    eval(got$call[["weights"]],
         list(DEAD = c(0, 0), C3 = c(3, 0), WT = c(2, 2))),
    c(6, 1)
  )
})

test_that("a WEIGHT job's censored rows are not weighted by the WEIGHT variable", {
  # #158's decisive shape at the emitted-call level: the translator used to
  # emit a bare WT here, which over-weights every censored row and deletes
  # them entirely when WT is 0 there.
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; WEIGHT WT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(
    eval(got$call[["weights"]], list(DEAD = c(1, 0), WT = c(4, 0))),
    c(4, 1)
  )
})

test_that("a plain 0/1 EVENT job's emitted weights are unit weights, as before", {
  # #154 non-regression: nothing about a job whose counts are all 0/1 may
  # change. Asserted on the evaluated weights, not on the call's shape.
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(eval(got$call[["weights"]], list(DEAD = c(1, 0, 1))),
               c(1, 1, 1))
  # And an ICENSOR job with a 0/1 EVENT and no WEIGHT still weights interval
  # rows by C3 and everything else by 1 -- what #154 shipped.
  txt2 <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; ICENSOR C3 = CT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got2 <- .hzr_parse_hazard(.hzr_sas_blocks(txt2)[[1L]])
  env <- list(DEAD = c(1, 0, 0), C3 = c(0, 3, 0))
  expect_equal(eval(got2$call[["weights"]], env), c(1, 3, 1))
})

test_that("the refusal reaches both the untranslated callout and a stop() chunk", {
  # Two places, and the second is the one that matters: a reader who scrolls
  # past a callout still gets whatever the fit chunk produces.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(c(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; LCENSOR ST; ICENSOR C3 = CT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  expect_warning(hzr_translate_sas(f), "LCENSOR")
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_true("LCENSOR + ICENSOR" %in% job$untranslated$construct)
  qmd <- paste(.hzr_render_qmd(job), collapse = "\n")
  expect_match(qmd, "UNTRANSLATED: LCENSOR \\+ ICENSOR")
  expect_match(qmd, "stop\\(")
  # No fit is emitted at all -- not a fit guarded by a warning.
  expect_false(grepl("fit <- hazard(", qmd, fixed = TRUE))
})
