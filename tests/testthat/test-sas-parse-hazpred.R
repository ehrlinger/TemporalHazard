test_that("a HAZPRED block becomes a predict() call with se.fit", {
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=EX.HZD OUT=PREDICT; TIME MONTHS; );"
  )
  b <- .hzr_sas_blocks(txt)[[1L]]
  got <- .hzr_parse_hazpred(b, txt)
  expect_equal(got$inhaz, "EX.HZD")
  expect_equal(got$call[["se.fit"]], TRUE)
  expect_equal(got$call[["newdata"]], as.name("PREDICT"))
})

test_that("a log-spaced DO grid becomes an exp(seq(...)) call", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; MAX=180; LN_MAX=LOG(MAX); INC=(5+LN_MAX)/99.9;",
    "DO LN_TIME=-5 TO LN_MAX BY INC, LN_MAX; MONTHS=EXP(LN_TIME); OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  # The span is log(hi) - lo, not the 5 + log(hi) that only coincides with it
  # because this job starts at -5, and the /99.9 denominator is read from the
  # job's own INC=, not assumed. The DO's `, LN_MAX` trailing element is SAS's
  # own extra list value -- `a TO b BY c, d` runs the loop and then takes
  # `d` -- so the emitted grid appends the exact bound, 180, after the loop.
  expect_equal(got$grid,
               quote(data.frame(
                 time = c(exp(-5 + seq(0, 99) * ((log(180) - -5) / 99.9)), 180)
               )))
})

test_that("a grid built by SET is untranslated, not guessed at", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; SET COHORT; ",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_null(got$grid)
  expect_true(any(grepl("grid", got$untranslated$reason)))
})

test_that("the default HAZPRED block emits both call and call_haz with se.fit", {
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=EX.HZD OUT=PREDICT; TIME MONTHS; );"
  )
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_equal(got$call[["type"]], "survival")
  expect_equal(got$call[["se.fit"]], TRUE)
  expect_equal(got$call_haz[["type"]], "hazard")
  expect_equal(got$call_haz[["se.fit"]], TRUE)
})

test_that("NOHAZ makes call_haz NULL", {
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=EX.HZD OUT=PREDICT NOHAZ; TIME MONTHS; );"
  )
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_equal(got$call[["type"]], "survival")
  expect_null(got$call_haz)
})

test_that("NOSURV makes call the hazard prediction", {
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=EX.HZD OUT=PREDICT NOSURV; TIME MONTHS; );"
  )
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_equal(got$call[["type"]], "hazard")
  expect_null(got$call_haz)
})

test_that("an explicit DO grid resolves DATA-step constants, incl. 1*DTY", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; DIGITAL=0;",
    "DTY=12/365.2425;",
    "DO MONTHS=1*DTY,2*DTY,24 TO 180 BY 12;",
    "OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_false(is.null(got$grid))

  dty <- 12 / 365.2425
  expected <- c(1 * dty, 2 * dty, seq(24, 180, by = 12))
  expect_equal(eval(got$grid)$time, expected)
})

test_that("a DO list referencing an unknown name still refuses and records", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; DO MONTHS=1*FOO,2*FOO; OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_null(got$grid)
  expect_true(any(grepl("grid", got$untranslated$reason)))
})

test_that("a constant defined in terms of an earlier constant resolves", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; A=10; B=A*2; DO MONTHS=1*B,2*B; OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_false(is.null(got$grid))
  expect_equal(eval(got$grid)$time, c(20, 40))
})

test_that("a DO list constant defined via a function call refuses, not evaluates", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; A=SQRT(4); DO MONTHS=1*A,2*A; OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_null(got$grid)
  expect_true(any(grepl("grid", got$untranslated$reason)))
})

test_that("NOSURV and NOHAZ together yield no predict() call at all", {
  # The ternary that picks `call` tests only want_surv, so without a guard
  # this degenerate input would silently produce a hazard predict() nobody
  # asked for. Both call and call_haz must be NULL, and the suppression must
  # be recorded, not dropped.
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=EX.HZD OUT=PREDICT NOSURV NOHAZ; TIME MONTHS; );"
  )
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_null(got$call)
  expect_null(got$call_haz)
  expect_true(any(grepl("NOSURV", got$untranslated$construct)))
})

test_that("the survival predict call sets conf.type = logit for SAS parity", {
  # SAS HAZPRED's survival confidence limits are logit-scale
  # (hzp_calc_srv_CL.c); predict.hazard() defaults to "log-log". Emitting the
  # default reproduces the job with silently different bounds.
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=E.H OUT=P; TIME MONTHS; );"
  )
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_equal(got$call[["conf.type"]], "logit")
  # Hazard CLs are log scale in both engines, so that call must NOT be steered.
  expect_null(got$call_haz[["conf.type"]])
})

test_that("conf.type is omitted when NOCL suppresses confidence limits", {
  txt <- .hzr_sas_normalise(
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=E.H OUT=P NOCL; TIME MONTHS; );"
  )
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_equal(got$call[["se.fit"]], FALSE)
  expect_null(got$call[["conf.type"]])
})

test_that("the log grid matches SAS's INC = (5 + LN_MAX)/99.9 step", {
  # hp.death.AVC.sas's own DO reads `DO LN_TIME=-5 TO LN_MAX BY INC,LN_MAX;`
  # -- the trailing `, LN_MAX` is SAS's own extra DO-list value, not
  # decoration: `a TO b BY c, d` runs the loop and then takes `d`. An
  # earlier version of this test omitted the trailing element and asserted
  # 100 points ending at 164.207 -- one point short of SAS's actual 101,
  # because it wasn't exercising the trailing element at all. Do not "fix"
  # this back to 100 points; that was the bug.
  txt <- .hzr_sas_normalise(paste(
    "DATA PGRID; MAX = 180; LN_MAX = LOG(MAX);",
    "INC = (5 + LN_MAX)/99.9;",
    "DO LN_TIME = -5 TO LN_MAX BY INC, LN_MAX; MONTHS = EXP(LN_TIME);",
    "OUTPUT; END; RUN;"
  ))
  cl <- .hzr_parse_grid(txt, "PGRID")
  grid <- eval(cl)
  expect_equal(nrow(grid), 101L)
  # The 100th (last loop) point is exp(-5 + 99*(log(180)+5)/99.9) = 164.207,
  # still correct -- it is simply no longer the last row.
  expect_equal(grid$time[100L], exp(-5 + 99 * (log(180) + 5) / 99.9),
               tolerance = 1e-8)
  # The 101st (trailing DO-list) point is the job's own MAX, exactly.
  expect_equal(grid$time[101L], 180)
})

test_that("a directly-assigned log bound is not mistaken for MAX", {
  # No separate MAX variable is ever assigned here -- LN_MAX is set directly,
  # as #153 describes. The unrelated LOG(2) is only there so the body still
  # trips the log-grid branch's own " DO "/"LOG("/EXP() detection heuristic.
  txt <- .hzr_sas_normalise(paste(
    "DATA PGRID; LN_MAX = 5.2; X = LOG(2); INC = (5 + LN_MAX)/99.9;",
    "DO LN_TIME = -5 TO LN_MAX BY INC; MONTHS = EXP(LN_TIME);",
    "OUTPUT; END; RUN;"
  ))
  cl <- .hzr_parse_grid(txt, "PGRID")
  grid <- if (is.null(cl)) NULL else eval(cl)
  # Must not silently produce a grid ending at t = 5.2 by reading LN_MAX=5.2
  # as MAX=5.2 and taking its log (#153). There is no MAX to log here, so
  # refusing to translate (NULL) is the correct outcome.
  expect_true(is.null(grid) || max(grid$time) > 100)
})

# The log-grid step is the job's own INC=, not a hardcoded /99.9. Three
# denominators appear across the public corpus -- /49.9, /99.9 and /999.9 --
# and the denominator sets both the step and the point count, because SAS's
# DO lo TO hi BY INC runs floor((hi - lo)/INC) + 1 times.

test_that("the log grid reads its step denominator from INC=, not /99.9", {
  # hp.dthip.PAIVS.time.sas and hmdeadp.sas both use /999.9. Read as /99.9
  # they produced 100 points on a step ten times too large -- every time
  # wrong, with `untranslated` empty and coverage reported as full. Both
  # jobs' DO also carries a `, LN_MAX` trailing element (SAS's own DO-list
  # extra value), which adds one more point beyond the loop: /999.9 lands
  # 1000 loop points + 1 trailing = 1001 rows.
  txt <- .hzr_sas_normalise(paste(
    "DATA PGRID; MAX = 84; LN_MAX = LOG(MAX);",
    "INC = (5 + LN_MAX)/999.9;",
    "DO LN_TIME = -5 TO LN_MAX BY INC, LN_MAX; MONTHS = EXP(LN_TIME);",
    "OUTPUT; END; RUN;"
  ))
  grid <- eval(.hzr_parse_grid(txt, "PGRID"))
  expect_equal(nrow(grid), 1001L)
  # The trailing point is the job's own MAX, exactly -- same value regardless
  # of the step, so it is the *loop's* last point (row 1000, not 1001) that
  # has to reflect the /999.9 step.
  expect_equal(grid$time[1000L], exp(-5 + 999 * (log(84) + 5) / 999.9),
               tolerance = 1e-8)
  expect_equal(grid$time[1001L], 84)

  # And the same job read at /99.9 is a different grid: if the loop's last
  # points agreed, the denominator would not be doing anything. The trailing
  # point is unaffected by the step -- both still land on 84 -- so the
  # comparison has to be on the loop's own last point, not on max().
  txt99 <- sub("/999.9", "/99.9", txt, fixed = TRUE)
  grid99 <- eval(.hzr_parse_grid(txt99, "PGRID"))
  expect_equal(nrow(grid99), 101L)
  expect_equal(grid99$time[101L], 84)
  expect_false(isTRUE(all.equal(grid99$time[100L], grid$time[1000L])))
})

test_that("a /49.9 step gives the 50 loop points SAS's DO loop lands, plus its trailing point", {
  # hs.dthar.TGA.setup.sas's form, with its bound named MAX rather than MAX0.
  # Its DO also carries the corpus's `, LN_MAX` trailing element, so the
  # emitted grid is 50 loop points + 1 trailing point = 51 rows.
  txt <- .hzr_sas_normalise(paste(
    "DATA PGRID; MAX = 96; LN_MAX = LOG(MAX);",
    "INC = (5 + LN_MAX)/49.9;",
    "DO LN_TIME = -5 TO LN_MAX BY INC, LN_MAX; MONTHS = EXP(LN_TIME);",
    "OUTPUT; END; RUN;"
  ))
  grid <- eval(.hzr_parse_grid(txt, "PGRID"))
  expect_equal(nrow(grid), 51L)
  # The 50th (last loop) point reflects the /49.9 step; the 51st (trailing
  # DO-list) point is the job's own MAX, exactly.
  expect_equal(grid$time[50L], exp(-5 + 49 * (log(96) + 5) / 49.9),
               tolerance = 1e-8)
  expect_equal(grid$time[51L], 96)
})

test_that("a log grid starting somewhere other than -5 uses its own span", {
  # Every corpus job starts at -5, which is the only reason a numerator of
  # 5 + log(hi) ever matched: the general span is log(hi) - lo. A job
  # starting at 0 must step by log(180)/99.9, not (5 + log(180))/99.9.
  txt <- .hzr_sas_normalise(paste(
    "DATA PGRID; MAX = 180; LN_MAX = LOG(MAX);",
    "INC = (LN_MAX - 0)/99.9;",
    "DO LN_TIME = 0 TO LN_MAX BY INC; MONTHS = EXP(LN_TIME);",
    "OUTPUT; END; RUN;"
  ))
  grid <- eval(.hzr_parse_grid(txt, "PGRID"))
  expect_equal(nrow(grid), 100L)
  expect_equal(min(grid$time), 1)
  expect_equal(max(grid$time), exp(99 * log(180) / 99.9), tolerance = 1e-8)
  # The old hardcoded (5 + log(180))/99.9 step would have landed here.
  expect_false(isTRUE(all.equal(max(grid$time),
                                exp(99 * (5 + log(180)) / 99.9))))
})

test_that("an INC= this cannot parse is refused, not stepped by a guess", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; MAX = 180; LN_MAX = LOG(MAX);",
    "INC = SQRT(LN_MAX)/99.9;",
    "DO LN_TIME = -5 TO LN_MAX BY INC; MONTHS = EXP(LN_TIME); OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_null(got$grid)
  expect_true(any(grepl("grid", got$untranslated$reason)))
})

test_that("an INC= numerator that is not the loop's span is refused", {
  # A numerator of (5 + LN_MAX) when the loop starts at 0 is not the span,
  # so the emitted (log(hi) - lo)/denominator would not be this job's step.
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; MAX = 180; LN_MAX = LOG(MAX);",
    "INC = (5 + LN_MAX)/99.9;",
    "DO LN_TIME = 0 TO LN_MAX BY INC; MONTHS = EXP(LN_TIME); OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_null(got$grid)
  expect_true(any(grepl("grid", got$untranslated$reason)))
})

test_that("a DO with no BY clause is refused, not given the /99.9 step", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PREDICT; MAX = 180; LN_MAX = LOG(MAX);",
    "DO LN_TIME = -5 TO LN_MAX; MONTHS = EXP(LN_TIME); OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; TIME MONTHS; );"
  ))
  got <- .hzr_parse_hazpred(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_null(got$grid)
  expect_true(any(grepl("grid", got$untranslated$reason)))
})
