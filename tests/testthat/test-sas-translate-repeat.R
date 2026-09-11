# test-sas-translate-repeat.R -- hzr_translate_sas() on jobs that call the
# SAS macro %repeat. Every test evaluates the emitted chunks; expected values
# are derived by hand from ~/Documents/macro.library/repeat.sas (see
# inst/dev/REPEAT-TRANSLATOR-PLAN.md, "The synthetic fixture"), never from
# hzr_repeated_events() itself.

translate_lines <- function(lines) {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(lines, f)
  suppressWarnings(hzr_translate_sas(f))
}

test_that("a job with %repeat but no PROC HAZARD or HAZPRED is still refused", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines("%repeat(in=bd, out=events);", f)
  expect_error(hzr_translate_sas(f), "no HAZARD or HAZPRED block found")
})

# Subject A: events at 1 and 3, a non-event at 2; B: no event, missing time;
# C: events at 1 and 6, where 6 is also end of follow-up.
bd_card <- function() {
  data.frame(
    CCFID = c("A", "A", "A", "B", "C", "C"),
    IV_EVENT = c(1, 2, 3, NA, 1, 6),
    IV_END = c(5, 5, 5, 4, 6, 6),
    CE_CARD = c(1, 0, 1, 0, 1, 1),
    stringsAsFactors = FALSE
  )
}

# Derived by hand from repeat.sas; see the plan's fixture table.
expected_events <- function() {
  data.frame(
    CCFID = c("A", "A", "A", "B", "C", "C"),
    IV_EVENT = c(1, 3, 5, 4, 1, 6),
    IV_END = c(5, 5, 5, 4, 6, 6),
    CE_CARD = c(1, 1, 0, 0, 1, 1),
    CN_CARD = c(0, 0, 1, 1, 0, 1),
    FIRST = c(1, 0, 0, 1, 1, 0),
    LAST = c(0, 1, 1, 1, 0, 1),
    EV_CARD = c(1, 1, 0, 0, 1, 1),
    IV_START = c(0, 1, 3, 0, 0, 1),
    EVENT_NO = c(1, 2, 2, 0, 1, 2),
    IV_SEG = c(1, 2, 2, 4, 1, 5),
    RENEWAL = c(1, 2, 3, 1, 1, 2),
    stringsAsFactors = FALSE
  )
}

parse_repeat <- function(args) {
  txt <- .hzr_sas_normalise(paste0("%repeat(", args, ");"))
  .hzr_parse_repeat(.hzr_sas_blocks(txt)[[1L]])
}

cardio_args <- "in=bd_card, out=events, id=ccfid, eventype=ce_card, iv_end=iv_end, rcensor=cn_card, event=ev_card"

test_that("the emitted %repeat call builds the macro's output under the job's names", {
  r <- parse_repeat(cardio_args)
  expect_equal(r$in_name, "BD_CARD")
  expect_equal(r$out_name, "EVENTS")
  expect_equal(nrow(r$untranslated), 0L)
  env <- new.env(parent = globalenv())
  env$BD_CARD <- bd_card()
  expect_no_warning(eval(r$call, env))
  expect_equal(dim(env$EVENTS), c(6L, 12L))
  expect_equal(env$EVENTS, expected_events())
})

test_that("an input column the macro would overwrite is dropped, with a warning", {
  r <- parse_repeat(cardio_args)
  env <- new.env(parent = globalenv())
  env$BD_CARD <- transform(bd_card(), EV_CARD = 99)
  expect_warning(eval(r$call, env), "EV_CARD")
  # Without the drop, names<- would leave two EV_CARD columns and $EV_CARD
  # would return the stale 99s.
  expect_equal(sum(names(env$EVENTS) == "EV_CARD"), 1L)
  expect_equal(env$EVENTS$EV_CARD, c(1, 1, 0, 0, 1, 1))
})

test_that("a NUMBER column in the input warns that SAS's event counts are not comparable", {
  r <- parse_repeat(cardio_args)
  env <- new.env(parent = globalenv())
  env$BD_CARD <- transform(bd_card(), NUMBER = 7)
  expect_warning(eval(r$call, env), "NUMBER")
  # R's result is the macro's intent: the column is neither dropped nor read.
  expect_equal(env$EVENTS$EVENT_NO, c(1, 2, 2, 0, 1, 2))
  expect_equal(env$EVENTS$NUMBER, rep(7, 6))
})

test_that("defaults are the macro's own, upper-cased", {
  r <- parse_repeat("")
  expect_equal(r$in_name, "BUILT")
  expect_equal(r$out_name, "EVENTS")
  env <- new.env(parent = globalenv())
  d <- bd_card()
  names(d) <- c("ID", "IV_EVENT", "IV_END", "EVENTYPE")
  env$BUILT <- d
  eval(r$call, env)
  expect_equal(
    names(env$EVENTS),
    c("ID", "IV_EVENT", "IV_END", "EVENTYPE", "RCENSOR", "FIRST", "LAST",
      "EVENT", "IV_START", "EVENT_NO", "IV_SEG", "RENEWAL")
  )
  expect_equal(env$EVENTS$IV_START, c(0, 1, 3, 0, 0, 1))
})

test_that("uncallable %repeat arguments are refused, each with its reason", {
  cases <- list(
    list("in=bd_card, eventype=rcensor", "RCENSOR is also an output"),
    list("in=bd_card, foo=bar", "`FOO=` is not a %repeat parameter"),
    list("in=&dsn", "`IN=&DSN` is not a plain SAS name"),
    list("event=x, rcensor=x", "two outputs are both named X"),
    list("bd_card", "`BD_CARD` is not a KEY=VALUE argument"),
    list("in=bd_card, event_no=number", "output NUMBER has the name of a %repeat loop counter"),
    list("in=bd_card, id=a, id=b", "`ID=` is given more than once"),
    list("in=bd_card, iv_start=lag_iv", "output LAG_IV has the name of a %repeat loop counter")
  )
  for (cs in cases) {
    r <- parse_repeat(cs[[1]])
    expect_null(r$in_name)
    expect_null(r$out_name)
    expect_identical(r$call[[1L]], as.name("stop"))
    expect_error(eval(r$call, new.env()), cs[[2]], fixed = TRUE)
    expect_equal(r$untranslated$construct, "%repeat")
    expect_length(r$untranslated$reason, 1L)
    expect_true(grepl(cs[[2]], r$untranslated$reason, fixed = TRUE))
    expect_equal(r$tokens_mapped, 0L)
  }
})

repeat_call <- paste0("%repeat(", cardio_args, ");")
hz_block <- c("%hazard( proc hazard data=events condition=14;",
              "lcensor iv_start; event ce_card; time iv_event; parms muc=0.5; );")

eval_upto <- function(job, env, last) {
  for (nm in names(job$calls)[seq_len(match(last, names(job$calls)))]) {
    eval(job$calls[[nm]], env)
  }
  env
}

test_that("a %repeat job emits guard, macro, status and fit, in that order", {
  job <- translate_lines(c(repeat_call, hz_block))
  expect_equal(names(job$calls), c("data", "repeated", "status", "fit"))
  # The guard is on the macro's input; EVENTS is built by a chunk, so it gets none.
  expect_error(eval(job$calls$data, new.env(parent = globalenv())), "Assign BD_CARD")
  env <- new.env(parent = globalenv())
  env$BD_CARD <- bd_card()
  eval_upto(job, env, "repeated")
  expect_equal(env$EVENTS, expected_events())
})

test_that("the translated fit reads the renamed IV_START: mu is 4 events over 15 units", {
  skip_on_cran()
  job <- translate_lines(c(repeat_call, hz_block))
  env <- new.env(parent = globalenv())
  env$BD_CARD <- bd_card()
  eval_upto(job, env, "fit")
  expect_true(env$fit$fit$converged)
  # 4/15 needs time_lower = IV_START; without it the exposure is sum(IV_EVENT) = 20.
  expect_equal(unname(exp(coef(env$fit))) / (4 / 15), 1, tolerance = 1e-4)
})

test_that("a non-default output name is what the job's fit reads", {
  job <- translate_lines(c(
    "%repeat(in=bd_card, id=ccfid, eventype=ce_card, event=ev_x);",
    "%hazard( proc hazard data=events; event ev_x; time iv_event; parms muc=0.5; );"
  ))
  env <- new.env(parent = globalenv())
  env$BD_CARD <- bd_card()
  eval_upto(job, env, "repeated")
  expect_false("event" %in% names(env$EVENTS))
  expect_equal(env$EVENTS$EV_X, c(1, 1, 0, 0, 1, 1))
})

test_that("a DATA step rewriting OUT= stops the document before the fit, quoted", {
  job <- translate_lines(c(
    repeat_call,
    "data events; set events; if iv_start ge iv_event then iv_event=iv_event+0.0001;",
    hz_block
  ))
  expect_equal(names(job$calls), c("data", "repeated", "rewrite", "status", "fit"))
  env <- new.env(parent = globalenv())
  env$BD_CARD <- bd_card()
  eval_upto(job, env, "repeated")
  expect_error(
    eval(job$calls$rewrite, env),
    "DATA EVENTS; SET EVENTS; IF IV_START GE IV_EVENT THEN IV_EVENT=IV_EVENT+0.0001;",
    fixed = TRUE
  )
  expect_true("EVENTS changed after %repeat" %in% job$untranslated$construct)
})

test_that("steps that do not rewrite OUT= emit no rewrite stop", {
  job <- translate_lines(c(
    repeat_call,
    "data other; set events; x=1;",
    "proc sort data=events; by ccfid;",
    hz_block
  ))
  expect_equal(names(job$calls), c("data", "repeated", "status", "fit"))
})

test_that("CREATE TABLE OUT is a rewrite too", {
  job <- translate_lines(c(
    # `select ccfid`, not `select *`: a `*` risks the normaliser's comment stripping.
    repeat_call, "proc sql; create table events as select ccfid from events; quit;", hz_block
  ))
  expect_equal(names(job$calls), c("data", "repeated", "rewrite", "status", "fit"))
})

test_that("one rewrite stops once, however many fits read OUT= after it", {
  job <- translate_lines(c(
    repeat_call, "data events; set events; x=1;", hz_block, hz_block
  ))
  expect_equal(sum(grepl("^rewrite", names(job$calls))), 1L)
  expect_equal(which(names(job$calls) == "rewrite"), 3L)
})

test_that("a refused %repeat still leaves a guard on the fit's DATA=", {
  job <- translate_lines(c("%repeat(in=bd_card, eventype=rcensor);", hz_block))
  expect_equal(names(job$calls), c("repeated", "data", "status", "fit"))
  expect_error(eval(job$calls[["repeated"]], new.env()), "RCENSOR is also an output", fixed = TRUE)
})

test_that("the rewrite scan fails closed: any step naming OUT= stops, bar a plain sort or a read", {
  stops <- c(
    "data work.events; set events; x=1; run;",
    "proc sort data=events nodupkey; by ccfid; run;",
    "proc sort data=bd_card out=events; by ccfid;",
    "proc sort data=events(where=(iv_seg>0)); by ccfid;",
    "proc sql; delete from events where iv_seg=0; quit;",
    "proc sql; create table work.events(drop=x) as select ccfid from events; quit;",
    "proc append base=events data=extra; run;",
    "proc datasets; modify events; rename a=b; quit;",
    "%vars(in=bd, out=events);",
    "proc sort data=events; by ccfid; where iv_seg > 0;",
    "proc sort data=events; by ccfid; %fix(data=events);",
    "data other; set events; x=1; %fix(data=events);",
    "data &dsn; set &dsn; x=1;",
    "proc sort data=&dsn nodupkey; by ccfid;"
  )
  for (s in stops) {
    job <- translate_lines(c(repeat_call, s, hz_block))
    expect_equal(names(job$calls), c("data", "repeated", "rewrite", "status", "fit"), info = s)
  }
  passes <- c(
    "data other; set events; x=1;",
    "data _null_; set events; put ccfid;",
    "proc sort data=events; by ccfid; run;",
    "proc sort data=work.events; by ccfid;",
    "proc sort data=bd_card; by ccfid; data x; set bd_card;"
  )
  for (s in passes) {
    job <- translate_lines(c(repeat_call, s, hz_block))
    expect_equal(names(job$calls), c("data", "repeated", "status", "fit"), info = s)
  }
})

test_that("a fit whose block encloses the %repeat call cannot render a fit", {
  # The fit's parse reads PROC HAZARD options from the block's first
  # statement, which here is the %repeat call, so DATA= is lost and the
  # rewrite scan never runs. That is loud, not silent: evaluation stops before
  # a fit exists. If a parser change ever keeps DATA= here, the first
  # assertion fails, and the rewrite scan must then handle an enclosing block.
  job <- translate_lines(c(
    "%run( %repeat(in=bd_card, id=ccfid, eventype=ce_card);",
    "data events; set events; x=1;",
    "proc hazard data=events; event ce_card; time iv_event; parms muc=0.5; );"
  ))
  expect_null(job$calls$fit[[3L]][["data"]])
  env <- new.env(parent = globalenv())
  env$BD_CARD <- bd_card()
  expect_error(for (nm in names(job$calls)) eval(job$calls[[nm]], env), "CE_CARD")
  expect_false(exists("fit", envir = env, inherits = FALSE))
})

test_that("a step between chained %repeat calls that changes the first OUT= stops", {
  first <- "%repeat(in=bd_card, out=ev1, id=ccfid, eventype=ce_card);"
  second <- "%repeat(in=ev1, out=events, id=ccfid, eventype=ce_card);"
  job <- translate_lines(c(first, "data ev1; set ev1; iv_event=iv_event+1;", second, hz_block))
  expect_equal(names(job$calls), c("data", "repeated", "rewrite", "repeated_2", "status", "fit"))
  expect_true("EV1 changed after %repeat" %in% job$untranslated$construct)
  # Without a step between them, the chain emits no stop.
  job <- translate_lines(c(first, second, hz_block))
  expect_equal(names(job$calls), c("data", "repeated", "repeated_2", "status", "fit"))
})

test_that("a WORK. libref on IN= and OUT= names the same dataset as the bare name", {
  job <- translate_lines(c(
    paste("%repeat(in=work.bd_card, out=work.events, id=ccfid, eventype=ce_card,",
          "iv_end=iv_end, rcensor=cn_card, event=ev_card);"),
    hz_block
  ))
  expect_equal(names(job$calls), c("data", "repeated", "status", "fit"))
  env <- new.env(parent = globalenv())
  env$BD_CARD <- bd_card()
  eval_upto(job, env, "repeated")
  expect_equal(env$EVENTS, expected_events())
})
