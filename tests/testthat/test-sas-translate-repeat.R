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
    list("bd_card", "`BD_CARD` is not a KEY=VALUE argument")
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
