# SAS macro control flow decides whether a fit runs at all (#687).
#
# The translator reads a job's text in order, so it translated the fits the
# text shows rather than the fits SAS runs, and reported a clean translation:
#
# - a `%HAZARD(...)` call in a `* ... ;` statement comment vanished with the
#   comment, though SAS's macro processor runs a call there;
# - a `%HAZARD(...)` inside a `%MACRO ... %MEND` definition, or under an
#   open-code `%IF`/`%DO`, was emitted as a plain fit, though SAS runs it
#   once per call of the macro, or only when the condition holds.
#
# Each is now recorded in `$untranslated`; a block under macro control is
# still emitted, with a warning chunk above it, since the usual case is a
# macro called once.

macro_job <- function(...) {
  f <- tempfile(fileext = ".sas")
  on.exit(unlink(f))
  writeLines(c(...), f)
  # The translator warns once per job that has untranslated rows; kept, so a
  # case can assert the user is told.
  warned <- character()
  job <- withCallingHandlers(hzr_translate_sas(f), warning = function(w) {
    warned <<- c(warned, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  attr(job, "warned") <- warned
  job
}

macro_live <- paste(
  "%HAZARD(PROC HAZARD DATA=AVCS; TIME INT_DEAD; EVENT DEAD;",
  "PARMS MUE=0.3 THALF=0.2 NU=1 M=1 FIXTHALF FIXNU FIXM;);"
)

macro_fit_slots <- function(job) grep("^fit(_[0-9]+)?$", names(job$calls), value = TRUE)

macro_rows <- function(job, pattern) {
  u <- job$untranslated
  u[grepl(pattern, paste(u$construct, u$reason), fixed = TRUE), , drop = FALSE]
}

test_that("a live %HAZARD call alone translates cleanly", {
  job <- macro_job(macro_live)
  expect_length(macro_fit_slots(job), 1L)
  expect_equal(nrow(job$untranslated), 0L)
  expect_false(any(grepl("^refusal", names(job$calls))))
  expect_length(attr(job, "warned"), 0L)
})

test_that("a %HAZARD call in a * comment is recorded, not dropped", {
  commented <- paste("*", sub("MUE=0.3", "MUE=0.2", macro_live, fixed = TRUE))
  job <- macro_job(commented, macro_live)
  expect_length(macro_fit_slots(job), 1L)
  rows <- macro_rows(job, "comment")
  expect_equal(nrow(rows), 1L)
  expect_match(rows$construct, "%HAZARD", fixed = TRUE)
  expect_true(any(grepl("untranslated construct", attr(job, "warned"), fixed = TRUE)))
})

test_that("a %HAZARD inside a %MACRO definition is recorded and warned", {
  job <- macro_job("%MACRO RUNIT;", macro_live, "%MEND RUNIT;")
  expect_length(macro_fit_slots(job), 1L)
  expect_equal(nrow(macro_rows(job, "%MACRO")), 1L)
  warn <- job$calls[grep("^refusal", names(job$calls))]
  expect_length(warn, 1L)
  expect_match(deparse1(warn[[1L]]), "%MACRO", fixed = TRUE)
  expect_true(any(grepl("untranslated construct", attr(job, "warned"), fixed = TRUE)))
})

test_that("a %HAZARD under an open-code %IF %DO block is recorded and warned", {
  job <- macro_job("%IF &RUN = 1 %THEN %DO;", macro_live, "%END;")
  expect_equal(nrow(macro_rows(job, "%IF")), 1L)
  expect_length(grep("^refusal", names(job$calls)), 1L)
})

test_that("a %HAZARD as the %THEN of a one-statement %IF is recorded", {
  job <- macro_job(paste("%IF &RUN = 1 %THEN", macro_live))
  expect_equal(nrow(macro_rows(job, "%IF")), 1L)
})

test_that("a %HAZARD after a %MACRO definition ends is not flagged", {
  # The boundary: a call after %MEND, and after a closed %DO block, runs once
  # in open code like any other.
  job <- macro_job("%MACRO NOTE; %PUT HELLO; %MEND NOTE;",
                   "%IF &RUN = 1 %THEN %DO; %PUT X; %END;", macro_live)
  expect_length(macro_fit_slots(job), 1L)
  expect_equal(nrow(macro_rows(job, "%MACRO")), 0L)
  expect_equal(nrow(macro_rows(job, "%IF")), 0L)
  expect_false(any(grepl("^refusal", names(job$calls))))
})

test_that("%THEN%DO with no space opens a block too", {
  # SAS's macro tokenizer needs no blank before a `%`.
  job <- macro_job("%IF &RUN = 1 %THEN%DO;", macro_live, "%END;")
  expect_equal(nrow(macro_rows(job, "%IF around PROC HAZARD")), 1L)
})

macro_repeat <- paste0("%repeat(in=bd_card, out=events, id=ccfid, eventype=ce_card, ",
                       "iv_end=iv_end, rcensor=cn_card, event=ev_card);")
macro_repeat_fit <- c("%hazard( proc hazard data=events condition=14;",
                      "lcensor iv_start; event ce_card; time iv_event; parms muc=0.5; );")

test_that("a %REPEAT call inside a %MACRO definition is recorded and warned", {
  # If the macro is never called, SAS fits EVENTS as built elsewhere, while
  # the document builds it with %repeat: different data under the same fit.
  open_code <- macro_job(macro_repeat, macro_repeat_fit)
  expect_equal(nrow(macro_rows(open_code, "around %REPEAT")), 0L)

  job <- macro_job("%MACRO R;", macro_repeat, "%MEND R;", macro_repeat_fit)
  expect_equal(nrow(macro_rows(job, "%MACRO around %REPEAT call")), 1L)
  expect_equal(nrow(macro_rows(job, "around PROC HAZARD")), 0L)
  warn <- job$calls[grep("^refusal", names(job$calls))]
  expect_length(warn, 1L)
  expect_match(deparse1(warn[[1L]]), "%REPEAT call", fixed = TRUE)
  # The warning chunk comes before the chunk that builds the dataset.
  expect_lt(match(names(warn), names(job$calls)),
            grep("^repeated", names(job$calls))[[1L]])
})

test_that("a PROC HAZPRED under an open-code %DO is recorded and warned", {
  pred <- "%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.HZD OUT=P; TIME MONTHS; );"
  open_code <- macro_job(macro_live, pred)
  expect_equal(nrow(macro_rows(open_code, "around PROC HAZPRED")), 0L)

  job <- macro_job(macro_live, "%DO I = 1 %TO 3;", pred, "%END;")
  expect_equal(nrow(macro_rows(job, "%IF around PROC HAZPRED")), 1L)
  expect_equal(nrow(macro_rows(job, "around PROC HAZARD")), 0L)
  expect_length(grep("^refusal", names(job$calls)), 1L)
})

test_that("compact macro branches are flagged: a name ends at `%` or `(`", {
  # Each spelling is what SAS's tokenizer reads as the control keyword; each
  # was silent while the check required a blank after it.
  jobs <- list(
    else_branch = paste("%IF &RUN = 1 %THEN %PUT OK; %ELSE", macro_live, sep = ""),
    if_paren = paste0("%IF(&RUN=1)%THEN", macro_live),
    do_while = c("%DO%WHILE(&I < 3);", macro_live, "%END;"),
    then_do_while = c("%IF &RUN = 1 %THEN %DO%UNTIL(&I > 3);", macro_live, "%END;")
  )
  for (name in names(jobs)) {
    job <- do.call(macro_job, as.list(jobs[[name]]))
    expect_equal(nrow(macro_rows(job, "%IF around PROC HAZARD")), 1L, info = name)
    expect_length(grep("^refusal", names(job$calls)), 1L)
    expect_length(macro_fit_slots(job), 1L)
  }
})

test_that("a closed compact %DO block does not flag what follows it", {
  job <- macro_job("%DO%WHILE(&I < 3); %PUT X; %END;", macro_live)
  expect_equal(nrow(macro_rows(job, "%IF")), 0L)
  expect_false(any(grepl("^refusal", names(job$calls))))
})

test_that("a macro name that starts with a control keyword is not control", {
  # The other side of the boundary: these are calls, and the %ENDIT must not
  # close the open %DO.
  for (call in c("%IFFY;", "%ELSEWHERE;", "%DOIT;")) {
    job <- macro_job(call, macro_live)
    expect_equal(nrow(macro_rows(job, "around PROC HAZARD")), 0L, info = call)
  }
  # In the statement that holds the fit, which is the one the check reads.
  for (call in c("%IFFY(1)", "%ELSEWHERE(1)")) {
    job <- macro_job(paste(call, macro_live))
    expect_equal(nrow(macro_rows(job, "around PROC HAZARD")), 0L, info = call)
  }
  job <- macro_job("%DO;", "%ENDIT;", macro_live, "%END;")
  expect_equal(nrow(macro_rows(job, "%IF around PROC HAZARD")), 1L)
})

# A %MACRO the file never calls never runs (#708). Its blocks are still emitted,
# so the model stays visible, but in chunks the document does not evaluate.

macro_run_slots <- function(job) setdiff(macro_fit_slots(job), names(job$not_run))

test_that("a fit inside a macro the file never calls is not run", {
  job <- macro_job(macro_live, "%MACRO FINAL;", sub("0.3", "0.4", macro_live), "%MEND FINAL;")
  expect_length(macro_fit_slots(job), 2L)
  expect_identical(macro_run_slots(job), "fit")
  expect_identical(unname(job$not_run[["fit_2"]]), "FINAL")
  # Every chunk the block made is not run, its warning chunk included.
  made <- names(job$calls)[seq(match("fit", names(job$calls)) + 1L, length(job$calls))]
  expect_setequal(names(job$not_run), made)
  row <- macro_rows(job, "%MACRO around PROC HAZARD")
  expect_equal(nrow(row), 1L)
  expect_match(row$reason, "FINAL, which this file never calls", fixed = TRUE)
})

test_that("a called macro's fit still runs", {
  for (call in c("%FINAL;", "%IF &X = 1 %THEN %FINAL;", "DATA _NULL_; CALL EXECUTE('%FINAL'); RUN;")) {
    job <- macro_job(macro_live, "%MACRO FINAL;", sub("0.3", "0.4", macro_live),
                     "%MEND FINAL;", call)
    expect_length(job$not_run, 0L)
    expect_length(macro_run_slots(job), 2L)
  }
})

test_that("a call in a * comment counts as a call", {
  # SAS's macro processor runs a call inside a statement comment (#687).
  job <- macro_job(macro_live, "%MACRO FINAL;", sub("0.3", "0.4", macro_live),
                   "%MEND FINAL;", "* %FINAL;")
  expect_length(job$not_run, 0L)
})

test_that("a bare PROC HAZARD in an uncalled macro is not run", {
  job <- macro_job(macro_live, "%MACRO V;",
                   "PROC HAZARD DATA=AVCS; TIME INT_DEAD; EVENT DEAD; PARMS MUE=0.4 THALF=0.2 NU=1 M=1;",
                   "%MEND V;")
  expect_identical(macro_run_slots(job), "fit")
})

test_that("a macro that calls only itself is never called", {
  job <- macro_job(macro_live, "%MACRO LOOP;", sub("0.3", "0.4", macro_live),
                   "%LOOP;", "%MEND LOOP;")
  expect_identical(macro_run_slots(job), "fit")
})

test_that("a nested macro runs only when it and its enclosing macro are both called", {
  inner <- c("%MACRO OUTER;", "%MACRO INNER;", sub("0.3", "0.4", macro_live),
             "%MEND INNER;", "%INNER;", "%MEND OUTER;")
  expect_length(do.call(macro_job, as.list(c(macro_live, inner, "%OUTER;")))$not_run, 0L)
  job <- do.call(macro_job, as.list(c(macro_live, inner)))
  expect_identical(macro_run_slots(job), "fit")
  expect_identical(unname(job$not_run[["fit_2"]]), "OUTER")
})

test_that("a skipped block does not take a later fit's data guard", {
  # Both fits read AVCS. The guard belongs to the fit that runs.
  job <- macro_job("%MACRO OLD;", sub("0.3", "0.4", macro_live), "%MEND OLD;", macro_live)
  guards <- grep("^data", names(job$calls), value = TRUE)
  expect_length(setdiff(guards, names(job$not_run)), 1L)
  expect_lt(match(setdiff(guards, names(job$not_run)), names(job$calls)),
            match(macro_run_slots(job), names(job$calls)))
})

test_that("the rendered document does not evaluate a skipped chunk", {
  f <- tempfile(fileext = ".sas")
  writeLines(c(macro_live, "%MACRO FINAL;", sub("0.3", "0.4", macro_live), "%MEND FINAL;"), f)
  out <- tempfile()
  suppressWarnings(hzr_translate_sas(f, out_dir = out))
  qmd <- readLines(file.path(out, sub("[.]sas$", ".qmd", basename(f))))
  at <- grep("^#\\| label: fit(_2)?$", qmd)
  expect_length(at, 2L)
  expect_false(identical(qmd[at[[1L]] + 1L], "#| eval: false"))
  expect_identical(qmd[at[[2L]] + 1L], "#| eval: false")
  expect_match(qmd[at[[2L]] + 2L], "FINAL", fixed = TRUE)
})

test_that("a prediction after a skipped fit predicts from the fit that ran", {
  # Both fits write OUTHAZ=OUTEST; only the first one runs in SAS.
  live <- sub("DATA=AVCS;", "DATA=AVCS OUTHAZ=OUTEST;", macro_live, fixed = TRUE)
  pred <- "%HAZPRED( PROC HAZPRED DATA=P INHAZ=OUTEST OUT=P; TIME MONTHS; );"
  grid <- "DATA P; DO MONTHS = 1 TO 3; OUTPUT; END; RUN;"
  job <- macro_job(live, "%MACRO OLD;", sub("0.3", "0.4", live), "%MEND OLD;", grid, pred)
  expect_identical(macro_run_slots(job), "fit")
  expect_identical(deparse1(job$calls$pred[[1L]]), "predict")
  expect_identical(deparse1(job$calls$pred[[2L]]), "fit")
})

test_that("a skipped prediction does not replace the job's grid", {
  live <- sub("DATA=AVCS;", "DATA=AVCS OUTHAZ=OUTEST;", macro_live, fixed = TRUE)
  pred <- "%HAZPRED( PROC HAZPRED DATA=P INHAZ=OUTEST OUT=P; TIME MONTHS; );"
  ran <- c(live, "DATA P; DO MONTHS = 1 TO 3; OUTPUT; END; RUN;", pred)
  old <- c("%MACRO OLD;", "DATA P; DO MONTHS = 10 TO 30; OUTPUT; END; RUN;", pred, "%MEND OLD;")
  alone <- do.call(macro_job, as.list(ran))
  job <- do.call(macro_job, as.list(c(ran, old)))
  expect_false(is.null(alone$grid))
  expect_true(any(startsWith(names(job$not_run), "pred")))
  expect_identical(job$grid, alone$grid)
})

test_that("a skipped prediction's INHAZ neither stops the document nor loads its fit", {
  pred <- "%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.HZD OUT=P; TIME MONTHS; );"
  grid <- "DATA P; DO MONTHS = 1 TO 3; OUTPUT; END; RUN;"
  # Unresolved: only the skipped block reads EX.HZD.
  job <- macro_job(macro_live, "%MACRO OLD;", grid, pred, "%MEND OLD;")
  expect_null(job$inhaz)
  expect_true(job$inhaz_resolved)
  # Resolved through librefs: the live prediction loads its own fit.
  f <- tempfile(fileext = ".sas")
  writeLines(c(macro_live, "%MACRO OLD;", grid, pred, "%MEND OLD;", grid, pred), f)
  job <- suppressWarnings(hzr_translate_sas(f, librefs = c(EX = "estimates")))
  live_pred <- setdiff(grep("^pred", names(job$calls), value = TRUE), names(job$not_run))
  loaded <- deparse1(job$calls[[live_pred[[1L]]]][[2L]])
  expect_false(loaded %in% names(job$not_run))
  expect_true(loaded %in% names(job$calls))
})

test_that("a skipped %repeat does not stop a later fit on its dataset", {
  # SAS fits EVENTS as the job's own steps built it; the %repeat never ran, so a
  # step that changes EVENTS is not a rewrite of the macro's output.
  job <- macro_job("%MACRO R;", macro_repeat, "%MEND R;",
                   "DATA EVENTS; SET EVENTS; X = 1; RUN;", macro_repeat_fit)
  expect_true(any(startsWith(names(job$not_run), "repeated")))
  live <- setdiff(names(job$calls), names(job$not_run))
  expect_false(any(startsWith(live, "rewrite")))
  expect_true(any(startsWith(live, "data")))
})

test_that("a file whose definitions do not close skips nothing", {
  # A call with no semicolon before %MEND, or an argument list that never
  # closes, leaves the definition open to the end of the file, so the real
  # call reads as inside the macro (#708 review). Both macros are called.
  for (body_end in c("%CLEANUP", "%NOTE(fit done :-( );")) {
    job <- macro_job("%MACRO RUNIT;", macro_live, body_end, "%MEND RUNIT;", "%RUNIT;")
    expect_length(job$not_run, 0L)
    expect_equal(nrow(macro_rows(job, "%MACRO around PROC HAZARD")), 1L)
  }
})

test_that("an unresolved prediction after a skipped first fit names the live fit", {
  pred <- "%HAZPRED( PROC HAZPRED DATA=P INHAZ=OUTEST OUT=P; TIME MONTHS; );"
  grid <- "DATA P; DO MONTHS = 1 TO 3; OUTPUT; END; RUN;"
  job <- macro_job("%MACRO OLD;", sub("0.3", "0.4", macro_live), "%MEND OLD;",
                   macro_live, grid, pred)
  expect_identical(macro_run_slots(job), "fit_2")
  expect_identical(deparse1(job$calls$pred[[2L]]), "fit_2")
})

test_that("print() names the chunks SAS never ran", {
  job <- macro_job(macro_live, "%MACRO FINAL;", sub("0.3", "0.4", macro_live), "%MEND FINAL;")
  out <- capture.output(print(job))
  expect_true(any(grepl("not run *: .*fit_2", out)))
  expect_false(any(grepl("not run", capture.output(print(macro_job(macro_live))))))
})

test_that("a rewrite inside a skipped macro does not stop a later live fit", {
  # A live %repeat builds EVENTS; a never-called macro changes it and fits it;
  # the live fit after it reads EVENTS as the %repeat left it (#709 review).
  job <- macro_job(macro_repeat,
                   "%MACRO OLD;", "DATA EVENTS; SET EVENTS; X = 1; RUN;",
                   macro_repeat_fit, "%MEND OLD;",
                   macro_repeat_fit)
  live <- setdiff(names(job$calls), names(job$not_run))
  expect_false(any(startsWith(live, "rewrite")))
  expect_length(setdiff(macro_fit_slots(job), names(job$not_run)), 1L)
  # The same step outside any macro still stops the fit.
  open_code <- macro_job(macro_repeat, "DATA EVENTS; SET EVENTS; X = 1; RUN;",
                         macro_repeat_fit)
  expect_true(any(startsWith(names(open_code$calls), "rewrite")))
})

test_that("a rewrite inside a skipped macro does not stop a chained %repeat", {
  second <- sub("in=bd_card, out=events", "in=events, out=ev2", macro_repeat, fixed = TRUE)
  fit2 <- sub("data=events", "data=ev2", macro_repeat_fit, fixed = TRUE)
  job <- macro_job(macro_repeat,
                   "%MACRO OLD;", "DATA EVENTS; SET EVENTS; X = 1; RUN;", "%MEND OLD;",
                   second, fit2)
  expect_false(any(startsWith(names(job$calls), "rewrite")))
  open_code <- macro_job(macro_repeat, "DATA EVENTS; SET EVENTS; X = 1; RUN;", second, fit2)
  expect_true(any(startsWith(names(open_code$calls), "rewrite")))
})
