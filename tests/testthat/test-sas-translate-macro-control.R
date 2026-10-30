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
