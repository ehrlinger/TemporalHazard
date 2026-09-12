test_that("block comments are stripped across lines", {
  src <- c("PROC HAZARD /* keep", "going */ DATA=A;")
  expect_equal(.hzr_sas_normalise(src), "PROC HAZARD DATA=A;")
})

test_that("statement comments are stripped, including commented-out PARMS", {
  # Commented-out PARMS lines are common in these jobs. Left in, they inflate
  # every count and can inject a second parameter set into the model.
  src <- c("PARMS MUE=0.2;", "*   PARMS MUE=0.9;", "EVENT DEAD;")
  expect_equal(.hzr_sas_normalise(src), "PARMS MUE=0.2; EVENT DEAD;")
})

test_that("an apostrophe inside a comment does not swallow later code", {
  src <- c("* patient's status ;", "EVENT DEAD;")
  expect_equal(.hzr_sas_normalise(src), "EVENT DEAD;")
})

test_that("inline comments after a semicolon are stripped", {
  # HAZARD's own lexer defines this: <STMT>\*[^;]*; in hazard_l.l
  src <- "TIME T; * a note ; EVENT D;"
  expect_equal(.hzr_sas_normalise(src), "TIME T; EVENT D;")
})

test_that("an apostrophe inside an inline comment does not swallow the following statement", {
  # HAZARD's lexer rule (<STMT>\*[^;]*;) is quote-agnostic: the comment ends
  # at the first literal `;`, apostrophes included. A quote-aware search here
  # would see an unbalanced quote in "patient's" and read on past the real
  # terminator, silently dropping "EVENT D;".
  src <- "TIME T; * patient's note ; EVENT D;"
  expect_equal(.hzr_sas_normalise(src), "TIME T; EVENT D;")
})

test_that("a block is bounded by balanced parens, not the first close", {
  txt <- .hzr_sas_normalise(
    "%HAZARD( PROC HAZARD DATA=A; PARMS MUE=EXP(1); ); DATA NEXT;"
  )
  b <- .hzr_sas_blocks(txt)
  expect_length(b, 1L)
  expect_equal(b[[1]]$proc, "HAZARD")
  expect_equal(b[[1]]$terminator, "paren")
  # The nested EXP( ... ) must not terminate the block early.
  expect_true(grepl("MUE=EXP(1)", b[[1]]$text, fixed = TRUE))
  expect_false(grepl("DATA NEXT", b[[1]]$text, fixed = TRUE))
})

test_that("HAZARD and HAZPRED blocks are both found, in order", {
  txt <- .hzr_sas_normalise(
    "%HAZARD( PROC HAZARD DATA=A; ); %HAZPRED( PROC HAZPRED DATA=P; );"
  )
  b <- .hzr_sas_blocks(txt)
  expect_equal(vapply(b, `[[`, "", "proc"), c("HAZARD", "HAZPRED"))
})

test_that("an unbalanced block is reported, never silently extended", {
  # Running to end of file is how comment prose reaches the token tables.
  txt <- .hzr_sas_normalise("%HAZARD( PROC HAZARD DATA=A; PARMS MUE=1;")
  b <- .hzr_sas_blocks(txt)
  expect_equal(b[[1]]$terminator, "none")
})

test_that("a bare parenthesised fragment is extracted", {
  txt <- .hzr_sas_normalise("(  PROC HAZARD DATA=A; EVENT D; TIME T; )")
  b <- .hzr_sas_blocks(txt)
  expect_length(b, 1L)
  expect_equal(b[[1]]$proc, "HAZARD")
  expect_equal(b[[1]]$terminator, "paren")
  expect_true(grepl("EVENT D;", b[[1]]$text, fixed = TRUE))
})

test_that("a PROC with no enclosing paren is returned, never dropped", {
  txt <- .hzr_sas_normalise("PROC HAZARD DATA=A; EVENT D; DATA NEXT;")
  b <- .hzr_sas_blocks(txt)
  expect_length(b, 1L)
  expect_equal(b[[1]]$terminator, "none")
  expect_false(grepl("DATA NEXT", b[[1]]$text, fixed = TRUE))
})

test_that("an unenclosed PROC with no following boundary is bounded and reported", {
  txt <- .hzr_sas_normalise("PROC HAZARD DATA=A; EVENT D; TIME T; RUN;")
  b <- .hzr_sas_blocks(txt)
  expect_length(b, 1L)
  expect_equal(b[[1]]$terminator, "none")
  expect_false(grepl("RUN;", b[[1]]$text, fixed = TRUE))
})

test_that("an unenclosed PROC with nothing following extends to end of text", {
  # Documented behaviour, not an accident: comments are already stripped by
  # .hzr_sas_normalise() before this function sees the text.
  txt <- .hzr_sas_normalise("PROC HAZARD DATA=A; EVENT D;")
  b <- .hzr_sas_blocks(txt)
  expect_length(b, 1L)
  expect_equal(b[[1]]$terminator, "none")
  expect_true(grepl("EVENT D;", b[[1]]$text, fixed = TRUE))
})

test_that("%REPEAT calls are blocks too, in file order, with offsets", {
  txt <- .hzr_sas_normalise(c(
    "%repeat(in=bd, out=events, id=ccfid);",
    "%hazard( proc hazard data=events; time t; event e; parms muc=0.1; );",
    "%hazpred( proc hazpred data=g inhaz=outest; time t; );"
  ))
  b <- .hzr_sas_blocks(txt)
  expect_equal(vapply(b, `[[`, "", "proc"), c("REPEAT", "HAZARD", "HAZPRED"))
  expect_equal(b[[1]]$text, "IN=BD, OUT=EVENTS, ID=CCFID")
  expect_equal(b[[1]]$terminator, "paren")
  expect_equal(substring(txt, b[[1]]$start, b[[1]]$end), "%REPEAT(IN=BD, OUT=EVENTS, ID=CCFID)")
  # A paren-bounded PROC block starts at its own opening paren.
  expect_equal(substr(txt, b[[2]]$start, b[[2]]$start), "(")
  expect_equal(substr(txt, b[[2]]$end, b[[2]]$end), ")")
  expect_lt(b[[1]]$end, b[[2]]$start)
  expect_lt(b[[2]]$end, b[[3]]$start)
})

test_that("a macro definition or a longer macro name is not a %repeat call", {
  txt <- .hzr_sas_normalise(c(
    "%macro repeat(in=built, out=events);",
    "%repeated(in=x);",
    "%hazard( proc hazard data=a; time t; event e; );"
  ))
  expect_equal(vapply(.hzr_sas_blocks(txt), `[[`, "", "proc"), "HAZARD")
})

test_that("an unenclosed PROC records where its bounded body ends", {
  txt <- .hzr_sas_normalise("PROC HAZARD DATA=A; EVENT D; DATA NEXT;")
  b <- .hzr_sas_blocks(txt)
  expect_equal(b[[1]]$start, 1L)
  expect_equal(substring(txt, b[[1]]$start, b[[1]]$end), "PROC HAZARD DATA=A; EVENT D; ")
})
