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
