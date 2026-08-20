# Guards on .hzr_sas_grammar, the keyword table generated from HAZARD's own
# lex sources by data-raw/hazard-grammar.R.

test_that("the grammar table is populated", {
  # A silently empty or truncated table would make every downstream coverage
  # figure a confident zero. Assert the size, not merely that it exists.
  expect_s3_class(.hzr_sas_grammar, "data.frame")
  expect_gt(nrow(.hzr_sas_grammar), 100L)
  expect_setequal(unique(.hzr_sas_grammar$proc), c("HAZARD", "HAZPRED"))
})

test_that("every keyword resolves to a parser token and a status", {
  g <- .hzr_sas_grammar
  expect_true(all(nzchar(g$keyword)))
  expect_true(all(nzchar(g$token)))
  expect_setequal(
    unique(g$implementation_status),
    c("mapped", "unmapped", "no_r_equivalent")
  )
  expect_false(anyNA(g$implementation_status))
})

test_that("a keyword is unique within its proc and lexer context", {
  g <- .hzr_sas_grammar
  key <- paste(g$proc, g$context, g$keyword, sep = "|")
  expect_equal(anyDuplicated(key), 0L)
})

test_that("every alias run ends in a canonical spelling of its token", {
  # Aliases fall through to the next rule's action in lex. If an alias pointed
  # at a token that never appears as a canonical keyword, the run was misparsed.
  g <- .hzr_sas_grammar
  canonical <- g$token[!g$is_alias]
  expect_true(all(g$token[g$is_alias] %in% canonical))
})

test_that("context disambiguates spellings that collide", {
  # M is the early-phase shape parameter in PARM context and MOVE in PHOP/STEP.
  # NOS is NOPRINTS in STEP context and NOSURV in HAZPRED's HZPP context.
  # A context-free lookup would translate both wrongly, so assert the collision
  # is real and that context separates it.
  g <- .hzr_sas_grammar

  # The general guarantee: within one proc and one lexer context, a spelling
  # resolves to exactly one token. This is what makes a context-aware lookup
  # total and a context-free one unsound.
  per_context <- tapply(
    g$token, paste(g$proc, g$context, g$keyword), function(x) length(unique(x))
  )
  expect_true(all(per_context == 1L))

  # And the collisions themselves are real, so the guarantee is load-bearing
  # rather than vacuous.
  m <- g[g$keyword == "M", ]
  expect_setequal(m$token, c("M", "MOVE"))
  expect_setequal(m$context[m$token == "M"], "PARM")

  nos <- g[g$keyword == "NOS", ]
  expect_setequal(nos$token, c("NOPRINTS", "NOSURV"))
})

test_that("MI is an alias for MAXITER, not a distinct option", {
  g <- .hzr_sas_grammar
  expect_equal(g$token[g$keyword == "MI"], "MAXITER")
  expect_true(all(g$is_alias[g$keyword == "MI"]))
})
