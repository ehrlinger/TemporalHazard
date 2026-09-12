# `%repeat` → `hzr_repeated_events()` in `hzr_translate_sas()`: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A SAS job that calls `%repeat(...)` is translated into an emitted chunk that calls
`hzr_repeated_events()` and renames the function's outputs to the job's names, so that the
job's `PROC HAZARD` fit reads the right columns.

**Architecture:** `.hzr_sas_blocks()` gains a third block kind, `REPEAT`, and records every
block's `start`/`end` offset. That keeps blocks in file order. `.hzr_parse_repeat()` turns a
`REPEAT` block into one emitted call, or into a `stop()` when the call cannot be translated.
`hzr_translate_sas()` dispatches the new kind: it guards `IN=`, marks `OUT=` as built, and
before each fit that reads `OUT=` it scans the text since the macro for DATA steps that
rewrite `OUT=`, emitting a quoted `stop()` for each.

**Tech Stack:** base R only; testthat 3e; withr (already used in tests); haven (Suggests,
gated test only).

**Spec:** `inst/dev/SAS-JOB-TRANSLATOR-DESIGN.md` §5.5. Branch `feat/translate-repeat`.

## Global Constraints

- No new dependency. `Imports:` stays `survival`; `haven` and `withr` are already Suggests.
- `stats::`-prefix any stats generic in `R/`. No `browser()`, bare `print()`, `cat()` or
  `library()` in `R/`.
- Version stays **1.2.11**. No bump. The NEWS bullet goes under the existing
  `# TemporalHazard 1.2.11` heading (several PRs share one patch version; that is deliberate).
- `man/` and `NAMESPACE` are generated. Run `devtools::document()`; never hand-edit them.
- Every emitted name is upper case (the translator runs `.hzr_sas_normalise()`, which calls
  `toupper()`).
- The chunk that calls the function is labelled **`repeated`**, never `repeat`: `repeat` is an
  R reserved word, so `job$calls$repeat` is a syntax error (verified 2026-09-10).
- Tests must be able to fail. Execute the emitted calls with `eval()`. Compare against
  vectors derived by hand from the macro, never from the function. Assert the shape before
  the values. For small numbers compare a ratio to 1, not a difference (waldo goes absolute
  when `expected < tolerance`).
- Anything that fits a model gets `skip_on_cran()`.
- PHI: nothing from `/Volumes/qhsstudies` enters the repository. The gated test reads the
  volume at run time and asserts only aggregates.
- Gates, in this order: `devtools::document()`, `lintr::lint_package()` (0 lints),
  `devtools::test()` (0 failures), `spelling::spell_check_package(use_wordlist = TRUE)`, then
  once per PR `R CMD check --as-cran` **with** the manual, built from a `git archive` of a
  **committed** tree under the session scratchpad.
- Never push to `main`. Branch, PR, stop; the maintainer merges.

## The synthetic fixture (used by Tasks 2–3)

Four columns, three subjects. The expected output below was derived **by hand** from
`repeat.sas`, stage by stage, and then confirmed to match `hzr_repeated_events()` on
2026-09-10:

- Subject A has events at 1 and 3 and a non-event row at 2, which stage 2 drops. Its
  follow-up ends at 5, so stage 4 appends a censored row at 5.
- Subject B has no event and a missing time. Stage 3 pads the row to `IV_END` = 4 with
  `rcensor` 1.
- Subject C has events at 1 and 6, with 6 equal to `IV_END`. Stage 4 appends a row at 6.
  Stage 6 sets `rcensor` to 1 on the event row at 6, so it is both an event and censored.
  Stage 7 drops the appended row, because it has zero length.

| CCFID | IV_EVENT | IV_END | CE_CARD | CN_CARD | FIRST | LAST | EV_CARD | IV_START | EVENT_NO | IV_SEG | RENEWAL |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 1 | 5 | 1 | 0 | 1 | 0 | 1 | 0 | 1 | 1 | 1 |
| A | 3 | 5 | 1 | 0 | 0 | 1 | 1 | 1 | 2 | 2 | 2 |
| A | 5 | 5 | 0 | 1 | 0 | 1 | 0 | 3 | 2 | 2 | 3 |
| B | 4 | 4 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 4 | 1 |
| C | 1 | 6 | 1 | 0 | 1 | 0 | 1 | 0 | 1 | 1 | 1 |
| C | 6 | 6 | 1 | 1 | 0 | 1 | 1 | 1 | 2 | 5 | 2 |

A constant-hazard fit of this frame with `LCENSOR IV_START` has a closed-form MLE:
4 events over 15 units of exposure (the sum of `IV_SEG`), so `mu = 4/15`. The prototype gave
`exp(log_mu) = 0.266666555`. If `time_lower` were ignored, the exposure would be the sum of
`IV_EVENT`, 20, giving 4/20. That is why the end-to-end test fails if the rename is broken.

---

### Task 1: `.hzr_sas_blocks()` recognises `%REPEAT(` and records offsets

**Files:**
- Modify: `R/sas-lex.R` (roxygen for `.hzr_sas_blocks()` and the function itself, lines 166–259)
- Modify: `R/translate-sas.R:157-161` (refuse a job with only `REPEAT` blocks) and the top of the `for (b in blocks)` loop at line 175 (skip `REPEAT` for now; Task 3 replaces this)
- Modify: `tests/testthat/test-sas-lex.R` (append)
- Create: `tests/testthat/test-sas-translate-repeat.R`

**Interfaces:**
- Produces: `.hzr_sas_blocks(txt)` returns a list of `list(proc, text, terminator, start, end)`.
  `proc` is one of `"HAZARD"`, `"HAZPRED"` or `"REPEAT"`. For `REPEAT`, `text` is the
  argument list inside the parentheses, `start` is the offset of `%`, and `end` is the offset
  of the closing `)`. For `HAZARD`/`HAZPRED`, `start` is the opening `(`, or the `P` of `PROC`
  when there is no enclosing paren, and `end` is the closing `)`, or the last character of
  the bounded body.
- Produces: `.hzr_sas_close_paren(txt, open_at)` returns the integer offset of the `)` that
  balances the `(` at `open_at`, or `NA_integer_`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-sas-lex.R`:

```r
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
```

Create `tests/testthat/test-sas-translate-repeat.R`:

```r
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
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "sas-lex|sas-translate-repeat")'`
Expected: FAIL. The first new lex test gets `c("HAZARD", "HAZPRED")` where it expects three
kinds; the offset tests fail on `NULL` `start`; the refusal test fails because translation
proceeds past the check.

- [ ] **Step 3: Implement**

In `R/sas-lex.R`, add the helper directly above `.hzr_sas_blocks()`'s roxygen block:

```r
#' Offset of the `)` that balances the `(` at `open_at`, or NA if none does.
#' @noRd
.hzr_sas_close_paren <- function(txt, open_at) {
  depth <- 0L
  for (i in seq(open_at, nchar(txt))) {
    ch <- substr(txt, i, i)
    if (ch == "(") depth <- depth + 1L
    if (ch == ")") {
      depth <- depth - 1L
      if (depth == 0L) return(i)
    }
  }
  NA_integer_
}
```

Change the first line of `.hzr_sas_blocks()`'s roxygen title and add a paragraph before
`#' @noRd`:

```r
#' Extract PROC HAZARD / PROC HAZPRED blocks and %repeat calls from normalised source.
```

```r
#'
#' A `%REPEAT(` call is returned as a third kind, `proc = "REPEAT"`, whose text
#' is its argument list. It is an ordinary macro call, so its group is the one
#' that opens right after the name and no backwards search is needed. Every
#' block also carries `start` and `end`, offsets into `txt`, so the caller can
#' read the source text between two blocks.
```

Replace the body of `.hzr_sas_blocks()` with:

```r
.hzr_sas_blocks <- function(txt) {
  out <- list()
  procs <- gregexpr("PROC HAZ(ARD|PRED) |%REPEAT *[(]", txt)[[1L]]
  if (procs[1L] == -1L) return(out)
  proc_lens <- attr(procs, "match.length")

  for (k in seq_along(procs)) {
    proc_at <- procs[k]

    if (identical(substring(txt, proc_at, proc_at + 6L), "%REPEAT")) {
      open_at <- proc_at + proc_lens[k] - 1L
      close_at <- .hzr_sas_close_paren(txt, open_at)
      end_at <- if (is.na(close_at)) nchar(txt) else close_at
      body_end <- if (is.na(close_at)) end_at else close_at - 1L
      out[[length(out) + 1L]] <- list(
        proc = "REPEAT", text = trimws(substring(txt, open_at + 1L, body_end)),
        terminator = if (is.na(close_at)) "none" else "paren",
        start = proc_at, end = end_at
      )
      next
    }

    proc <- if (identical(substring(txt, proc_at, proc_at + 10L), "PROC HAZARD")) {
      "HAZARD"
    } else {
      "HAZPRED"
    }

    # Scan backwards from the PROC keyword for the nearest unmatched `(` --
    # the paren that opens the group containing this PROC statement.
    open_at <- NA_integer_
    balance <- 0L
    if (proc_at > 1L) {
      for (i in seq(proc_at - 1L, 1L)) {
        ch <- substr(txt, i, i)
        if (ch == ")") {
          balance <- balance + 1L
        } else if (ch == "(") {
          if (balance == 0L) {
            open_at <- i
            break
          }
          balance <- balance - 1L
        }
      }
    }

    if (is.na(open_at)) {
      # No enclosing paren anywhere before this PROC. Bound the text at the
      # next PROC / DATA / RUN; boundary, never dropping the block. If none
      # of those follow either, the block genuinely extends to the end of
      # txt -- that is safe here because .hzr_sas_blocks() only ever sees
      # output from .hzr_sas_normalise(), which has already stripped all
      # comments, so the end-of-file-prose hazard this bounding rule exists
      # to prevent cannot arise at this point.
      search_from <- proc_at + proc_lens[k]
      rest <- substring(txt, search_from)
      b <- regexpr("PROC |DATA |RUN;", rest)
      end_at <- if (b == -1L) nchar(txt) else search_from + b - 2L
      body <- substring(txt, proc_at, end_at)
      term <- "none"
      start_at <- proc_at
    } else {
      close_at <- .hzr_sas_close_paren(txt, open_at)
      term <- if (is.na(close_at)) "none" else "paren"
      body <- substring(txt, open_at + 1L,
                        if (is.na(close_at)) nchar(txt) else close_at - 1L)
      start_at <- open_at
      end_at <- if (is.na(close_at)) nchar(txt) else close_at
    }

    out[[length(out) + 1L]] <- list(proc = proc, text = trimws(body),
                                    terminator = term,
                                    start = start_at, end = end_at)
  }
  out
}
```

In `R/translate-sas.R`, replace

```r
  blocks <- .hzr_sas_blocks(txt)
  if (!length(blocks)) {
    stop("no HAZARD or HAZPRED block found in ", path, call. = FALSE)
  }
```

with

```r
  blocks <- .hzr_sas_blocks(txt)
  # A %repeat call alone is a dataset-building job (the tp.bd.* templates),
  # not a HAZARD job, and there is nothing to fit.
  if (!any(vapply(blocks, function(b) b$proc != "REPEAT", logical(1L)))) {
    stop("no HAZARD or HAZPRED block found in ", path, call. = FALSE)
  }
```

and make the first statement inside `for (b in blocks) {`:

```r
    if (identical(b$proc, "REPEAT")) next # translated in Task 3
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `Rscript -e 'devtools::test(filter = "sas-lex|sas-translate-repeat|translate-sas|sas-parse")'`
Expected: 0 failures. The existing lex, parse and translate tests still pass: the new fields
are additive and their fixtures contain no `%REPEAT`.

- [ ] **Step 5: Commit**

```bash
git add R/sas-lex.R R/translate-sas.R tests/testthat/test-sas-lex.R tests/testthat/test-sas-translate-repeat.R
git commit -m "feat: recognise %repeat calls as SAS job blocks, with source offsets

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: `.hzr_parse_repeat()` builds the call, or refuses

**Files:**
- Modify: `R/sas-parse-job.R` (append at the end of the file)
- Modify: `tests/testthat/test-sas-translate-repeat.R` (append)

**Interfaces:**
- Consumes: the `REPEAT` block from Task 1, `list(proc = "REPEAT", text = <chr>, ...)`.
- Produces: `.hzr_parse_repeat(block)` returns `list(call, untranslated, tokens_seen,
  tokens_mapped, in_name, out_name)`.
  - `call` is `<OUT> <- local({...})`, or `stop("<msg>")` when refused.
  - `untranslated` is a `.hzr_untranslated_frame()` with construct `"%repeat"`, one row per
    problem, and no rows when accepted.
  - `in_name` and `out_name` are `character(1)` (upper case, possibly `LIB.MEMBER`) when
    accepted, and `NULL` when refused.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-sas-translate-repeat.R`:

```r
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
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "sas-translate-repeat")'`
Expected: FAIL with `could not find function ".hzr_parse_repeat"`.

- [ ] **Step 3: Implement**

Append to `R/sas-parse-job.R`:

```r
# ---------------------------------------------------------------------------
# %repeat -> hzr_repeated_events()
# ---------------------------------------------------------------------------

# %repeat's keyword parameters and defaults, as the macro declares them
# (~/Documents/macro.library/repeat.sas), upper-cased as
# .hzr_sas_normalise() leaves every name.
.hzr_repeat_defaults <- c(
  IN = "BUILT", OUT = "EVENTS", EVENTYPE = "EVENTYPE", IV_EVENT = "IV_EVENT",
  IV_END = "IV_END", ID = "ID", EVENT = "EVENT", EVENT_NO = "EVENT_NO",
  RCENSOR = "RCENSOR", IV_START = "IV_START", IV_SEG = "IV_SEG", RENEWAL = "RENEWAL"
)

# hzr_repeated_events()'s fixed output names, each mapped to the macro
# parameter that renames it. first and last have no parameter: the macro
# always writes them under those names.
.hzr_repeat_outputs <- c(
  event = "EVENT", event_no = "EVENT_NO", rcensor = "RCENSOR",
  iv_start = "IV_START", iv_seg = "IV_SEG", renewal = "RENEWAL"
)

#' Translate one `%repeat(...)` call into an `hzr_repeated_events()` chunk.
#'
#' The emitted chunk renames the function's lower-case outputs to the names the
#' job uses, upper-cased like every name this translator emits; without that the
#' fit reads a column that does not exist. Before the call it drops any input
#' column named like one of those outputs. SAS overwrites such a column, and the
#' macro assigns each output on every row before reading it, so dropping is
#' exactly what SAS does. Left in place, the rename would produce two columns of
#' one name, and `$` would return the stale one.
#'
#' A call this cannot express is refused rather than guessed at: the chunk is a
#' `stop()`, and each problem is an `$untranslated` row. That covers an unknown
#' keyword, a positional argument, a value that is not a plain SAS name (a
#' `&macro` reference, say), an input column that is also an output
#' (`EVENTYPE=RCENSOR`: SAS zeroes that indicator before reading it, so the job's
#' own answer is not the model it describes), and two outputs with one name. The
#' body is split on every comma; a value holding parentheses, where a comma could
#' be nested, is refused as not a plain name either way.
#' @noRd
.hzr_parse_repeat <- function(block) {
  parts <- if (nzchar(block$text)) trimws(strsplit(block$text, ",", fixed = TRUE)[[1L]]) else character(0)
  args <- .hzr_repeat_defaults
  problems <- character(0)

  for (p in parts) {
    kv <- regmatches(p, regexec("^([A-Z_][A-Z0-9_]*) ?= ?(.*)$", p))[[1L]]
    if (!length(kv)) {
      problems <- c(problems, sprintf("`%s` is not a KEY=VALUE argument", p))
      next
    }
    key <- kv[2L]
    val <- trimws(kv[3L])
    if (!key %in% names(args)) {
      problems <- c(problems, sprintf("`%s=` is not a %%repeat parameter", key))
      next
    }
    name_re <- if (key %in% c("IN", "OUT")) {
      "^[A-Z_][A-Z0-9_]*([.][A-Z_][A-Z0-9_]*)?$"
    } else {
      "^[A-Z_][A-Z0-9_]*$"
    }
    if (!grepl(name_re, val)) {
      problems <- c(problems, sprintf("`%s=%s` is not a plain SAS name", key, val))
      next
    }
    args[[key]] <- val
  }

  inputs <- unname(args[c("ID", "EVENTYPE", "IV_EVENT", "IV_END")])
  targets <- c(unname(args[.hzr_repeat_outputs]), "FIRST", "LAST")
  for (hit in intersect(inputs, targets)) {
    problems <- c(problems, sprintf(paste(
      "input column %s is also an output %%repeat writes; SAS overwrites it before reading it,",
      "so the job's own result is not the model it describes"
    ), hit))
  }
  dup <- unique(targets[duplicated(targets)])
  if (length(dup)) {
    problems <- c(problems, sprintf("two outputs are both named %s", paste(dup, collapse = ", ")))
  }

  if (length(problems)) {
    msg <- paste0("hzr_translate_sas() did not translate this job's %repeat call: ",
                  paste(problems, collapse = "; "), ".")
    return(list(
      call = bquote(stop(.(msg))),
      untranslated = .hzr_untranslated_frame(rep(NA_integer_, length(problems)),
                                             rep("%repeat", length(problems)), problems),
      tokens_seen = length(parts), tokens_mapped = 0L, in_name = NULL, out_name = NULL
    ))
  }

  in_name <- args[["IN"]]
  out_name <- args[["OUT"]]
  from <- c(names(.hzr_repeat_outputs), "first", "last")
  overwrite_msg <- paste0(" in ", in_name, ": %repeat writes columns of these names, so they were dropped ",
                          "before the call, as SAS overwrites them.")
  loop_msg <- paste0(" in ", in_name, ": SAS's %repeat reads a LAG_IV or NUMBER column in place of its ",
                     "own loop counters, so for this job the SAS listing's segment starts and event ",
                     "counts are not comparable with these.")
  call <- bquote(.(as.name(out_name)) <- local({
    d <- .(as.name(in_name))
    drop <- intersect(names(d), .(targets))
    if (length(drop)) {
      warning(paste(drop, collapse = ", "), .(overwrite_msg), call. = FALSE)
      d <- d[setdiff(names(d), drop)]
    }
    loop <- intersect(names(d), c("LAG_IV", "NUMBER"))
    if (length(loop)) warning(paste(loop, collapse = ", "), .(loop_msg), call. = FALSE)
    out <- hzr_repeated_events(d, id = .(args[["ID"]]), time = .(args[["IV_EVENT"]]),
                               followup = .(args[["IV_END"]]), indicator = .(args[["EVENTYPE"]]))
    names(out)[match(.(from), names(out))] <- .(targets)
    out
  }))

  list(call = call, untranslated = .hzr_untranslated_frame(),
       tokens_seen = length(parts), tokens_mapped = length(parts),
       in_name = in_name, out_name = out_name)
}
```

`targets` is in the same order as `from`: the six outputs in `.hzr_repeat_outputs` order,
then `FIRST` and `LAST`. That is what makes the `names<-` rename correct, and the
`expected_events()` test fails if the two orders ever diverge.

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `Rscript -e 'devtools::test(filter = "sas-translate-repeat")'`
Expected: 0 failures.

Then **mutate the implementation once to prove the tests can fail**, and revert each change:
- Delete the `d <- d[setdiff(names(d), drop)]` line. The clash test should now fail on
  `sum(names(...) == "EV_CARD")`.
- Swap `"FIRST", "LAST"` to `"LAST", "FIRST"` in `targets`. `expected_events()` should now
  fail.

- [ ] **Step 5: Commit**

```bash
git add R/sas-parse-job.R tests/testthat/test-sas-translate-repeat.R
git commit -m "feat: translate a %repeat call into an hzr_repeated_events() chunk

Renames outputs to the job's upper-case names, drops input columns SAS would
overwrite (with a warning), warns on LAG_IV/NUMBER, and refuses calls it
cannot express.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 3: wire `REPEAT` blocks into `hzr_translate_sas()` and stop on post-macro rewrites

**Files:**
- Modify: `R/sas-parse-job.R` (append `.hzr_repeat_rewrites()`)
- Modify: `R/translate-sas.R` (loop body at line 175 onward: replace Task 1's `next` and add the rewrite check inside the `HAZARD` branch)
- Modify: `tests/testthat/test-sas-translate-repeat.R` (append)

**Interfaces:**
- Consumes: `.hzr_sas_blocks()` `start`/`end` (Task 1); `.hzr_parse_repeat()` (Task 2).
- Produces: `.hzr_repeat_rewrites(segment, out)` returns a `character` vector, one element
  per step that writes `out`, each the step's normalised statements joined by `"; "` and
  ending in `;`.
- Produces: new chunk kinds in `job$calls`: `data` (the `IN=` guard, via
  `.hzr_next_call_name(calls, "data")`), `repeated`, and `rewrite`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-sas-translate-repeat.R`:

```r
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
  expect_true("DATA EVENTS after %repeat" %in% job$untranslated$construct)
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
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `Rscript -e 'devtools::test(filter = "sas-translate-repeat")'`
Expected: FAIL. Chunk names come out as `c("data", "status", "fit")`, because Task 1's `next`
skips the macro and the guard is on `EVENTS`. The call to `.hzr_repeat_rewrites` errors as
not found once it is referenced.

- [ ] **Step 3: Implement**

Append to `R/sas-parse-job.R`:

```r
#' Steps between a `%repeat` call and a fit that rewrite the macro's `OUT=`.
#'
#' `segment` is the normalised source between the two. A hit is a `DATA`
#' statement naming `out` among its output datasets (options in parentheses
#' ignored) or a `PROC SQL` `CREATE TABLE out`. Each is returned quoted, from
#' that statement up to the next `DATA`, `PROC`, `%HAZ`, `%REPEAT` or `RUN`.
#' `PROC SORT` is not a rewrite: reordering rows does not change the
#' likelihood.
#' @noRd
.hzr_repeat_rewrites <- function(segment, out) {
  stmts <- trimws(strsplit(segment, ";", fixed = TRUE)[[1L]])
  stmts <- stmts[nzchar(stmts)]
  boundary <- "^(DATA |PROC |%HAZ|%REPEAT|RUN$)"
  hits <- character(0)
  for (i in seq_along(stmts)) {
    s <- stmts[i]
    tok <- strsplit(s, " ", fixed = TRUE)[[1L]]
    writes <- if (startsWith(s, "DATA ")) {
      out %in% strsplit(trimws(gsub("[(][^)]*[)]", " ", substring(s, 6L))), " +")[[1L]]
    } else {
      j <- which(tok == "CREATE")
      any(tok[j + 1L] %in% "TABLE" & tok[j + 2L] %in% out)
    }
    if (writes) {
      last <- i
      while (last < length(stmts) && !grepl(boundary, stmts[last + 1L])) last <- last + 1L
      hits <- c(hits, paste0(paste(stmts[i:last], collapse = "; "), ";"))
    }
  }
  hits
}
```

In `R/translate-sas.R`, after `first_unresolved_inhaz <- NULL` / `n_unresolved_inhaz <- 0L`,
add:

```r
  repeat_scan <- list()      # %repeat OUT= -> offset its rewrite scan resumes from
```

Replace Task 1's `if (identical(b$proc, "REPEAT")) next # translated in Task 3` with:

```r
    if (identical(b$proc, "REPEAT")) {
      r <- .hzr_parse_repeat(b)
      # The macro's input is built by the job's own DATA steps, which this
      # translator does not translate: the same loud guard a PROC HAZARD
      # DATA= gets, once per name.
      if (!is.null(r$in_name) && !(r$in_name %in% guarded_data)) {
        calls[[.hzr_next_call_name(calls, "data")]] <- bquote(
          if (!exists(.(r$in_name))) {
            stop("This job built ", .(r$in_name), " in SAS DATA steps, which ",
                 "hzr_translate_sas() does not translate. Assign ", .(r$in_name),
                 " as it stood at the job's %repeat call, with its columns named ",
                 "as the job spells them, in upper case, before rendering.")
          }
        )
        guarded_data <- c(guarded_data, r$in_name)
      }
      calls[[.hzr_next_call_name(calls, "repeated")]] <- r$call
      # OUT= is now built by a chunk, so a fit reading it needs no guard.
      if (!is.null(r$out_name)) {
        guarded_data <- c(guarded_data, r$out_name)
        repeat_scan[[r$out_name]] <- b$end
      }
    } else if (identical(b$proc, "HAZARD")) {
```

and remove the original `if (identical(b$proc, "HAZARD")) {` line that this now replaces.
The `} else {` HAZPRED branch and the trailing `untr <- rbind(untr, r$untranslated)` /
`seen` / `mapped` lines are unchanged. They already accumulate the `REPEAT` branch's
fields.

Inside the `HAZARD` branch, immediately after the `tryCatch(.hzr_parse_hazard(b), ...)`
assignment and **before** the existing `if (!is.null(r$call[["data"]]))` guard block, add:

```r
      # A DATA step between %repeat and this fit that rewrites the macro's
      # OUT= is job code this translator does not fold in. Fitting without it
      # fits data SAS did not fit, so stop here -- ahead of the status chunk,
      # which writes into the same data frame. The scan resumes where the
      # last one ended, so one rewrite stops once however many fits follow.
      dname <- if (is.null(r$call[["data"]])) NULL else as.character(r$call[["data"]])
      if (!is.null(dname) && !is.null(repeat_scan[[dname]])) {
        seg <- substring(txt, repeat_scan[[dname]] + 1L, b$start - 1L)
        for (step in .hzr_repeat_rewrites(seg, dname)) {
          calls[[.hzr_next_call_name(calls, "rewrite")]] <- bquote(stop(.(paste0(
            "This job changes ", dname, " after %repeat, in a SAS DATA step that ",
            "hzr_translate_sas() does not translate: ", step, " Replace this chunk ",
            "with R code that makes the same change to ", dname, "."
          ))))
          untr <- rbind(untr, .hzr_untranslated_frame(
            NA_integer_, paste0("DATA ", dname, " after %repeat"), step
          ))
        }
        repeat_scan[[dname]] <- b$end
      }
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `Rscript -e 'devtools::test(filter = "sas-translate-repeat|translate-sas|sas-translate")'`
Expected: 0 failures, and the fit test reports no skip. `devtools::test()` sets `NOT_CRAN`,
so check the `[ FAIL 0 | WARN 0 | SKIP n | PASS m ]` line and confirm the 4/15 test is not
among the skips.

Mutate once, then revert. Change `.hzr_repeat_rewrites(seg, dname)` to
`.hzr_repeat_rewrites("", dname)`, so the scan sees nothing. The rewrite test should now fail
on the chunk names. If it passes, it cannot detect a missed rewrite.

- [ ] **Step 5: Commit**

```bash
git add R/sas-parse-job.R R/translate-sas.R tests/testthat/test-sas-translate-repeat.R
git commit -m "feat: hzr_translate_sas() translates %repeat and stops on post-macro rewrites

Guards the macro's IN=, drops the guard on its OUT=, and emits a quoted stop()
before any fit that reads an OUT= a DATA step rewrote after the macro.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4: volume-gated parity on the real cardioversion job

**Files:**
- Modify: `tests/testthat/test-repeated-events-parity.R` (append)

**Interfaces:**
- Consumes: `skip_if_no_maze_datasets()` (defined at the top of this file);
  `.hzr_derive_bd_card(dir)` and `.hzr_maze_datasets_dir()` (`inst/sas-parity/helper-sas-parity.R`,
  sourced by `tests/testthat/helper-sas-parity.R`); the chunk names from Task 3.

- [ ] **Step 1: Write the test**

```r
test_that("hzr_translate_sas() on the cardioversion job builds EVENTS, then stops at line 65", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("haven")
  dir <- skip_if_no_maze_datasets()
  sas <- file.path(dirname(dir), "distributions", "hz.ce_cardioversion_repeated.ehb.sas")
  testthat::skip_if_not(file.exists(sas), "hz.ce_cardioversion_repeated.ehb.sas not on the volume")

  job <- suppressWarnings(hzr_translate_sas(sas))
  expect_equal(names(job$calls)[1:5], c("data", "repeated", "rewrite", "status", "fit"))

  bd_card <- .hzr_derive_bd_card(dir)
  names(bd_card) <- toupper(names(bd_card))
  env <- new.env(parent = globalenv())
  env$BD_CARD <- bd_card
  eval(job$calls$data, env)
  expect_no_warning(eval(job$calls[["repeated"]], env))
  # Shape first, then values: the .log's 962 rows, and 357 input columns plus
  # the function's eight.
  expect_equal(dim(env$EVENTS), c(962L, 365L))
  expect_equal(sum(env$EVENTS$CE_CARD == 1), 388L)
  # The job's line 65 is job code, not macro code; the document stops on it.
  expect_error(eval(job$calls$rewrite, env), "IV_EVENT=IV_EVENT+0.0001141553", fixed = TRUE)
})
```

- [ ] **Step 2: Run it with the volume mounted**

Run: `Rscript -e 'devtools::test(filter = "repeated-events-parity")'`
Expected: PASS, and **not** SKIP. If it skips, the volume is not mounted and the test is not
evidence. Say so in the PR rather than counting it green.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-repeated-events-parity.R
git commit -m "test: translated cardioversion job builds EVENTS and stops at its line 65

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 5: documentation, NEWS, spec pointers

**Files:**
- Modify: `R/translate-sas.R` (roxygen: add a section after `@section Experimental:`)
- Modify: `NEWS.md` (bullet under `# TemporalHazard 1.2.11` → `## New features`)
- Modify: `inst/dev/REPEATED-EVENTS-DESIGN.md` ("Out of scope", the `hzr_translate_sas()` bullet)
- Modify: `inst/dev/SAS-JOB-TRANSLATOR-DESIGN.md` §5.5 "Evidence" (final census count)
- Modify: `inst/WORDLIST` (only the words spelling flags)
- Regenerate: `man/hzr_translate_sas.Rd`

- [ ] **Step 1: Roxygen.** In `R/translate-sas.R`, insert before `#' @section Comparing a translated fit against a SAS listing:`:

```r
#' @section Repeated events:
#' A job that builds its fit input with the SAS macro `%repeat` has that call
#' translated to [hzr_repeated_events()]. The emitted chunk renames the
#' function's output columns to the names the job gives them, upper-cased like
#' every name the translator emits, so the fit reads them. An input column named
#' like one of those outputs is dropped first, with a warning, because the macro
#' overwrites it. The macro's input is built by the job's own DATA steps, which
#' are not translated, so the document stops until that input is assigned. A DATA
#' step that changes the macro's output before the fit is not translated either.
#' Its chunk stops and quotes the step, for the reader to replace with R code.
#'
```

- [ ] **Step 2: Regenerate and check the Rd escaping**

Run: `Rscript -e 'devtools::document()'` and then `grep -n 'repeat' man/hzr_translate_sas.Rd`
Expected: the section is present, and each backticked `%repeat` appears as
`\verb{\%repeat}`. Roxygen escapes it inside a code span, as `man/hzr_repeated_events.Rd`
line 5 already shows (checked 2026-09-10). An unescaped `%` in Rd starts a comment and
silently truncates the line, so keep every `%repeat` in the roxygen text inside backticks.

- [ ] **Step 3: NEWS.** In `NEWS.md`, insert immediately before the line `# TemporalHazard 1.2.10` (keep one blank line on each side):

```markdown
* **`hzr_translate_sas()` translates a `%repeat` call** into
  `hzr_repeated_events()` (#241), renaming its outputs to the names the job
  gives them so the job's fit reads them. The macro's input is still the
  reader's to supply. A DATA step that changes the macro's output before the
  fit now stops the document with the step quoted; the translation no longer
  fits data that SAS did not fit.
```

- [ ] **Step 4: Spec pointers.** In `inst/dev/REPEATED-EVENTS-DESIGN.md`, replace the "Out of scope" bullet that begins `- Wiring \`hzr_repeated_events()\` into \`hzr_translate_sas()\`.` with:

```markdown
- Wiring `hzr_repeated_events()` into `hzr_translate_sas()`. Done separately; see
  `SAS-JOB-TRANSLATOR-DESIGN.md` §5.5.
```

In `SAS-JOB-TRANSLATOR-DESIGN.md` §5.5 "Evidence", replace "found **at least 526**" and "138
of them" with the census's final counts. Read them from the scratchpad
`repeat-census.txt` once its background scan has printed `exit=`. If it has not finished, keep
"at least" and give the latest count and the time it was read.

- [ ] **Step 5: Commit**

```bash
git add R/translate-sas.R man/hzr_translate_sas.Rd NEWS.md inst/dev/REPEATED-EVENTS-DESIGN.md inst/dev/SAS-JOB-TRANSLATOR-DESIGN.md
git commit -m "docs: document %repeat translation in hzr_translate_sas() and NEWS

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 6: gates, review, PR

**Files:** `inst/WORDLIST` only if spelling flags words.

- [ ] **Step 1: The definition of done, in order.** Each must pass before the next:

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
Rscript -e 'spelling::spell_check_package(use_wordlist = TRUE)'
```

Expected: `document()` leaves `git status` clean (commit it if not). Lint prints no lints.
Tests: `FAIL 0`. Record the `SKIP`/`PASS` counts for the PR. Spelling: no words, or add
each flagged SAS identifier to `inst/WORDLIST` (sorted) and commit.

- [ ] **Step 2: `R CMD check --as-cran` from a committed tree, under the scratchpad**

Run the full suite in the **foreground** with a 600000 ms timeout; never background it and
return early.

```bash
SCRATCH=/private/tmp/claude-504/-Users-ehrlinj-Documents-GitHub-TemporalHazard--claude-worktrees-charming-mirzakhani-331055/831a925a-4e5d-4c94-a0fd-c82969ba4fef/scratchpad
rm -rf "$SCRATCH/tree" "$SCRATCH/check" && mkdir -p "$SCRATCH/tree" "$SCRATCH/check"
git archive HEAD | tar -x -C "$SCRATCH/tree"
cd "$SCRATCH/check" && R CMD build "$SCRATCH/tree" && R CMD check --as-cran TemporalHazard_1.2.11.tar.gz
tar tzf TemporalHazard_1.2.11.tar.gz | grep -iE 'CLAUDE|AGENTS|REPEAT-TRANSLATOR-PLAN' || echo "no dev files in tarball"
```

Expected: `Status: OK`, or only the NOTEs already present on `main`; no WARNING or ERROR.
Record the overall time against the 3m 41s (221s) baseline in `AGENTS.md`.

- [ ] **Step 3: Advisory review.** Dispatch the `r-reviewer` agent over `git diff origin/main...HEAD`. Verify each finding against the code before acting on it; a clean report does not replace Steps 1–2.

- [ ] **Step 4: Push the branch and open the PR.** Never push `main`.

```bash
git push -u origin feat/translate-repeat
gh pr create --base main --title "hzr_translate_sas(): translate %repeat into hzr_repeated_events()" --body-file <body>
```

The body states:
- what is translated;
- the three decisions (clash drop, rewrite stop, `IN=` guard) and the `LAG_IV`/`NUMBER` warning;
- the gate results with counts, and whether the volume-gated test **ran or skipped**;
- that fitting the cardioversion job to LL −267.885 is still the next piece of work.

End the body with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`. Then
stop: the maintainer merges.
