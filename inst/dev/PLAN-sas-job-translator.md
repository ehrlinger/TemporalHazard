# SAS Job Translator Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Translate a SAS `PROC HAZARD` / `PROC HAZPRED` job into a `.qmd` that
reproduces it with `TemporalHazard`, or fails loudly rather than emitting a
document over a model it could not build.

**Architecture:** A context-aware lexer and parser drive off `.hzr_sas_grammar`
— a keyword table generated from HAZARD's own lex sources (already landed in
`59a68fd`). The parser produces an internal `hzr_sas_job` holding **unevaluated
R calls**, so rendering is `deparse()` rather than string templating and the
parity harness can `eval()` the fit call directly.

**Tech Stack:** Base R plus `survival` (the sole `Imports`). `testthat` 3e.
Quarto for rendered output, already in `Suggests`.

**Spec:** `inst/dev/SAS-JOB-TRANSLATOR-DESIGN.md`. Read it before Task 1.

## Global Constraints

- **Branch:** `feat/sas-job-translator`. Never push to `main` or `dev`.
- **No new dependencies.** `Imports:` stays `survival` alone. `librefs` is a
  named character vector, never a config-file path.
- **One new export only:** `hzr_translate_sas()`. `hzr_sas_job` stays internal.
- **`stats::` prefixing is house style.** An unqualified `stats` generic is an
  `R CMD check` NOTE that no test will catch.
- **Censoring codes are `-1` left, `0` right, `1` event, `2` interval.**
  `survival::Surv()` uses different integers for the same meanings.
- **No `browser()`, no bare `print()`, no `library()` in `R/`.** `cat()` only in
  a `print.*` method.
- **Definition of done, every task:** `devtools::document()` → `lintr::lint_package()`
  (must be 0) → `devtools::test()` (0 failures). Run them; reading the code is
  not evidence.
- **Assertions must be able to fail.** `expect_equal(status, c(1, -1))`, never
  `expect_true(all(status %in% <every possible code>))`.
- **No production job content in the repo.** Fixtures are synthesised or come
  from the public `examples/` corpus.

---

## File Structure

| Path | Responsibility |
|---|---|
| `R/sas-lex.R` | comment stripping, normalisation, balanced-paren block extraction |
| `R/sas-grammar-lookup.R` | context-aware keyword → token resolution against `.hzr_sas_grammar` |
| `R/sas-parse-job.R` | statement dispatch; builds the `hazard()` / `predict()` calls |
| `R/sas-parse-parms.R` | `PARMS` → `phases` list + `theta`; the G3/Weibull logic |
| `R/sas-job.R` | `hzr_sas_job` constructor, validator, `print` method |
| `R/sas-render-qmd.R` | `hzr_sas_job` → `.qmd` text |
| `R/translate-sas.R` | `hzr_translate_sas()`, `librefs`, `INHAZ` resolution |

Split by responsibility rather than layer. `sas-parse-parms.R` is separate from
`sas-parse-job.R` because the parameter/phase mapping is the densest logic in the
feature and benefits from being readable on its own.

---

### Task 1: Lexer

`tools/sas-hazard-profile.R` PART 1 is a validated prototype of this — it ran
clean over 990 files across four corpora. Port it, do not rewrite it.

**Files:**
- Create: `R/sas-lex.R`
- Test: `tests/testthat/test-sas-lex.R`

**Interfaces:**
- Consumes: nothing.
- Produces: `.hzr_sas_strip_comments(lines)` → character vector;
  `.hzr_sas_normalise(lines)` → length-1 uppercase single-spaced string.

- [ ] **Step 1: Write the failing test**

```r
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-lex")'`
Expected: FAIL, `could not find function ".hzr_sas_normalise"`.

- [ ] **Step 3: Port the prototype**

Copy `.idx`, `.re`, `.first_semi`, `sas_strip_comments`, `sas_strip_inline_comments`
and `sas_normalise` from `tools/sas-hazard-profile.R` into `R/sas-lex.R`, renaming
the exported-looking ones to the dotted internal form: `.hzr_sas_strip_comments()`,
`.hzr_sas_strip_inline_comments()`, `.hzr_sas_normalise()`. Keep `.idx`, `.re` and
`.first_semi` as-is. Add roxygen `@noRd` to each.

Do not change the logic. It matches HAZARD's own comment rule and has been
validated against 990 files.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-lex")'`
Expected: PASS, 4 tests.

- [ ] **Step 5: Full gate, then commit**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
git add R/sas-lex.R tests/testthat/test-sas-lex.R NAMESPACE man
git commit -m "feat(sas): SAS comment stripping and normalisation"
```

---

### Task 2: Balanced-paren block extraction

**Files:**
- Modify: `R/sas-lex.R`
- Test: `tests/testthat/test-sas-lex.R`

**Interfaces:**
- Consumes: `.hzr_sas_normalise()`.
- Produces: `.hzr_sas_blocks(txt)` → list of
  `list(proc = "HAZARD"|"HAZPRED", text = <chr>, terminator = "paren"|"none")`.

Rationale, from the spec §3.3: `hazard_l.l` puts `)` in its **whitespace class**,
so block delimitation is the `%HAZARD(...)` macro's job. Scanning forward for the
first `);` failed on **38 of 82 blocks** in the `lv_function` study because the
block contained a nested paren. Match balanced parens from `%HAZARD(` instead.

- [ ] **Step 1: Write the failing test**

```r
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-lex")'`
Expected: FAIL, `could not find function ".hzr_sas_blocks"`.

- [ ] **Step 3: Implement**

```r
#' Extract PROC HAZARD / PROC HAZPRED blocks from normalised source.
#'
#' Blocks are delimited by the `%HAZARD(...)` macro's parentheses, not by the
#' PROC statement: HAZARD's own lexer treats `)` as whitespace, so the macro
#' owns delimitation. Bounding on the first `);` fails whenever the block
#' contains a nested paren.
#' @noRd
.hzr_sas_blocks <- function(txt) {
  out <- list()
  starts <- gregexpr("%HAZ(ARD|PRED) *\\(", txt)[[1L]]
  if (starts[1L] == -1L) return(out)

  for (s in starts) {
    open_at <- s + attr(starts, "match.length")[which(starts == s)] - 1L
    depth <- 0L
    close_at <- NA_integer_
    for (i in seq(open_at, nchar(txt))) {
      ch <- substr(txt, i, i)
      if (ch == "(") depth <- depth + 1L
      if (ch == ")") {
        depth <- depth - 1L
        if (depth == 0L) { close_at <- i; break }
      }
    }
    term <- if (is.na(close_at)) "none" else "paren"
    body <- substring(txt, open_at + 1L,
                      if (is.na(close_at)) nchar(txt) else close_at - 1L)
    proc <- if (grepl("PROC HAZARD", body)) "HAZARD" else "HAZPRED"
    out[[length(out) + 1L]] <- list(proc = proc, text = trimws(body),
                                    terminator = term)
  }
  out
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-lex")'`
Expected: PASS, 7 tests.

- [ ] **Step 5: Regression-check against the public corpus**

This is the measurement that justified the rewrite, so verify the fix rather
than assume it.

```bash
Rscript -e '
devtools::load_all(quiet = TRUE)
fs <- list.files("~/Documents/GitHub/hazard", pattern = "[.]sas$",
                 recursive = TRUE, full.names = TRUE)
n <- c(paren = 0L, none = 0L)
for (f in fs) {
  b <- .hzr_sas_blocks(.hzr_sas_normalise(readLines(f, warn = FALSE)))
  for (x in b) n[x$terminator] <- n[x$terminator] + 1L
}
print(n)'
```

Expected: `none` is 0 or very close to it. If `none` is large, the paren matcher
is wrong — do not proceed.

- [ ] **Step 6: Full gate, then commit**

```bash
Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && Rscript -e 'devtools::test()'
git add R/sas-lex.R tests/testthat/test-sas-lex.R
git commit -m "feat(sas): bound HAZ blocks on balanced macro parens"
```

---

### Task 3: Context-aware grammar lookup

**Files:**
- Create: `R/sas-grammar-lookup.R`
- Test: `tests/testthat/test-sas-grammar-lookup.R`

**Interfaces:**
- Consumes: `.hzr_sas_grammar` (in `R/sysdata.rda`, already landed).
- Produces: `.hzr_sas_token(keyword, proc, context)` → length-1 character token,
  or `NA_character_` if the keyword is not in that context.

Spec §3.2: `M` is the early-phase shape parameter in `PARM` context but `MOVE`
in `PHOP`/`STEP`; `NOS` is `NOPRINTS` in `STEP` but `NOSURV` in `HZPP`. `M`
occurs 76–115 times **per study** as a `PARMS` parameter, so a context-free
lookup produces a wrong model with no error.

- [ ] **Step 1: Write the failing test**

```r
test_that("context decides which token a colliding spelling resolves to", {
  expect_equal(.hzr_sas_token("M", "HAZARD", "PARM"), "M")
  expect_equal(.hzr_sas_token("M", "HAZARD", "PHOP"), "MOVE")
  expect_equal(.hzr_sas_token("NOS", "HAZPRED", "HZPP"), "NOSURV")
})

test_that("aliases resolve to their canonical token", {
  expect_equal(.hzr_sas_token("MI", "HAZARD", "HZRP"), "MAXITER")
  expect_equal(.hzr_sas_token("QUASI", "HAZARD", "HZRP"), "QUASINEWTON")
  expect_equal(.hzr_sas_token("PARMS", "HAZARD", "STMT"), "PARAMETERS")
})

test_that("an unknown keyword returns NA rather than guessing", {
  expect_true(is.na(.hzr_sas_token("NOTAKEYWORD", "HAZARD", "HZRP")))
  # Right spelling, wrong context is still unknown.
  expect_true(is.na(.hzr_sas_token("SLENTRY", "HAZARD", "HZRP")))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-grammar-lookup")'`
Expected: FAIL, `could not find function ".hzr_sas_token"`.

- [ ] **Step 3: Implement**

```r
#' Resolve a SAS keyword to its parser token within a proc and lexer context.
#'
#' Returns NA_character_ when the keyword is unknown in that context. Callers
#' must treat NA as untranslated rather than falling back to a context-free
#' match: `M` is the early-phase shape parameter in PARM context and MOVE in
#' PHOP, and guessing between them silently changes the model.
#' @noRd
.hzr_sas_token <- function(keyword, proc, context) {
  g <- .hzr_sas_grammar
  hit <- g$keyword == keyword & g$proc == proc & g$context == context
  if (!any(hit)) return(NA_character_)
  unique(g$token[hit])[[1L]]
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-grammar-lookup")'`
Expected: PASS, 8 tests.

- [ ] **Step 5: Full gate, then commit**

```bash
Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && Rscript -e 'devtools::test()'
git add R/sas-grammar-lookup.R tests/testthat/test-sas-grammar-lookup.R
git commit -m "feat(sas): context-aware grammar lookup"
```

---

### Task 4: `hzr_sas_job` container

**Files:**
- Create: `R/sas-job.R`
- Test: `tests/testthat/test-sas-job.R`

**Interfaces:**
- Produces: `.hzr_sas_job(source, calls, grid, inhaz, outhaz, untranslated, coverage)`
  where `source` is `list(path =, checksum =)`. **`checksum` is an md5 from
  `tools::md5sum()`** — base R has no SHA-256, and adding `digest` for a
  provenance stamp would breach the no-new-dependency constraint.
  → object of class `hzr_sas_job`; `print.hzr_sas_job(x, ...)`;
  `.hzr_validate_sas_job(x)` → `x` invisibly, or `stop()`.

`untranslated` is a data frame (`line`, `construct`, `reason`) so the report can
name the construct, never just count it. `coverage` is
`list(tokens_seen =, tokens_mapped =)`.

- [ ] **Step 1: Write the failing test**

```r
test_that("a job with no statements seen is rejected, not returned empty", {
  # Three shapes of absent, ascending danger: NULL, empty container, and a
  # hollow object -- right shape, empty inside -- which defeats is.null() and
  # a length check both. This is the hollow case.
  j <- .hzr_sas_job(
    source = list(path = "x.sas", checksum = "abc"),
    calls = list(), grid = NULL, inhaz = NULL, outhaz = NULL,
    untranslated = .hzr_untranslated_frame(),
    coverage = list(tokens_seen = 0L, tokens_mapped = 0L)
  )
  expect_error(.hzr_validate_sas_job(j), "no SAS statements")
})

test_that("a populated job validates and prints its coverage", {
  j <- .hzr_sas_job(
    source = list(path = "x.sas", checksum = "abc"),
    calls = list(fit = quote(hazard(time = t, status = s))),
    grid = NULL, inhaz = NULL, outhaz = "LIB.D",
    untranslated = .hzr_untranslated_frame(),
    coverage = list(tokens_seen = 4L, tokens_mapped = 4L)
  )
  expect_silent(.hzr_validate_sas_job(j))
  expect_output(print(j), "4/4")
})

test_that("the untranslated frame keeps construct and line, not a count", {
  f <- .hzr_untranslated_frame(line = 12L, construct = "RESTRICT",
                               reason = "no R equivalent")
  expect_equal(f$construct, "RESTRICT")
  expect_equal(f$line, 12L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-job")'`
Expected: FAIL, `could not find function ".hzr_sas_job"`.

- [ ] **Step 3: Implement**

```r
#' Empty (or single-row) untranslated-construct frame.
#' @noRd
.hzr_untranslated_frame <- function(line = integer(0),
                                    construct = character(0),
                                    reason = character(0)) {
  data.frame(line = as.integer(line), construct = as.character(construct),
             reason = as.character(reason), stringsAsFactors = FALSE)
}

#' Construct an hzr_sas_job.
#' @noRd
.hzr_sas_job <- function(source, calls, grid, inhaz, outhaz,
                         untranslated, coverage) {
  structure(list(source = source, calls = calls, grid = grid,
                 inhaz = inhaz, outhaz = outhaz,
                 untranslated = untranslated, coverage = coverage),
            class = "hzr_sas_job")
}

#' Validate an hzr_sas_job.
#'
#' A job where nothing was seen is a parse failure wearing a green badge: the
#' object has the right shape and is empty inside. Reject it here rather than
#' letting a renderer emit a document over nothing.
#' @noRd
.hzr_validate_sas_job <- function(x) {
  stopifnot(inherits(x, "hzr_sas_job"))
  if (!is.numeric(x$coverage$tokens_seen) || x$coverage$tokens_seen < 1L) {
    stop("no SAS statements were recognised in ", x$source$path,
         "; this is a parse failure, not an empty job.", call. = FALSE)
  }
  invisible(x)
}

#' @export
print.hzr_sas_job <- function(x, ...) {
  cat("<hzr_sas_job>", basename(x$source$path), "\n")
  cat("  calls       :", paste(names(x$calls), collapse = ", "), "\n")
  cat("  coverage    : ", x$coverage$tokens_mapped, "/",
      x$coverage$tokens_seen, " tokens mapped\n", sep = "")
  if (nrow(x$untranslated)) {
    cat("  untranslated:", nrow(x$untranslated), "construct(s)\n")
    for (i in seq_len(nrow(x$untranslated))) {
      cat("    -", x$untranslated$construct[i], "-",
          x$untranslated$reason[i], "\n")
    }
  }
  invisible(x)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-job")'`
Expected: PASS, 5 tests.

- [ ] **Step 5: Full gate, then commit**

```bash
Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && Rscript -e 'devtools::test()'
git add R/sas-job.R tests/testthat/test-sas-job.R NAMESPACE man
git commit -m "feat(sas): hzr_sas_job container with hollow-object guard"
```

---

### Task 5: `PARMS` → phases and theta

The densest task. Read spec §5 and §7.1 first.

**Files:**
- Create: `R/sas-parse-parms.R`
- Test: `tests/testthat/test-sas-parse-parms.R`

**Interfaces:**
- Consumes: `.hzr_sas_token()`.
- Produces: `.hzr_parse_parms(operands, covars = list())` →
  `list(phases = <call>, theta = <call>, untranslated = <data.frame>)`
  where `operands` is a character vector of `PARMS` tokens such as
  `c("MUE=0.2", "THALF=0.15", "NU=1.4", "M=1", "FIXM", "MUC=0.0005")`.

  `covars` is an optional named list — `list(early = c("X1", "X2"), constant = ,
  late = )` — carrying the operands of the `EARLY` / `CONSTANT` / `LATE`
  statements. When a phase has covariates, its `hzr_phase()` call gains
  `formula = ~X1 + X2`, built with
  `str2lang(paste("~", paste(v, collapse = " + ")))`. These statements appear
  16–44 times per production study, so this is v1 scope, not an extension.

Mapping (spec §3.1 and §7.1):

| SAS | R |
|---|---|
| `MUE=`, `MUC=`, `MUL=` | `theta` entries, `log()` of the value |
| `THALF=`, `NU=`, `M=` | `hzr_phase("cdf", t_half=, nu=, m=)` |
| `TAU=`, `GAMMA=`, `ALPHA=`, `ETA=` | `hzr_phase("g3", tau=, gamma=, alpha=, eta=)` |
| `FIX<param>` | appended to that phase's `fixed=` |
| `WEIBULL` | G3 with `alpha = 1, eta = 1`, both fixed |
| presence of `MUC` | a `hzr_phase("constant")` between them |

- [ ] **Step 1: Write the failing test**

```r
test_that("an early-plus-constant PARMS maps to two phases and a theta", {
  ops <- c("MUE=0.2361727", "THALF=0.1512095", "NU=1.438652", "M=1", "FIXM",
           "MUC=0.0005436977")
  got <- .hzr_parse_parms(ops)
  expect_equal(
    got$phases,
    quote(list(
      hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                fixed = "m"),
      hzr_phase("constant")
    ))
  )
  expect_equal(got$theta, quote(c(log(0.2361727), log(0.0005436977))))
})

test_that("WEIBULL becomes a G3 constrained at alpha = 1, eta = 1", {
  # PARMS ... WEIBULL is setopt(6) -> SETG3_weibull, g3flag += 2. The R general
  # form (((t/tau)^gamma + 1)^(1/alpha) - 1)^eta collapses at alpha = eta = 1
  # to (t/tau)^gamma, a Weibull cumulative hazard. See spec 7.1.
  ops <- c("MUL=0.01", "TAU=2", "GAMMA=1.5", "WEIBULL", "FIXTAU", "FIXGAMMA")
  got <- .hzr_parse_parms(ops)
  expect_equal(
    got$phases,
    quote(list(
      hzr_phase("g3", tau = 2, gamma = 1.5, alpha = 1, eta = 1,
                fixed = c("tau", "gamma", "alpha", "eta"))
    ))
  )
})

test_that("phase covariates become a formula on the owning phase", {
  ops <- c("MUE=0.2", "THALF=1", "NU=1", "MUC=0.001")
  got <- .hzr_parse_parms(ops, covars = list(early = c("AGE", "SEX"),
                                             constant = "AGE"))
  expect_equal(
    got$phases,
    quote(list(
      hzr_phase("cdf", t_half = 1, nu = 1, formula = ~AGE + SEX),
      hzr_phase("constant", formula = ~AGE)
    ))
  )
})

test_that("an unrecognised PARMS token is recorded, never dropped", {
  got <- .hzr_parse_parms(c("MUE=0.2", "THALF=1", "NU=1", "FIXGAE2"))
  expect_equal(got$untranslated$construct, "FIXGAE2")
  expect_equal(nrow(got$untranslated), 1L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-parse-parms")'`
Expected: FAIL, `could not find function ".hzr_parse_parms"`.

- [ ] **Step 3: Implement**

Build three accumulators — `early` (`t_half`, `nu`, `m`), `late`
(`tau`, `gamma`, `alpha`, `eta`), and `mu` (`MUE`, `MUC`, `MUL` in that order) —
plus a `fixed` character vector per phase and an untranslated frame. Walk the
operands once:

- `KEY=VALUE` → resolve `KEY` via `.hzr_sas_token(KEY, "HAZARD", "PARM")`;
  route `MUE`/`MUC`/`MUL` to `mu`, `THALF`/`NU`/`M` to `early`, and
  `TAU`/`GAMMA`/`ALPHA`/`ETA` to `late`. `NA` token → untranslated.
- Bare `FIX<X>` → append the **R parameter name** to the owning phase's `fixed`,
  via an explicit map. `tolower(X)` is WRONG: `FIXTHALF`'s R parameter is
  `t_half`, not `thalf`, and `hzr_phase()` accepts only `"t_half"`, `"nu"`,
  `"m"` (cdf/hazard) or `"tau"`, `"gamma"`, `"alpha"`, `"eta"` (g3). A name
  outside that set leaves the parameter free — a different model, no error.
  `FIXDELTA`, `FIXMNU1`, `FIXGE2` and `FIXGAE2` have no owning `hzr_phase()`
  parameter and go to `untranslated`.
- Bare `WEIBULL` → set `late$alpha <- 1`, `late$eta <- 1`, and append
  `"alpha"` and `"eta"` to the late phase's `fixed`.
- Any other bare token with an `NA` token → untranslated.

Then assemble with `bquote()` so the result is a call, not a string:

```r
phase_calls <- list()
if (length(early)) {
  phase_calls[[length(phase_calls) + 1L]] <- as.call(c(
    quote(hzr_phase), "cdf", early,
    if (length(fixed_early)) list(fixed = fixed_early)
  ))
}
if (has_muc) phase_calls[[length(phase_calls) + 1L]] <- quote(hzr_phase("constant"))
if (length(late)) {
  phase_calls[[length(phase_calls) + 1L]] <- as.call(c(
    quote(hzr_phase), "g3", late,
    if (length(fixed_late)) list(fixed = fixed_late)
  ))
}
phases <- as.call(c(quote(list), phase_calls))
theta  <- as.call(c(quote(c), lapply(mu, function(v) bquote(log(.(v))))))
```

Order matters: early, then constant, then late. `fixed` is a length-1 character
when there is one entry and a `c(...)` call otherwise — matching the tests above.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-parse-parms")'`
Expected: PASS, 3 tests.

- [ ] **Step 5: Full gate, then commit**

```bash
Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && Rscript -e 'devtools::test()'
git add R/sas-parse-parms.R tests/testthat/test-sas-parse-parms.R
git commit -m "feat(sas): map PARMS to phases and theta, incl. the Weibull G3"
```

---

### Task 6: Censoring statements → status codes

**The highest-risk mapping in the feature.** Read spec §5.1 before starting.

**Files:**
- Modify: `R/sas-parse-job.R` (created here)
- Test: `tests/testthat/test-sas-censoring.R`

**Interfaces:**
- Produces: `.hzr_censor_spec(statements)` → `list(status_expr = <call>, time_lower = <name|NULL>, time_upper = <name|NULL>, untranslated = <df>)`
  where `statements` is a named list such as
  `list(EVENT = "DEAD", TIME = "T", ICENSOR = c("LO", "HI"))`.

This package codes `-1` left, `0` right, `1` event, `2` interval.
`survival::Surv(type = "interval")` uses `0`/`1`/`2`/`3` for right/event/left/interval.
Carrying one coding into the other is a wrong answer with no error and has
shipped here before.

- [ ] **Step 1: Write the failing test**

```r
test_that("ICENSOR produces interval code 2 and both bounds", {
  st <- list(EVENT = "DEAD", TIME = "T", ICENSOR = c("LO", "HI"))
  got <- .hzr_censor_spec(st)
  expect_equal(got$time_lower, as.name("LO"))
  expect_equal(got$time_upper, as.name("HI"))
  # Assert the exact codes. `all(status %in% c(0,1,2,3))` cannot fail and is
  # the assertion shape that hid this class of bug for years.
  expect_equal(eval(got$status_expr, list(DEAD = c(1, 0), LO = c(1, 1), HI = c(2, 2))),
               c(1, 2))
})

test_that("LCENSOR produces -1, not survival's 2", {
  st <- list(EVENT = "DEAD", TIME = "T", LCENSOR = "LFLAG")
  got <- .hzr_censor_spec(st)
  expect_equal(eval(got$status_expr, list(DEAD = c(0, 0), LFLAG = c(1, 0))),
               c(-1, 0))
})

test_that("RCENSOR produces 0", {
  st <- list(EVENT = "DEAD", TIME = "T", RCENSOR = "RFLAG")
  got <- .hzr_censor_spec(st)
  expect_equal(eval(got$status_expr, list(DEAD = c(0, 1), RFLAG = c(1, 0))),
               c(0, 1))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-censoring")'`
Expected: FAIL, `could not find function ".hzr_censor_spec"`.

- [ ] **Step 3: Implement**

Build a nested `ifelse()` call, most specific first, defaulting to the `EVENT`
variable itself:

```r
#' Translate SAS censoring statements to this package's status coding.
#'
#' Codes are -1 left, 0 right, 1 event, 2 interval. `survival::Surv()` uses
#' different integers for the same meanings under `type = "interval"`; never
#' carry one coding into the other.
#' @noRd
.hzr_censor_spec <- function(statements) {
  ev <- as.name(statements$EVENT)
  expr <- ev
  if (!is.null(statements$LCENSOR)) {
    expr <- bquote(ifelse(.(as.name(statements$LCENSOR)) == 1, -1, .(expr)))
  }
  if (!is.null(statements$ICENSOR)) {
    lo <- as.name(statements$ICENSOR[[1L]])
    expr <- bquote(ifelse(.(ev) == 0 & !is.na(.(lo)), 2, .(expr)))
  }
  list(
    status_expr = expr,
    time_lower = if (!is.null(statements$ICENSOR)) as.name(statements$ICENSOR[[1L]]),
    time_upper = if (!is.null(statements$ICENSOR)) as.name(statements$ICENSOR[[2L]]),
    untranslated = .hzr_untranslated_frame()
  )
}
```

`RCENSOR` needs no branch: right-censoring is code `0`, which is what the
`EVENT` variable already carries when the event did not occur. Add a comment
saying so, because the absence of a branch otherwise looks like an omission.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-censoring")'`
Expected: PASS, 3 tests.

- [ ] **Step 5: Full gate, then commit**

```bash
Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && Rscript -e 'devtools::test()'
git add R/sas-parse-job.R tests/testthat/test-sas-censoring.R
git commit -m "feat(sas): translate ICENSOR/LCENSOR/RCENSOR to package status codes"
```

---

### Task 7: Parse a `PROC HAZARD` block into a `hazard()` call

**Files:**
- Modify: `R/sas-parse-job.R`
- Test: `tests/testthat/test-sas-parse-hazard.R`

**Interfaces:**
- Consumes: `.hzr_sas_token()`, `.hzr_parse_parms()`, `.hzr_censor_spec()`.
- Produces: `.hzr_parse_hazard(block)` → `list(call = <call>, outhaz = <chr|NULL>, untranslated = <df>, tokens_seen = <int>, tokens_mapped = <int>)`.

Control mapping (spec §7): `MI`/`MAXITER` → `control$maxit`; `CONDITION` →
`control$condition`; `CONSERVE`/`NOCONSERVE` → `control$conserve` **always
emitted, never defaulted**; `QUASINEWTON` → `control$method = "bfgs"`;
`STEEPEST` → untranslated with the reason naming issue #145;
`PRINTIT`/`NOCOR`/`NOCOV`/`NOLOG`/`NONOTES`/`NOPRINT` → no R effect, not recorded
as untranslated.

- [ ] **Step 1: Write the failing test**

```r
test_that("a canonical AVC-style block becomes a hazard() call", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONSERVE P OUTHAZ=EX.HZD",
    "STEEPEST QUASI CONDITION=14 MI=200;",
    "EVENT DEAD; TIME INT_DEAD;",
    "PARMS MUE=0.2361727 THALF=0.1512095 NU=1.438652 M=1 FIXM MUC=0.0005436977; );"
  ))
  b <- .hzr_sas_blocks(txt)[[1L]]
  got <- .hzr_parse_hazard(b)

  expect_equal(got$outhaz, "EX.HZD")
  expect_equal(got$call[["data"]], as.name("AVCS"))
  expect_equal(got$call[["time"]], as.name("INT_DEAD"))
  expect_equal(got$call[["control"]],
               quote(list(maxit = 200, condition = 14, conserve = TRUE,
                          method = "bfgs")))
  # STEEPEST has no R equivalent and must be surfaced, not dropped.
  expect_true("STEEPEST" %in% got$untranslated$construct)
})

test_that("NOCONSERVE is emitted explicitly rather than defaulted", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A NOCONSERVE CONDITION=14;",
    "EVENT D; TIME T; PARMS MUE=1 THALF=1 NU=1; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(got$call[["control"]][["conserve"]], FALSE)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-parse-hazard")'`
Expected: FAIL, `could not find function ".hzr_parse_hazard"`.

- [ ] **Step 3: Implement**

```r
#' Parse a PROC HAZARD block into a hazard() or hzr_stepwise() call.
#'
#' Every keyword is resolved through .hzr_sas_token() in its lexer context.
#' A keyword that does not resolve is recorded in `untranslated`, never
#' skipped: a dropped option changes the model silently.
#' @noRd
.hzr_parse_hazard <- function(block) {
  st <- strsplit(block$text, ";", fixed = TRUE)[[1L]]
  untr <- .hzr_untranslated_frame()
  seen <- 0L
  mapped <- 0L
  note <- function(kw, reason) {
    untr <<- rbind(untr, .hzr_untranslated_frame(NA_integer_, kw, reason))
  }

  # --- statement 1: the PROC line and its options -------------------------
  toks <- strsplit(trimws(st[[1L]]), " ", fixed = TRUE)[[1L]]
  toks <- toks[nzchar(toks)]
  ctl <- list()
  data_name <- NULL
  outhaz <- NULL

  for (tok in toks) {
    eqp <- .idx(tok, "=")
    key <- if (eqp > 0L) substring(tok, 1L, eqp - 1L) else tok
    val <- if (eqp > 0L) substring(tok, eqp + 1L) else ""
    token <- .hzr_sas_token(key, "HAZARD", "HZRP")
    if (identical(token, "PROC") || identical(token, "HAZARD")) next
    seen <- seen + 1L
    if (is.na(token)) {
      note(key, "unknown PROC HAZARD option")
      next
    }
    mapped <- mapped + 1L
    switch(token,
      DATA        = data_name <- val,
      OUTHAZ      = outhaz <- val,
      MAXITER     = ctl$maxit <- as.numeric(val),
      CONDITION   = ctl$condition <- as.numeric(val),
      CONSERVE    = ctl$conserve <- TRUE,
      NOCONSERVE  = ctl$conserve <- FALSE,
      QUASINEWTON = ctl$method <- "bfgs",
      STEEPEST    = {
        mapped <- mapped - 1L
        note("STEEPEST", "no R equivalent for steepest descent (issue #145)")
      },
      # Print-only options. Mapped, because "no R effect" is the correct
      # translation, not a gap.
      PRINTIT = NULL, NOPRINT = NULL, NOCOR = NULL, NOCOV = NULL,
      NOLOG = NULL, NONOTES = NULL,
      {
        mapped <- mapped - 1L
        note(key, "no hazard() equivalent")
      }
    )
  }

  # --- statements 2..n ----------------------------------------------------
  statements <- list()
  parms_ops <- character(0)
  sel_ops <- NULL
  covars <- list()

  for (i in seq_along(st)[-1L]) {
    w <- strsplit(trimws(st[[i]]), " ", fixed = TRUE)[[1L]]
    w <- w[nzchar(w)]
    if (!length(w)) next
    kw <- w[[1L]]
    ops <- w[-1L]
    token <- .hzr_sas_token(kw, "HAZARD", "STMT")
    seen <- seen + 1L
    if (is.na(token)) {
      note(kw, "unknown HAZARD statement")
      next
    }
    mapped <- mapped + 1L
    switch(token,
      TIME       = statements$TIME <- ops[[1L]],
      EVENT      = statements$EVENT <- ops[[1L]],
      ICENSOR    = statements$ICENSOR <- ops,
      LCENSOR    = statements$LCENSOR <- ops[[1L]],
      RCENSOR    = statements$RCENSOR <- ops[[1L]],
      WEIGHT     = statements$WEIGHT <- ops[[1L]],
      PARAMETERS = parms_ops <- ops,
      STEPWISE   = sel_ops <- ops,
      EARLY      = covars$early <- ops,
      CONSTANT   = covars$constant <- ops,
      LATE       = covars$late <- ops,
      {
        mapped <- mapped - 1L
        note(kw, "no R equivalent")
      }
    )
  }

  # --- assemble -----------------------------------------------------------
  parms <- .hzr_parse_parms(parms_ops, covars = covars)
  untr <- rbind(untr, parms$untranslated)
  cens <- .hzr_censor_spec(statements)
  untr <- rbind(untr, cens$untranslated)

  args <- list()
  if (!is.null(data_name)) args$data <- as.name(data_name)
  args$time <- as.name(statements$TIME)
  args$status <- cens$status_expr
  if (!is.null(cens$time_lower)) args$time_lower <- cens$time_lower
  if (!is.null(cens$time_upper)) args$time_upper <- cens$time_upper
  args$phases <- parms$phases
  args$theta <- parms$theta
  if (!is.null(statements$WEIGHT)) args$weights <- as.name(statements$WEIGHT)

  # Canonical control order, so the emitted call does not depend on the order
  # the options happened to appear in the SAS text.
  ctl <- ctl[intersect(c("maxit", "condition", "conserve", "method"),
                       names(ctl))]
  if (length(ctl)) args$control <- as.call(c(quote(list), ctl))

  head <- quote(hazard)
  if (!is.null(sel_ops)) {
    sel <- .hzr_selection_spec(sel_ops)
    untr <- rbind(untr, sel$untranslated)
    if (!is.null(sel$direction)) {
      head <- quote(hzr_stepwise)
      args$direction <- sel$direction
      if (!is.null(sel$slentry)) args$slentry <- sel$slentry
      if (!is.null(sel$slstay)) args$slstay <- sel$slstay
    }
  }

  list(call = as.call(c(head, args)), outhaz = outhaz,
       untranslated = untr, tokens_seen = seen, tokens_mapped = mapped)
}
```

**Dependency note:** this calls `.hzr_selection_spec()` from Task 9. Implement
the `sel_ops` branch as written; the Task 9 tests exercise it, and until Task 9
lands `sel_ops` is always `NULL` for the Task 7 fixtures.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-parse-hazard")'`
Expected: PASS, 2 tests.

- [ ] **Step 5: Full gate, then commit**

```bash
Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && Rscript -e 'devtools::test()'
git add R/sas-parse-job.R tests/testthat/test-sas-parse-hazard.R
git commit -m "feat(sas): parse PROC HAZARD into a hazard() call"
```

---

### Task 8: Parse `PROC HAZPRED` and its prediction grid

**Files:**
- Modify: `R/sas-parse-job.R`
- Test: `tests/testthat/test-sas-parse-hazpred.R`

**Interfaces:**
- Produces: `.hzr_parse_hazpred(block, txt)` → `list(call = <call>, inhaz = <chr>, grid = <call|NULL>, untranslated = <df>, ...)`.
  `txt` is the whole normalised source, needed to find the grid `DATA` step.

Spec §5.4: 92% of grids are `log_grid` or `explicit_do`; a `derived_set` grid
emits `UNTRANSLATED`. `HAZPRED` emits `_SURVIV`/`_CLLSURV`/`_CLUSURV` and
`_HAZARD`/`_CLLHAZ`/`_CLUHAZ`, so the call needs `se.fit = TRUE`.

**Note:** `conf.type` is deliberately **not** set — spec §8 open question 2
records that whether SAS uses `log-log` or `logit` limits is unresolved, and
guessing produces silently wrong bounds. Emit the package default and add a
`.qmd` callout naming the uncertainty.

- [ ] **Step 1: Write the failing test**

```r
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
  expect_equal(got$grid,
               quote(data.frame(MONTHS = exp(seq(-5, log(180), length.out = 100)))))
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-parse-hazpred")'`
Expected: FAIL, `could not find function ".hzr_parse_hazpred"`.

- [ ] **Step 3: Implement the grid helper**

```r
#' Translate the DATA step that builds a HAZPRED prediction grid.
#'
#' Returns an unevaluated data.frame() call, or NULL when the step is not one
#' of the two stereotyped forms. NULL means untranslated, never "no grid":
#' a predict() with no newdata is a hollow result, so the caller must record it.
#' @noRd
.hzr_parse_grid <- function(txt, name) {
  if (is.null(name) || !nzchar(name)) return(NULL)
  p <- .idx(txt, paste0("DATA ", name, ";"))
  if (p == 0L) return(NULL)

  rest <- substring(txt, p + 1L)
  ends <- c(.idx(rest, "DATA "), .idx(rest, "PROC "), .idx(rest, "%HAZ"))
  ends <- ends[ends > 0L]
  body <- if (length(ends)) substring(rest, 1L, min(ends)) else rest

  time_var <- local({
    m <- regmatches(body, regexec("([A-Z_][A-Z0-9_]*) *= *EXP\\(", body))[[1L]]
    if (length(m) > 1L) m[[2L]] else NA_character_
  })

  # --- log-spaced grid: DO LN_TIME=-5 TO LN_MAX BY INC; t = EXP(LN_TIME) ----
  if (grepl(" DO ", body) && grepl("LOG(", body, fixed = TRUE) &&
      !is.na(time_var)) {
    lo <- local({
      m <- regmatches(body,
             regexec("DO [A-Z_][A-Z0-9_]* *= *(-?[0-9.]+) TO", body))[[1L]]
      if (length(m) > 1L) as.numeric(m[[2L]]) else NA_real_
    })
    hi <- local({
      m <- regmatches(body, regexec("MAX *= *([0-9.]+)", body))[[1L]]
      if (length(m) > 1L) as.numeric(m[[2L]]) else NA_real_
    })
    if (is.na(lo) || is.na(hi)) return(NULL)
    # SAS writes INC=(5+LN_MAX)/99.9, i.e. 100 points inclusive.
    inner <- bquote(exp(seq(.(lo), log(.(hi)), length.out = 100)))
    cl <- as.call(list(quote(data.frame), inner))
    names(cl) <- c("", time_var)
    return(cl)
  }

  # --- explicit DO list: DO MONTHS=1,2,3,6,12,24 TO 180 BY 12; --------------
  if (grepl(" DO ", body)) {
    m <- regmatches(body,
           regexec("DO ([A-Z_][A-Z0-9_]*) *= *([^;]+);", body))[[1L]]
    if (length(m) < 3L) return(NULL)
    var <- m[[2L]]
    parts <- trimws(strsplit(m[[3L]], ",", fixed = TRUE)[[1L]])
    elems <- list()
    for (part in parts) {
      rng <- regmatches(part,
               regexec("^([0-9.]+) +TO +([0-9.]+)( +BY +([0-9.]+))?$",
                       part))[[1L]]
      if (length(rng) >= 3L) {
        by <- if (length(rng) >= 5L && nzchar(rng[[5L]])) as.numeric(rng[[5L]]) else 1
        elems[[length(elems) + 1L]] <- bquote(
          seq(.(as.numeric(rng[[2L]])), .(as.numeric(rng[[3L]])), by = .(by))
        )
      } else if (grepl("^[0-9.]+$", part)) {
        elems[[length(elems) + 1L]] <- as.numeric(part)
      } else {
        # An element we cannot read (a SAS expression such as 1*DTY). Refuse
        # the whole grid rather than emit a partial one.
        return(NULL)
      }
    }
    inner <- as.call(c(quote(c), elems))
    cl <- as.call(list(quote(data.frame), inner))
    names(cl) <- c("", var)
    return(cl)
  }

  # SET-derived or anything else: not translatable.
  NULL
}
```

- [ ] **Step 3b: Implement the HAZPRED parser**

```r
#' Parse a PROC HAZPRED block into predict() call(s).
#'
#' HAZPRED emits _SURVIV/_CLLSURV/_CLUSURV and _HAZARD/_CLLHAZ/_CLUHAZ, so it
#' maps to two predict() calls with se.fit. `conf.type` is deliberately NOT set:
#' whether SAS uses log-log or logit limits is unresolved (spec section 8), and
#' guessing produces silently wrong bounds.
#' @noRd
.hzr_parse_hazpred <- function(block, txt) {
  st <- strsplit(block$text, ";", fixed = TRUE)[[1L]]
  untr <- .hzr_untranslated_frame()
  seen <- 0L
  mapped <- 0L
  note <- function(kw, reason) {
    untr <<- rbind(untr, .hzr_untranslated_frame(NA_integer_, kw, reason))
  }

  toks <- strsplit(trimws(st[[1L]]), " ", fixed = TRUE)[[1L]]
  toks <- toks[nzchar(toks)]
  data_name <- NULL
  inhaz <- NULL
  want_surv <- TRUE
  want_haz <- TRUE
  want_cl <- TRUE

  for (tok in toks) {
    eqp <- .idx(tok, "=")
    key <- if (eqp > 0L) substring(tok, 1L, eqp - 1L) else tok
    val <- if (eqp > 0L) substring(tok, eqp + 1L) else ""
    token <- .hzr_sas_token(key, "HAZPRED", "HZPP")
    if (identical(token, "PROC") || identical(token, "HAZPRED")) next
    seen <- seen + 1L
    if (is.na(token)) {
      note(key, "unknown PROC HAZPRED option")
      next
    }
    mapped <- mapped + 1L
    switch(token,
      DATA    = data_name <- val,
      INHAZ   = inhaz <- val,
      OUT     = NULL,
      NOSURV  = want_surv <- FALSE,
      NOHAZ   = want_haz <- FALSE,
      NOCL    = want_cl <- FALSE,
      CLIMITS = want_cl <- TRUE,
      NOLOG = NULL, NONOTES = NULL,
      {
        mapped <- mapped - 1L
        note(key, "no predict() equivalent")
      }
    )
  }

  for (i in seq_along(st)[-1L]) {
    w <- strsplit(trimws(st[[i]]), " ", fixed = TRUE)[[1L]]
    w <- w[nzchar(w)]
    if (!length(w)) next
    kw <- w[[1L]]
    token <- .hzr_sas_token(kw, "HAZPRED", "STMT")
    seen <- seen + 1L
    if (is.na(token)) {
      note(kw, "unknown HAZPRED statement")
      next
    }
    mapped <- mapped + 1L
    if (!(token %in% c("TIME", "ID"))) {
      mapped <- mapped - 1L
      note(kw, "SAS listing control; no R effect")
    }
  }

  grid <- .hzr_parse_grid(txt, data_name)
  if (is.null(grid) && !is.null(data_name)) {
    note(paste0("DATA=", data_name),
         "prediction grid DATA step is not one of the translatable forms")
  }

  mk <- function(type) {
    args <- list(quote(fit))
    if (!is.null(data_name)) args$newdata <- as.name(data_name)
    args$type <- type
    args$se.fit <- want_cl
    as.call(c(quote(predict), args))
  }

  list(
    call = if (want_surv) mk("survival") else mk("hazard"),
    call_haz = if (want_surv && want_haz) mk("hazard") else NULL,
    inhaz = inhaz, grid = grid, untranslated = untr,
    tokens_seen = seen, tokens_mapped = mapped
  )
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-parse-hazpred")'`
Expected: PASS, 3 tests.

- [ ] **Step 5: Full gate, then commit**

```bash
Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && Rscript -e 'devtools::test()'
git add R/sas-parse-job.R tests/testthat/test-sas-parse-hazpred.R
git commit -m "feat(sas): parse PROC HAZPRED and its prediction grid"
```

---

### Task 9: `SELECTION` forks to `hzr_stepwise()`

Read spec §5.3 before starting. This task exists because of a mapping that is
wrong in the obvious direction.

**Files:**
- Modify: `R/sas-parse-job.R`
- Test: `tests/testthat/test-sas-selection.R`

**Interfaces:**
- Consumes: `.hzr_sas_token()`, `.hzr_parse_parms()`.
- Produces: `.hzr_selection_spec(operands)` →
  `list(direction = <chr|NULL>, slentry = <num|NULL>, slstay = <num|NULL>, untranslated = <df>)`.
  `.hzr_parse_hazard()` gains: when a `SELECTION` statement is present, the
  emitted call is `hzr_stepwise(...)` rather than `hazard(...)`.

**`SELECTION FORWARD` maps to `direction = "both"`, not `"forward"`.** In
`hazard_y.y`, `STEPWISE` sets option 21, `BACKWARD` 22, `ONEWAY` 34 — and
`FORWARD`, `FW`, `SW`, `SELECT` and `STEPWISE` all lex to the same `STEPWISE`
token, so HAZARD does not distinguish them. Option 21 is two-way. Corroborated
by `slstay`, which would be meaningless under a forward-only search.

Phase statements carry up to **320** covariates in the `lv_function` study;
those are candidate pools, not fitted covariates.

- [ ] **Step 1: Write the failing test**

```r
test_that("SELECTION FORWARD maps to direction = both, not forward", {
  # The obvious mapping is wrong and fails silently on every stepwise job.
  got <- .hzr_selection_spec(c("FORWARD", "SLENTRY=0.05", "SLSTAY=0.10"))
  expect_equal(got$direction, "both")
  expect_equal(got$slentry, 0.05)
  expect_equal(got$slstay, 0.10)
})

test_that("every forward/stepwise spelling reaches the same direction", {
  for (kw in c("FORWARD", "FW", "SW", "SELECT", "STEPWISE")) {
    expect_equal(.hzr_selection_spec(kw)$direction, "both",
                 info = paste("spelling:", kw))
  }
})

test_that("BACKWARD and ONEWAY are distinct from stepwise", {
  expect_equal(.hzr_selection_spec("BACKWARD")$direction, "backward")
  expect_equal(.hzr_selection_spec("BW")$direction, "backward")
  # ONEWAY means no stepwise at all: the caller emits a plain hazard() fit.
  expect_null(.hzr_selection_spec("NOSTEPWISE")$direction)
})

test_that("a SELECTION block emits hzr_stepwise, not hazard", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT D; TIME T; EARLY X1 X2 X3;",
    "PARMS MUE=1 THALF=1 NU=1;",
    "SELECTION FORWARD SLENTRY=0.05 SLSTAY=0.1; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(got$call[[1L]], as.name("hzr_stepwise"))
  expect_equal(got$call[["direction"]], "both")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-selection")'`
Expected: FAIL, `could not find function ".hzr_selection_spec"`.

- [ ] **Step 3: Implement**

```r
#' Translate a SELECTION statement to hzr_stepwise() arguments.
#'
#' HAZARD's lexer collapses FORWARD, FW, SW, SELECT and STEPWISE into one
#' token, which hazard_y.y maps to option 21; BACKWARD is option 22 and ONEWAY
#' is 34. Option 21 is therefore two-way, so FORWARD maps to "both". Mapping it
#' to "forward" would change the search on every stepwise job with no error.
#' @noRd
.hzr_selection_spec <- function(operands) {
  out <- list(direction = NULL, slentry = NULL, slstay = NULL,
              untranslated = .hzr_untranslated_frame())
  for (op in operands) {
    eqp <- .idx(op, "=")
    key <- if (eqp > 0L) substring(op, 1L, eqp - 1L) else op
    val <- if (eqp > 0L) as.numeric(substring(op, eqp + 1L)) else NA_real_
    token <- .hzr_sas_token(key, "HAZARD", "STEP")
    if (is.na(token)) {
      out$untranslated <- rbind(out$untranslated, .hzr_untranslated_frame(
        NA_integer_, key, "unknown SELECTION option"
      ))
      next
    }
    switch(token,
      STEPWISE = { out$direction <- "both" },
      BACKWARD = { out$direction <- "backward" },
      ONEWAY   = { out$direction <- NULL },
      SLENTRY  = { out$slentry <- val },
      SLSTAY   = { out$slstay <- val },
      {
        out$untranslated <- rbind(out$untranslated, .hzr_untranslated_frame(
          NA_integer_, key, "no hzr_stepwise() equivalent"
        ))
      }
    )
  }
  out
}
```

In `.hzr_parse_hazard()`, after the statement loop: if a `SELECTION` statement
was seen **and** `spec$direction` is non-`NULL`, build the call with
`quote(hzr_stepwise)` as the head instead of `quote(hazard)`, appending
`direction`, `slentry` and `slstay`. Phase-statement covariates become the
candidate pool rather than fixed model terms.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-selection")'`
Expected: PASS, 8 tests.

- [ ] **Step 5: Full gate, then commit**

```bash
Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && Rscript -e 'devtools::test()'
git add R/sas-parse-job.R tests/testthat/test-sas-selection.R
git commit -m "feat(sas): fork SELECTION jobs to hzr_stepwise() with direction=both"
```

---
### Task 10: Renderer

**Files:**
- Create: `R/sas-render-qmd.R`
- Test: `tests/testthat/test-sas-render-qmd.R`

**Interfaces:**
- Consumes: `hzr_sas_job`.
- Produces: `.hzr_render_qmd(job)` → character vector of `.qmd` lines.

Because calls are stored unevaluated, rendering is `deparse()`, not string
templating.

- [ ] **Step 1: Write the failing test**

```r
test_that("calls render as R chunks", {
  j <- .hzr_sas_job(
    source = list(path = "hz.death.AVC.sas", checksum = "abc"),
    calls = list(fit = quote(hazard(time = T, status = D))),
    grid = NULL, inhaz = NULL, outhaz = "EX.HZD",
    untranslated = .hzr_untranslated_frame(),
    coverage = list(tokens_seen = 6L, tokens_mapped = 6L)
  )
  out <- .hzr_render_qmd(j)
  expect_true(any(grepl("^```\\{r\\}", out)))
  expect_true(any(grepl("hazard(time = T, status = D)", out, fixed = TRUE)))
})

test_that("an untranslated construct becomes a visible callout", {
  j <- .hzr_sas_job(
    source = list(path = "x.sas", checksum = "abc"),
    calls = list(fit = quote(hazard(time = T, status = D))),
    grid = NULL, inhaz = NULL, outhaz = NULL,
    untranslated = .hzr_untranslated_frame(12L, "STEEPEST",
                                           "no R equivalent (see #145)"),
    coverage = list(tokens_seen = 7L, tokens_mapped = 6L)
  )
  out <- .hzr_render_qmd(j)
  expect_true(any(grepl("UNTRANSLATED", out)))
  expect_true(any(grepl("STEEPEST", out)))
})

test_that("an unresolved INHAZ renders a stop(), not a comment", {
  # The document must FAIL to render rather than produce a report over a model
  # it never loaded.
  j <- .hzr_sas_job(
    source = list(path = "hp.sas", checksum = "abc"),
    calls = list(pred = quote(predict(fit, newdata = PREDICT))),
    grid = NULL, inhaz = "CABGKUL.HMDEADP", outhaz = NULL,
    untranslated = .hzr_untranslated_frame(),
    coverage = list(tokens_seen = 4L, tokens_mapped = 4L)
  )
  out <- .hzr_render_qmd(j)
  expect_true(any(grepl("stop(", out, fixed = TRUE)))
  expect_true(any(grepl("CABGKUL.HMDEADP", out, fixed = TRUE)))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-render-qmd")'`
Expected: FAIL, `could not find function ".hzr_render_qmd"`.

- [ ] **Step 3: Implement**

```r
#' Render an hzr_sas_job as Quarto source.
#'
#' Calls are stored unevaluated, so rendering is deparse() rather than string
#' templating. An unresolved INHAZ emits stop(), not a comment: the document
#' must fail to render rather than produce a report over a model it never
#' loaded.
#' @noRd
.hzr_render_qmd <- function(job) {
  out <- character(0)
  add <- function(...) out <<- c(out, ...)
  base <- basename(job$source$path)

  add("---",
      sprintf('title: "%s"', sub("[.]sas$", "", base)),
      "format: html",
      "---",
      "")
  add(sprintf("<!-- Translated from %s (md5 %s) by hzr_translate_sas(). -->",
              base, substr(job$source$checksum, 1L, 12L)),
      "")
  add("```{r}", "#| label: setup", "library(TemporalHazard)", "```", "")

  if (!is.null(job$inhaz) && !isTRUE(job$inhaz_resolved)) {
    lib <- sub("[.].*$", "", job$inhaz)
    add("```{r}",
        "#| label: inhaz-unresolved",
        sprintf(
          "stop('unresolved INHAZ=%s -- pass librefs = c(%s = \"<path>\") to hzr_translate_sas()')",
          job$inhaz, lib
        ),
        "```", "")
  }

  for (nm in names(job$calls)) {
    add("```{r}",
        sprintf("#| label: %s", nm),
        deparse(job$calls[[nm]], width.cutoff = 60L),
        "```", "")
  }

  for (i in seq_len(nrow(job$untranslated))) {
    add("::: {.callout-warning}",
        sprintf("## UNTRANSLATED: %s", job$untranslated$construct[i]),
        job$untranslated$reason[i],
        ":::", "")
  }

  add(sprintf("<!-- coverage: %d/%d tokens mapped -->",
              job$coverage$tokens_mapped, job$coverage$tokens_seen))
  out
}
```

`job$inhaz_resolved` is set by `hzr_translate_sas()` in Task 11; on a job that
lacks it, `job$inhaz_resolved` is `NULL` and `!isTRUE(NULL)` is `TRUE`, so an
unresolved `INHAZ` fails closed.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-render-qmd")'`
Expected: PASS, 3 tests.

- [ ] **Step 5: Full gate, then commit**

```bash
Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && Rscript -e 'devtools::test()'
git add R/sas-render-qmd.R tests/testthat/test-sas-render-qmd.R
git commit -m "feat(sas): render an hzr_sas_job to Quarto"
```

---

### Task 11: `hzr_translate_sas()` and `INHAZ` resolution

The only export. Read spec §4 and §4.1.

**Files:**
- Create: `R/translate-sas.R`
- Test: `tests/testthat/test-translate-sas.R`

**Interfaces:**
- Consumes: everything above.
- Produces: exported `hzr_translate_sas(path, out_dir = NULL, librefs = NULL)`
  → `hzr_sas_job` invisibly, writing `<out_dir>/<basename>.qmd`.

`librefs` is a **named character vector**, `c(EX = "estimates/")`. Never a path
to a config file: that would put `yaml` in the dependency graph, and `Imports:`
is deliberately `survival` alone.

- [ ] **Step 1: Write the failing test**

```r
test_that("a job with no recognisable block errors and writes nothing", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(c("DATA X; SET Y;", "PROC PRINT;"), f)
  out <- withr::local_tempdir()
  expect_error(hzr_translate_sas(f, out_dir = out), "no SAS statements|no HAZARD")
  expect_length(list.files(out, pattern = "[.]qmd$"), 0L)
})

test_that("librefs resolves an INHAZ that no translated job wrote", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.HZD OUT=P; TIME MONTHS; );"
  ), f)
  out <- withr::local_tempdir()
  job <- hzr_translate_sas(f, out_dir = out, librefs = c(EX = "estimates"))
  qmd <- readLines(file.path(out, sub("[.]sas$", ".qmd", basename(f))))
  expect_true(any(grepl("hzr_read_outhaz", qmd, fixed = TRUE)))
  expect_false(any(grepl("stop(", qmd, fixed = TRUE)))
})

test_that("without librefs an unresolved INHAZ warns and emits stop()", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines("%HAZPRED( PROC HAZPRED DATA=P INHAZ=EX.HZD OUT=P; TIME MONTHS; );", f)
  out <- withr::local_tempdir()
  expect_warning(hzr_translate_sas(f, out_dir = out), "unresolved")
  qmd <- readLines(file.path(out, sub("[.]sas$", ".qmd", basename(f))))
  expect_true(any(grepl("stop(", qmd, fixed = TRUE)))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="translate-sas")'`
Expected: FAIL, `could not find function "hzr_translate_sas"`.

- [ ] **Step 3: Implement**

```r
#' Translate a SAS HAZARD job into a Quarto document
#'
#' Reads a SAS program containing `PROC HAZARD` and/or `PROC HAZPRED` blocks and
#' emits a Quarto document that reproduces the analysis with [hazard()] and
#' [predict.hazard()].
#'
#' Constructs the translator does not cover are recorded on the returned object
#' and rendered as visible callouts, never dropped. A `PROC HAZPRED` job whose
#' `INHAZ=` fitted model cannot be located emits a `stop()`, so the document
#' fails to render rather than reporting over a model it did not load.
#'
#' @param path Path to a `.sas` file.
#' @param out_dir Directory to write the `.qmd` into. `NULL` (default) parses
#'   without writing.
#' @param librefs Optional named character vector mapping SAS librefs to
#'   directories, e.g. `c(EX = "estimates")`, used to resolve `INHAZ=`.
#' @return An `hzr_sas_job` object, invisibly.
#' @examples
#' \donttest{
#' job <- hzr_translate_sas(
#'   system.file("extdata", "hz-example.sas", package = "TemporalHazard")
#' )
#' }
#' @export
hzr_translate_sas <- function(path, out_dir = NULL, librefs = NULL) {
  stopifnot(is.character(path), length(path) == 1L)
  if (!file.exists(path)) stop("no such file: ", path, call. = FALSE)
  if (!is.null(librefs) &&
      (!is.character(librefs) || is.null(names(librefs)))) {
    stop('librefs must be a named character vector, e.g. c(EX = "estimates").',
         call. = FALSE)
  }

  txt <- .hzr_sas_normalise(readLines(path, warn = FALSE))
  blocks <- .hzr_sas_blocks(txt)
  if (!length(blocks)) {
    stop("no HAZARD or HAZPRED block found in ", path, call. = FALSE)
  }

  calls <- list()
  untr <- .hzr_untranslated_frame()
  seen <- 0L
  mapped <- 0L
  inhaz <- NULL
  outhaz <- NULL
  grid <- NULL

  for (b in blocks) {
    if (identical(b$proc, "HAZARD")) {
      r <- .hzr_parse_hazard(b)
      calls$fit <- r$call
      outhaz <- r$outhaz
    } else {
      r <- .hzr_parse_hazpred(b, txt)
      inhaz <- r$inhaz
      grid <- r$grid
      calls$pred <- r$call
      if (!is.null(r$call_haz)) calls$pred_haz <- r$call_haz
    }
    untr <- rbind(untr, r$untranslated)
    seen <- seen + r$tokens_seen
    mapped <- mapped + r$tokens_mapped
  }

  # The grid is an assignment, placed before the predict() chunks that use it.
  if (!is.null(grid) && !is.null(calls$pred[["newdata"]])) {
    nm <- as.character(calls$pred[["newdata"]])
    calls <- c(list(grid = call("<-", as.name(nm), grid)), calls)
  }

  # --- INHAZ resolution: this job's own OUTHAZ, then librefs, then fail ----
  resolved <- FALSE
  if (!is.null(inhaz)) {
    if (!is.null(outhaz) && identical(inhaz, outhaz)) {
      resolved <- TRUE
    } else if (!is.null(librefs)) {
      lib <- sub("[.].*$", "", inhaz)
      mem <- tolower(sub("^[^.]*[.]", "", inhaz))
      if (lib %in% names(librefs)) {
        calls <- c(
          list(fit = call("<-", as.name("fit"),
                          bquote(hzr_read_outhaz(
                            file.path(.(unname(librefs[[lib]])), .(mem)))))),
          calls
        )
        resolved <- TRUE
      }
    }
  }

  job <- .hzr_sas_job(
    source = list(path = path,
                  checksum = unname(tools::md5sum(path))),
    calls = calls, grid = grid, inhaz = inhaz, outhaz = outhaz,
    untranslated = untr,
    coverage = list(tokens_seen = seen, tokens_mapped = mapped)
  )
  job$inhaz_resolved <- resolved
  .hzr_validate_sas_job(job)

  if (!is.null(inhaz) && !resolved) {
    warning("unresolved INHAZ=", inhaz, " in ", basename(path),
            "; the emitted document will stop() rather than render.",
            call. = FALSE)
  }
  if (nrow(untr)) {
    warning(nrow(untr), " untranslated construct(s) in ", basename(path), ": ",
            paste(unique(untr$construct), collapse = ", "), call. = FALSE)
  }

  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    writeLines(
      .hzr_render_qmd(job),
      file.path(out_dir, sub("[.]sas$", ".qmd", basename(path)))
    )
  }
  invisible(job)
}
```

`tools::md5sum()` is base R, so this adds no dependency. The `@examples` block
references `inst/extdata/hz-example.sas`; create it in this task by copying a
de-identified `PROC HAZARD` block from the public `examples/` corpus.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="translate-sas")'`
Expected: PASS, 3 tests.

- [ ] **Step 5: Full gate, then commit**

```bash
Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && Rscript -e 'devtools::test()'
git add R/translate-sas.R tests/testthat/test-translate-sas.R NAMESPACE man
git commit -m "feat(sas): add hzr_translate_sas()"
```

---

### Task 12: End-to-end corpus smoke test

**Files:**
- Test: `tests/testthat/test-sas-translate-corpus.R`

Guards the property the whole design rests on: over a real corpus the translator
either produces a job or fails loudly, and never returns a hollow one.

- [ ] **Step 1: Write the test**

```r
test_that("every public-corpus job either translates or fails loudly", {
  skip_on_cran()
  repo <- Sys.getenv("HAZARD_REPO", "~/Documents/GitHub/hazard")
  skip_if_not(dir.exists(path.expand(repo)), "hazard checkout not available")

  fs <- list.files(path.expand(repo), pattern = "[.]sas$",
                   recursive = TRUE, full.names = TRUE)
  skip_if(length(fs) == 0L, "no .sas files found")

  out <- withr::local_tempdir()
  n_ok <- 0L
  for (f in fs) {
    job <- tryCatch(
      suppressWarnings(hzr_translate_sas(f, out_dir = out)),
      error = function(e) NULL
    )
    if (is.null(job)) next
    # A translated job must have seen something. A hollow job is the failure
    # this whole design exists to prevent.
    expect_gt(job$coverage$tokens_seen, 0L)
    expect_gte(job$coverage$tokens_seen, job$coverage$tokens_mapped)
    n_ok <- n_ok + 1L
  }
  # Warn loudly if nothing translated: a pass over zero jobs is not a pass.
  expect_gt(n_ok, 20L)
})
```

- [ ] **Step 2: Run it**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="sas-translate-corpus")'`
Expected: PASS, with `n_ok` above 20. If it passes with `n_ok` at exactly the
threshold, investigate — the corpus has 59 files with HAZ blocks.

- [ ] **Step 3: Full gate, then commit**

```bash
Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && Rscript -e 'devtools::test()'
git add tests/testthat/test-sas-translate-corpus.R
git commit -m "test(sas): corpus-wide smoke test for the translator"
```

---

### Task 13: G3 Weibull correspondence fixture

Spec §7.1. The CoE half of this task was **split out to
[#146](https://github.com/ehrlinger/temporal_hazard/issues/146)** so this PR does
not touch `R/likelihood-multiphase.R` — AGENTS.md calls the multiphase path the
place the bugs are, and a translator PR has no business editing it.

**Files:**
- Test: `tests/testthat/test-g3-weibull-correspondence.R`

The identity was checked by hand on 2026-08-20 and holds (8.6e-16 at moderate
times, 5.0e-13 across 1e-6 to 1e6, hazard 1.9e-16). This task makes it a
permanent test so it cannot regress unnoticed.

- [ ] **Step 1: Write the test**

```r
test_that("a G3 constrained at alpha = eta = 1 is a Weibull cumulative hazard", {
  # PARMS ... WEIBULL is setopt(6) -> SETG3_weibull, g3flag += 2 in the
  # reference implementation. The R general form
  #   (((t/tau)^gamma + 1)^(1/alpha) - 1)^eta
  # collapses at alpha = eta = 1 to (t/tau)^gamma. This is the most common
  # production configuration: 76-88 blocks in every study profiled.
  t <- c(0.01, 0.5, 1, 5, 50, 180)
  tau <- 2.5
  gamma <- 1.4
  got <- hzr_decompos_g3(t, tau = tau, gamma = gamma, alpha = 1, eta = 1)
  expect_equal(got$G3, (t / tau)^gamma, tolerance = 1e-12)
})

test_that("the constrained hazard matches the Weibull hazard", {
  t <- c(0.5, 5, 50)
  tau <- 2.5
  gamma <- 1.4
  got <- hzr_decompos_g3(t, tau = tau, gamma = gamma, alpha = 1, eta = 1)
  expect_equal(got$g3, (gamma / tau) * (t / tau)^(gamma - 1),
               tolerance = 1e-12)
})

test_that("the constrained path survives extreme times", {
  # HAZARD branches on g3flag partly for numerical reasons: at the constraint
  # the general path computes expm1(log1p(exp(y))) where the answer is exp(y).
  # Measured worst case is 5e-13 relative, so the tolerance is set just wide
  # enough to admit it and no wider -- a looser one would stop detecting drift.
  t <- c(1e-6, 1e-3, 1e3, 1e6)
  got <- hzr_decompos_g3(t, tau = 1, gamma = 2, alpha = 1, eta = 1)
  expect_true(all(is.finite(got$G3)))
  expect_equal(got$G3, t^2, tolerance = 1e-11)
})

test_that("the comparison is not trivially satisfied", {
  # A max discrepancy of exactly zero across every case would more likely mean
  # nothing was compared than that the identity is exact. Assert the general
  # form differs from the Weibull form away from the constraint.
  t <- c(0.5, 5, 50)
  general <- hzr_decompos_g3(t, tau = 2.5, gamma = 1.4, alpha = 2, eta = 1.5)
  expect_false(isTRUE(all.equal(general$G3, (t / 2.5)^1.4)))
})
```

- [ ] **Step 2: Run it**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(filter="g3-weibull")'`
Expected: PASS, 4 tests.

**If any of the first three fail, stop and report.** That would mean the Weibull
late phase is *not* reachable as a constrained G3, turning spec §7.1 from a
verification gap into a feature gap and changing v1 scope.

- [ ] **Step 3: Full gate, then commit**

```bash
Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && Rscript -e 'devtools::test()'
git add tests/testthat/test-g3-weibull-correspondence.R
git commit -m "test: pin the G3 Weibull correspondence at alpha = eta = 1"
```

**Note:** this verifies the identity *inside R*. It does not compare against SAS.
The SAS parity leg still needs a reference run, so §7.1 stays open.

---

### Task 14: Documentation and PR

**Files:**
- Modify: `NEWS.md`, `DESCRIPTION` (patch digit only), `_pkgdown.yml`
- Modify: `vignettes/sas-to-r-migration.qmd`

- [ ] **Step 1: Bump the patch digit**

Edit `DESCRIPTION` `Version:` — patch digit only. **Never** roll MINOR or MAJOR;
that is the maintainer's call. Plain three digits, no `.9000`, no fourth digit.

- [ ] **Step 2: Update NEWS.md by hand**

There is no test here that greps `NEWS.md` for the `DESCRIPTION` version, so
nothing will catch a mismatch. Keep the top heading and `DESCRIPTION` in sync
manually.

- [ ] **Step 3: Add `hzr_translate_sas()` to `_pkgdown.yml`**

- [ ] **Step 4: Add a translator section to the migration vignette**

Written for persona (c), the bilingual SAS-to-R migrant, per `.claude/house-style.md`.
Watch the vignette budget: rebuild is already 50s of a 3m44s check against
CRAN's ~10 minute ceiling. Keep new chunks cheap or `eval: false`.

- [ ] **Step 5: Full release gate**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
git commit -am "docs(sas): document the SAS job translator"
mkdir -p "$TMPDIR/tree"
git archive HEAD | tar -x -C "$TMPDIR/tree"
R CMD build "$TMPDIR/tree"
R CMD check --as-cran TemporalHazard_*.tar.gz
```

Commit **before** `git archive` — it exports the committed tree, so an
uncommitted fix is silently absent. Check **with** the manual. Confirm developer
files stayed out:

```bash
tar tzf TemporalHazard_*.tar.gz | grep -iE 'CLAUDE|AGENTS|data-raw' || echo "clean"
```

- [ ] **Step 6: Run the r-reviewer agent over the diff, then open the PR**

`hzr_translate_sas()` is a new export, so this is required. Verify each finding
against the code rather than acting on it, and do not let a clean report stand
in for the commands above.

```bash
git push -u origin feat/sas-job-translator
gh pr create --base dev --title "feat(sas): SAS PROC HAZARD/HAZPRED job translator"
```

Then stop. The maintainer merges. `Closes #NNN` does not fire on merges to
`dev` — close issues by hand.
