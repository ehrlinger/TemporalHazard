# What a Fit Did Not Do — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Every `hazard` object records what it did not do (`$degraded`) and
why (`$degraded_causes`). `print()` and `summary()` always show that record as
a "Not done in this run" block, "none" included, and a validator keeps the
record, the object and the printed text in agreement.

**Architecture:** Reasons are carried up as fields from where they arise:
`.hzr_safe_solve()` → `.hzr_optim_generic()` → `.hzr_optim_multiphase()`, and
the weak-direction check through a new `_impl`. The object's own state decides
*which* entries appear; the reasons decide only *why*. One pure builder stamps
the record at the two places a `hazard` object is built, a validator `stop()`s
if the record disagrees with the object, and one formatter feeds both printers.
An output checker parses the printed text back into entries and compares them
with the raw record, never through the formatter.

**Tech Stack:** Base R; testthat 3e (`local_mocked_bindings()`, installed
3.3.2); roxygen2; lintr.

**Spec:** `inst/dev/DEGRADED-RECORD-DESIGN.md` (including its Amendments
section). **Issue:** #242. **Branch:** `feat/degraded-record` (already cut from
`origin/main` at `84931a7`, upstream unset).

> Note (2026-09-10): the `"no starting values (theta = NULL)..."` cause
> quoted in Tasks 1, 4 and 8 was dropped after #243 made that call `stop()`;
> see spec Amendment 1.

## Global Constraints

- Base R only. **No new dependency** (numDeriv stays a Suggests).
- Version stays **1.2.10**: no bump. `DESCRIPTION` and the `NEWS.md` top heading must both still read 1.2.10. PR #241 also folds into 1.2.10, so expect a `NEWS.md` rebase conflict with it.
- Capability names and canonical order, verbatim: `fitting`, `standard_errors`, `conserved_phase_variance`, `weak_direction_check`, `conservation_of_events`.
- Cause strings are exactly those in the spec's "Entries and causes" table. Tests hard-code them from the spec, never from the implementation's constants.
- Heading text, verbatim: `Not done in this run:`. Clean: `  Not done in this run: none`. Legacy: `  Not done in this run: not recorded (object built before this record existed)`. Entries: four spaces, then `name: cause`, one per line, never `strwrap()`ped.
- `fit$fit$weak` keeps its list / `NULL` / `NA` meaning. `.hzr_weak_direction()`'s return value must not change.
- `stats::` prefix on stats functions; no `print()`, `browser()` or `library()` in `R/`; `cat()` only inside `print.*` methods; lines ≤ 120 characters (`.lintr`).
- No PHI anywhere.
- Never push to `main`. Branch, PR, stop; the maintainer merges.
- Run tests with `NOT_CRAN=true`, or `skip_on_cran()` tests are silently skipped.
- Gate order, every task: `devtools::document()` → `lintr::lint_package()` (0 lints) → the task's tests. The full gate (plus `spelling`, plus `R CMD check`) runs in Task 8.

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `R/degraded.R` | create | constants, builder, validator, formatter, output checker |
| `R/hessian-invert.R` | modify | `.hzr_safe_solve()` gains `reason`; `.hzr_weak_direction_impl()` + wrapper |
| `R/optimizer.R` | modify | `.hzr_numderiv_available()`; `se_unavailable_reason` in the result |
| `R/likelihood-multiphase.R` | modify | `conserved_variance_reason`; reset the SE reason when the recompute succeeds |
| `R/hazard_api.R` | modify | `fit_ran`, reasons, stamping, validation; `summary()` carries the record; both printers |
| `R/read-outhaz.R` | modify | stamp the SAS import |
| `tests/testthat/test-degraded.R` | create | builder, validator, formatter, output checker (pure) |
| `tests/testthat/test-degraded-reasons.R` | create | reasons from the Hessian path and the weak check |
| `tests/testthat/test-degraded-fits.R` | create | the record on real fits and imports |
| `tests/testthat/test-not-done-print.R` | create | printed block, legacy objects, falsifiability |
| `tests/testthat/test-weak-direction.R` | modify | the removed NA note → the block |
| `tests/testthat/test-hessian-stability.R` | modify | the hand-built object gets a record |
| `NEWS.md` | modify | one bullet under 1.2.10 |

---

### Task 1: The record: builder, validator, formatter, output checker

**Files:**
- Create: `R/degraded.R`
- Test: `tests/testthat/test-degraded.R`

**Interfaces:**
- Consumes: `%||%` (already defined in `R/stepwise-step.R`).
- Produces:
  - `.hzr_capabilities` (character, 5)
  - `.hzr_cause_not_recorded` (`"cause not recorded"`)
  - `.hzr_coe_causes` (named character)
  - `.hzr_is_na_scalar(x)` → logical(1)
  - `.hzr_degraded_record(vcov, weak, control, dist, fitted = TRUE, imported = FALSE, not_fitted_cause = "not requested (fit = FALSE)", reasons = list())` → `list(degraded = <chr>, degraded_causes = <named chr>)`. `reasons` may carry `se`, `weak` and `conserved_variance`, each a character(1) or `NA`.
  - `.hzr_validate_degraded(object, fitted, imported = FALSE)` → `invisible(TRUE)` or `stop()`
  - `.hzr_format_not_done(degraded, causes)` → character lines
  - `.hzr_check_not_done_output(lines, object)` → character problems; `character(0)` means clean

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-degraded.R`:

```r
# The record of what a fit did not do (#242, following #197): builder,
# validator, formatter and output check. Every expected cause string below is
# copied from inst/dev/DEGRADED-RECORD-DESIGN.md, never from the package's own
# constants -- a test that reads the table it is checking cannot fail.

vc <- diag(2)

fake_fit <- function(degraded = character(0), causes = character(0),
                     vcov = vc, weak = NULL, dist = "weibull",
                     control = list()) {
  structure(list(spec = list(dist = dist, control = control),
                 fit = list(vcov = vcov, weak = weak),
                 degraded = degraded, degraded_causes = causes),
            class = "hazard")
}

record <- function(...) {
  .hzr_degraded_record(control = list(), dist = "weibull", ...)
}

test_that("a clean single-distribution fit records nothing", {
  rec <- record(vcov = vc, weak = NULL)
  expect_identical(rec$degraded, character(0))
  expect_length(rec$degraded_causes, 0L)
})

test_that("each capability is recorded with the cause it was given", {
  cases <- list(
    list(args = list(vcov = NULL, weak = NA, fitted = FALSE),
         want = c(fitting = "not requested (fit = FALSE)",
                  standard_errors = "model not fitted",
                  weak_direction_check = "model not fitted")),
    list(args = list(vcov = NULL, weak = NA, fitted = FALSE,
                     not_fitted_cause = paste0("no starting values (theta = NULL); ",
                                               "the optimizer did not run")),
         want = c(fitting = "no starting values (theta = NULL); the optimizer did not run",
                  standard_errors = "model not fitted",
                  weak_direction_check = "model not fitted")),
    list(args = list(vcov = vc, weak = NA, imported = TRUE),
         want = c(fitting = "imported from SAS output; not fitted in R",
                  weak_direction_check = "imported from SAS output; no R Hessian")),
    list(args = list(vcov = NULL, weak = NA, imported = TRUE),
         want = c(fitting = "imported from SAS output; not fitted in R",
                  standard_errors = "covariance not imported",
                  weak_direction_check = "imported from SAS output; no R Hessian")),
    list(args = list(vcov = NA, weak = NA,
                     reasons = list(se = "Hessian not invertible",
                                    weak = "standard errors unavailable")),
         want = c(standard_errors = "Hessian not invertible",
                  weak_direction_check = "standard errors unavailable")),
    # State forces the entries; a missing reason must not remove them.
    list(args = list(vcov = NA, weak = NA),
         want = c(standard_errors = "cause not recorded",
                  weak_direction_check = "cause not recorded"))
  )
  for (i in seq_along(cases)) {
    rec <- do.call(record, cases[[i]]$args)
    expect_identical(rec$degraded, names(cases[[i]]$want), info = i)
    expect_identical(rec$degraded_causes, cases[[i]]$want, info = i)
  }
})

test_that("every conserve_disabled_reason becomes the spec's cause", {
  want <- c(
    not_requested         = "not requested (conserve = FALSE)",
    unsupported_censoring = "status outside {0, 1} (left or interval censoring)",
    single_phase          = "fewer than two phases",
    no_events             = "no events",
    setup_failed          = "setup did not complete"
  )
  for (r in names(want)) {
    rec <- .hzr_degraded_record(
      vcov = vc, weak = NULL, dist = "multiphase",
      control = list(conserve_applied = FALSE, conserve_disabled_reason = r))
    expect_identical(rec$degraded_causes,
                     c(conservation_of_events = want[[r]]), info = r)
  }
  # A reason the optimizer adds later still produces an entry.
  rec <- .hzr_degraded_record(
    vcov = vc, weak = NULL, dist = "multiphase",
    control = list(conserve_applied = FALSE, conserve_disabled_reason = "new"))
  expect_identical(rec$degraded_causes[["conservation_of_events"]],
                   "cause not recorded")
  # Applied CoE, and CoE on a non-multiphase fit, record nothing.
  expect_length(.hzr_degraded_record(
    vcov = vc, weak = NULL, dist = "multiphase",
    control = list(conserve_applied = TRUE))$degraded, 0L)
  expect_length(.hzr_degraded_record(
    vcov = vc, weak = NULL, dist = "weibull",
    control = list(conserve_applied = FALSE))$degraded, 0L)
})

test_that("the conserved phase's variance is recorded only when CoE ran", {
  rec <- .hzr_degraded_record(
    vcov = vc, weak = NULL, dist = "multiphase",
    control = list(conserve_applied = TRUE),
    reasons = list(conserved_variance = "numDeriv not installed"))
  expect_identical(rec$degraded_causes,
                   c(conserved_phase_variance = "numDeriv not installed"))
  rec <- .hzr_degraded_record(
    vcov = vc, weak = NULL, dist = "multiphase",
    control = list(conserve_applied = TRUE),
    reasons = list(conserved_variance = NA_character_))
  expect_length(rec$degraded, 0L)
})

test_that("entries come out in canonical order", {
  rec <- .hzr_degraded_record(
    vcov = NA, weak = NA, dist = "multiphase",
    control = list(conserve_applied = FALSE,
                   conserve_disabled_reason = "unsupported_censoring"),
    reasons = list(weak = "standard errors unavailable",
                   se = "Hessian not invertible"))
  expect_identical(rec$degraded, c("standard_errors", "weak_direction_check",
                                   "conservation_of_events"))
})

test_that("the validator accepts a record that matches the object", {
  expect_true(.hzr_validate_degraded(fake_fit(), fitted = TRUE))
  ok <- fake_fit(c("standard_errors", "weak_direction_check"),
                 c(standard_errors = "x", weak_direction_check = "y"),
                 vcov = NA, weak = NA)
  expect_true(.hzr_validate_degraded(ok, fitted = TRUE))
})

test_that("the validator rejects each way the record can go wrong", {
  bad <- list(
    list(fake_fit(c("standard_errors", "standard_errors"),
                  c(standard_errors = "x", standard_errors = "x"), vcov = NA),
         "duplicate"),
    list(fake_fit("coffee", c(coffee = "x")), "unknown capability"),
    list(fake_fit(c("weak_direction_check", "standard_errors"),
                  c(weak_direction_check = "y", standard_errors = "x"),
                  vcov = NA, weak = NA),
         "canonical order"),
    list(fake_fit("standard_errors", c(weak_direction_check = "x"), vcov = NA),
         "names\\(degraded_causes\\)"),
    list(fake_fit("standard_errors", c(standard_errors = ""), vcov = NA),
         "empty cause"),
    list(fake_fit(vcov = NA), "'standard_errors' is absent"),
    list(fake_fit("standard_errors", c(standard_errors = "x")),
         "'standard_errors' is listed"),
    list(fake_fit(weak = NA), "'weak_direction_check' is absent"),
    list(fake_fit(dist = "multiphase", control = list(conserve_applied = FALSE)),
         "'conservation_of_events' is absent"),
    list(fake_fit("conserved_phase_variance", c(conserved_phase_variance = "x")),
         "was not applied")
  )
  for (b in bad) {
    expect_error(.hzr_validate_degraded(b[[1]], fitted = TRUE), b[[2]],
                 info = b[[2]])
  }
  expect_error(.hzr_validate_degraded(fake_fit(), fitted = FALSE),
               "'fitting' is absent")
  expect_error(.hzr_validate_degraded(fake_fit(), fitted = TRUE, imported = TRUE),
               "'fitting' is absent")
})

test_that("the formatter prints none, the entries, or not recorded", {
  expect_identical(.hzr_format_not_done(character(0), character(0)),
                   "  Not done in this run: none")
  expect_identical(
    .hzr_format_not_done(c("standard_errors", "weak_direction_check"),
                         c(standard_errors = "a", weak_direction_check = "b")),
    c("  Not done in this run:", "    standard_errors: a",
      "    weak_direction_check: b"))
  expect_identical(
    .hzr_format_not_done(NULL, NULL),
    "  Not done in this run: not recorded (object built before this record existed)")
})

test_that("the output check passes a faithful block and names each fault", {
  obj <- fake_fit(c("standard_errors", "weak_direction_check"),
                  c(standard_errors = "a", weak_direction_check = "b"),
                  vcov = NA, weak = NA)
  good <- c("hazard object", "  Not done in this run:",
            "    standard_errors: a", "    weak_direction_check: b",
            "  other: 1")
  chk <- function(lines, o = obj) .hzr_check_not_done_output(lines, o)

  expect_identical(chk(good), character(0))
  expect_identical(chk(good[-(2:4)]), "heading missing")
  expect_identical(chk(c(good, "  Not done in this run: none")),
                   "heading printed more than once")
  expect_true("entry missing: weak_direction_check" %in% chk(good[-4]))
  expect_true("'none' printed over a non-empty record" %in%
                chk("  Not done in this run: none"))
  expect_true("cause differs for standard_errors" %in%
                chk(sub(": a$", ": z", good)))
  expect_true("entry not in the record: fitting" %in%
                chk(append(good, "    fitting: q", after = 2)))
  expect_identical(chk(good[c(1, 2, 4, 3, 5)]), "entries out of order")

  expect_identical(chk("  Not done in this run: none", fake_fit()), character(0))

  legacy <- fake_fit()
  legacy$degraded <- NULL
  legacy$degraded_causes <- NULL
  expect_identical(
    chk("  Not done in this run: not recorded (object built before this record existed)",
        legacy),
    character(0))
  expect_identical(chk("  Not done in this run: none", legacy),
                   "an object with no record was not reported as unrecorded")
})
```

- [ ] **Step 2: Run the tests to confirm they fail**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-degraded.R")'`
Expected: FAIL with `could not find function ".hzr_degraded_record"`.

- [ ] **Step 3: Write `R/degraded.R`**

```r
# degraded.R -- what a fit did not do (#242, following #197)
#
# Every hazard object carries $degraded, the capabilities this run did not
# deliver in canonical order, and $degraded_causes, the same names each with
# the reason. print() and summary() always show the record, "none" included:
# a line that appears only on bad news cannot be told from a line its author
# forgot to write. Design: inst/dev/DEGRADED-RECORD-DESIGN.md.
#
# The object's state decides WHICH entries appear; the reason fields carried
# up from the optimizer supply only WHY. A lost step whose reason did not
# reach the builder still gets its entry, with the cause "cause not recorded",
# so a path the reasons miss is reported rather than hidden.

.hzr_capabilities <- c(
  "fitting", "standard_errors", "conserved_phase_variance",
  "weak_direction_check", "conservation_of_events"
)

.hzr_cause_not_recorded <- "cause not recorded"

.hzr_not_done_heading <- "Not done in this run:"

.hzr_not_done_legacy <- "not recorded (object built before this record existed)"

# Readable causes for the codes .hzr_optim_multiphase() stores in
# spec$control$conserve_disabled_reason.
.hzr_coe_causes <- c(
  not_requested         = "not requested (conserve = FALSE)",
  unsupported_censoring = "status outside {0, 1} (left or interval censoring)",
  single_phase          = "fewer than two phases",
  no_events             = "no events",
  setup_failed          = "setup did not complete"
)

# TRUE for the weak = NA "not examined" value only: a list (ridge found) and
# NULL (examined, clean) are both FALSE.
.hzr_is_na_scalar <- function(x) {
  !is.list(x) && length(x) == 1L && is.na(x)
}

.hzr_reason_or_unrecorded <- function(r) {
  if (is.null(r) || length(r) != 1L || is.na(r) || !nzchar(r)) {
    .hzr_cause_not_recorded
  } else {
    r
  }
}

#' Build the record of what a fit did not do
#'
#' @param vcov,weak The fit's covariance and weak-direction result.
#' @param control The fit's control list (reads `conserve_applied` and
#'   `conserve_disabled_reason`).
#' @param dist Distribution name.
#' @param fitted Whether the optimizer ran.
#' @param imported Whether the object was read from SAS output.
#' @param not_fitted_cause The `fitting` cause when `fitted` is `FALSE`.
#' @param reasons List of carried-up reasons: `se`, `weak`,
#'   `conserved_variance`, each a single string or `NA`.
#' @return `list(degraded, degraded_causes)`, in canonical order.
#' @noRd
.hzr_degraded_record <- function(vcov, weak, control, dist,
                                 fitted = TRUE, imported = FALSE,
                                 not_fitted_cause = "not requested (fit = FALSE)",
                                 reasons = list()) {
  # Entries are added in canonical order, so names(causes) is already ordered.
  causes <- character(0)

  if (!fitted) {
    causes["fitting"] <- not_fitted_cause
  } else if (imported) {
    causes["fitting"] <- "imported from SAS output; not fitted in R"
  }

  if (!is.matrix(vcov)) {
    causes["standard_errors"] <- if (!fitted) {
      "model not fitted"
    } else if (imported) {
      "covariance not imported"
    } else {
      .hzr_reason_or_unrecorded(reasons$se)
    }
  }

  # No state marks a failed recompute (the conserved log_mu simply keeps an NA
  # variance), so this one entry is keyed on its reason.
  cv <- reasons$conserved_variance
  if (isTRUE(control$conserve_applied) && length(cv) == 1L && !is.na(cv)) {
    causes["conserved_phase_variance"] <- cv
  }

  if (.hzr_is_na_scalar(weak)) {
    causes["weak_direction_check"] <- if (!fitted) {
      "model not fitted"
    } else if (imported) {
      "imported from SAS output; no R Hessian"
    } else {
      .hzr_reason_or_unrecorded(reasons$weak)
    }
  }

  if (identical(dist, "multiphase") && isFALSE(control$conserve_applied)) {
    r <- control$conserve_disabled_reason
    known <- length(r) == 1L && !is.na(r) && r %in% names(.hzr_coe_causes)
    causes["conservation_of_events"] <- if (known) {
      .hzr_coe_causes[[r]]
    } else {
      .hzr_cause_not_recorded
    }
  }

  list(degraded = names(causes) %||% character(0), degraded_causes = causes)
}

#' Check the record against the object it describes
#'
#' A disagreement is a package bug, not a user error, so this stops.
#'
#' @param object A `hazard` object carrying `degraded` / `degraded_causes`.
#' @param fitted Whether the optimizer ran.
#' @param imported Whether the object was read from SAS output.
#' @return `TRUE`, invisibly; otherwise an error.
#' @noRd
.hzr_validate_degraded <- function(object, fitted, imported = FALSE) {
  d <- object$degraded
  cz <- object$degraded_causes
  fail <- function(...) {
    stop("Internal error: the 'Not done in this run' record disagrees with ",
         "the object it describes (", ..., "). This is a TemporalHazard ",
         "bug; please report it at ",
         "https://github.com/ehrlinger/TemporalHazard/issues.", call. = FALSE)
  }

  if (!is.character(d) || !is.character(cz)) fail("the record is not character")
  if (anyDuplicated(d)) fail("duplicate entry")
  unknown <- setdiff(d, .hzr_capabilities)
  if (length(unknown)) fail("unknown capability '", unknown[1], "'")
  if (!identical(d, intersect(.hzr_capabilities, d))) {
    fail("entries out of canonical order")
  }
  if (!identical(as.character(names(cz)), d)) {
    fail("names(degraded_causes) differ from degraded")
  }
  if (anyNA(cz) || !all(nzchar(cz))) fail("empty cause")

  ctl <- object$spec$control
  agree <- function(cap, state) {
    if (!identical(cap %in% d, state)) {
      fail("'", cap, "' is ", if (state) "absent" else "listed",
           " but the object's state says the opposite")
    }
  }
  agree("fitting", !fitted || imported)
  agree("standard_errors", !is.matrix(object$fit$vcov))
  agree("weak_direction_check", .hzr_is_na_scalar(object$fit$weak))
  agree("conservation_of_events",
        identical(object$spec$dist, "multiphase") &&
          isFALSE(ctl$conserve_applied))
  if ("conserved_phase_variance" %in% d && !isTRUE(ctl$conserve_applied)) {
    fail("'conserved_phase_variance' is listed but Conservation of Events ",
         "was not applied")
  }
  invisible(TRUE)
}

#' Lines of the "Not done in this run" block
#'
#' The one formatter behind print.hazard() and print.summary.hazard().
#' `degraded = NULL` means the object predates the record.
#'
#' @param degraded,causes The record.
#' @return Character vector, one element per printed line.
#' @noRd
.hzr_format_not_done <- function(degraded, causes) {
  if (is.null(degraded)) {
    return(paste0("  ", .hzr_not_done_heading, " ", .hzr_not_done_legacy))
  }
  if (!length(degraded)) {
    return(paste0("  ", .hzr_not_done_heading, " none"))
  }
  c(paste0("  ", .hzr_not_done_heading),
    paste0("    ", degraded, ": ", unname(causes[degraded])))
}

#' Check printed output against the record
#'
#' Parses the printed block back into entries and compares them with the raw
#' record. Deliberately independent of .hzr_format_not_done(): it uses its
#' own literal strings and never calls the formatter, so a printer that drops
#' or rewrites an entry cannot also rewrite what it is checked against.
#'
#' @param lines Printed output, as from capture.output().
#' @param object The `hazard` (or `summary.hazard`) object that was printed.
#' @return Character vector of problems; `character(0)` when faithful.
#' @noRd
.hzr_check_not_done_output <- function(lines, object) {
  lines <- sub("[[:space:]]+$", "", lines)
  at <- grep("^[[:space:]]*Not done in this run:", lines)
  if (length(at) == 0L) return("heading missing")
  if (length(at) > 1L) return("heading printed more than once")
  rest <- sub("^[[:space:]]*Not done in this run:[[:space:]]*", "", lines[at])

  if (is.null(object$degraded)) {
    legacy <- "not recorded (object built before this record existed)"
    if (identical(rest, legacy)) return(character(0))
    return("an object with no record was not reported as unrecorded")
  }

  record <- object$degraded_causes
  rec_names <- as.character(names(record))
  entry_re <- "^    ([a-z_]+): (.+)$"
  problems <- character(0)

  if (identical(rest, "none")) {
    printed <- character(0)
    if (length(record)) {
      problems <- c(problems, "'none' printed over a non-empty record")
    }
  } else if (nzchar(rest)) {
    return(paste0("unexpected text after the heading: '", rest, "'"))
  } else {
    j <- at + 1L
    while (j <= length(lines) && grepl(entry_re, lines[j])) j <- j + 1L
    ent <- lines[seq.int(at + 1L, length.out = j - at - 1L)]
    printed <- sub(entry_re, "\\2", ent)
    names(printed) <- sub(entry_re, "\\1", ent)
    if (!length(printed)) {
      problems <- c(problems, "heading with neither entries nor 'none'")
    }
  }

  pr_names <- as.character(names(printed))
  if (anyDuplicated(pr_names)) problems <- c(problems, "an entry printed twice")
  # sprintf(), not paste0(): paste0("x", character(0)) is "x", not character(0).
  problems <- c(problems,
                sprintf("entry missing: %s", setdiff(rec_names, pr_names)),
                sprintf("entry not in the record: %s", setdiff(pr_names, rec_names)))
  shared <- intersect(rec_names, pr_names)
  differ <- shared[record[shared] != printed[shared]]
  problems <- c(problems, sprintf("cause differs for %s", differ))
  if (!length(problems) && !identical(pr_names, rec_names)) {
    problems <- "entries out of order"
  }
  problems
}
```

- [ ] **Step 4: Run the tests to confirm they pass**

Run: `Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-degraded.R")'`
Expected: 0 lints; all tests in `test-degraded.R` PASS, with 0 failures and 0 skips.

- [ ] **Step 5: Commit**

```bash
git add R/degraded.R tests/testthat/test-degraded.R
git commit -m "feat(degraded): the record of what a fit did not do, with its validator (#242)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: Carry the standard-error reason up from the Hessian path

**Files:**
- Modify: `R/hessian-invert.R` (`.hzr_safe_solve()`, around lines 20-76)
- Modify: `R/optimizer.R` (a new helper before `.hzr_optim_generic()` around line 51; the Hessian block around lines 163-221)
- Test: `tests/testthat/test-degraded-reasons.R`

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces:
  - `.hzr_safe_solve()` result gains `reason`: `NA_character_`, `"Hessian has non-finite entries"`, or `"Hessian not invertible"`.
  - `.hzr_numderiv_available()` → logical(1). Mockable; used in Task 5 too.
  - The `.hzr_optim_generic()` result gains `se_unavailable_reason`: `NA_character_` when `vcov` is a matrix, otherwise one of `"numDeriv not installed and no analytic Hessian"`, `"numDeriv::hessian() failed"`, `"Hessian has non-finite entries"`, `"Hessian not invertible"`.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-degraded-reasons.R`:

```r
# Reasons carried up to the record (#242): the Hessian path's reason for losing
# standard errors, and the weak-direction check's reason for not looking.
# Each cause string is the spec's, verbatim.

# A one-parameter exponential model, log-rate scale, small enough to fit in
# milliseconds and with a Hessian known in closed form: exp(theta) * sum(time).
exp_logl <- function(theta, time, status, time_lower = NULL, time_upper = NULL,
                     x = NULL, weights = NULL, return_gradient = FALSE) {
  sum(status * theta[1] - exp(theta[1]) * time)
}
exp_grad <- function(theta, time, status, ...) {
  sum(status) - exp(theta[1]) * sum(time)
}
fit_exp <- function(hess) {
  set.seed(3)
  time <- stats::rexp(40, 0.5)
  suppressWarnings(.hzr_optim_generic(
    logl_fn = exp_logl, gradient_fn = exp_grad,
    time = time, status = rep(1, 40), theta_start = 0,
    hessian_fn = function(p) hess(p, time)
  ))
}

test_that(".hzr_safe_solve() says why it returned no covariance", {
  expect_identical(suppressWarnings(.hzr_safe_solve(matrix(NaN, 2, 2)))$reason,
                   "Hessian has non-finite entries")
  expect_identical(suppressWarnings(.hzr_safe_solve(matrix(0, 2, 2)))$reason,
                   "Hessian not invertible")
  expect_identical(.hzr_safe_solve(diag(2))$reason, NA_character_)
})

test_that("an analytic Hessian that inverts leaves no reason", {
  res <- fit_exp(function(p, time) matrix(exp(p) * sum(time), 1, 1))
  expect_true(is.matrix(res$vcov))                 # guard: SEs exist
  expect_identical(res$se_unavailable_reason, NA_character_)
})

test_that("each way the Hessian path loses standard errors names its cause", {
  res <- fit_exp(function(p, time) matrix(NaN, 1, 1))
  expect_false(is.matrix(res$vcov))
  expect_identical(res$se_unavailable_reason, "Hessian has non-finite entries")

  res <- fit_exp(function(p, time) matrix(0, 1, 1))
  expect_false(is.matrix(res$vcov))
  expect_identical(res$se_unavailable_reason, "Hessian not invertible")

  local_mocked_bindings(.hzr_numderiv_available = function() FALSE)
  res <- fit_exp(function(p, time) NULL)
  expect_false(is.matrix(res$vcov))
  expect_identical(res$se_unavailable_reason,
                   "numDeriv not installed and no analytic Hessian")
})

test_that("a numDeriv failure is named, not folded into 'not installed'", {
  skip_if_not_installed("numDeriv")
  local_mocked_bindings(hessian = function(...) stop("boom"),
                        .package = "numDeriv")
  res <- fit_exp(function(p, time) NULL)
  expect_false(is.matrix(res$vcov))
  expect_identical(res$se_unavailable_reason, "numDeriv::hessian() failed")
})
```

- [ ] **Step 2: Run the tests to confirm they fail**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-degraded-reasons.R")'`
Expected: FAIL. `$reason` and `$se_unavailable_reason` are `NULL`, so `expect_identical` reports "actual is NULL", and `.hzr_numderiv_available` cannot be mocked because it does not exist yet.

- [ ] **Step 3: Add `reason` to `.hzr_safe_solve()`**

In `R/hessian-invert.R`, make these three replacements:

```r
    warning("Hessian contains non-finite entries; standard errors unavailable")
    return(list(vcov = NA, rcond = NA_real_, pd = NA))
```
→
```r
    warning("Hessian contains non-finite entries; standard errors unavailable")
    return(list(vcov = NA, rcond = NA_real_, pd = NA,
                reason = "Hessian has non-finite entries"))
```

```r
      warning("Hessian not invertible; standard errors unavailable")
      return(list(vcov = NA, rcond = rc, pd = NA))
```
→
```r
      warning("Hessian not invertible; standard errors unavailable")
      return(list(vcov = NA, rcond = rc, pd = NA,
                  reason = "Hessian not invertible"))
```

```r
  list(vcov = vcov, rcond = rc, pd = pd)
}
```
→
```r
  list(vcov = vcov, rcond = rc, pd = pd, reason = NA_character_)
}
```

In the roxygen `@return` of `.hzr_safe_solve()`, change `#'   if not invertible).` to:

```r
#'   if not invertible); \code{reason} (why \code{vcov} is \code{NA}:
#'   \code{"Hessian has non-finite entries"} or \code{"Hessian not
#'   invertible"}; \code{NA_character_} when a matrix was returned).
```

- [ ] **Step 4: Add `.hzr_numderiv_available()` and `se_unavailable_reason` to the optimizer**

In `R/optimizer.R`, insert immediately before the roxygen block of `.hzr_optim_generic()`, after the file-header comment:

```r
# numDeriv is a Suggests. Wrapped so tests can simulate an install without it:
# requireNamespace() is base R and cannot be mocked directly.
.hzr_numderiv_available <- function() {
  requireNamespace("numDeriv", quietly = TRUE)
}

```

Then, still in `.hzr_optim_generic()`, replace:

```r
  if (is.null(hess_result)) {
    # Reached when no analytic Hessian was supplied, or the hook declined by
```
with
```r
  # Why standard errors are unavailable, when they are, for the record that
  # print() and summary() show (#242). Set beside each warning below, so the
  # warning and the record name the same cause.
  se_reason <- NA_character_
  if (is.null(hess_result)) {
    # Reached when no analytic Hessian was supplied, or the hook declined by
```

Replace:
```r
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
```
with
```r
    if (!.hzr_numderiv_available()) {
```

Replace:
```r
              "not pull it by default.", call. = FALSE)
    } else {
```
with
```r
              "not pull it by default.", call. = FALSE)
      se_reason <- "numDeriv not installed and no analytic Hessian"
    } else {
```

Replace:
```r
          warning("numDeriv::hessian() failed, so standard errors are ",
                  "unavailable: ", conditionMessage(e), call. = FALSE)
          NULL
        }
      )
    }
  }
```
with
```r
          warning("numDeriv::hessian() failed, so standard errors are ",
                  "unavailable: ", conditionMessage(e), call. = FALSE)
          NULL
        }
      )
      if (is.null(hess_result)) se_reason <- "numDeriv::hessian() failed"
    }
  }
```

Replace:
```r
    list(vcov = NA, rcond = NA_real_, pd = NA)
  }
```
with
```r
    list(vcov = NA, rcond = NA_real_, pd = NA, reason = se_reason)
  }
```

Replace:
```r
    rcond = inv$rcond,
    pd = inv$pd
  )
}
```
with
```r
    rcond = inv$rcond,
    pd = inv$pd,
    se_unavailable_reason = if (is.matrix(inv$vcov)) NA_character_ else inv$reason
  )
}
```

Add to the `.hzr_optim_generic()` roxygen `@return`, if it lists fields, the line
`#'   \code{se_unavailable_reason} (why \code{vcov} is \code{NA}, else \code{NA_character_}),`.
If `@return` is a single sentence, append ` Includes \code{se_unavailable_reason}.` to it.

- [ ] **Step 5: Run the new tests and the existing Hessian tests**

Run: `Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && NOT_CRAN=true Rscript -e 'devtools::load_all(); for (f in c("test-degraded-reasons.R","test-hessian-invert.R","test-hessian-stability.R")) testthat::test_file(file.path("tests/testthat", f))'`
Expected: 0 lints; `test-degraded-reasons.R` all PASS; `test-hessian-invert.R` and `test-hessian-stability.R` unchanged, 0 failures.

If the `numDeriv` mock errors with "Can't find binding for `hessian`", the mock is not reaching the namespace. Report that; do not delete the test.

- [ ] **Step 6: Commit**

```bash
git add R/hessian-invert.R R/optimizer.R tests/testthat/test-degraded-reasons.R
git commit -m "feat(degraded): the Hessian path says why standard errors were lost (#242)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 3: The weak-direction check says why it did not look

**Files:**
- Modify: `R/hessian-invert.R` (`.hzr_weak_direction()`, around lines 142-248)
- Test: `tests/testthat/test-degraded-reasons.R` (append)

**Interfaces:**
- Consumes: nothing.
- Produces: `.hzr_weak_direction_impl(vcov, rcond, param_names = NULL, tol = .hzr_rcond_tol, cor_tol = .hzr_ridge_cor_tol, share = 0.9)` → `list(weak = <list|NULL|NA>, reason = <chr(1)>)`. `reason` is `NA_character_` unless `weak` is `NA`. `.hzr_weak_direction()` keeps its exact signature and return value.

- [ ] **Step 1: Append the failing tests**

Append to `tests/testthat/test-degraded-reasons.R`:

```r
test_that("every 'could not look' exit of the weak-direction check names why", {
  why <- function(...) .hzr_weak_direction_impl(...)$reason
  expect_identical(why(NULL, NA_real_), "standard errors unavailable")
  expect_identical(why(diag(2), NA_real_), "Hessian condition number unavailable")
  expect_identical(why(NA, 1e-10), "standard errors unavailable")
  # Finite, positive diagonal, but a non-finite entry among estimated params.
  expect_identical(why(matrix(c(1, Inf, Inf, 1), 2), 1e-10),
                   "covariance has non-finite entries")
  # Finite covariance whose correlation overflows: 1e300 / (1e-150)^2.
  expect_identical(why(matrix(c(1e-300, 1e300, 1e300, 1e-300), 2), 1e-10),
                   "covariance has non-finite entries")
})

test_that("the split leaves .hzr_weak_direction() returning what it did", {
  expect_identical(.hzr_weak_direction(NULL, NA_real_), NA)
  expect_null(.hzr_weak_direction(diag(2), 0.5))
  res <- .hzr_weak_direction_impl(diag(2), 0.5)
  expect_true("weak" %in% names(res))     # a NULL result keeps its slot
  expect_null(res$weak)
  expect_identical(res$reason, NA_character_)
})
```

- [ ] **Step 2: Run to confirm they fail**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-degraded-reasons.R")'`
Expected: FAIL with `could not find function ".hzr_weak_direction_impl"`.

- [ ] **Step 3: Split the function**

In `R/hessian-invert.R`:

Rename the definition line `.hzr_weak_direction <- function(vcov, rcond, param_names = NULL,` to `.hzr_weak_direction_impl <- function(vcov, rcond, param_names = NULL,` (the other arguments are unchanged). Its roxygen block stays attached to it. Add this line to that block, just before `#' @noRd`:

```r
#'   The \code{_impl} returns \code{list(weak, reason)}, where \code{reason}
#'   names why \code{weak} is \code{NA}; \code{.hzr_weak_direction()} returns
#'   \code{weak} alone, unchanged, for every existing caller.
```

As the first lines of the function body, before the `# (1) Gate on ...` comment, insert:

```r
  na_because <- function(reason) list(weak = NA, reason = reason)
  looked <- function(weak) list(weak = weak, reason = NA_character_)

```

Replace each exit exactly:

| Old line | New line |
|---|---|
| `  if (length(rcond) != 1L \|\| is.na(rcond)) return(NA)` | `  if (length(rcond) != 1L \|\| is.na(rcond)) {`<br>`    return(na_because(if (is.matrix(vcov)) "Hessian condition number unavailable"`<br>`                      else "standard errors unavailable"))`<br>`  }` |
| `  if (rcond >= tol) return(NULL)` | `  if (rcond >= tol) return(looked(NULL))` |
| `  if (is.null(vcov) \|\| !is.matrix(vcov)) return(NA)` | `  if (is.null(vcov) \|\| !is.matrix(vcov)) return(na_because("standard errors unavailable"))` |
| `  if (nrow(vcov) < 2L) return(NULL)` | `  if (nrow(vcov) < 2L) return(looked(NULL))` |
| `  if (length(keep) < 2L) return(NULL)` | `  if (length(keep) < 2L) return(looked(NULL))` |
| `  if (anyNA(V) \|\| any(!is.finite(V))) return(NA)` | `  if (anyNA(V) \|\| any(!is.finite(V))) return(na_because("covariance has non-finite entries"))` |
| `  if (anyNA(R) \|\| any(!is.finite(R))) return(NA)` | `  if (anyNA(R) \|\| any(!is.finite(R))) return(na_because("covariance has non-finite entries"))` |
| `  if (is.null(e)) return(NA)` | `  if (is.null(e)) return(na_because("eigendecomposition failed"))` |
| `  if (is.null(found)) return(NULL)` | `  if (is.null(found)) return(looked(NULL))` |
| `  found$n_directions <- length(seen)`<br>`  found` | `  found$n_directions <- length(seen)`<br>`  looked(found)` |

(The `\|` above is a literal `|` escaped for the table.) Then run
`grep -n "return(NA)\|return(NULL)" R/hessian-invert.R`. It must print **no line
between the `_impl` definition and its closing brace**.

Immediately after the `_impl` function's closing brace, add the wrapper:

```r

# The shape every existing caller relies on: list / NULL / NA, unchanged.
.hzr_weak_direction <- function(vcov, rcond, param_names = NULL,
                                tol = .hzr_rcond_tol,
                                cor_tol = .hzr_ridge_cor_tol,
                                share = 0.9) {
  .hzr_weak_direction_impl(vcov, rcond, param_names = param_names,
                           tol = tol, cor_tol = cor_tol, share = share)$weak
}
```

- [ ] **Step 4: Run the new tests and the whole weak-direction suite**

Run: `Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && NOT_CRAN=true Rscript -e 'devtools::load_all(); for (f in c("test-degraded-reasons.R","test-weak-direction.R")) testthat::test_file(file.path("tests/testthat", f))'`
Expected: 0 lints; both files 0 failures. `test-weak-direction.R` must pass unmodified. It exercises the wrapper, and that is the proof its return is unchanged.

- [ ] **Step 5: Commit**

```bash
git add R/hessian-invert.R tests/testthat/test-degraded-reasons.R
git commit -m "feat(degraded): the weak-direction check names why it did not look (#242)

The eigen() failure exit has a cause string but no test: eigen() is base R
and cannot be mocked, and no finite symmetric input is known to fail it.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4: Stamp and validate the record in `hazard()`

**Files:**
- Modify: `R/hazard_api.R` (the dispatch around line 754; the multiphase branch around 803; the single-distribution branch around 831; the weak check around 857; the assembly around 910; `summary.hazard()` around 1599)
- Test: `tests/testthat/test-degraded-fits.R`

**Interfaces:**
- Consumes: `.hzr_degraded_record()`, `.hzr_validate_degraded()` (Task 1); `optim_result$se_unavailable_reason` (Task 2); `.hzr_weak_direction_impl()` (Task 3); `.hzr_numderiv_available()` (Task 2, mocked in tests).
- Produces: every object from `hazard()` has `$degraded` and `$degraded_causes`; `summary.hazard()` objects carry both. Task 5 fills `degraded_reasons$conserved_variance` from `optim_result$conserved_variance_reason`, and this task already reads it.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-degraded-fits.R`:

```r
# The record on real fits (#242). Every test asserts a guard first -- that the
# fit really is in the state the record is about -- and that no entry fell
# back to "cause not recorded", which exists for paths the reasons miss and
# which no path this package knows about should reach.

weib_data <- function(n = 80) {
  set.seed(7)
  data.frame(t = stats::rweibull(n, 1.4, 3), d = rep(1L, n))
}

# Copied from test-coe-applied.R's coe_fit(): explicit starts, n_starts = 1.
mp_fit <- function(status_pattern, conserve = TRUE) {
  set.seed(42)
  n <- 60
  time <- stats::rexp(n, 0.2) + 0.05
  status <- rep(status_pattern, length.out = n)
  args <- list(time = time, status = status)
  if (any(status == 2)) {
    args$time_lower <- time
    args$time_upper <- ifelse(status == 2, time * 1.5, time)
  }
  suppressWarnings(do.call(hazard, c(args, list(
    phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                  const = hzr_phase("constant")),
    dist = "multiphase", fit = TRUE,
    control = list(conserve = conserve, n_starts = 1)
  ))))
}

test_that("a clean fit records nothing", {
  fit <- hazard(survival::Surv(t, d) ~ 1, data = weib_data(), dist = "weibull",
                fit = TRUE, theta = c(mu = 1, nu = 1))
  expect_true(is.matrix(fit$fit$vcov))             # guard: SEs exist
  expect_null(fit$fit$weak)                        # guard: examined, clean
  expect_identical(fit$degraded, character(0))
  expect_length(fit$degraded_causes, 0L)
})

test_that("an unfitted object records that it was not fitted, and why", {
  dat <- weib_data()
  fit0 <- hazard(time = dat$t, status = dat$d, dist = "weibull", fit = FALSE)
  expect_identical(fit0$degraded,
                   c("fitting", "standard_errors", "weak_direction_check"))
  expect_identical(fit0$degraded_causes[["fitting"]],
                   "not requested (fit = FALSE)")
  expect_false("cause not recorded" %in% fit0$degraded_causes)
})

test_that("fit = TRUE without theta says the optimizer did not run", {
  # hazard() runs the single-distribution optimizer only when theta is given,
  # and has never said so. If this call errors or fits, the premise has
  # changed: stop and report, do not adjust the test.
  dat <- weib_data()
  fit1 <- hazard(time = dat$t, status = dat$d, dist = "weibull", fit = TRUE)
  expect_null(fit1$fit$counts)                     # guard: no optimizer ran
  expect_identical(fit1$degraded_causes[["fitting"]],
                   "no starting values (theta = NULL); the optimizer did not run")
})

test_that("CoE switched off by the user is listed as not requested", {
  skip_on_cran()
  fit <- mp_fit(c(1, 0), conserve = FALSE)
  expect_false(fit$spec$control$conserve_applied)  # guard
  expect_identical(fit$degraded_causes[["conservation_of_events"]],
                   "not requested (conserve = FALSE)")
  expect_false("cause not recorded" %in% fit$degraded_causes)
})

test_that("interval rows disable CoE, and the record says why", {
  skip_on_cran()
  fit <- mp_fit(c(1, 0, 2))
  expect_false(is.null(fit$fit$theta))             # guard: the fit ran
  expect_identical(fit$degraded_causes[["conservation_of_events"]],
                   "status outside {0, 1} (left or interval censoring)")
  expect_false("cause not recorded" %in% fit$degraded_causes)
})

test_that("no analytic Hessian and no numDeriv: SEs and the weak check are listed", {
  skip_on_cran()
  # Interval rows make the analytic multiphase Hessian decline by design, so
  # the fit needs numDeriv; mocked absent, no Hessian exists at all.
  local_mocked_bindings(.hzr_numderiv_available = function() FALSE)
  fit <- mp_fit(c(1, 0, 2))
  expect_false(is.matrix(fit$fit$vcov))            # guard: SEs really were lost
  expect_identical(
    fit$degraded_causes[c("standard_errors", "weak_direction_check")],
    c(standard_errors = "numDeriv not installed and no analytic Hessian",
      weak_direction_check = "standard errors unavailable"))
  expect_false("cause not recorded" %in% fit$degraded_causes)
})

test_that("summary() carries the record", {
  dat <- weib_data()
  fit0 <- hazard(time = dat$t, status = dat$d, dist = "weibull", fit = FALSE)
  s <- summary(fit0)
  expect_identical(s$degraded, fit0$degraded)
  expect_identical(s$degraded_causes, fit0$degraded_causes)
})
```

- [ ] **Step 2: Run to confirm they fail**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-degraded-fits.R")'`
Expected: FAIL, because `fit$degraded` is `NULL` in every test.

- [ ] **Step 3: Track whether the optimizer ran, and collect the reasons**

In `R/hazard_api.R`, replace:

```r
  # Distribution dispatch -- select the distribution-specific optimizer and fit.
  if (fit && dist == "multiphase") {
```
with
```r
  # Filled by whichever optimizer branch runs, and read by the record of what
  # this fit did not do (#242). fit_ran stays FALSE when fit = TRUE but the
  # single-distribution branch below is skipped because theta is NULL -- a
  # skip that raises no error, and that the record is now the only report of.
  fit_ran <- FALSE
  degraded_reasons <- list()

  # Distribution dispatch -- select the distribution-specific optimizer and fit.
  if (fit && dist == "multiphase") {
```

Replace:
```r
    control$conserve_disabled_reason <- optim_result$conserve_disabled_reason
```
with
```r
    control$conserve_disabled_reason <- optim_result$conserve_disabled_reason
    fit_ran <- TRUE
    degraded_reasons$se <- optim_result$se_unavailable_reason
    degraded_reasons$conserved_variance <- optim_result$conserved_variance_reason
```

Replace (the end of the single-distribution branch):
```r
    fit_state$counts <- optim_result$counts
    fit_state$message <- optim_result$message
  }

  # An ill-conditioned Hessian already warns that standard errors are
```
with
```r
    fit_state$counts <- optim_result$counts
    fit_state$message <- optim_result$message
    fit_ran <- TRUE
    degraded_reasons$se <- optim_result$se_unavailable_reason
  }

  # An ill-conditioned Hessian already warns that standard errors are
```

- [ ] **Step 4: Take the weak reason from the same exit that returned NA**

Replace:
```r
  fit_state$weak <- .hzr_weak_direction(fit_state$vcov, fit_state$rcond,
                                        weak_names)
```
with
```r
  weak_check <- .hzr_weak_direction_impl(fit_state$vcov, fit_state$rcond,
                                         weak_names)
  fit_state$weak <- weak_check$weak
  degraded_reasons$weak <- weak_check$reason
```

- [ ] **Step 5: Stamp and validate**

Replace:
```r
  class(obj) <- "hazard"
  obj
}

#' Predict from a hazard model object
```
with
```r
  class(obj) <- "hazard"

  # What this fit did not do, and why (#242). The object's state decides which
  # entries appear; the reasons carried up from the optimizer say why. The
  # validator stops if the two disagree, because that is a package bug.
  record <- .hzr_degraded_record(
    vcov = fit_state$vcov, weak = fit_state$weak, control = control,
    dist = dist, fitted = fit_ran,
    not_fitted_cause = if (fit) {
      "no starting values (theta = NULL); the optimizer did not run"
    } else {
      "not requested (fit = FALSE)"
    },
    reasons = degraded_reasons
  )
  obj$degraded <- record$degraded
  obj$degraded_causes <- record$degraded_causes
  .hzr_validate_degraded(obj, fitted = fit_ran)
  obj
}

#' Predict from a hazard model object
```

- [ ] **Step 6: Let `summary()` carry the record**

In `summary.hazard()`, replace:
```r
    weak = object$fit$weak,
    phases = object$spec$phases
```
with
```r
    weak = object$fit$weak,
    degraded = object$degraded,
    degraded_causes = object$degraded_causes,
    phases = object$spec$phases
```

- [ ] **Step 7: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-degraded-fits.R")'`
Expected: 0 lints; all PASS.

Then run: `NOT_CRAN=true Rscript -e 'devtools::test()'`
Expected: 0 failures. **Any `Internal error: the 'Not done in this run' record disagrees` is a real finding**: some fit path produces state the builder does not model. Report the call and the message. Do not loosen the validator.

- [ ] **Step 8: Commit**

```bash
git add R/hazard_api.R tests/testthat/test-degraded-fits.R
git commit -m "feat(degraded): hazard() stamps and validates the record (#242)

Also records, rather than changes, a silent skip: a single-distribution
hazard(fit = TRUE) with no theta never reaches the optimizer.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 5: Record a failed CoE full-information recompute

**Files:**
- Modify: `R/likelihood-multiphase.R` (the recompute block around lines 2138-2197)
- Test: `tests/testthat/test-degraded-fits.R` (append)

**Interfaces:**
- Consumes: `.hzr_numderiv_available()` (Task 2); the `.hzr_safe_solve()` `reason` (Task 2); `degraded_reasons$conserved_variance` wiring (Task 4).
- Produces: the `.hzr_optim_multiphase()` result gains `conserved_variance_reason`: `"numDeriv not installed"`, `"Hessian could not be computed"`, `"Hessian has non-finite entries"` or `"Hessian not invertible"`, set only when CoE was applied and the recompute failed. When the recompute succeeds, `se_unavailable_reason` resets to `NA_character_`.

- [ ] **Step 1: Append the failing test**

Append to `tests/testthat/test-degraded-fits.R`:

```r
test_that("a failed full-information recompute under CoE is recorded", {
  skip_on_cran()
  # With the analytic Hessian declining and numDeriv absent, both the fit's
  # own vcov and the CoE recompute have no Hessian to work from.
  local_mocked_bindings(
    .hzr_hessian_multiphase = function(...) NULL,
    .hzr_numderiv_available = function() FALSE
  )
  fit <- mp_fit(c(1, 0))
  expect_true(fit$spec$control$conserve_applied)   # guard: the recompute ran
  expect_identical(fit$degraded_causes[["conserved_phase_variance"]],
                   "numDeriv not installed")
  expect_false("cause not recorded" %in% fit$degraded_causes)
})
```

- [ ] **Step 2: Run to confirm it fails**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-degraded-fits.R")'`
Expected: this test FAILS with `subscript out of bounds`, because `conserved_phase_variance` is not in `degraded_causes`. The Task 4 tests still pass.

- [ ] **Step 3: Carry the reason**

In `R/likelihood-multiphase.R`, inside `if (use_conserve && !is.null(fixmu_pos)) {`, replace:

```r
      # numDeriv fallback only when analytic declines (left/interval rows)
      if (is.null(H_unc) && requireNamespace("numDeriv", quietly = TRUE)) {
```
with
```r
      # numDeriv fallback only when analytic declines (left/interval rows)
      if (is.null(H_unc) && .hzr_numderiv_available()) {
```

Directly after that `if` block's closing `}` (just before `if (is.matrix(H_unc)) {` and `inv_unc <- .hzr_safe_solve(H_unc)`), insert:

```r
      # Why the recompute failed, if it does, for the record print() shows
      # (#242). Refined below when there is a Hessian that will not invert.
      cv_reason <- if (is.matrix(H_unc)) {
        NA_character_
      } else if (.hzr_numderiv_available()) {
        "Hessian could not be computed"
      } else {
        "numDeriv not installed"
      }
```

Replace:
```r
          best_result$fixed_mask <- !free_unc
          recomputed <- TRUE
        }
      }
```
with
```r
          best_result$fixed_mask <- !free_unc
          recomputed <- TRUE
          # The full-information vcov replaces the reduced one, so a reason the
          # reduced one carried for missing SEs no longer applies.
          best_result$se_unavailable_reason <- NA_character_
        } else {
          cv_reason <- inv_unc$reason
        }
      }
```

Replace the warning block:
```r
      if (!recomputed) {
        # Fall back to the reduced (search-only) vcov, in which the conserved
        # log_mu has no variance. Warn loudly so the understated uncertainty is
        # visible rather than silent.
        warning(
          "Could not compute the full-information variance for the ",
          "Conservation-of-Events-conserved phase 'log_mu' (numDeriv ",
          "unavailable or the Hessian was not invertible). Its standard error ",
          "stays NA and downstream standard errors / confidence limits for ",
          "that phase may be understated.", call. = FALSE)
      }
```
with
```r
      if (!recomputed) {
        # Fall back to the reduced (search-only) vcov, in which the conserved
        # log_mu has no variance. Warn loudly so the understated uncertainty is
        # visible rather than silent, and record it, so it is not only a
        # warning that scrolled past.
        warning(
          "Could not compute the full-information variance for the ",
          "Conservation-of-Events-conserved phase 'log_mu' (", cv_reason,
          "). Its standard error stays NA and downstream standard errors / ",
          "confidence limits for that phase may be understated.", call. = FALSE)
        best_result$conserved_variance_reason <- cv_reason
      }
```

- [ ] **Step 4: Run the tests and the CoE suite**

Run: `Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && NOT_CRAN=true Rscript -e 'devtools::load_all(); for (f in c("test-degraded-fits.R","test-coe-applied.R")) testthat::test_file(file.path("tests/testthat", f))'`
Expected: 0 lints; both files 0 failures.

If the mocked fit errors before the record is built (for example inside phase-identifiability code that needs a Hessian), report the error. Do not remove the mock of `.hzr_hessian_multiphase`: it is the only known way to reach this branch.

- [ ] **Step 5: Commit**

```bash
git add R/likelihood-multiphase.R tests/testthat/test-degraded-fits.R
git commit -m "feat(degraded): record a failed CoE full-information recompute (#242)

Until now a warning was the only trace that the conserved phase's standard
error stayed NA and its confidence limits were understated.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 6: Stamp the SAS import

**Files:**
- Modify: `R/read-outhaz.R` (the `structure(... class = "hazard")` at the end of `.hzr_outhaz_to_spec()`, around lines 259-272)
- Test: `tests/testthat/test-degraded-fits.R` (append)

**Interfaces:**
- Consumes: `.hzr_degraded_record()`, `.hzr_validate_degraded()` (Task 1).
- Produces: every `hazard` object built from SAS output carries the record.

- [ ] **Step 1: Append the failing tests**

```r
outhaz_fixture <- function() {
  system.file("extdata", "outhaz-fixture.rds", package = "TemporalHazard")
}

test_that("an imported SAS fit records that R did not fit or examine it", {
  obj <- .hzr_outhaz_to_spec(hzr_read_outhaz(outhaz_fixture()), need_vcov = TRUE)
  expect_true(is.matrix(obj$fit$vcov))             # guard: covariance was read
  expect_identical(
    obj$degraded_causes,
    c(fitting = "imported from SAS output; not fitted in R",
      weak_direction_check = "imported from SAS output; no R Hessian"))
})

test_that("an import read without its covariance also lists standard errors", {
  obj <- .hzr_outhaz_to_spec(hzr_read_outhaz(outhaz_fixture()), need_vcov = FALSE)
  expect_false(is.matrix(obj$fit$vcov))            # guard
  expect_identical(obj$degraded,
                   c("fitting", "standard_errors", "weak_direction_check"))
  expect_identical(obj$degraded_causes[["standard_errors"]],
                   "covariance not imported")
})
```

- [ ] **Step 2: Run to confirm they fail**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-degraded-fits.R")'`
Expected: the two new tests FAIL, because `obj$degraded_causes` is `NULL`.

- [ ] **Step 3: Stamp**

In `R/read-outhaz.R`, replace:

```r
  structure(
    list(
      call = NULL, call_env = NULL,
```
with
```r
  obj <- structure(
    list(
      call = NULL, call_env = NULL,
```

and replace the end of that call:
```r
      engine = "sas-outhaz"
    ),
    class = "hazard"
  )
}
```
with
```r
      engine = "sas-outhaz"
    ),
    class = "hazard"
  )

  # What was not done in R, and why (#242). An import was never fitted or
  # examined here, so the record says so rather than printing "none" over a
  # SAS fit.
  record <- .hzr_degraded_record(
    vcov = fit$vcov, weak = fit$weak, control = obj$spec$control,
    dist = "multiphase", fitted = TRUE, imported = TRUE
  )
  obj$degraded <- record$degraded
  obj$degraded_causes <- record$degraded_causes
  .hzr_validate_degraded(obj, fitted = TRUE, imported = TRUE)
  obj
}
```

- [ ] **Step 4: Run the new tests and the SAS-import suites**

Run: `Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && NOT_CRAN=true Rscript -e 'devtools::load_all(); for (f in c("test-degraded-fits.R","test-read-outhaz.R","test-predict-outhaz.R")) testthat::test_file(file.path("tests/testthat", f))'`
Expected: 0 lints; all three files 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/read-outhaz.R tests/testthat/test-degraded-fits.R
git commit -m "feat(degraded): an imported SAS fit records what R did not do (#242)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 7: Always print the block; prove the check can fail

**Files:**
- Modify: `R/hazard_api.R` (`print.hazard()` around 1481-1501; `print.summary.hazard()` around 1660-1684)
- Modify: `tests/testthat/test-weak-direction.R:357-367`
- Modify: `tests/testthat/test-hessian-stability.R:112-125`
- Test: `tests/testthat/test-not-done-print.R`

**Interfaces:**
- Consumes: `.hzr_format_not_done()`, `.hzr_check_not_done_output()` (Task 1); records from Tasks 4-6.
- Produces: the printed block in both printers.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-not-done-print.R`:

```r
# The printed block is checked against the record by parsing the printed text
# back into entries (#242). The check never calls the printer's formatter, so
# a printer that drops or rewrites an entry cannot also rewrite what it is
# checked against. The last three tests make the printer wrong on purpose and
# assert that the check says so: a check that cannot fail proves nothing.

clean_fit <- function() {
  set.seed(7)
  dat <- data.frame(t = stats::rweibull(80, 1.4, 3), d = rep(1L, 80))
  hazard(survival::Surv(t, d) ~ 1, data = dat, dist = "weibull",
         fit = TRUE, theta = c(mu = 1, nu = 1))
}

unfitted <- function() {
  set.seed(7)
  hazard(time = stats::rweibull(80, 1.4, 3), status = rep(1L, 80),
         dist = "weibull", fit = FALSE)
}

printed <- function(x) {
  list(print = capture.output(print(x)),
       summary = capture.output(print(summary(x))))
}

test_that("print() and summary() say 'none' for a clean fit", {
  fit <- clean_fit()
  expect_identical(fit$degraded, character(0))     # guard
  for (out in printed(fit)) {
    expect_true("  Not done in this run: none" %in% out)
    expect_identical(.hzr_check_not_done_output(out, fit), character(0))
  }
})

test_that("print() and summary() show every entry of a degraded record", {
  fit0 <- unfitted()
  expect_length(fit0$degraded, 3L)                 # guard: something to print
  for (out in printed(fit0)) {
    expect_identical(.hzr_check_not_done_output(out, fit0), character(0))
  }
})

test_that("an object from before the record existed says so, not 'none'", {
  fit <- clean_fit()
  fit$degraded <- NULL
  fit$degraded_causes <- NULL
  for (out in printed(fit)) {
    expect_true(paste0("  Not done in this run: not recorded ",
                       "(object built before this record existed)") %in% out)
    expect_identical(.hzr_check_not_done_output(out, fit), character(0))
  }
})

test_that("the check fails when the printer drops an entry", {
  fit0 <- unfitted()
  local_mocked_bindings(.hzr_format_not_done = function(degraded, causes) {
    keep <- degraded[-1]
    c("  Not done in this run:", paste0("    ", keep, ": ", unname(causes[keep])))
  })
  out <- capture.output(print(fit0))
  expect_true("entry missing: fitting" %in% .hzr_check_not_done_output(out, fit0))
})

test_that("the check fails when the printer always says 'none'", {
  fit0 <- unfitted()
  local_mocked_bindings(.hzr_format_not_done = function(degraded, causes) {
    "  Not done in this run: none"
  })
  out <- capture.output(print(summary(fit0)))
  expect_true("'none' printed over a non-empty record" %in%
                .hzr_check_not_done_output(out, fit0))
})

test_that("the check fails when the printer leaves the block out", {
  fit <- clean_fit()
  local_mocked_bindings(.hzr_format_not_done = function(degraded, causes) {
    character(0)
  })
  expect_identical(.hzr_check_not_done_output(capture.output(print(fit)), fit),
                   "heading missing")
})
```

- [ ] **Step 2: Run to confirm they fail**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-not-done-print.R")'`
Expected: the first three tests FAIL, with `.hzr_check_not_done_output()` returning `"heading missing"`. The last three pass already, and that is fine: they test the checker, and they only become meaningful once the printer calls the formatter.

- [ ] **Step 3: Print the block from `print.hazard()`**

In `print.hazard()`, replace:
```r
    cat("  converged:   ", x$fit$converged, "\n")
  }
  invisible(x)
}
```
with
```r
    cat("  converged:   ", x$fit$converged, "\n")
  }
  # Always printed, "none" included (#242).
  cat(.hzr_format_not_done(x$degraded, x$degraded_causes), sep = "\n")
  invisible(x)
}
```

- [ ] **Step 4: Print the block from `print.summary.hazard()`, replacing the two notes**

Remove the `NA` branch of the weak note. Replace:
```r
    cat("\n")
  } else if (length(x$weak) == 1L && is.na(x$weak)) {
    # Say that the ridge check did not run, rather than leaving its silence to
    # be read as a clean bill of health.
    cat(strwrap(paste0(
          "Note: the fit was not examined for a weakly identified ",
          "direction; no usable Hessian was available."),
                width = 76, indent = 2, exdent = 8),
        sep = "\n")
    cat("\n")
  }
```
with
```r
    cat("\n")
  }
```

Replace:
```r
  if (!is.null(x$converged) && !is.na(x$converged) && isFALSE(x$has_vcov)) {
    cat("  Note: standard errors unavailable; the Hessian could not be ",
        "inverted.\n", sep = "")
  }
```
with
```r
  # Always printed, "none" included (#242): a line that appears only on bad
  # news cannot be told from a line its author forgot to write. It replaces
  # the notes that said "not examined for a weakly identified direction" and
  # "standard errors unavailable; the Hessian could not be inverted", both of
  # which could name the wrong cause.
  cat(.hzr_format_not_done(x$degraded, x$degraded_causes), sep = "\n")
```

- [ ] **Step 5: Update the two tests that asserted the removed notes**

In `tests/testthat/test-weak-direction.R`, replace the body of `test_that("summary() says so when the ridge check could not run", {...})`:
```r
  out <- paste(capture.output(print(summary(fit0))), collapse = " ")
  out <- gsub("\\s+", " ", out)
  expect_match(out, "not examined for a weakly identified direction")
```
with
```r
  out <- capture.output(print(summary(fit0)))
  expect_true("    weak_direction_check: model not fitted" %in% out)
  expect_identical(.hzr_check_not_done_output(out, fit0), character(0))
```

In `tests/testthat/test-hessian-stability.R`, in `test_that("summary reports when standard errors are unavailable", {...})`, give the hand-built object a record. Replace:
```r
         data = list(time = 1:5, x = NULL),
         call = quote(hazard())),
    class = "hazard"))
  out <- capture.output(print(s))
  expect_true(any(grepl("standard errors unavailable", out)))
```
with
```r
         data = list(time = 1:5, x = NULL),
         call = quote(hazard()),
         degraded = c("standard_errors", "weak_direction_check"),
         degraded_causes = c(standard_errors = "Hessian not invertible",
                             weak_direction_check = "standard errors unavailable")),
    class = "hazard"))
  out <- capture.output(print(s))
  expect_true("    standard_errors: Hessian not invertible" %in% out)
```

Before editing, confirm that the `data = list(time = 1:5, x = NULL),` / `call = quote(hazard())),` lines are the ones at `test-hessian-stability.R:121-123`. The same pattern appears at lines 49 and 63; edit only the block inside this `test_that()`.

- [ ] **Step 6: Run the print tests, then the full suite**

Run: `Rscript -e 'devtools::document()' && Rscript -e 'lintr::lint_package()' && NOT_CRAN=true Rscript -e 'devtools::load_all(); for (f in c("test-not-done-print.R","test-weak-direction.R","test-hessian-stability.R")) testthat::test_file(file.path("tests/testthat", f))'`
Expected: 0 lints; all three files 0 failures.

Then run: `NOT_CRAN=true Rscript -e 'devtools::test()'`
Expected: 0 failures. Any other test that matched the old note text will fail here; update it the same way as Step 5 and list it in the PR.

- [ ] **Step 7: Commit, then hand-mutate and count the kills**

```bash
git add R/hazard_api.R tests/testthat/test-not-done-print.R tests/testthat/test-weak-direction.R tests/testthat/test-hessian-stability.R
git commit -m "feat(degraded): print() and summary() always say what was not done (#242)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

Do each mutation below only after that commit. For each one: apply it, run
`NOT_CRAN=true Rscript -e 'devtools::test(filter = "degraded|not-done")'`,
record the failure count, then revert it with `git checkout -- <file>`. A
mutant **survives** if the count is 0. A survivor is a finding: report it and
add a test that kills it before moving on.

| # | File | Mutation | Must fail |
|---|---|---|---|
| M1 | `R/degraded.R` | in `.hzr_degraded_record()`, change `if (.hzr_is_na_scalar(weak)) {` to `if (FALSE) {` | builder and fit tests (validator stops) |
| M2 | `R/hazard_api.R` | delete the `cat(.hzr_format_not_done(...))` line in `print.hazard()` | `test-not-done-print.R` |
| M3 | `R/hazard_api.R` | same, in `print.summary.hazard()` | `test-not-done-print.R` |
| M4 | `R/degraded.R` | in `.hzr_format_not_done()`, change `degraded` to `degraded[-1]` in the entries line | `test-degraded.R`, `test-not-done-print.R` |
| M5 | `R/degraded.R` | swap `"model not fitted"` and `"covariance not imported"` in the `standard_errors` branch | `test-degraded.R` |
| M6 | `R/optimizer.R` | delete `se_reason <- "numDeriv not installed and no analytic Hessian"` | `test-degraded-reasons.R`, `test-degraded-fits.R` |

Confirm after the last one: `git status --short` prints nothing.

---

### Task 8: Documentation, NEWS, the full gate, review, PR

**Files:**
- Modify: `R/hazard_api.R` (the `hazard()` `@return` around line 445; the `print.hazard()` and `print.summary.hazard()` roxygen)
- Modify: `NEWS.md`
- Regenerated: `man/hazard.Rd`, `man/print.hazard.Rd`, `man/print.summary.hazard.Rd`

**Interfaces:**
- Consumes: everything above.
- Produces: the PR.

- [ ] **Step 1: `hazard()` `@return`**

Replace:
```r
#'   and \code{engine} (implementation tag, \code{"native-r-m2"}).
```
with
```r
#'   \code{engine} (implementation tag, \code{"native-r-m2"}), and two
#'   fields recording what the fit did not do: \code{degraded}, a character
#'   vector of the steps not performed, in the fixed order
#'   \code{"fitting"}, \code{"standard_errors"},
#'   \code{"conserved_phase_variance"}, \code{"weak_direction_check"},
#'   \code{"conservation_of_events"}, and empty when nothing was lost; and
#'   \code{degraded_causes}, a character vector with the same names giving
#'   the reason for each. \code{print()} and \code{summary()} always show
#'   them as a "Not done in this run" block, which reads "none" when nothing
#'   was lost. \code{fit$fit$weak} is \code{NA} exactly when
#'   \code{"weak_direction_check"} is listed.
```

- [ ] **Step 2: The printers' roxygen**

In `print.hazard()`'s roxygen, replace
`#' number of predictors, distribution, theta vector, and log-likelihood.`
with
```r
#' number of predictors, distribution, theta vector, and log-likelihood,
#' followed by the "Not done in this run" block described in [hazard()].
```

In `print.summary.hazard()`'s roxygen, replace:
```r
#' positive-definite, a note warns that the standard errors may be unreliable;
#' when the Hessian could not be inverted at all, a note reports that standard
#' errors are unavailable.  A further note names the parameters spanning a
#' weakly identified direction when one was found, or records that the check
#' could not run when no Hessian was available.  S3 dispatch only -- users
```
with
```r
#' positive-definite, a note warns that the standard errors may be unreliable,
#' and a further note names the parameters spanning a weakly identified
#' direction when one was found.  A "Not done in this run" block is always
#' printed: it lists each step this fit did not perform, with the reason, and
#' reads "none" when nothing was lost.  S3 dispatch only -- users
```

- [ ] **Step 3: `NEWS.md`**

Under `# TemporalHazard 1.2.10`, if there is no `## New features` subsection, add one directly under the version heading, above `## Bug fixes`. Under it add:

```markdown
* **Every fit now says what it did not do** (#242, following #197). A
  `hazard` object carries `degraded`, the steps the fit did not perform, and
  `degraded_causes`, the reason for each. `print()` and `summary()` always
  show them as a "Not done in this run" block, and the block reads "none"
  when nothing was lost: a line that appears only on bad news cannot be told
  from one that was never written. Five steps are recorded:
  - fitting itself: `fit = FALSE`; a single-distribution `fit = TRUE` call
    with no `theta`, which skips the optimizer without an error or warning;
    or a fit imported from SAS output;
  - standard errors, naming whether numDeriv was missing, `numDeriv::hessian()`
    failed, or the Hessian was non-finite or singular;
  - the variance of the phase that Conservation of Events conserves;
  - the weak-direction check;
  - Conservation of Events itself.

  The block replaces two notes that could name the wrong cause: "not examined
  for a weakly identified direction" and "standard errors unavailable; the
  Hessian could not be inverted". `fit$fit$weak` keeps its meaning, and is
  `NA` exactly when `"weak_direction_check"` is listed. An object saved by an
  earlier version prints "not recorded" rather than "none".
```

`DESCRIPTION` stays at `Version: 1.2.10`. Confirm with `sed -n 4p DESCRIPTION; sed -n 1p NEWS.md`.
Expected: `Version: 1.2.10` and `# TemporalHazard 1.2.10`.

- [ ] **Step 4: The full local gate, in order**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
NOT_CRAN=true Rscript -e 'devtools::test()'
Rscript -e 'spelling::spell_check_package()'
```
Expected: `document()` leaves `git status` showing only `man/` changes from this task. There are 0 lints and 0 failures; record the PASS and SKIP counts. `spell_check_package()` reports no words from this change; add any genuine term to `inst/WORDLIST`. Then commit:

```bash
git add R/hazard_api.R NEWS.md man/ inst/WORDLIST
git commit -m "docs(degraded): document the record and the Not-done block (#242)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

- [ ] **Step 5: `R CMD check --as-cran` from a clean archive of the committed tree, with the manual**

Build under the session scratchpad, not `$TMPDIR/tree`, which collides across sessions. Let `S` be `<scratchpad>/check`:

```bash
rm -rf "$S" && mkdir -p "$S/tree"
git archive HEAD | tar -x -C "$S/tree"
(cd "$S" && R CMD build tree && R CMD check --as-cran TemporalHazard_1.2.10.tar.gz)
tar tzf "$S/TemporalHazard_1.2.10.tar.gz" | grep -iE 'CLAUDE|AGENTS|inst/dev' || echo "clean tarball"
```
Expected: `Status: OK`, and "clean tarball". Record the overall time and the tests and vignettes `[Ns]` timings from `00check.log`, together with the commit.

- [ ] **Step 6: Review**

Dispatch the `r-reviewer` agent over `git diff origin/main...HEAD`. Verify every finding against the code before acting on it. The report is advisory, and a clean report does not stand in for Steps 4 and 5. Fix what holds up, re-run Step 4 (and Step 5 if `R/` changed), and commit.

- [ ] **Step 7: Push the branch and open the PR, then stop**

```bash
git push -u origin feat/degraded-record
gh pr create --base main --title "Record what a fit did not do: degraded / degraded_causes and an always-printed 'Not done in this run' block" --body-file <scratchpad>/pr-body.md
```

The PR body must contain:
- `Closes #242`, and what changed;
- the gate results with counts;
- the check time and timings;
- the mutation table from Task 7 with a kill count for each mutant;
- the two acknowledged gaps: the untested `eigen()` exit, and the unchanged theta-NULL silent skip, now recorded but not fixed;
- the expected `NEWS.md` conflict with #241;
- a closing line `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.

Do not merge. The maintainer merges.
