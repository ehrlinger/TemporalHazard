# `hzr_repeated_events()` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an exported `hzr_repeated_events()` that reproduces the SAS `%repeat` macro in R, so repeated-events
HAZARD jobs can be reproduced without SAS.

**Architecture:** One new file, `R/repeated-events.R`, holding seven small internal stage functions that mirror the
macro's seven DATA steps one-for-one, plus the exported wrapper that chains them. Each stage takes a data frame and
returns a data frame, so the row ladder is observable between stages and each stage is testable alone. Base R only.

**Tech Stack:** R, roxygen2, testthat 3rd edition, lintr. No new dependency.

**Design doc:** `inst/dev/REPEATED-EVENTS-DESIGN.md`. Read it before starting — it records two traps this plan is
shaped around.

## Global Constraints

- **Never push to `main`.** Branch, commit, push the branch, open a PR, stop. Branch is `feat/hzr-repeated-events`,
  already cut from `main` at `f93d00c`.
- **No new dependency.** Base R only. `numDeriv`-style Suggests additions are not in scope here.
- **Version:** `DESCRIPTION` is at `1.2.10`. Bump the **patch** digit only, to `1.2.11`, and update the `NEWS.md`
  top heading to match **by hand** — nothing in this repo greps for the mismatch. Never a fourth digit, never
  `.9000`, never roll minor or major.
- **`stats::` prefixing is the house style** for any `stats` generic.
- **Line length 120** (`.lintr`). `object_usage_linter` and `indentation_linter` are off; `assignment_linter`
  permits `<<-`.
- **No `browser()`, no bare `print()`, no `library()` inside `R/`.**
- **Generated files are never hand-edited:** `man/`, `NAMESPACE` come from `devtools::document()`.
- **Gate order, every time:** `devtools::document()` → `lintr::lint_package()` (0 lints) → `devtools::test()`
  (0 failures) → `spelling::spell_check_package(use_wordlist = TRUE)`. Then once per PR, `R CMD check --as-cran`
  **with** the manual from a clean `git archive` export of a **committed** tree, built under the session scratchpad
  (not `$TMPDIR/tree` — that path collides across concurrent sessions and fabricates an ERROR).
- **Set `NOT_CRAN=true`** when running tests, or `skip_on_cran()` cases hide.

---

## File structure

| Path | Responsibility |
|---|---|
| `R/repeated-events.R` (create) | Validation, the sort/flag helpers, seven stage functions, exported `hzr_repeated_events()` |
| `tests/testthat/test-repeated-events.R` (create) | Unconditional synthetic tests, one fixture per macro branch |
| `NAMESPACE`, `man/` (generated) | via `devtools::document()` |
| `inst/WORDLIST` (modify) | SAS identifiers introduced into roxygen prose |
| `DESCRIPTION`, `NEWS.md` (modify) | patch bump to 1.2.11 |

## Stage-to-macro map

The macro is seven DATA steps. This plan's stage numbering matches the design doc's:

| Stage | Macro source | What it does |
|---|---|---|
| 1 | `data &out; set &in; &rcensor=0;` | add `rcensor`, all zero |
| 2 | `if first.&id=0 and (&eventype=0 or .) then delete` | drop non-first non-events |
| 3 | solo-nonevent pad; `if first.&id>last.&id and nonevent then delete` | pad no-first-event subjects, drop leading non-events |
| 4 | `first=first.&id; last=last.&id; output; ... output;` | record flags, append terminal censored row |
| 5 | `if &rcensor=1 or &eventype=1; ... &event=` | subset, derive `event` |
| 6 | `retain lag_iv 0 number 0; ...` | `iv_start`, `event_no`, `iv_seg`, `rcensor` update, `renewal` |
| 7 | `if &iv_seg=0 and &eventype=0 and (first.&id NE 1) then delete` | drop zero-duration rows |
| 8 | *(the job's own guard, not the macro)* | **BLOCKED — see Task 7** |

## Three SAS subtleties that must survive translation

Get these wrong and the function returns a plausible, wrong data frame with no error.

1. **"Non-event" and "event" are not complements.** Deletion tests `(&eventype=0 or &eventype=.)`; the `event`
   assignment tests `&eventype=1`. An `eventype` of `2` is neither. Two separate predicates, never one negated.
2. **Stage 7 tests `&eventype=0` exactly**, *not* "0 or missing". A missing `eventype` row is not deleted there.
3. **SAS sorts missing numerics first.** `order()` puts `NA` last. Stage sorts must place `NA` times first.

---

### Task 1: Validation, sorting and group flags

**Files:**
- Create: `R/repeated-events.R`
- Create: `tests/testthat/test-repeated-events.R`

**Interfaces:**
- Consumes: nothing.
- Produces: `.hzr_re_order(data, id, time)` returning an integer permutation; `.hzr_re_flags(ids)` returning
  `list(first = logical, last = logical)`; `.hzr_re_nonevent(x)` returning logical;
  `.hzr_re_validate(data, id, time, followup, indicator)` returning `invisible(NULL)` or stopping.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-repeated-events.R`:

```r
test_that(".hzr_re_flags marks group boundaries", {
  f <- .hzr_re_flags(c("a", "a", "b", "c", "c", "c"))
  expect_equal(f$first, c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(f$last, c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
})

test_that(".hzr_re_flags handles a single row and zero rows", {
  expect_equal(.hzr_re_flags("a"), list(first = TRUE, last = TRUE))
  expect_equal(.hzr_re_flags(character(0)), list(first = logical(0), last = logical(0)))
})

test_that(".hzr_re_order sorts by id then time, with missing times first, as SAS does", {
  d <- data.frame(id = c(2, 1, 1, 1), t = c(5, 3, NA, 1))
  expect_equal(.hzr_re_order(d, "id", "t"), c(3L, 4L, 2L, 1L))
})

test_that(".hzr_re_order is stable at tied times", {
  d <- data.frame(id = c("a", "a", "a"), t = c(1, 1, 1), tag = 1:3)
  expect_equal(d$tag[.hzr_re_order(d, "id", "t")], 1:3)
})

test_that(".hzr_re_nonevent follows SAS (&eventype=0 or &eventype=.) and excludes other codes", {
  expect_equal(.hzr_re_nonevent(c(0, 1, NA, 2)), c(TRUE, FALSE, TRUE, FALSE))
})

test_that(".hzr_re_validate rejects a missing column", {
  d <- data.frame(id = 1, t = 1, fu = 1)
  expect_error(.hzr_re_validate(d, "id", "t", "fu", "ev"), "not found")
})

test_that(".hzr_re_validate rejects an input column that would be overwritten by an output column", {
  d <- data.frame(id = 1, t = 1, fu = 1, ev = 1, rcensor = 0)
  expect_error(.hzr_re_validate(d, "id", "t", "fu", "ev"), "rcensor")
})

test_that(".hzr_re_validate rejects a missing id", {
  d <- data.frame(id = c(1, NA), t = 1, fu = 1, ev = 1)
  expect_error(.hzr_re_validate(d, "id", "t", "fu", "ev"), "missing")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-repeated-events.R")'
```

Expected: every test errors with `could not find function ".hzr_re_flags"` (and the other three names).

- [ ] **Step 3: Write the implementation**

Create `R/repeated-events.R`:

```r
# repeated-events.R -- R reimplementation of the SAS %repeat macro
#
# PURPOSE
# -------
# Repeated-events HAZARD jobs build their fit input with the SAS macro
# %repeat (~/Documents/macro.library/repeat.sas, 161 lines, 2003-10-17).
# That input was never saved, so those jobs cannot be reproduced without
# rebuilding it.  This file rebuilds it in R.
#
# The macro is long-in, long-out: its input is already one row per
# candidate event per subject.  The transform is gap-filling and
# segmentation -- keep the event rows, guarantee at least one row per
# subject, append a terminal censored row at end of follow-up when the
# last event is earlier, then lag the event time within subject to turn
# a set of event *times* into a set of left-truncated *intervals*.
#
# SAS/C BRIDGE
# ------------
# Seven internal stage functions mirror the macro's seven DATA steps
# one-for-one, so the row ladder in the reference log is observable
# between stages.  See inst/dev/REPEATED-EVENTS-DESIGN.md.
#
# Three SAS subtleties are load bearing and are commented at their sites:
# "non-event" and "event" are not complements; stage 7 tests eventype=0
# exactly; and SAS sorts missing numerics first where order() sorts them
# last.

# Output columns the macro creates.  An input column of any of these
# names would be silently overwritten, so validation refuses it.
.hzr_re_output_cols <- c(
  "rcensor", "first", "last", "event", "event_no", "iv_start", "iv_seg", "renewal"
)

# SAS: (&eventype=0 or &eventype=.)
#
# NOT the complement of the event test, which is &eventype=1.  An
# eventype of 2 is neither a non-event nor an event.  Keep the two
# predicates separate.
.hzr_re_nonevent <- function(x) {
  is.na(x) | x == 0
}

# SAS: proc sort data=&out; by &id &iv_event;
#
# Missing numerics sort FIRST in SAS and LAST under order(), so missing
# times are mapped to -Inf for the sort key.  method = "radix" is forced
# because character-vector sort order for other methods depends on the
# locale's collating sequence, which would change which row is a subject's
# first row between machines.
.hzr_re_order <- function(data, id, time) {
  key <- data[[time]]
  key[is.na(key)] <- -Inf
  order(data[[id]], key, method = "radix")
}

# SAS: first.&id / last.&id, for data already sorted by id.
.hzr_re_flags <- function(ids) {
  n <- length(ids)
  if (n == 0L) {
    return(list(first = logical(0), last = logical(0)))
  }
  if (n == 1L) {
    return(list(first = TRUE, last = TRUE))
  }
  boundary <- ids[-1L] != ids[-n]
  list(first = c(TRUE, boundary), last = c(boundary, TRUE))
}

.hzr_re_validate <- function(data, id, time, followup, indicator) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  cols <- c(id = id, time = time, followup = followup, indicator = indicator)
  for (arg in names(cols)) {
    value <- cols[[arg]]
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("`%s` must be a single column name.", arg), call. = FALSE)
    }
    if (!value %in% names(data)) {
      stop(sprintf("Column \"%s\" (argument `%s`) not found in `data`.", value, arg), call. = FALSE)
    }
  }
  clash <- intersect(names(data), .hzr_re_output_cols)
  if (length(clash) > 0L) {
    stop(
      sprintf(
        "`data` already has column(s) %s, which this function creates and would overwrite.",
        paste(sprintf("\"%s\"", clash), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (anyNA(data[[id]])) {
    stop(sprintf("Column \"%s\" (argument `id`) has missing values.", id), call. = FALSE)
  }
  invisible(NULL)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-repeated-events.R")'
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 12 ]` or more.

- [ ] **Step 5: Lint**

```bash
Rscript -e 'lintr::lint("R/repeated-events.R")'
```

Expected: no output.

- [ ] **Step 6: Commit**

```bash
git add R/repeated-events.R tests/testthat/test-repeated-events.R
git commit -m "feat(repeated-events): validation, SAS-faithful sort and group flags"
```

---

### Task 2: Stages 1 to 3 — selection and padding

**Files:**
- Modify: `R/repeated-events.R` (append)
- Modify: `tests/testthat/test-repeated-events.R` (append)

**Interfaces:**
- Consumes: `.hzr_re_order()`, `.hzr_re_flags()`, `.hzr_re_nonevent()`.
- Produces: `.hzr_re_stage1(data)`, `.hzr_re_stage2(data, id, time, indicator)`,
  `.hzr_re_stage3(data, id, time, followup, indicator)`. Each takes and returns a data frame; stages 2 and 3 return
  data sorted by `(id, time)`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-repeated-events.R`:

```r
# A four-subject fixture covering the branches stages 2-3 discriminate:
#   s1 -- two events, no leading non-event
#   s2 -- a leading non-event row then one event
#   s3 -- a single non-event row (no first event at all)
#   s4 -- one event followed by a trailing non-event row
re_fixture <- function() {
  data.frame(
    id = c("s1", "s1", "s2", "s2", "s3", "s4", "s4"),
    t  = c(1, 3, 0, 2, NA, 4, 6),
    fu = c(10, 10, 10, 10, 10, 10, 10),
    ev = c(1, 1, 0, 1, 0, 1, 0),
    stringsAsFactors = FALSE
  )
}

test_that("stage 1 adds rcensor as all zero and changes nothing else", {
  d <- re_fixture()
  out <- .hzr_re_stage1(d)
  expect_equal(out$rcensor, rep(0, 7))
  expect_equal(out[names(d)], d)
})

test_that("stage 2 keeps every first row and every event, dropping non-first non-events", {
  out <- .hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev")
  # s2's leading non-event is at t=0 so it sorts first and is KEPT as the
  # group's first row; s4's trailing non-event is dropped.
  expect_equal(out$id, c("s1", "s1", "s2", "s2", "s3", "s4"))
  expect_equal(out$t, c(1, 3, 0, 2, NA, 4))
})

test_that("stage 2 does not treat a non-1 non-0 indicator code as a non-event", {
  d <- data.frame(id = c("a", "a"), t = c(1, 2), fu = c(9, 9), ev = c(1, 2), stringsAsFactors = FALSE)
  out <- .hzr_re_stage2(.hzr_re_stage1(d), "id", "t", "ev")
  expect_equal(nrow(out), 2L)
})

test_that("stage 3 pads a subject that never had a first event", {
  out <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  s3 <- out[out$id == "s3", ]
  expect_equal(nrow(s3), 1L)
  expect_equal(s3$t, 10)
  expect_equal(s3$rcensor, 1)
})

test_that("stage 3 drops a leading non-event when the subject has later rows", {
  out <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  expect_equal(out$id, c("s1", "s1", "s2", "s3", "s4"))
  expect_equal(out$t, c(1, 3, 2, 10, 4))
})
```

- [ ] **Step 2: Run to verify they fail**

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-repeated-events.R")'
```

Expected: the five new tests error with `could not find function ".hzr_re_stage1"` etc.

- [ ] **Step 3: Write the implementation**

Append to `R/repeated-events.R`:

```r
# Stage 1 -- SAS: data &out; set &in; &rcensor=0;
.hzr_re_stage1 <- function(data) {
  data$rcensor <- 0
  data
}

# Stage 2 -- SAS: if first.&id=0 and (&eventype=0 or &eventype=.) then delete;
#
# Keep the group's first row unconditionally, plus every row that is not a
# non-event.  Note the double negative is deliberate: "not a non-event" is
# wider than "is an event", and SAS's predicate is the non-event one.
.hzr_re_stage2 <- function(data, id, time, indicator) {
  data <- data[.hzr_re_order(data, id, time), , drop = FALSE]
  flags <- .hzr_re_flags(data[[id]])
  keep <- flags$first | !.hzr_re_nonevent(data[[indicator]])
  data[keep, , drop = FALSE]
}

# Stage 3 -- SAS:
#   if first.&id=1 and last.&id=1 and nonevent then do; &iv_event=&iv_end; &rcensor=1; end;
#   if first.&id>last.&id and nonevent then delete;
#
# Both statements read the SAME first./last. values, so the flags are
# computed once, before the padding mutation.  The two conditions are
# disjoint (first & last against first & !last), so order does not matter
# between them -- but recomputing the flags after the mutation would be
# wrong in principle and is avoided here on purpose.
.hzr_re_stage3 <- function(data, id, time, followup, indicator) {
  data <- data[.hzr_re_order(data, id, time), , drop = FALSE]
  flags <- .hzr_re_flags(data[[id]])
  nonevent <- .hzr_re_nonevent(data[[indicator]])

  solo <- flags$first & flags$last & nonevent
  data[[time]][solo] <- data[[followup]][solo]
  data$rcensor[solo] <- 1

  drop <- flags$first & !flags$last & nonevent
  data[!drop, , drop = FALSE]
}
```

- [ ] **Step 4: Run to verify they pass**

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-repeated-events.R")'
```

Expected: `FAIL 0`, pass count up by at least 10 assertions.

- [ ] **Step 5: Lint and commit**

```bash
Rscript -e 'lintr::lint("R/repeated-events.R")'
git add R/repeated-events.R tests/testthat/test-repeated-events.R
git commit -m "feat(repeated-events): stages 1-3, row selection and no-first-event padding"
```

---

### Task 3: Stage 4 — group flags and the terminal censored row

**Files:**
- Modify: `R/repeated-events.R` (append)
- Modify: `tests/testthat/test-repeated-events.R` (append)

**Interfaces:**
- Consumes: `.hzr_re_order()`, `.hzr_re_flags()`.
- Produces: `.hzr_re_stage4(data, id, time, followup, indicator)`, returning a data frame with numeric `first` and
  `last` columns added and terminal rows appended. The appended row is a copy of the subject's last row with
  `rcensor = 1`, `time = followup`, `indicator = 0`, and `first`/`last` inherited from the row it was copied from.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-repeated-events.R`:

```r
test_that("stage 4 records first and last as numeric data columns", {
  d <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  out <- .hzr_re_stage4(d, "id", "t", "fu", "ev")
  expect_type(out$first, "double")
  expect_type(out$last, "double")
})

test_that("stage 4 appends a terminal censored row only where the last row is not already censored", {
  d <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  out <- .hzr_re_stage4(d, "id", "t", "fu", "ev")
  # s1, s2 and s4 gain a terminal row; s3 was already rcensor=1 at stage 3.
  expect_equal(as.vector(table(out$id)[c("s1", "s2", "s3", "s4")]), c(3L, 2L, 1L, 2L))
  expect_equal(nrow(out), 8L)
})

test_that("stage 4's appended row carries followup as the time and a zero indicator", {
  d <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  out <- .hzr_re_stage4(d, "id", "t", "fu", "ev")
  added <- out[out$id == "s1" & out$rcensor == 1, ]
  expect_equal(nrow(added), 1L)
  expect_equal(added$t, 10)
  expect_equal(added$ev, 0)
})

test_that("stage 4 puts the appended row after the row it was copied from", {
  d <- .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(re_fixture()), "id", "t", "ev"), "id", "t", "fu", "ev")
  out <- .hzr_re_stage4(d, "id", "t", "fu", "ev")
  s1 <- out[out$id == "s1", ]
  expect_equal(s1$t, c(1, 3, 10))
})
```

- [ ] **Step 2: Run to verify they fail**

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-repeated-events.R")'
```

Expected: `could not find function ".hzr_re_stage4"`.

- [ ] **Step 3: Write the implementation**

Append to `R/repeated-events.R`:

```r
# Stage 4 -- SAS:
#   first=first.&id; last=last.&id;
#   output;
#   if last.&id and &rcensor=0 then do;
#     &rcensor=1; &iv_event=&iv_end; &eventype=0; output;
#   end;
#
# `first` and `last` become ordinary data columns here.  Stage 6's renewal
# assignment reads THESE values, not the automatic variables of its own
# step, and stage 5 subsets the data in between.  They must be carried
# forward, never recomputed.  See REPEATED-EVENTS-DESIGN.md.
.hzr_re_stage4 <- function(data, id, time, followup, indicator) {
  data <- data[.hzr_re_order(data, id, time), , drop = FALSE]
  flags <- .hzr_re_flags(data[[id]])
  data$first <- as.numeric(flags$first)
  data$last <- as.numeric(flags$last)

  extra_rows <- which(flags$last & data$rcensor == 0)
  if (length(extra_rows) == 0L) {
    return(data)
  }
  extra <- data[extra_rows, , drop = FALSE]
  extra$rcensor <- 1
  extra[[time]] <- extra[[followup]]
  extra[[indicator]] <- 0

  # The appended row must land immediately after its source row, as the
  # second SAS `output` does.  Interleave rather than rbind-and-resort:
  # at a tie (last event exactly at &iv_end) a re-sort could place the
  # copy first, which would change which row stage 6 lags from.
  combined <- rbind(data, extra)
  place <- c(seq_len(nrow(data)), extra_rows + 0.5)
  combined[order(place, method = "radix"), , drop = FALSE]
}
```

- [ ] **Step 4: Run to verify they pass**

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-repeated-events.R")'
```

Expected: `FAIL 0`.

- [ ] **Step 5: Lint and commit**

```bash
Rscript -e 'lintr::lint("R/repeated-events.R")'
git add R/repeated-events.R tests/testthat/test-repeated-events.R
git commit -m "feat(repeated-events): stage 4, group flags and the terminal censored row"
```

---

### Task 4: Stages 5 and 6 — the segment derivation

This is the core task and the one where a wrong answer looks right. Read the design doc's "stale `first`/`last`
carry-forward" section before starting.

**Files:**
- Modify: `R/repeated-events.R` (append)
- Modify: `tests/testthat/test-repeated-events.R` (append)

**Interfaces:**
- Consumes: `.hzr_re_order()`, `.hzr_re_flags()`.
- Produces: `.hzr_re_stage5(data, indicator)` and `.hzr_re_stage6(data, id, time, followup)`. After stage 6 the
  frame carries numeric `event`, `event_no`, `iv_start`, `iv_seg`, `renewal`, and an updated `rcensor`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-repeated-events.R`:

```r
re_through_stage4 <- function(d = re_fixture()) {
  .hzr_re_stage4(
    .hzr_re_stage3(.hzr_re_stage2(.hzr_re_stage1(d), "id", "t", "ev"), "id", "t", "fu", "ev"),
    "id", "t", "fu", "ev"
  )
}

test_that("stage 5 derives event from eventype=1 exactly, not from the non-event complement", {
  d <- data.frame(
    id = c("a", "a", "a"), t = c(1, 2, 3), fu = c(9, 9, 9), ev = c(1, 2, 0),
    rcensor = c(0, 1, 1), first = c(1, 0, 0), last = c(0, 0, 1), stringsAsFactors = FALSE
  )
  out <- .hzr_re_stage5(d, "ev")
  expect_equal(out$ev, c(1, 2, 0))
  expect_equal(out$event, c(1, 0, 0))
})

test_that("stage 5 keeps only censored rows and events", {
  d <- data.frame(
    id = c("a", "a"), t = c(1, 2), fu = c(9, 9), ev = c(0, 1),
    rcensor = c(0, 0), first = c(1, 0), last = c(0, 1), stringsAsFactors = FALSE
  )
  out <- .hzr_re_stage5(d, "ev")
  expect_equal(nrow(out), 1L)
  expect_equal(out$t, 2)
})

test_that("stage 6 lags the event time within subject to give iv_start and iv_seg", {
  out <- .hzr_re_stage6(.hzr_re_stage5(re_through_stage4(), "ev"), "id", "t", "fu")
  s1 <- out[out$id == "s1", ]
  expect_equal(s1$t, c(1, 3, 10))
  expect_equal(s1$iv_start, c(0, 1, 3))
  expect_equal(s1$iv_seg, c(1, 2, 7))
})

test_that("stage 6 counts event repeats within subject and restarts at each subject", {
  out <- .hzr_re_stage6(.hzr_re_stage5(re_through_stage4(), "ev"), "id", "t", "fu")
  expect_equal(out$event_no[out$id == "s1"], c(1, 2, 2))
  expect_equal(out$event_no[out$id == "s2"], c(1, 1))
  expect_equal(out$event_no[out$id == "s3"], 0)
})

test_that("stage 6 sets rcensor where the event time equals end of follow-up", {
  d <- data.frame(
    id = c("a"), t = 9, fu = 9, ev = 1, rcensor = 0, first = 1, last = 1, stringsAsFactors = FALSE
  )
  out <- .hzr_re_stage6(.hzr_re_stage5(d, "ev"), "id", "t", "fu")
  expect_equal(out$rcensor, 1)
})

test_that("stage 6 bumps renewal for a subject with no events and for a trailing censored row", {
  out <- .hzr_re_stage6(.hzr_re_stage5(re_through_stage4(), "ev"), "id", "t", "fu")
  # s3 never had an event: first == 1 and event_no == 0, so renewal is bumped to 1.
  expect_equal(out$renewal[out$id == "s3"], 1)
  # s1's trailing censored row: last == 1, rcensor == 1, event == 0, so 2 -> 3.
  expect_equal(out$renewal[out$id == "s1"], c(1, 2, 3))
})

test_that("stage 6 reads the STALE first/last from stage 4, not flags recomputed after stage 5", {
  # Subject "a" has three rows at stage 4, with last == 1 on the third.
  # Stage 5 removes that third row (rcensor 0 and not an event), so the
  # carried flags leave NO row with last == 1, while flags recomputed
  # after the subset would put last == 1 on the second row.  That second
  # row satisfies the rest of the renewal bump (rcensor == 1, event == 0),
  # so the two readings give different renewal values -- which is what
  # makes this test able to fail.
  #
  #   carried    -> last c(0, 0) -> no bump -> renewal c(1, 1)
  #   recomputed -> last c(0, 1) ->    bump -> renewal c(1, 2)
  d <- data.frame(
    id      = c("a", "a", "a"),
    t       = c(1, 5, 7),
    fu      = c(9, 9, 9),
    ev      = c(1, 0, 0),
    rcensor = c(0, 1, 0),
    first   = c(1, 0, 0),
    last    = c(0, 0, 1),
    stringsAsFactors = FALSE
  )
  out <- .hzr_re_stage6(.hzr_re_stage5(d, "ev"), "id", "t", "fu")
  expect_equal(nrow(out), 2L)
  expect_equal(out$last, c(0, 0))
  expect_equal(out$renewal, c(1, 1))
})
```

- [ ] **Step 2: Run to verify they fail**

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-repeated-events.R")'
```

Expected: `could not find function ".hzr_re_stage5"`.

- [ ] **Step 3: Write the implementation**

Append to `R/repeated-events.R`:

```r
# Stage 5 -- SAS:
#   if &rcensor=1 or &eventype=1;
#   if &eventype=1 then &event=1; else &event=0;
#
# A subsetting IF, then the event flag.  &eventype=1 is an exact test:
# a missing eventype fails it, and so does any other code.
.hzr_re_stage5 <- function(data, indicator) {
  is_event <- !is.na(data[[indicator]]) & data[[indicator]] == 1
  data <- data[data$rcensor == 1 | is_event, , drop = FALSE]
  data$event <- as.numeric(!is.na(data[[indicator]]) & data[[indicator]] == 1)
  data
}

# Stage 6 -- SAS:
#   retain lag_iv 0 number 0;
#   &iv_start=0; &event_no=0;
#   if first.&id=1 then do; number=&event; &event_no=&event; end;
#   else do; &iv_start=lag_iv; &event_no=number + &event; number=&event_no; end;
#   lag_iv=&iv_event;
#   if &iv_event=&iv_end then &rcensor=1;
#   &iv_seg=&iv_event-&iv_start;
#   &renewal=&event_no;
#   if (first and &event_no=0) or (last and &rcensor=1 and &event=0)
#      then &renewal=&event_no+1;
#
# The loop is row-sequential because SAS's RETAIN is.  lag_iv and number
# are DATA-step loop state and are deliberately NOT returned -- see
# REPEATED-EVENTS-DESIGN.md, "Column-count discrepancy".
.hzr_re_stage6 <- function(data, id, time, followup) {
  data <- data[.hzr_re_order(data, id, time), , drop = FALSE]
  flags <- .hzr_re_flags(data[[id]])
  n <- nrow(data)

  event_time <- data[[time]]
  event <- data$event
  iv_start <- numeric(n)
  event_no <- numeric(n)
  lag_iv <- 0
  number <- 0

  for (i in seq_len(n)) {
    if (flags$first[i]) {
      number <- event[i]
      event_no[i] <- event[i]
    } else {
      iv_start[i] <- lag_iv
      event_no[i] <- number + event[i]
      number <- event_no[i]
    }
    lag_iv <- event_time[i]
  }

  data$iv_start <- iv_start
  data$event_no <- event_no

  at_end <- !is.na(event_time) & !is.na(data[[followup]]) & event_time == data[[followup]]
  data$rcensor[at_end] <- 1

  data$iv_seg <- event_time - iv_start

  # `first` and `last` here are the stage-4 columns, deliberately stale.
  # `rcensor` is the value just updated on the line above, not stage 5's.
  data$renewal <- event_no
  bump <- (data$first == 1 & event_no == 0) |
    (data$last == 1 & data$rcensor == 1 & data$event == 0)
  data$renewal[bump] <- event_no[bump] + 1

  data
}
```

- [ ] **Step 4: Run to verify they pass**

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-repeated-events.R")'
```

Expected: `FAIL 0`.

- [ ] **Step 5: Kill a mutant — prove the stale-flag test can fail**

Temporarily insert `flags <- .hzr_re_flags(data[[id]]); data$first <- as.numeric(flags$first); data$last <-
as.numeric(flags$last)` immediately before the `data$renewal <- event_no` line, so the flags are recomputed rather
than carried. Re-run the tests.

Expected: the test `"stage 6 reads the STALE first/last from stage 4"` **FAILS**, reporting `renewal` as
`c(1, 2)` against the expected `c(1, 1)`. If it passes, the test is not testing what it claims and must be
strengthened before you continue — this is the single most important check in the plan, because the carried and
recomputed flags agree on most fixtures and a test built on one of those constrains nothing. Remove the temporary
lines afterwards and re-run to confirm green.

- [ ] **Step 6: Lint and commit**

```bash
Rscript -e 'lintr::lint("R/repeated-events.R")'
git add R/repeated-events.R tests/testthat/test-repeated-events.R
git commit -m "feat(repeated-events): stages 5-6, segment derivation and modulated renewal"
```

---

### Task 5: Stage 7, the exported wrapper, and documentation

**Files:**
- Modify: `R/repeated-events.R` (append)
- Modify: `tests/testthat/test-repeated-events.R` (append)
- Generated: `NAMESPACE`, `man/hzr_repeated_events.Rd`

**Interfaces:**
- Consumes: every stage function.
- Produces: exported `hzr_repeated_events(data, id, time, followup, indicator)` returning a data frame.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-repeated-events.R`:

```r
test_that("stage 7 drops a zero-duration non-event row that is not the subject's first", {
  d <- data.frame(
    id = c("a", "a"), t = c(0, 5), fu = c(9, 9), ev = c(0, 0),
    rcensor = c(1, 1), first = c(1, 0), last = c(0, 1),
    event = c(0, 0), event_no = c(0, 0), iv_start = c(0, 5), iv_seg = c(0, 0),
    renewal = c(1, 1), stringsAsFactors = FALSE
  )
  out <- .hzr_re_stage7(d, "id", "t", "ev")
  expect_equal(nrow(out), 1L)
  expect_equal(out$iv_start, 0)
})

test_that("stage 7 tests eventype=0 exactly and keeps a zero-duration row with a missing indicator", {
  d <- data.frame(
    id = c("a", "a"), t = c(0, 5), fu = c(9, 9), ev = c(0, NA),
    rcensor = c(1, 1), first = c(1, 0), last = c(0, 1),
    event = c(0, 0), event_no = c(0, 0), iv_start = c(0, 5), iv_seg = c(0, 0),
    renewal = c(1, 1), stringsAsFactors = FALSE
  )
  expect_equal(nrow(.hzr_re_stage7(d, "id", "t", "ev")), 2L)
})

test_that("hzr_repeated_events returns the documented columns and drops SAS loop state", {
  out <- hzr_repeated_events(re_fixture(), "id", "t", "fu", "ev")
  expect_true(all(c("event", "event_no", "rcensor", "iv_start", "iv_seg", "renewal", "first", "last")
                  %in% names(out)))
  expect_false(any(c("lag_iv", "number") %in% names(out)))
  expect_equal(names(out)[1:4], c("id", "t", "fu", "ev"))
})

test_that("hzr_repeated_events gives every subject at least one row", {
  out <- hzr_repeated_events(re_fixture(), "id", "t", "fu", "ev")
  expect_setequal(unique(out$id), c("s1", "s2", "s3", "s4"))
})

test_that("hzr_repeated_events produces segments that tile each subject's follow-up without gaps", {
  out <- hzr_repeated_events(re_fixture(), "id", "t", "fu", "ev")
  for (subject in unique(out$id)) {
    rows <- out[out$id == subject, ]
    expect_equal(rows$iv_start, c(0, utils::head(rows$t, -1)), info = subject)
    expect_equal(rows$iv_seg, rows$t - rows$iv_start, info = subject)
    expect_equal(max(rows$t), 10, info = subject)
  }
})

test_that("hzr_repeated_events rejects a column name clash rather than overwriting", {
  d <- re_fixture()
  d$renewal <- 0
  expect_error(hzr_repeated_events(d, "id", "t", "fu", "ev"), "renewal")
})
```

- [ ] **Step 2: Run to verify they fail**

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-repeated-events.R")'
```

Expected: `could not find function ".hzr_re_stage7"` and `could not find function "hzr_repeated_events"`.

- [ ] **Step 3: Write the implementation**

Append to `R/repeated-events.R`:

```r
# Stage 7 -- SAS:
#   if &iv_seg=0 and &eventype=0 and (first.&id NE 1) then delete;
#
# Note &eventype=0 EXACTLY.  Unlike stages 2 and 3 this is not the
# "0 or missing" non-event test, so a zero-duration row with a missing
# indicator survives here.
.hzr_re_stage7 <- function(data, id, time, indicator) {
  data <- data[.hzr_re_order(data, id, time), , drop = FALSE]
  flags <- .hzr_re_flags(data[[id]])
  drop <- data$iv_seg == 0 &
    !is.na(data[[indicator]]) & data[[indicator]] == 0 &
    !flags$first
  data[!drop, , drop = FALSE]
}

#' Build repeated-event segments, reproducing the SAS `%repeat` macro
#'
#' Converts a long data set of candidate event times into one row per
#' inter-event segment, ready to fit as a repeated-events model. This is a
#' native R implementation of the SAS macro `%repeat` used by the
#' repeated-events HAZARD jobs, whose fit input was never saved.
#'
#' @details
#' The input holds one row per candidate event per subject, with an indicator
#' naming the event of interest: a value of 1 marks an event, and 0 or `NA`
#' marks its absence. For each subject the function keeps the event rows,
#' guarantees at least one row, appends a censored row at the end of follow-up
#' when the last event happened earlier, and then lags the event time within
#' subject so that each row describes the interval since the previous event.
#'
#' Two details are inherited from the macro and are not R conventions. An
#' indicator code that is neither 0, 1 nor missing is treated as neither an
#' event nor an absence, and rows are ordered with missing times first, as SAS
#' sorts them.
#'
#' @param data A data frame with one row per candidate event per subject.
#' @param id Name of the column identifying the subject.
#' @param time Name of the column giving the interval from time zero to each
#'   event.
#' @param followup Name of the column giving the interval from time zero to the
#'   end of follow-up.
#' @param indicator Name of the event indicator column: 1 marks an event, 0 or
#'   `NA` marks its absence.
#'
#' @return A data frame with one row per inter-event segment. Every input
#'   column is retained, with `time` and `indicator` altered where the macro
#'   alters them, and the following columns added:
#'   \describe{
#'     \item{`event`}{1 when the row is an event of interest, otherwise 0.}
#'     \item{`event_no`}{Running count of event repeats within the subject.}
#'     \item{`rcensor`}{1 when the row is right censored, otherwise 0.}
#'     \item{`iv_start`}{Interval from time zero to the start of the segment.}
#'     \item{`iv_seg`}{Duration of the segment, `time` minus `iv_start`.}
#'     \item{`renewal`}{Segment number under the modulated renewal
#'       formulation.}
#'     \item{`first`, `last`}{Whether the row was the subject's first or last
#'       before censored rows were appended.}
#'   }
#'
#' @examples
#' events <- data.frame(
#'   id = c("s1", "s1", "s2", "s3"),
#'   t = c(1, 3, 2, NA),
#'   fu = c(10, 10, 10, 10),
#'   ev = c(1, 1, 1, 0)
#' )
#' hzr_repeated_events(events, id = "id", time = "t", followup = "fu", indicator = "ev")
#'
#' @export
hzr_repeated_events <- function(data, id, time, followup, indicator) {
  .hzr_re_validate(data, id, time, followup, indicator)

  data <- .hzr_re_stage1(data)
  data <- .hzr_re_stage2(data, id, time, indicator)
  data <- .hzr_re_stage3(data, id, time, followup, indicator)
  data <- .hzr_re_stage4(data, id, time, followup, indicator)
  data <- .hzr_re_stage5(data, indicator)
  data <- .hzr_re_stage6(data, id, time, followup)
  data <- .hzr_re_stage7(data, id, time, indicator)

  row.names(data) <- NULL
  data
}
```

- [ ] **Step 4: Document, then run the tests**

```bash
Rscript -e 'devtools::document()'
NOT_CRAN=true Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-repeated-events.R")'
```

Expected: `NAMESPACE` gains `export(hzr_repeated_events)`, `man/hzr_repeated_events.Rd` is created, and `FAIL 0`.

- [ ] **Step 5: Run the example, do not just read it**

```bash
Rscript -e 'devtools::load_all(quiet = TRUE); print(example("hzr_repeated_events", package = "TemporalHazard", character.only = TRUE, give.lines = FALSE))'
```

Expected: a data frame, four subjects represented, no error and no warning. A roxygen example that has never been
run is not evidence that it works.

- [ ] **Step 6: Lint and commit**

```bash
Rscript -e 'lintr::lint_package()'
git add R/repeated-events.R tests/testthat/test-repeated-events.R NAMESPACE man/
git commit -m "feat(repeated-events): stage 7 and the exported hzr_repeated_events()"
```

---

### Task 6: Version bump, spelling, and the full gate

**Files:**
- Modify: `DESCRIPTION:3`
- Modify: `NEWS.md:1`
- Modify: `inst/WORDLIST`

- [ ] **Step 1: Bump the patch version in both places**

Edit `DESCRIPTION` line 3 to `Version: 1.2.11`. Edit the `NEWS.md` top heading to `# TemporalHazard 1.2.11` and add,
under a `## New features` section:

```markdown
* `hzr_repeated_events()` builds repeated-event segments from a long data set
  of candidate event times, reproducing the SAS `%repeat` macro used by the
  repeated-events HAZARD jobs.
```

Nothing in this repo greps for a `DESCRIPTION`/`NEWS.md` mismatch, so this is a hand check.

- [ ] **Step 2: Run the spelling check**

```bash
Rscript -e 'spelling::spell_check_package(use_wordlist = TRUE)'
```

Expected: flags SAS identifiers appearing in the new roxygen prose. `spelling.yaml` is a **required blocking check**
on `main` even though `AGENTS.md`'s definition-of-done omits it.

- [ ] **Step 3: Add only the genuinely-flagged words to the wordlist**

Append the words the previous step actually reported to `inst/WORDLIST`, keeping the file sorted. Do not add words
speculatively — an unused wordlist entry is dead weight and hides a later real hit.

```bash
Rscript -e 'spelling::spell_check_package(use_wordlist = TRUE)'
```

Expected: no output.

- [ ] **Step 4: Run the full local gate in order**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
NOT_CRAN=true Rscript -e 'devtools::test()'
Rscript -e 'spelling::spell_check_package(use_wordlist = TRUE)'
```

Expected: 0 lints; 0 failures; skips only from SAS fixture availability; no spelling output. Record the actual pass
and skip counts — do not report "tests pass" without them.

- [ ] **Step 5: Commit, then build and check from a committed tree**

```bash
git add DESCRIPTION NEWS.md inst/WORDLIST
git commit -m "chore(repeated-events): bump to 1.2.11, add SAS identifiers to the wordlist"

TREE="$(mktemp -d)/tree"
mkdir -p "$TREE"
git archive HEAD | tar -x -C "$TREE"
R CMD build "$TREE"
R CMD check --as-cran TemporalHazard_1.2.11.tar.gz
```

Expected: `Status: OK`. Build from the archive, not the working tree — an empty `inst/doc` fabricates two vignette
WARNINGs and a worktree's `.git` file lands in the tarball as a spurious hidden-files NOTE. `mktemp -d` is used
instead of `AGENTS.md`'s `$TMPDIR/tree` because that fixed path collides across concurrent sessions in this tree.

- [ ] **Step 6: Confirm developer files stayed out of the tarball**

```bash
tar tzf TemporalHazard_1.2.11.tar.gz | grep -iE 'CLAUDE|AGENTS|inst/dev' || echo "clean"
```

Expected: `clean`.

- [ ] **Step 7: Run the r-reviewer agent over the diff**

This change adds an exported function, so `AGENTS.md` requires it. Treat the findings as advisory — verify each one
against the code, and do not let a clean report substitute for the commands above.

- [ ] **Step 8: Push the branch and open a PR**

```bash
git push -u origin feat/hzr-repeated-events
gh pr create --title "feat: hzr_repeated_events(), an R reimplementation of the SAS %repeat macro" --body "..."
```

Then **stop**. The maintainer merges. Note that eight required checks green still leaves the PR `BLOCKED` until
someone approves it, and Copilot's review comes back `COMMENTED`, never `APPROVED`.

---

### Task 7: BLOCKED — stage 8 and the parity test

**Do not attempt this task while `/Volumes/qhsstudies` is unmounted.** It cannot be done correctly without the study
volume, and guessing at it would produce exactly the failure mode this package is prone to: a parity test that looks
like evidence and is not.

Three things must be read from the volume first, per `inst/dev/REPEATED-EVENTS-DESIGN.md` Open Items:

1. **Confirm which `repeat.sas` ran.** Diff the log's echoed source against
   `~/Documents/macro.library/repeat.sas`. This settles identity only — if it is a different macro, Tasks 2 to 5 are
   built on the wrong stage list and must be revisited.
2. **Read the guard predicate** from `hz.ce_cardioversion_repeated.ehb.sas`. It drops one row of 963, and which row
   depends on whether the test is `iv_start < iv_end`, a zero test, or a missing test. The maintainer chose to
   include this guard inside `hzr_repeated_events()` rather than leave it to the caller, so it becomes stage 8 and
   the function returns 962 rather than 963.
3. **Derive the ladder-to-stage mapping** from the log's step boundaries. The log records seven shapes and the macro
   has seven DATA steps, but they do not line up: stages 2 and 3 both leave `709 x 358`, and stages 6 and 7 both
   leave `963 x 367`, so at least two shapes are repeats. Until the mapping is derived from the log itself, only the
   endpoints are anchored.

Only then write the gated parity test, `skip_if` on the presence of the volume, asserting **row coverage before any
value comparison**: the ladder, then 962 rows, 388 events, 574 right censored, 387 left censored, `IV_EVENT` within
`[0.0001140795, 12.99137]` and `IV_START` within `[0.0001140795, 9.716832]`, and `ncol` equal to 367 minus the two
named loop-state columns. Also assert whether any tied `(id, iv_event)` pairs exist in `bd_card`, so the sort
stability caveat is settled with evidence rather than assumed inert.

A gated test prints green when it never ran, so it is not evidence on its own. The unconditional tests from Tasks 1
to 5 carry the burden.

---

## Self-review notes

- **Spec coverage.** Contract → Task 5. Fixed output names → Task 5 (no rename arguments). `renewal` in v1 → Task 4.
  Column-count discrepancy → Task 5 (`lag_iv`/`number` dropped, asserted) and Task 7 (the 367 assertion). Stages
  1–7 → Tasks 1–5. Stage 8 → Task 7, blocked. Stale-flag trap → Task 4 Steps 1 and 5. Sort stability → Task 1
  (`method = "radix"`) and Task 7 (the tie assertion). Synthetic tests → Tasks 1–5. Gated parity test → Task 7.
  Gates and mechanics → Task 6.
- **Not covered, by design.** Wiring into `hzr_translate_sas()` and running the cardioversion fit are out of scope
  per the spec.
- **Naming consistency.** `.hzr_re_stage1` … `.hzr_re_stage7`, `.hzr_re_order`, `.hzr_re_flags`,
  `.hzr_re_nonevent`, `.hzr_re_validate`, `.hzr_re_output_cols` used identically in every task.
