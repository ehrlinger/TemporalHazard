# bh.dead SAS Parity Harness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a development-only harness that compares `TemporalHazard`'s multiphase fit and `hzr_bootstrap(scope=)` embedded stepwise screen against the SAS/C HAZARD originals (`%HAZARD`, `%HAZBOOT`) on a clinical analysis held on a secure volume.

**Architecture:** Two SAS listing files are parsed into a single fixture list (`shape`, `early`, `constant`, `meta`). The fixture is written to the secure volume — never the repo. A loader resolves it via environment variables. A testthat parity test refits both SAS specifications in R and asserts the deterministic shape fit tightly and the stochastic bootstrap screen statistically, deriving every expectation from the fixture at run time. A Quarto notebook on the volume is the human-facing debug layer.

**Tech Stack:** R, testthat 3e, `haven` (Suggests, for `read_sas`), base R parsing (no new dependencies).

Design document: `inst/dev/BHDEAD-SAS-PARITY-DESIGN.md`

## Global Constraints

- **This repository is PUBLIC.** `.Rbuildignore` excludes files from the CRAN tarball but does nothing about git. No study values (estimates, selection frequencies, cohort size) may appear in any committed file — including code, comments, tests, and commit messages.
- **No PHI in the repo.** `bhblt.sas7bdat` must never be copied into `temporal_hazard`.
- **The built fixture stays off-repo**, beside the data on the secure volume.
- **Every expectation is derived from the fixture at run time.** No SAS number is hard-coded.
- **Test fixtures must use invented numbers, not renamed real ones.** Reproduce the SAS listing *format*, never its values. Renaming a variable to `VAR_A` while keeping its verbatim counts and percentages anonymises nothing — the values themselves are the study result and are trivially matched back against the listing. Use round, obviously fake numbers (`500`, `50.0`, `1.000`).
- **Harness code lives in `inst/dev/bhdead-parity/`** (`.Rbuildignore`d via the existing `^inst/dev$` rule). The only exception is the test, which must live in `tests/testthat/` to run.
- **The test must be inert for everyone else:** skip unless data AND fixture are readable, and unless `haven` is installed. It must never run on CRAN or CI.
- **No new hard dependencies.** `haven` goes in `Suggests`, guarded by `skip_if_not_installed()`.
- **Style:** follow existing repo conventions — internal helpers are `.hzr_`-prefixed with `@noRd`; two-space indent; `<-` assignment; lintr line limit 120.
- **Versioning:** do not touch the MINOR/MAJOR digit. This work is test-only and needs no version bump.
- Branch: `test/bhdead-sas-parity-harness`. Never push to `main`/`dev` directly; open a PR.

## Environment variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `TEMPORALHAZARD_BHBLT` | *(none — must be set)* | Row-level SAS dataset |
| `TEMPORALHAZARD_BHDEAD_FIXTURE` | *(none — must be set)* | Built SAS reference |

There is **no default path**. A study-identifying directory must not appear in
this public repository, so an unset variable yields `""` and the loader returns
`NULL`, which is how callers skip. Set both (e.g. in your `.Rprofile`).

---

### Task 1: Fixture schema, validator, and loader

**Files:**
- Create: `inst/dev/bhdead-parity/bhdead-fixture.R`
- Test: `tests/testthat/test-bhdead-fixture.R`

**Interfaces:**
- Consumes: nothing (first task).
- Produces:
  - `.hzr_bhdead_fixture_schema()` → named list: `shape = c("specified","converged")`, `phase = c("name","n","pct","min","max","mean","sd")`, `meta = c("resampl","sle","sls","n_obs","captured_on","source")`
  - `.hzr_validate_bhdead_fixture(fix)` → returns `fix` unchanged if valid; `stop()`s with a message listing every problem otherwise.
  - `.hzr_bhdead_fixture_path()` → `Sys.getenv("TEMPORALHAZARD_BHDEAD_FIXTURE", unset = "")`. No default; `""` when unset.
  - `.hzr_bhblt_path()` → `Sys.getenv("TEMPORALHAZARD_BHBLT", unset = "")`. No default; `""` when unset.
  - `.hzr_load_bhdead_fixture(path = .hzr_bhdead_fixture_path())` → validated fixture list, or `NULL` when the file does not exist.
  - `.hzr_bhdead_candidates(fix)` → character vector of candidate variable names (lowercased, intercept row removed).

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-bhdead-fixture.R`:

```r
# Unit tests for the bh.dead parity fixture schema/loader.
# These construct synthetic fixtures in-memory and never touch the secure
# volume, so they run everywhere (including CRAN) and contain no study values.

source(system.file("dev", "bhdead-parity", "bhdead-fixture.R",
                   package = "TemporalHazard"), local = TRUE)

.fake_phase <- function(names_vec) {
  data.frame(
    name = names_vec,
    n    = rep(10L, length(names_vec)),
    pct  = rep(50, length(names_vec)),
    min  = rep(-1, length(names_vec)),
    max  = rep(1, length(names_vec)),
    mean = rep(0, length(names_vec)),
    sd   = rep(1, length(names_vec)),
    stringsAsFactors = FALSE
  )
}

.fake_fixture <- function() {
  list(
    shape = list(
      specified = c(thalf = 1, nu = -1, m = 0, mue = 1, muc = 1),
      converged = c(thalf = 1, nu = -1, m = 0, mue = 1, muc = 1)
    ),
    early    = .fake_phase(c("E0", "VAR_A", "VAR_B")),
    constant = .fake_phase(c("C0", "VAR_A", "VAR_B")),
    meta = list(resampl = 1000L, sle = 0.12, sls = 0.10, n_obs = 100L,
                captured_on = "2026-07-16", source = "synthetic")
  )
}

test_that(".hzr_validate_bhdead_fixture accepts a well-formed fixture", {
  fix <- .fake_fixture()
  expect_identical(.hzr_validate_bhdead_fixture(fix), fix)
})

test_that(".hzr_validate_bhdead_fixture reports every missing field at once", {
  fix <- .fake_fixture()
  fix$meta$sle <- NULL
  fix$shape$converged <- NULL
  expect_error(.hzr_validate_bhdead_fixture(fix), "sle")
  expect_error(.hzr_validate_bhdead_fixture(fix), "converged")
})

test_that(".hzr_validate_bhdead_fixture rejects a phase missing a column", {
  fix <- .fake_fixture()
  fix$early$pct <- NULL
  expect_error(.hzr_validate_bhdead_fixture(fix), "early")
})

test_that(".hzr_bhdead_candidates drops the intercept row and lowercases", {
  expect_identical(.hzr_bhdead_candidates(.fake_fixture()), c("var_a", "var_b"))
})

test_that(".hzr_load_bhdead_fixture returns NULL when the file is absent", {
  expect_null(.hzr_load_bhdead_fixture(file.path(tempdir(), "nope.rds")))
})

test_that(".hzr_load_bhdead_fixture round-trips a valid fixture", {
  p <- file.path(tempdir(), "bhdead-test.rds")
  on.exit(unlink(p), add = TRUE)
  saveRDS(.fake_fixture(), p)
  expect_identical(.hzr_load_bhdead_fixture(p), .fake_fixture())
})

test_that("path helpers honour their environment variables", {
  withr::with_envvar(c(TEMPORALHAZARD_BHBLT = "/tmp/x.sas7bdat"), {
    expect_identical(.hzr_bhblt_path(), "/tmp/x.sas7bdat")
  })
  withr::with_envvar(c(TEMPORALHAZARD_BHDEAD_FIXTURE = "/tmp/y.rds"), {
    expect_identical(.hzr_bhdead_fixture_path(), "/tmp/y.rds")
  })
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-bhdead-fixture.R")'`
Expected: FAIL — `cannot open file` / the sourced script does not exist.

- [ ] **Step 3: Write the implementation**

Create `inst/dev/bhdead-parity/bhdead-fixture.R`:

```r
# bhdead-fixture.R -- Resolve, load and validate the bh.dead SAS parity fixture.
#
# DEVELOPMENT ONLY. inst/dev/ is .Rbuildignore'd and does not ship.
#
# This repository is public. The fixture itself is NOT stored here -- it lives
# beside the source data on a secure volume, resolved via the environment
# variables below. This file contains no study values; every expectation used
# by the parity test is read out of the fixture at run time.
#
# Build the fixture with inst/dev/bhdead-parity/parse-bhdead-lst.R.
# See inst/dev/bhdead-parity/README.md.

#' Path to the row-level SAS dataset
#'
#' From `TEMPORALHAZARD_BHBLT`. No default: a study-identifying path must not
#' live in this public repo. Unset yields "", which callers treat as absent.
#' @noRd
.hzr_bhblt_path <- function() {
  Sys.getenv("TEMPORALHAZARD_BHBLT", unset = "")
}

#' Path to the built parity fixture
#'
#' From `TEMPORALHAZARD_BHDEAD_FIXTURE`. No default, as above.
#' @noRd
.hzr_bhdead_fixture_path <- function() {
  Sys.getenv("TEMPORALHAZARD_BHDEAD_FIXTURE", unset = "")
}

#' Required fields of a bh.dead parity fixture
#' @noRd
.hzr_bhdead_fixture_schema <- function() {
  list(
    shape = c("specified", "converged"),
    phase = c("name", "n", "pct", "min", "max", "mean", "sd"),
    meta  = c("resampl", "sle", "sls", "n_obs", "captured_on", "source")
  )
}

#' Validate a bh.dead parity fixture against the schema
#'
#' @param fix A list as produced by the parse script.
#' @return `fix`, unchanged, if valid. Throws otherwise, listing every problem.
#' @noRd
.hzr_validate_bhdead_fixture <- function(fix) {
  schema <- .hzr_bhdead_fixture_schema()
  problems <- character()

  if (!is.list(fix)) {
    stop("bh.dead fixture must be a list.", call. = FALSE)
  }

  for (group in c("shape", "meta")) {
    if (!group %in% names(fix)) {
      problems <- c(problems, paste0("missing top-level: ", group))
      next
    }
    missing_fields <- setdiff(schema[[group]], names(fix[[group]]))
    if (length(missing_fields) > 0L) {
      problems <- c(problems, sprintf("missing %s: %s", group,
                                      paste(missing_fields, collapse = ", ")))
    }
  }

  for (phase in c("early", "constant")) {
    if (!phase %in% names(fix)) {
      problems <- c(problems, paste0("missing top-level: ", phase))
      next
    }
    missing_cols <- setdiff(schema$phase, names(fix[[phase]]))
    if (length(missing_cols) > 0L) {
      problems <- c(problems, sprintf("missing %s columns: %s", phase,
                                      paste(missing_cols, collapse = ", ")))
    }
  }

  if (length(problems) > 0L) {
    stop("Invalid bh.dead fixture:\n  - ",
         paste(problems, collapse = "\n  - "), call. = FALSE)
  }
  fix
}

#' Load the bh.dead parity fixture
#'
#' @param path Path to the `.rds`. Defaults to `.hzr_bhdead_fixture_path()`.
#' @return The validated fixture, or `NULL` when the file does not exist.
#' @noRd
.hzr_load_bhdead_fixture <- function(path = .hzr_bhdead_fixture_path()) {
  if (!nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  .hzr_validate_bhdead_fixture(readRDS(path))
}

#' Candidate variable names from a fixture (intercept dropped, lowercased)
#'
#' The intercept rows are SAS's `E0` (early) and `C0` (constant). Everything
#' else in the phase table is a screened candidate.
#'
#' @noRd
.hzr_bhdead_candidates <- function(fix) {
  nm <- union(fix$early$name, fix$constant$name)
  nm <- setdiff(nm, c("E0", "C0"))
  tolower(nm)
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-bhdead-fixture.R")'`
Expected: PASS, 8 assertions, 0 failures.

- [ ] **Step 5: Add `withr` to Suggests if absent**

Check: `grep -n "withr" DESCRIPTION`
If absent, add `withr,` to the `Suggests:` block (alphabetical, after `testthat`).

- [ ] **Step 6: Commit**

```bash
git add inst/dev/bhdead-parity/bhdead-fixture.R tests/testthat/test-bhdead-fixture.R DESCRIPTION
git commit -m "test(dev): bh.dead parity fixture schema, validator and loader

Development-only harness (inst/dev is .Rbuildignore'd). The fixture is
resolved from a secure volume via TEMPORALHAZARD_BHDEAD_FIXTURE and is
never stored in this public repo. Unit tests build synthetic fixtures
in-memory, so they run everywhere and carry no study values."
```

---

### Task 2: Parse the SAS listings into a fixture

**Files:**
- Create: `inst/dev/bhdead-parity/parse-bhdead-lst.R`
- Create: `inst/dev/bhdead-parity/README.md`
- Test: `tests/testthat/test-bhdead-parse.R`

**Interfaces:**
- Consumes: `.hzr_validate_bhdead_fixture()` from Task 1.
- Produces:
  - `.hzr_parse_bhdead_phase(lines, phase_label)` → data.frame with columns `name, n, pct, min, max, mean, sd`. `phase_label` is `"Early Phase"` or `"Constant Phase"`.
  - `.hzr_parse_bhdead_shape(lines)` → list with `specified` and `converged`, each a named numeric of `thalf, nu, m, mue, muc`.
  - `.hzr_build_bhdead_fixture(hz_lst, bh_lst, meta)` → validated fixture list.

The phase tables have this **shape** (leading `Obs` index, then the fields).
Column positions are not stable, so split on whitespace rather than fixed
widths. **All numbers shown below are invented** — never paste real listing
values into committed files:

```
                           Obs    _NAME_         N     PCT          MIN        MAX        MEAN        STD

                             1    E0           500    100.0      -1.000      1.000      0.1000      1.000
                             2    VAR_A        250     50.0      -2.000      2.000      0.2000      2.000
```

The shape block in the other listing has this shape (again, invented numbers):

```
                                          Phase       Parameter     Fixed?     Estimate
                                                      THALF         No            1.200000
```

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-bhdead-parse.R`:

```r
# Unit tests for the bh.dead .lst parsers.
#
# Inputs below reproduce the SAS listing FORMAT with entirely INVENTED numbers.
# Never paste real listing values here: this repository is public, and a
# renamed variable with verbatim values is not anonymised. Round, obviously
# fake numbers make that impossible to get wrong by accident.

source(system.file("dev", "bhdead-parity", "parse-bhdead-lst.R",
                   package = "TemporalHazard"), local = TRUE)

.fake_bh_lines <- c(
  "                        Study Title Redacted                       3",
  "                                                  Early Phase",
  "",
  "                           Obs    _NAME_         N     PCT          MIN        MAX        MEAN        STD",
  "",
  "                             1    E0           500    100.0    -100.000    100.000     -1.0000     10.000",
  "                             2    VAR_A        250     50.0      -2.000      2.000      0.2000      2.000",
  "                             3    VAR_B        100     20.0      -3.000      3.000      0.3000      3.000",
  "                        Study Title Redacted                       4",
  "                                                  Constant Phase",
  "",
  "                           Obs    _NAME_         N     PCT          MIN        MAX        MEAN        STD",
  "",
  "                             1    C0           500    100.0    -200.000    200.000     -2.0000     20.000",
  "                             2    VAR_A        150     30.0      -4.000      4.000      0.4000      4.000"
)

.fake_hz_lines <- c(
  "          Phase        Parameter       Specified        Used       |  Theta",
  "                       THALF             1.100000      1.100000    |  E2",
  "                       NU                 -0.5000       -0.5000    |  E3",
  "                       M                        0             0    |  E4",
  "                       MUE                 0.8000        0.8000    |  E0",
  "          Constant:    MUC               0.050000      0.050000    |  C0",
  "                                          Estimates for Model Parameters",
  "                                          Phase       Parameter     Fixed?     Estimate",
  "                                                      THALF         No            1.200000",
  "                                                      NU            No           -0.600000",
  "                                                      M             Yes                  0",
  "                                                      MUE                        0.850000",
  "                                          Constant:   MUC                       0.055000"
)

test_that(".hzr_parse_bhdead_phase extracts the Early table", {
  d <- .hzr_parse_bhdead_phase(.fake_bh_lines, "Early Phase")
  expect_identical(d$name, c("E0", "VAR_A", "VAR_B"))
  expect_identical(d$n, c(500L, 250L, 100L))
  expect_equal(d$pct, c(100.0, 50.0, 20.0))
  expect_equal(d$mean[2], 0.2)
  expect_equal(d$sd[3], 3.0)
})

test_that(".hzr_parse_bhdead_phase does not bleed across phase boundaries", {
  d <- .hzr_parse_bhdead_phase(.fake_bh_lines, "Constant Phase")
  expect_identical(d$name, c("C0", "VAR_A"))
  expect_equal(d$pct, c(100.0, 30.0))
})

test_that(".hzr_parse_bhdead_shape reads specified and converged values", {
  s <- .hzr_parse_bhdead_shape(.fake_hz_lines)
  expect_equal(s$specified[["thalf"]], 1.1)
  expect_equal(s$specified[["mue"]], 0.8)
  expect_equal(s$converged[["thalf"]], 1.2)
  expect_equal(s$converged[["nu"]], -0.6)
  expect_equal(s$converged[["m"]], 0)
  expect_equal(s$converged[["muc"]], 0.055)
})

test_that(".hzr_build_bhdead_fixture returns a schema-valid fixture", {
  hz <- tempfile(fileext = ".lst"); bh <- tempfile(fileext = ".lst")
  on.exit(unlink(c(hz, bh)), add = TRUE)
  writeLines(.fake_hz_lines, hz); writeLines(.fake_bh_lines, bh)
  fix <- .hzr_build_bhdead_fixture(
    hz, bh,
    meta = list(resampl = 1000L, sle = 0.12, sls = 0.10, n_obs = 100L,
                captured_on = "2026-07-16", source = "synthetic")
  )
  expect_identical(.hzr_validate_bhdead_fixture(fix), fix)
  expect_identical(fix$early$name[1], "E0")
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-bhdead-parse.R")'`
Expected: FAIL — the sourced script does not exist.

- [ ] **Step 3: Write the implementation**

Create `inst/dev/bhdead-parity/parse-bhdead-lst.R`:

```r
# parse-bhdead-lst.R -- Build the bh.dead parity fixture from SAS listings.
#
# DEVELOPMENT ONLY. inst/dev/ is .Rbuildignore'd and does not ship.
#
# This repository is public: the .lst inputs and the .rds output both stay on
# the secure volume. This script contains parsing logic only -- no study values.
#
# Usage (paths default to the secure volume; override as needed):
#   source("inst/dev/bhdead-parity/parse-bhdead-lst.R")
#   .hzr_write_bhdead_fixture()
#
# See README.md in this directory.

# Companion helpers (path resolution, schema, validator). Resolved through the
# installed/loaded package rather than a relative path, so this script can be
# sourced from any working directory.
source(system.file("dev", "bhdead-parity", "bhdead-fixture.R",
                   package = "TemporalHazard"))

#' Extract one phase table from a bootstrap listing
#'
#' Rows look like: `<obs> <NAME> <N> <PCT> <MIN> <MAX> <MEAN> <STD>`.
#' Column positions are not stable across SAS listings, so split on
#' whitespace and take fields positionally.
#'
#' @param lines Character vector -- the whole listing.
#' @param phase_label "Early Phase" or "Constant Phase".
#' @return data.frame(name, n, pct, min, max, mean, sd)
#' @noRd
.hzr_parse_bhdead_phase <- function(lines, phase_label) {
  starts <- grep(phase_label, lines, fixed = TRUE)
  if (length(starts) == 0L) {
    stop("Phase label not found in listing: ", phase_label, call. = FALSE)
  }
  # A listing repeats the phase label on each page header; take from the first
  # occurrence to the first *other* phase label that follows it.
  other <- if (identical(phase_label, "Early Phase")) {
    "Constant Phase"
  } else {
    "Late Phase"
  }
  from <- starts[1]
  ends <- grep(other, lines, fixed = TRUE)
  ends <- ends[ends > from]
  to <- if (length(ends) > 0L) ends[1] - 1L else length(lines)

  block <- lines[from:to]
  # Data rows: obs index, then a NAME of uppercase/digits/underscore, then 6 numbers.
  rx <- paste0("^\\s*([0-9]+)\\s+([A-Z0-9_]+)\\s+([0-9]+)\\s+",
               "(-?[0-9.]+)\\s+(-?[0-9.]+)\\s+(-?[0-9.]+)\\s+",
               "(-?[0-9.]+)\\s+(-?[0-9.]+)\\s*$")
  hits <- regmatches(block, regexec(rx, block))
  hits <- Filter(function(x) length(x) == 9L, hits)
  if (length(hits) == 0L) {
    stop("No data rows parsed for ", phase_label, call. = FALSE)
  }

  data.frame(
    name = vapply(hits, `[[`, character(1), 3L),
    n    = as.integer(vapply(hits, `[[`, character(1), 4L)),
    pct  = as.numeric(vapply(hits, `[[`, character(1), 5L)),
    min  = as.numeric(vapply(hits, `[[`, character(1), 6L)),
    max  = as.numeric(vapply(hits, `[[`, character(1), 7L)),
    mean = as.numeric(vapply(hits, `[[`, character(1), 8L)),
    sd   = as.numeric(vapply(hits, `[[`, character(1), 9L)),
    stringsAsFactors = FALSE
  )
}

#' Extract specified and converged shape parameters from a shape-fit listing
#'
#' @param lines Character vector -- the whole listing.
#' @return list(specified = named numeric, converged = named numeric), each
#'   over thalf, nu, m, mue, muc.
#' @noRd
.hzr_parse_bhdead_shape <- function(lines) {
  # "Specified" block: `<PARAM> <specified> <used> | ...`
  spec_of <- function(param) {
    rx <- paste0("^\\s*(?:Constant:\\s+)?", param,
                 "\\s+(-?[0-9.]+)\\s+(-?[0-9.]+)\\s*\\|")
    hit <- regmatches(lines, regexec(rx, lines))
    hit <- Filter(function(x) length(x) == 3L, hit)
    if (length(hit) == 0L) {
      stop("Specified value not found for ", param, call. = FALSE)
    }
    as.numeric(hit[[1]][2])
  }

  # "Estimates for Model Parameters" block, below its header.
  est_from <- grep("Estimates for Model Parameters", lines, fixed = TRUE)
  if (length(est_from) == 0L) {
    stop("'Estimates for Model Parameters' block not found.", call. = FALSE)
  }
  est_block <- lines[est_from[1]:length(lines)]
  conv_of <- function(param) {
    # `<PARAM> [Fixed?] <estimate>` -- the estimate is the last number on the line.
    rx <- paste0("^\\s*(?:Constant:\\s+)?", param, "\\s+.*?(-?[0-9.]+)\\s*$")
    hit <- regmatches(est_block, regexec(rx, est_block))
    hit <- Filter(function(x) length(x) == 2L, hit)
    if (length(hit) == 0L) {
      stop("Converged value not found for ", param, call. = FALSE)
    }
    as.numeric(hit[[1]][2])
  }

  params <- c(thalf = "THALF", nu = "NU", m = "M", mue = "MUE", muc = "MUC")
  list(
    specified = vapply(params, spec_of, numeric(1)),
    converged = vapply(params, conv_of, numeric(1))
  )
}

#' Build the fixture from the two listings
#'
#' @param hz_lst Path to the shape-fit listing.
#' @param bh_lst Path to the bootstrap listing.
#' @param meta Named list: resampl, sle, sls, n_obs, captured_on, source.
#' @return A validated fixture list.
#' @noRd
.hzr_build_bhdead_fixture <- function(hz_lst, bh_lst, meta) {
  hz <- readLines(hz_lst, warn = FALSE)
  bh <- readLines(bh_lst, warn = FALSE)
  fix <- list(
    shape    = .hzr_parse_bhdead_shape(hz),
    early    = .hzr_parse_bhdead_phase(bh, "Early Phase"),
    constant = .hzr_parse_bhdead_phase(bh, "Constant Phase"),
    meta     = meta
  )
  .hzr_validate_bhdead_fixture(fix)
}

#' Build and write the fixture to the secure volume
#'
#' @noRd
.hzr_write_bhdead_fixture <- function(
    hz_lst = file.path(.hzr_bhdead_default_root, "distributions", "hz.dead.lst"),
    bh_lst = file.path(.hzr_bhdead_default_root, "analyses", "bh.dead.lst"),
    out    = .hzr_bhdead_fixture_path(),
    meta   = list(resampl = 1000L, sle = 0.12, sls = 0.10, n_obs = NA_integer_,
                  captured_on = as.character(Sys.Date()),
                  source = "hz.dead.lst + bh.dead.lst")) {
  fix <- .hzr_build_bhdead_fixture(hz_lst, bh_lst, meta)
  saveRDS(fix, out)
  message("Wrote fixture: ", out,
          " (early: ", nrow(fix$early), " rows, constant: ",
          nrow(fix$constant), " rows)")
  invisible(fix)
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-bhdead-parse.R")'`
Expected: PASS, 10 assertions, 0 failures.

- [ ] **Step 5: Write the README**

Create `inst/dev/bhdead-parity/README.md`:

```markdown
# bh.dead SAS parity harness (development only)

`inst/dev/` is `.Rbuildignore`d: nothing here ships to CRAN.

**This repository is public.** Neither the row-level data nor the study's
aggregate results live here. The fixture is built from SAS listings on a
secure volume and written back to that volume. Committed code contains no
study values; the parity test derives every expectation from the fixture at
run time and skips when the volume is not mounted.

## Paths

| Environment variable | Purpose |
|----------------------|---------|
| `TEMPORALHAZARD_BHBLT` | Row-level SAS dataset (`bhblt.sas7bdat`) |
| `TEMPORALHAZARD_BHDEAD_FIXTURE` | Built fixture (`bhdead.rds`) |

Both default to the secure-volume locations. Set them to point elsewhere.

## Rebuilding the fixture

With the volume mounted:

```r
devtools::load_all()
source("inst/dev/bhdead-parity/parse-bhdead-lst.R")
.hzr_write_bhdead_fixture()
```

This reads the shape-fit listing and the bootstrap listing, and writes
`bhdead.rds` next to the data.

## Running the parity test

```r
devtools::load_all()
testthat::test_file("tests/testthat/test-bhdead-sas-parity.R")
```

It skips unless the data and the fixture are both readable and `haven` is
installed.

## What the SAS side does

* Shape fit: `PROC HAZARD ... noconserve`, fixing only `M`.
* Bootstrap: `%hazboot(seed=-1, resampl=1000, sle=0.12, sls=0.1)`, base model
  fixing `nu` and `m`.

`seed=-1` means SAS seeds from the time of day, so its selection frequencies
are one random realisation and cannot be reproduced exactly. Bootstrap
assertions are therefore statistical, not exact.
```

- [ ] **Step 6: Commit**

```bash
git add inst/dev/bhdead-parity/parse-bhdead-lst.R inst/dev/bhdead-parity/README.md tests/testthat/test-bhdead-parse.R
git commit -m "test(dev): parse SAS listings into the bh.dead parity fixture

Parsers for the bootstrap phase tables and the shape-fit estimate block,
plus a writer that emits the fixture to the secure volume. Unit tests use
synthetic listing fragments in the SAS output format, so they run anywhere
and carry no study values."
```

---

### Task 3: Build the real fixture and calibrate

**Files:**
- None committed. This task produces the off-repo fixture and records observed behaviour.

**Interfaces:**
- Consumes: `.hzr_write_bhdead_fixture()` from Task 2.
- Produces: `bhdead.rds` on the secure volume; the numbers needed to choose Task 4's thresholds.

This task is a **checkpoint, not a code change**. Its purpose is to replace
guessed thresholds with observed ones.

- [ ] **Step 1: Build the fixture against the real listings**

```r
devtools::load_all()
source("inst/dev/bhdead-parity/parse-bhdead-lst.R")
fix <- .hzr_write_bhdead_fixture()
```

Expected: a message reporting `early: 93 rows, constant: 93 rows`.

- [ ] **Step 2: Sanity-check the parse**

```r
stopifnot(nrow(fix$early) == 93L, nrow(fix$constant) == 93L)
stopifnot(fix$early$name[1] == "E0", fix$constant$name[1] == "C0")
stopifnot(length(.hzr_bhdead_candidates(fix)) == 92L)
stopifnot(all(is.finite(fix$early$pct)), all(fix$early$pct <= 100))
stopifnot(all(c("thalf","nu","m","mue","muc") %in% names(fix$shape$converged)))
```

Expected: all pass silently. **If the candidate count is not 92, stop and
re-read the parser** — the scope is the whole point of the comparison.

- [ ] **Step 3: Confirm every candidate exists in the data**

```r
caa <- haven::read_sas(.hzr_bhblt_path())
names(caa) <- tolower(names(caa))
missing <- setdiff(.hzr_bhdead_candidates(fix), names(caa))
print(missing)
```

Expected: `character(0)`. Any name printed here must be resolved before
Task 4 — a missing column would otherwise degrade into 92 per-candidate
warnings.

- [ ] **Step 4: Run the shape fit and record the gap**

```r
sp <- fix$shape$specified
shape_fit <- hazard(
  survival::Surv(iv_dead, dead) ~ 1, data = caa, dist = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = sp[["thalf"]], nu = sp[["nu"]],
                         m = sp[["m"]], fixed = "m"),
    constant = hzr_phase("constant")
  ),
  theta   = c(early.log_mu = log(sp[["mue"]]),
              early.log_t_half = log(sp[["thalf"]]),
              early.nu = sp[["nu"]], early.m = sp[["m"]],
              constant.log_mu = log(sp[["muc"]])),
  control = list(conserve = FALSE),
  fit = TRUE
)
th <- coef(shape_fit)
data.frame(
  param = c("thalf", "nu", "mue", "muc"),
  r     = c(exp(th[["early.log_t_half"]]), th[["early.nu"]],
            exp(th[["early.log_mu"]]), exp(th[["constant.log_mu"]])),
  sas   = fix$shape$converged[c("thalf", "nu", "mue", "muc")]
)
```

Record the relative gaps. **If any exceeds 1e-3, do not loosen the tolerance
in Task 4 — investigate.** A gap here means R and SAS disagree on a
deterministic fit, which is a finding.

- [ ] **Step 5: Run the bootstrap screen at n_boot = 5, then 50**

```r
scope_vars <- .hzr_bhdead_candidates(fix)
scope <- list(
  early    = reformulate(scope_vars),
  constant = reformulate(scope_vars)
)
base_fit <- hazard(
  survival::Surv(iv_dead, dead) ~ 1, data = caa, dist = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = sp[["thalf"]], nu = sp[["nu"]],
                         m = sp[["m"]], fixed = c("nu", "m")),
    constant = hzr_phase("constant")
  ),
  theta   = c(early.log_mu = log(sp[["mue"]]),
              early.log_t_half = log(sp[["thalf"]]),
              early.nu = sp[["nu"]], early.m = sp[["m"]],
              constant.log_mu = log(sp[["muc"]])),
  control = list(conserve = FALSE),
  fit = TRUE
)
bs5 <- hzr_bootstrap(base_fit, n_boot = 5, seed = 123, scope = scope,
                     slentry = fix$meta$sle, slstay = fix$meta$sls,
                     control = list(n_starts = 1, conserve = FALSE),
                     verbose = TRUE)
bs5$n_success
```

Expected at `n_boot = 5`: `n_success > 0`. Then repeat at `n_boot = 50` and
time it — this sets expectations for the 1000 run.

- [ ] **Step 6: Compute the statistics Task 4 will assert**

```r
cmp <- merge(
  data.frame(name = toupper(sub("^early\\.", "", bs50$summary$parameter)),
             r_pct = bs50$summary$pct, stringsAsFactors = FALSE),
  fix$early[, c("name", "pct")], by = "name"
)
cor(cmp$r_pct, cmp$pct, method = "spearman")
top10 <- head(fix$early$name[order(-fix$early$pct)], 10)
sum(top10 %in% cmp$name[cmp$r_pct > 0])
```

Record both. These observed values choose Task 4's floors.

- [ ] **Step 7: Record findings in the design doc**

Append a "Calibration (observed)" section to
`inst/dev/BHDEAD-SAS-PARITY-DESIGN.md` describing, **without study values**:
the shape-fit gap magnitude (order only, e.g. "within 1e-4"), the observed
rank correlation band, the top-10 recovery count, and the chosen floors with
their rationale.

- [ ] **Step 8: Commit**

```bash
git add inst/dev/BHDEAD-SAS-PARITY-DESIGN.md
git commit -m "docs(dev): record bh.dead parity calibration and chosen thresholds

Observed behaviour from the first real run against the secure volume,
recorded without study values. Fixes the floors used by the parity test."
```

---

### Task 4: The parity test

**Files:**
- Create: `tests/testthat/test-bhdead-sas-parity.R`
- Modify: `DESCRIPTION` (add `haven` to `Suggests`)

**Interfaces:**
- Consumes: `.hzr_load_bhdead_fixture()`, `.hzr_bhblt_path()`, `.hzr_bhdead_candidates()` from Task 1; the floors recorded in Task 3.
- Produces: nothing consumed downstream.

**Scope change (2026-07-16, owner decision).** The statistical assertions are
DEFERRED — there are no `<SPEARMAN_FLOOR>` / `<TOP10_FLOOR>` values to fill in.
Task 3 measured one bootstrap replicate over the 92-variable pool at ~5.4
minutes, so SAS's `resampl=1000` is ~90 hours; the feasible `n_boot = 5` gives a
`pct` restricted to {0, 20, 40, 60, 80, 100}, and a rank correlation of that
against SAS's 1000-resample percentages is dominated by ties and noise. A floor
calibrated from it would encode noise while reading as verified parity.

Task 4 therefore asserts: the shape fit tightly (deterministic, measured at
3.4e-04 against SAS), and the bootstrap screen as a **smoke check only**.
Statistical parity is unblocked by `criterion = "score"` — see the design doc's
"Deferred: score-test selection". Do not invent a floor to fill the gap.

- [ ] **Step 1: Write the test**

Create `tests/testthat/test-bhdead-sas-parity.R`:

```r
# SAS parity for the bh.dead analysis: multiphase shape fit (%HAZARD) and the
# hzr_bootstrap(scope=) embedded stepwise screen (%HAZBOOT).
#
# THIS REPOSITORY IS PUBLIC. No study value appears in this file: the data and
# the SAS reference both live on a secure volume, and every expectation is read
# from the fixture at run time. The test skips unless both are readable, so it
# never runs on CRAN, on CI, or on any machine without the volume mounted.
#
# Enable by building the fixture first:
#   source("inst/dev/bhdead-parity/parse-bhdead-lst.R"); .hzr_write_bhdead_fixture()
# See inst/dev/bhdead-parity/README.md.

.bhdead_helpers <- system.file("dev", "bhdead-parity", "bhdead-fixture.R",
                                package = "TemporalHazard")

# Shape-fit tolerance is calibrated from observed output: R reproduces SAS's
# converged shape fit to 3.4e-04 max relative error (2026-07-16), so 1e-3 is a
# real guard with headroom, not a guess.
#
# There is deliberately NO statistical tolerance for the bootstrap screen. SAS
# seeds from the time of day (seed=-1), so its selection frequencies are one
# random realisation; and one R replicate over the 92-variable pool costs ~5.4
# minutes, making a SAS-scale n_boot infeasible until criterion = "score"
# lands. At the feasible n_boot the frequencies are too coarse to support a
# meaningful floor. See BHDEAD-SAS-PARITY-DESIGN.md, "Deferred: score-test
# selection". The bootstrap test below is a smoke check by design.
.bhdead_tolerance <- list(
  shape_rel = 1e-3
)

.bhdead_setup <- function() {
  testthat::skip_if_not_installed("haven")
  if (!nzchar(.bhdead_helpers) || !file.exists(.bhdead_helpers)) {
    testthat::skip("bh.dead parity helpers not installed")
  }
  source(.bhdead_helpers, local = parent.frame())
  fix <- .hzr_load_bhdead_fixture()
  testthat::skip_if(is.null(fix), "bh.dead SAS fixture not available")
  data_path <- .hzr_bhblt_path()
  testthat::skip_if_not(file.exists(data_path), "bhblt data not available")
  caa <- haven::read_sas(data_path)
  names(caa) <- tolower(names(caa))
  list(fix = fix, caa = caa)
}

# Build the phase list for a given `fixed` specification, seeding every start
# value from the fixture.
.bhdead_phases <- function(sp, fixed) {
  list(
    early    = hzr_phase("cdf", t_half = sp[["thalf"]], nu = sp[["nu"]],
                         m = sp[["m"]], fixed = fixed),
    constant = hzr_phase("constant")
  )
}

.bhdead_theta <- function(sp) {
  c(early.log_mu     = log(sp[["mue"]]),
    early.log_t_half = log(sp[["thalf"]]),
    early.nu         = sp[["nu"]],
    early.m          = sp[["m"]],
    constant.log_mu  = log(sp[["muc"]]))
}

test_that("multiphase shape fit reproduces the SAS shape fit", {
  env <- .bhdead_setup()
  sp <- env$fix$shape$specified
  cv <- env$fix$shape$converged

  # hz.dead fixes only M; THALF and NU are estimated. noconserve -> conserve = FALSE.
  fit <- hazard(
    survival::Surv(iv_dead, dead) ~ 1, data = env$caa, dist = "multiphase",
    phases  = .bhdead_phases(sp, fixed = "m"),
    theta   = .bhdead_theta(sp),
    control = list(conserve = FALSE),
    fit     = TRUE
  )
  expect_true(isTRUE(fit$fit$converged))

  th <- coef(fit)
  expect_equal(exp(th[["early.log_t_half"]]), cv[["thalf"]],
               tolerance = .bhdead_tolerance$shape_rel)
  expect_equal(th[["early.nu"]], cv[["nu"]],
               tolerance = .bhdead_tolerance$shape_rel)
  expect_equal(exp(th[["early.log_mu"]]), cv[["mue"]],
               tolerance = .bhdead_tolerance$shape_rel)
  expect_equal(exp(th[["constant.log_mu"]]), cv[["muc"]],
               tolerance = .bhdead_tolerance$shape_rel)
})

test_that("every SAS candidate exists in the data", {
  env <- .bhdead_setup()
  expect_identical(
    setdiff(.hzr_bhdead_candidates(env$fix), names(env$caa)),
    character(0)
  )
})

test_that("bootstrap screen runs end-to-end on the full SAS candidate pool", {
  env <- .bhdead_setup()
  sp <- env$fix$shape$specified
  scope_vars <- .hzr_bhdead_candidates(env$fix)

  # bh.dead fixes nu and m; THALF stays free.
  base_fit <- hazard(
    survival::Surv(iv_dead, dead) ~ 1, data = env$caa, dist = "multiphase",
    phases  = .bhdead_phases(sp, fixed = c("nu", "m")),
    theta   = .bhdead_theta(sp),
    control = list(conserve = FALSE),
    fit     = TRUE
  )

  bs <- hzr_bootstrap(
    base_fit, n_boot = 5, seed = 123,
    scope   = list(early = stats::reformulate(scope_vars),
                   constant = stats::reformulate(scope_vars)),
    slentry = env$fix$meta$sle, slstay = env$fix$meta$sls,
    control = list(n_starts = 1, conserve = FALSE)
  )

  expect_s3_class(bs, "hzr_bootstrap")
  expect_gt(bs$n_success, 0)
  expect_gt(nrow(bs$summary), 0)
  # The intercepts are never dropped by stepwise, so they appear in every
  # successful replicate.
  intercepts <- bs$summary[bs$summary$parameter %in%
                             c("early.log_mu", "constant.log_mu"), ]
  expect_true(all(intercepts$n == bs$n_success))
})
```

- [ ] **Step 2: Run the test with the volume mounted**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-bhdead-sas-parity.R")'`
Expected: PASS. If the shape-fit assertions fail, that is a **finding** —
report it, do not widen `shape_rel`.

- [ ] **Step 3: Verify it skips cleanly without the volume**

Run:
```bash
TEMPORALHAZARD_BHBLT=/nonexistent TEMPORALHAZARD_BHDEAD_FIXTURE=/nonexistent \
  Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-bhdead-sas-parity.R")'
```
Expected: all tests SKIP, 0 failures. This is what CRAN and CI will see.

- [ ] **Step 4: Add `haven` to Suggests**

In `DESCRIPTION`, add `haven,` to `Suggests:` (alphabetical, after `ggplot2`).

- [ ] **Step 5: Run the full suite**

Run: `Rscript -e 'devtools::test()'`
Expected: 0 failures; the three bh.dead parity tests skip or pass.

- [ ] **Step 6: Commit**

```bash
git add tests/testthat/test-bhdead-sas-parity.R DESCRIPTION
git commit -m "test(dev): bh.dead SAS parity test for shape fit and bootstrap screen

Refits both SAS specifications in R -- the shape fit (fixed = 'm',
conserve = FALSE) asserted tightly, and the bootstrap screen over the full
SAS candidate pool as a smoke check. Skips unless the data and fixture are
readable, so it never runs on CRAN or CI. No study values are committed:
every start value and expectation is read from the fixture at run time."
```

---

### Task 5: Rewrite the analysis notebook

**Files:**
- Modify: `<volume>/analyses/bh.dead_r_bootstrap.qmd` (on the secure volume — **not** committed)

**Interfaces:**
- Consumes: `.hzr_load_bhdead_fixture()`, `.hzr_bhdead_candidates()` from Task 1.
- Produces: the R-vs-SAS comparison table (the debugging deliverable).

This file lives on the secure volume and is not committed. It mirrors the SAS
step for step and fixes the three fidelity defects.

- [ ] **Step 1: Replace the notebook**

```markdown
---
title: "BH Death — R vs SAS bootstrap variable screening"
date: today
format: html
params:
  n_boot: 50
---

Mirrors `distributions/hz.dead.sas` (shape fit) and `analyses/bh.dead.sas`
(bootstrap screen). Escalate `n_boot` 5 → 50 → 1000; SAS used `resampl=1000`.

```{r}
#| label: setup
library(TemporalHazard)
library(haven)

source(system.file("dev", "bhdead-parity", "bhdead-fixture.R",
                   package = "TemporalHazard"))

fix <- .hzr_load_bhdead_fixture()
stopifnot(!is.null(fix))

caa <- read_sas(.hzr_bhblt_path())
names(caa) <- tolower(names(caa))

sp <- fix$shape$specified
```

## Step 1 — shape fit (`hz.dead.sas`)

SAS runs `PROC HAZARD ... noconserve` and fixes **only M**: THALF and NU are
estimated. Hence `fixed = "m"` and `conserve = FALSE`.

```{r}
#| label: hz-shape
theta0 <- c(early.log_mu = log(sp[["mue"]]),
            early.log_t_half = log(sp[["thalf"]]),
            early.nu = sp[["nu"]], early.m = sp[["m"]],
            constant.log_mu = log(sp[["muc"]]))

shape_fit <- hazard(
  survival::Surv(iv_dead, dead) ~ 1, data = caa, dist = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = sp[["thalf"]], nu = sp[["nu"]],
                         m = sp[["m"]], fixed = "m"),
    constant = hzr_phase("constant")
  ),
  theta = theta0, control = list(conserve = FALSE), fit = TRUE
)
summary(shape_fit)

th <- coef(shape_fit)
data.frame(
  param = c("thalf", "nu", "mue", "muc"),
  R     = c(exp(th[["early.log_t_half"]]), th[["early.nu"]],
            exp(th[["early.log_mu"]]), exp(th[["constant.log_mu"]])),
  SAS   = as.numeric(fix$shape$converged[c("thalf", "nu", "mue", "muc")])
) |> transform(rel_diff = abs(R - SAS) / abs(SAS))
```

## Step 2 — bootstrap base (`bh.dead.sas` `parms … fixnu fixm`)

SAS's bootstrap base fixes **nu and m**, leaving THALF free.

```{r}
#| label: bh-base
base_fit <- hazard(
  survival::Surv(iv_dead, dead) ~ 1, data = caa, dist = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = sp[["thalf"]], nu = sp[["nu"]],
                         m = sp[["m"]], fixed = c("nu", "m")),
    constant = hzr_phase("constant")
  ),
  theta = theta0, control = list(conserve = FALSE), fit = TRUE
)
summary(base_fit)
```

## Step 3 — candidate scope, straight from the SAS listing

The pool is read from the fixture, never hand-written: `bh.dead.sas`'s macro
comments out several covariates via a nested-comment trap (SAS `/* */` does
not nest), and a hand-copied list silently disagrees with what SAS ran.

```{r}
#| label: scope
scope_vars <- .hzr_bhdead_candidates(fix)
missing_vars <- setdiff(scope_vars, names(caa))
if (length(missing_vars) > 0) {
  stop("Candidates absent from the data: ", paste(missing_vars, collapse = ", "))
}
length(scope_vars)

scope <- list(early    = reformulate(scope_vars),
              constant = reformulate(scope_vars))
```

## Step 4 — bootstrap screen (`%hazboot`)

```{r}
#| label: bh-bootstrap-screen
bs <- hzr_bootstrap(
  base_fit, n_boot = params$n_boot, seed = 123,
  scope   = scope,
  slentry = fix$meta$sle, slstay = fix$meta$sls,
  control = list(n_starts = 1, conserve = FALSE),
  verbose = TRUE
)
bs$n_success
bs$n_failed
```

## Step 5 — R vs SAS, biggest gaps first

SAS seeds from the time of day (`seed=-1`), so its percentages are one random
realisation. Compare ranks and the high-frequency set, not exact values.

```{r}
#| label: compare
r_pct <- function(prefix) {
  rows <- grep(paste0("^", prefix, "\\."), bs$summary$parameter)
  data.frame(
    name  = toupper(sub(paste0("^", prefix, "\\."), "", bs$summary$parameter[rows])),
    r_pct = bs$summary$pct[rows],
    stringsAsFactors = FALSE
  )
}

compare_phase <- function(prefix, sas_tbl) {
  out <- merge(r_pct(prefix), sas_tbl[, c("name", "pct")], by = "name",
               all = TRUE)
  names(out)[names(out) == "pct"] <- "sas_pct"
  out$r_pct[is.na(out$r_pct)] <- 0
  out$gap <- out$r_pct - out$sas_pct
  out[order(-abs(out$gap)), ]
}

early_cmp <- compare_phase("early", fix$early)
head(early_cmp, 20)
```

```{r}
#| label: rank-agreement
cc <- early_cmp[is.finite(early_cmp$sas_pct), ]
cor(cc$r_pct, cc$sas_pct, method = "spearman")
```
```

- [ ] **Step 2: Render at n_boot = 5**

Run: `quarto render bh.dead_r_bootstrap.qmd -P n_boot:5`
Expected: renders; Step 3 reports 92 candidates; Step 4 reports `n_success > 0`.

- [ ] **Step 3: Render at n_boot = 50**

Run: `quarto render bh.dead_r_bootstrap.qmd -P n_boot:50`
Expected: renders; the Step 5 comparison table is populated.

- [ ] **Step 4: Render at n_boot = 1000**

Run: `quarto render bh.dead_r_bootstrap.qmd -P n_boot:1000`
Expected: matches SAS's `resampl=1000`. Long-running; the progress bar shows
movement. This is the comparison of record.

- [ ] **Step 5: No commit**

The notebook lives on the secure volume and is deliberately not committed.
Confirm the repo is clean: `git status --short` shows nothing from this task.

---

## Self-review

**Spec coverage:**

| Spec requirement | Task |
|------------------|------|
| Parse `hz.dead.lst` → shape reference | 2 |
| Parse `bh.dead.lst` → Early/Constant tables | 2 |
| Fixture written off-repo to the volume | 2 (`.hzr_write_bhdead_fixture`) |
| 92-var pool derived, never hand-transcribed | 1 (`.hzr_bhdead_candidates`), verified in 3 |
| Env-var resolution + skip gate | 1, 4 |
| `haven` → Suggests, guarded | 4 |
| Shape fit asserted tightly (~1e-3) | 4 |
| Bootstrap asserted (smoke only; statistical parity DEFERRED — infeasible n_boot, see Task 4 scope change) | 4 |
| Thresholds calibrated, not invented | 3 — shape_rel 1e-3 from measured 3.4e-04; no bootstrap floor invented |
| `n_boot` staged 5 → 50 → 1000 | 4 (5, smoke). 50/1000 PARKED: ~5.4 min/replicate → 1000 ≈ 90 h. Unblocked by `criterion = "score"`. |
| qmd fixes `fixed=`, `conserve=`, scope | 5 |
| Fail loud on missing columns | 4 (test), 5 (qmd `stop()`) |
| No study values committed | Global constraint; enforced in 1, 2, 4 |
| No PHI in repo | Global constraint; data read from volume only |

**Placeholder scan:** none. Task 4 originally carried `<SPEARMAN_FLOOR>` /
`<TOP10_FLOOR>` for Task 3 to fill; Task 3's measurements showed no honest
value exists at a feasible `n_boot` (see the scope change in Task 4), so the
statistical assertions were removed rather than filled with a guess. The only
surviving tolerance, `shape_rel = 1e-3`, is calibrated from observed output
(measured 3.4e-04).

**Type consistency:** `.hzr_bhdead_candidates()`, `.hzr_load_bhdead_fixture()`,
`.hzr_bhblt_path()`, `.hzr_bhdead_fixture_path()`, `.hzr_validate_bhdead_fixture()`
are defined in Task 1 and used with the same signatures in Tasks 2-5.
`.hzr_parse_bhdead_phase()` / `.hzr_parse_bhdead_shape()` /
`.hzr_build_bhdead_fixture()` / `.hzr_write_bhdead_fixture()` are defined in
Task 2 and used in Task 3. Fixture fields (`shape$specified`,
`shape$converged`, `early`, `constant`, `meta$sle`, `meta$sls`) are consistent
across all tasks.
