# bhdead-fixture.R -- Resolve, load and validate the bh.dead SAS parity fixture.
#
# DEVELOPMENT ONLY. inst/dev/ is .Rbuildignore'd and does not ship.
#
# This repository is public. The fixture itself is NOT stored here -- it lives
# beside the source data on a secure volume. Both paths below come entirely
# from environment variables; there is no default. When a variable is unset,
# the corresponding helper returns "", and .hzr_load_bhdead_fixture() treats
# that as "not configured" and returns NULL, so callers skip gracefully. This
# file contains no study values; every expectation used by the parity test is
# read out of the fixture at run time.
#
# Build the fixture with inst/dev/bhdead-parity/parse-bhdead-lst.R.
# See inst/dev/bhdead-parity/README.md.

#' Path to the row-level SAS dataset (env-only, no default)
#' @noRd
.hzr_bhblt_path <- function() {
  Sys.getenv("TEMPORALHAZARD_BHBLT", unset = "")
}

#' Path to the built parity fixture (env-only, no default)
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
