# parse-bhdead-lst.R -- Build the bh.dead parity fixture from SAS listings.
#
# DEVELOPMENT ONLY. inst/dev/ is .Rbuildignore'd and does not ship.
#
# This repository is public: the .lst inputs and the .rds output both stay on
# the secure volume. This script contains parsing logic only -- no study values.
#
# Usage (hz_lst/bh_lst are required -- no default paths are stored here):
#   source("inst/dev/bhdead-parity/parse-bhdead-lst.R")
#   .hzr_write_bhdead_fixture(hz_lst = "...", bh_lst = "...")
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
  # occurrence to the end of the listing. The Obs sequence below -- not a
  # following label -- is what bounds the table, so this works regardless of
  # what (if anything) follows.
  from <- starts[1]
  block <- lines[from:length(lines)]

  # Data rows: obs index, then a NAME of uppercase/digits/underscore, then 6 numbers.
  rx <- paste0("^\\s*([0-9]+)\\s+([A-Z0-9_]+)\\s+([0-9]+)\\s+",
               "(-?[0-9.]+)\\s+(-?[0-9.]+)\\s+(-?[0-9.]+)\\s+",
               "(-?[0-9.]+)\\s+(-?[0-9.]+)\\s*$")
  hits <- regmatches(block, regexec(rx, block))
  hits <- Filter(function(x) length(x) == 9L, hits)
  if (length(hits) == 0L) {
    stop("No data rows parsed for ", phase_label, call. = FALSE)
  }

  # The Obs column increments monotonically within one table and restarts at 1
  # in a different one. Truncate at the first break in that run so trailing
  # content that merely matches the data-row shape (a different table, or
  # nothing at all) is never silently absorbed.
  obs <- as.integer(vapply(hits, `[[`, character(1), 2L))
  if (obs[1] != 1L) {
    stop("First data row for ", phase_label, " has Obs = ", obs[1],
         " (expected 1); the table did not start where expected.",
         call. = FALSE)
  }
  breaks <- which(diff(obs) != 1L)
  last <- if (length(breaks) > 0L) breaks[1] else length(obs)
  hits <- hits[seq_len(last)]

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
#' `hz_lst` and `bh_lst` are required: no study-identifying path is stored in
#' this public repository. `out` defaults to the path named by the
#' `TEMPORALHAZARD_BHDEAD_FIXTURE` environment variable.
#'
#' @param hz_lst Path to the shape-fit listing (required).
#' @param bh_lst Path to the bootstrap listing (required).
#' @param out Output `.rds` path. Defaults to `.hzr_bhdead_fixture_path()`.
#' @param meta Named list: resampl, sle, sls, n_obs, captured_on, source.
#' @noRd
.hzr_write_bhdead_fixture <- function(
    hz_lst,
    bh_lst,
    out    = .hzr_bhdead_fixture_path(),
    meta   = list(resampl = 1000L, sle = 0.12, sls = 0.10, n_obs = NA_integer_,
                  captured_on = as.character(Sys.Date()),
                  source = "hz.dead.lst + bh.dead.lst")) {
  if (!nzchar(out)) {
    stop("Set TEMPORALHAZARD_BHDEAD_FIXTURE (or pass `out`) before writing ",
         "the fixture.", call. = FALSE)
  }
  fix <- .hzr_build_bhdead_fixture(hz_lst, bh_lst, meta)
  saveRDS(fix, out)
  message("Wrote fixture: ", out,
          " (early: ", nrow(fix$early), " rows, constant: ",
          nrow(fix$constant), " rows)")
  invisible(fix)
}
