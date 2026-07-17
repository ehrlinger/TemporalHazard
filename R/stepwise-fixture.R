# stepwise-fixture.R -- Load and validate SAS stepwise parity fixtures.
#
# Fixtures live under `inst/fixtures/stepwise-<scenario>.rds` (and so
# reside under `<pkg root>/fixtures/` at install time).  They are
# deliberately optional -- tests that consume them skip when the file
# is absent, so the package installs and tests pass without any SAS
# output.
#
# The workflow for adding a new fixture is documented in
# `inst/extdata/stepwise-fixtures/README.md`.  .hzr_build_stepwise_fixture
# packages three SAS capture files into the expected list shape;
# .hzr_load_stepwise_fixture is the lookup called from the parity test.

#' Required top-level fields of a stepwise parity fixture
#'
#' Used by both the loader and the builder to validate the list shape.
#' @keywords internal
#' @noRd
.hzr_stepwise_fixture_schema <- function() {
  list(
    meta = c(
      "dataset", "dist", "criterion", "direction",
      "slentry", "slstay", "sas_version", "captured_on",
      "proc_hazard"
    ),
    scope = c("candidates", "force_in", "force_out"),
    steps = c("step_num", "action", "variable", "phase",
              "stat", "df", "p_value"),
    final = c("coef", "logLik", "iterations")
  )
}


#' Validate a stepwise parity fixture against the schema
#'
#' @param fix A list as produced by `.hzr_build_stepwise_fixture()` or
#'   deserialised from `inst/fixtures/stepwise-*.rds`.
#' @return The fixture, unchanged, if valid.  Throws with a list of
#'   missing fields otherwise.
#'
#' @keywords internal
#' @noRd
.hzr_validate_stepwise_fixture <- function(fix) {
  schema <- .hzr_stepwise_fixture_schema()
  problems <- character()

  for (group in names(schema)) {
    if (!group %in% names(fix)) {
      problems <- c(problems, paste0("missing top-level: ", group))
      next
    }
    missing_fields <- setdiff(schema[[group]], names(fix[[group]]))
    if (length(missing_fields) > 0L) {
      problems <- c(problems,
                    sprintf("missing %s: %s", group,
                            paste(missing_fields, collapse = ", ")))
    }
  }

  if (is.data.frame(fix$steps) && nrow(fix$steps) > 0L) {
    bad_actions <- setdiff(fix$steps$action, c("enter", "drop"))
    if (length(bad_actions) > 0L) {
      problems <- c(problems,
                    sprintf("unknown steps$action: %s",
                            paste(bad_actions, collapse = ", ")))
    }
  } else if (!is.null(fix$steps)) {
    problems <- c(problems, "steps must be a data.frame")
  }

  if (length(problems) > 0L) {
    stop("Invalid stepwise parity fixture:\n  - ",
         paste(problems, collapse = "\n  - "),
         call. = FALSE)
  }
  invisible(fix)
}


#' Load a stepwise parity fixture by name
#'
#' Returns NULL (not an error) when the fixture file is absent, so
#' callers can `skip_if_null()`.
#'
#' @param name Short name like `"cabgkul-forward-wald"`.
#' @return The fixture list or NULL.
#'
#' @keywords internal
#' @noRd
.hzr_load_stepwise_fixture <- function(name) {
  file_name <- paste0("stepwise-", name, ".rds")
  path <- system.file("fixtures", file_name, package = "TemporalHazard")
  if (!nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  fix <- readRDS(path)
  .hzr_validate_stepwise_fixture(fix)
}


#' Build a stepwise parity fixture from SAS capture files
#'
#' Reads the three files emitted by
#' `inst/extdata/stepwise-fixtures/cabgkul-forward-wald.sas` and
#' assembles them into the fixture list (see `schema.md`).  Writes the
#' result to `inst/fixtures/` if `out_path` is supplied.
#'
#' @param trace_csv,final_csv,meta_txt Paths to the SAS outputs.
#' @param out_path Optional output `.rds` path.  When set, the
#'   validated fixture is written there before being returned.
#'
#' @return The fixture list.
#'
#' @keywords internal
#' @noRd
.hzr_build_stepwise_fixture <- function(trace_csv, final_csv, meta_txt,
                                         out_path = NULL) {
  if (!file.exists(trace_csv)) {
    stop("trace_csv not found: ", trace_csv, call. = FALSE)
  }
  if (!file.exists(final_csv)) {
    stop("final_csv not found: ", final_csv, call. = FALSE)
  }
  if (!file.exists(meta_txt)) {
    stop("meta_txt not found: ", meta_txt, call. = FALSE)
  }

  trace <- utils::read.csv(trace_csv, stringsAsFactors = FALSE)
  # SAS may emit "." for missing phase -- translate to NA for the R side.
  if ("phase" %in% names(trace)) {
    trace$phase[trace$phase == "."] <- NA_character_
  }

  final_tbl <- utils::read.csv(final_csv, stringsAsFactors = FALSE)
  if ("phase" %in% names(final_tbl)) {
    final_tbl$phase[final_tbl$phase == "."] <- NA_character_
  }

  meta_lines <- readLines(meta_txt)
  meta_kv    <- .hzr_parse_keyvals(meta_lines)

  loglik <- as.numeric(meta_kv[["logLik"]] %||% NA_real_)
  iters  <- as.integer(meta_kv[["iterations"]] %||% NA_integer_)

  fix <- list(
    meta = list(
      dataset     = meta_kv[["dataset"]]     %||% NA_character_,
      dist        = meta_kv[["dist"]]        %||% NA_character_,
      shaping     = meta_kv[["shaping"]],
      criterion   = meta_kv[["criterion"]]   %||% "wald",
      direction   = meta_kv[["direction"]]   %||% NA_character_,
      slentry     = as.numeric(meta_kv[["slentry"]] %||% NA_real_),
      slstay      = as.numeric(meta_kv[["slstay"]]  %||% NA_real_),
      sas_version = meta_kv[["sas_version"]] %||% NA_character_,
      captured_on = meta_kv[["captured_on"]] %||% NA_character_,
      proc_hazard = meta_kv[["proc_hazard"]] %||% NA_character_,
      notes       = meta_kv[["notes"]]
    ),
    scope = list(
      candidates = unique(trace$variable),
      force_in   = character(),
      force_out  = character()
    ),
    steps = trace,
    final = list(
      coef       = final_tbl,
      logLik     = loglik,
      iterations = iters
    )
  )

  .hzr_validate_stepwise_fixture(fix)

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(fix, out_path)
  }
  fix
}


#' Parse a `key=value` text file into a named character list
#'
#' @keywords internal
#' @noRd
.hzr_parse_keyvals <- function(lines) {
  nonblank <- lines[nzchar(lines) & !startsWith(trimws(lines), "#")]
  kv <- strsplit(nonblank, "=", fixed = TRUE)
  out <- vapply(kv, function(p) trimws(paste(p[-1], collapse = "=")),
                character(1L))
  names(out) <- vapply(kv, function(p) trimws(p[[1L]]), character(1L))
  as.list(out)
}
