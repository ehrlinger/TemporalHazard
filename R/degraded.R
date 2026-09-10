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

# Estimated parameters whose variance is missing from a covariance matrix
# that does exist. .hzr_safe_solve() masks a non-positive variance to NA and
# still returns a matrix, so "vcov is a matrix" does not mean every estimated
# parameter has a standard error. A parameter held fixed carries an NA row by
# design and is not counted. Returns character(0) when vcov is not a matrix:
# that case is recorded separately, with its own causes.
.hzr_params_missing_variance <- function(vcov, fixed_mask = NULL,
                                         param_names = NULL) {
  if (!is.matrix(vcov)) return(character(0))
  d <- diag(vcov)
  p <- length(d)
  estimated <- if (length(fixed_mask) == p) !as.logical(fixed_mask) else rep(TRUE, p)
  estimated[is.na(estimated)] <- TRUE
  lacking <- estimated & !(is.finite(d) & d > 0)
  nms <- if (length(param_names) == p) as.character(param_names) else paste0("par", seq_len(p))
  nms[lacking]
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
#' @param fixed_mask Logical, `TRUE` for a parameter held fixed; its missing
#'   variance is by design and not recorded.
#' @param param_names Names used in the partial-loss cause.
#' @return `list(degraded, degraded_causes)`, in canonical order.
#' @noRd
.hzr_degraded_record <- function(vcov, weak, control, dist,
                                 fitted = TRUE, imported = FALSE,
                                 not_fitted_cause = "not requested (fit = FALSE)",
                                 reasons = list(),
                                 fixed_mask = NULL, param_names = NULL) {
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
  } else {
    lacking <- .hzr_params_missing_variance(vcov, fixed_mask, param_names)
    if (length(lacking)) {
      causes["standard_errors"] <- paste0("no finite positive variance for: ",
                                          paste(lacking, collapse = ", "))
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
  agree("standard_errors",
        !is.matrix(object$fit$vcov) ||
          length(.hzr_params_missing_variance(object$fit$vcov,
                                              object$fit$fixed_mask)) > 0L)
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
