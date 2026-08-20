# sas-job.R -- the hzr_sas_job container that the SAS parser fills and the
# Quarto renderer consumes.
#
# The package's signature defect is an output that looks like a result and
# is not. A job where nothing was parsed has the right shape and is empty
# inside -- it defeats both is.null() and a length check. .hzr_validate_sas_job()
# exists to reject that hollow case rather than let a renderer emit a
# document over nothing.

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
