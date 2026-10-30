#!/usr/bin/env Rscript
# release-suite.R -- the full local test suite, run as a release gate.
#
# WHY
#   R CMD check runs with NOT_CRAN false, and no CI runner has the study
#   volume, so the two parity files that read the private study datasets run
#   nowhere but a local suite with the data in reach. Twice a gate passed with
#   them skipped: v1.2.13 with the volume unmounted, which shipped two red
#   cardioversion tests, and the #600 gate on 2026-10-12 with a stale mount
#   (/Volumes/qhsstudies present but unreadable, the share live at
#   /Volumes/qhsstudies-1). Both runs printed a skip count that a person had
#   to read and judge.
#
# WHAT IT DOES
#   Refuses before running if either data directory does not resolve. It asks
#   the tests' own resolvers (inst/sas-parity/helper-sas-parity.R), so the
#   gate and the tests cannot disagree about where the data are.
#   Refuses after running unless both files ran, skipped nothing, and the
#   whole suite has no failure and no error.
#   Only the final "RELEASE SUITE PASSED" line means it passed.
#
# USAGE (from the package root)
#   Rscript tools/release-suite.R [results.rds]
#   The resolvers have no default path, so point them at the data first,
#   once in ~/.Renviron or on the command line:
#   HAZARD_MAZE_DATASETS=<dir> HAZARD_ACHALASIA_DIR=<dir> \
#     Rscript tools/release-suite.R

volume_files <- c(
  "test-repeated-events-parity.R",
  "test-repeated-events-renewal-parity.R"
)

refuse <- function(...) {
  message("RELEASE SUITE REFUSED: ", ...)
  quit(save = "no", status = 1)
}

if (!file.exists("DESCRIPTION") ||
      read.dcf("DESCRIPTION", "Package")[1, 1] != "TemporalHazard") {
  refuse("run this from the TemporalHazard package root.")
}

out <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out)) out <- file.path(tempdir(), "release-suite.rds")

resolvers <- new.env()
sys.source(file.path("inst", "sas-parity", "helper-sas-parity.R"),
           envir = resolvers)
dirs <- c(HAZARD_MAZE_DATASETS = resolvers$.hzr_maze_datasets_dir(),
          HAZARD_ACHALASIA_DIR = resolvers$.hzr_achalasia_dir())
if (anyNA(dirs)) {
  refuse(
    "the study data do not resolve for ",
    paste(names(dirs)[is.na(dirs)], collapse = " and "),
    ". Set the variable to the data directory (once, in ~/.Renviron) ",
    "and run again. A stale mount reads as missing, so check for ",
    "/Volumes/qhsstudies-1."
  )
}
message("study data: ", paste(names(dirs), dirs, sep = " = ", collapse = "; "))

# What produced the run travels with it, for tools/feature-matrix.R (#682).
source(file.path("tools", "suite-provenance.R"))
Sys.setenv(NOT_CRAN = "true")
res <- tryCatch(
  suite_record(function() {
    as.data.frame(devtools::test(reporter = "summary", stop_on_failure = FALSE))
  }),
  suite_changed = function(e) refuse(conditionMessage(e))
)
saveRDS(res, out)

problems <- character()
for (f in volume_files) {
  r <- res[res$file == f, ]
  if (nrow(r) == 0L) {
    problems <- c(problems, paste(f, "did not run"))
    next
  }
  message(sprintf("%s: %d tests, %d passed, %d failed, %d skipped",
                  f, nrow(r), sum(r$passed), sum(r$failed), sum(r$skipped)))
  if (any(r$skipped)) {
    problems <- c(problems, sprintf("%s skipped %d of its %d tests",
                                    f, sum(r$skipped), nrow(r)))
  }
}
if (sum(res$failed) > 0L) {
  problems <- c(problems, sprintf("%d failed expectations", sum(res$failed)))
}
if (any(res$error)) {
  problems <- c(problems, sprintf("%d tests errored", sum(res$error)))
}

other <- res[res$skipped & !res$file %in% volume_files, ]
message(sprintf(
  "suite: %d files, %d tests, %d passed, %d failed, %d errored, %d skipped",
  length(unique(res$file)), nrow(res), sum(res$passed), sum(res$failed),
  sum(res$error), sum(res$skipped)
))
if (nrow(other)) {
  message("other skips, by file (not gated):")
  counts <- table(other$file)
  message(paste0("  ", names(counts), ": ", counts, collapse = "\n"))
}
message("results saved to ", out)

if (length(problems)) refuse(paste(problems, collapse = "; "), ".")
message("RELEASE SUITE PASSED")
