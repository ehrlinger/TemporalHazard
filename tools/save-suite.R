#!/usr/bin/env Rscript
# save-suite.R -- run the full test suite and save it with its provenance.
#
# The same run tools/release-suite.R makes, without its study-data gate, so a
# run made without the study volume can still feed tools/feature-matrix.R.
# The record of what produced it is attached as attr(results, "provenance").
#
# USAGE (from the package root, on a committed tree)
#   Rscript tools/save-suite.R results.rds

if (!file.exists("DESCRIPTION") ||
      read.dcf("DESCRIPTION", "Package")[1, 1] != "TemporalHazard") {
  stop("run this from the TemporalHazard package root.", call. = FALSE)
}
out <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out)) stop("give the .rds to write as the first argument.",
                     call. = FALSE)
source(file.path("tools", "suite-provenance.R"))

started <- Sys.time()
Sys.setenv(NOT_CRAN = "true")
res <- as.data.frame(devtools::test(reporter = "summary",
                                    stop_on_failure = FALSE))
attr(res, "provenance") <- suite_provenance(started)
saveRDS(res, out)
message(sprintf("%d tests, %d passed, %d failed, %d errored, %d skipped; ",
                nrow(res), sum(res$passed), sum(res$failed), sum(res$error),
                sum(res$skipped)),
        "saved to ", out)
