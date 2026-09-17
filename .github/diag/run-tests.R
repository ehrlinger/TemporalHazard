# Mirrors tests/testthat.R, but keeps the per-test results.
library(testthat)
library(TemporalHazard)
cat("MODE", Sys.getenv("DIAG_MODE"), "TESTTHAT_PARALLEL=", Sys.getenv("TESTTHAT_PARALLEL"),
    "Ncpus=", format(getOption("Ncpus")), "\n")
cat("SEARCH", search(), "\n")
w <- callr::r(function() search())
cat("CALLR_SEARCH", w, "\n")
cat("ENV R_DEFAULT_PACKAGES=", Sys.getenv("R_DEFAULT_PACKAGES", "<unset>"),
    " R_TESTS=", Sys.getenv("R_TESTS"), "\n")
res <- test_dir("testthat", package = "TemporalHazard", load_package = "installed",
                reporter = "check", stop_on_failure = FALSE)
df <- as.data.frame(res)
df <- df[, c("file", "context", "test", "nb", "passed", "failed", "skipped", "error", "warning", "real")]
write.csv(df, file.path(Sys.getenv("DIAG_OUT"), paste0(Sys.getenv("DIAG_MODE"), ".csv")), row.names = FALSE)
cat(sprintf("TOTALS %s tests=%d passed=%d failed=%d skipped=%d warning=%d error=%d\n",
            Sys.getenv("DIAG_MODE"), nrow(df), sum(df$passed), sum(df$failed),
            sum(df$skipped), sum(df$warning), sum(df$error)))
