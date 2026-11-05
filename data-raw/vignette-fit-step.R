# Regenerates inst/extdata/vignette-fit-step.rds, the stepwise fit the
# clinical-analysis-walkthrough vignette loads instead of refitting.
#
# The fit costs about 26 s on a fast Mac and several times that on CRAN's
# Windows machines, which pushed the overall check past CRAN's 10-minute
# limit. The vignette loads this file when it is present and refits when it
# is not; tests/testthat/test-vignette-fit-step.R refits and compares, so a
# stale file fails CI.
#
# Run from the package root, with the current package installed:
#   Rscript data-raw/vignette-fit-step.R
#
# It runs the vignette's own code, every chunk up to and including
# `stepwise-forward`, through the helper the test uses, so the stored fit
# cannot drift from the code the vignette shows.

library(TemporalHazard)
source("tests/testthat/helper-vignette-code.R")
env <- run_vignette_through(
  "vignettes/clinical-analysis-walkthrough.qmd", "stepwise-forward"
)
fit_step <- env$fit_step
stopifnot(inherits(fit_step, "hzr_stepwise"), nrow(fit_step$steps) > 0L)

attr(fit_step, "generated_with") <- as.character(
  utils::packageVersion("TemporalHazard")
)
saveRDS(fit_step, "inst/extdata/vignette-fit-step.rds", compress = "xz")
