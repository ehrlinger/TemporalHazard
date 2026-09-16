# Rebuild a multiphase fit as hazard() returned it before #299, when a phase
# formula given with no `data` was ignored: the phase took the global `x`, or
# no columns. hazard() refuses that call now, so such a fit exists only as a
# saved object. Fit the same call without the phase formulas, then write
# them back where hazard() kept them. Every element but `call` and
# `call_env` then equals the old object; test-phase-formula-needs-data.R
# checks that against hazard() with the refusal mocked away.
hzr_saved_before_299 <- function(phases, ...) {
  bare <- lapply(phases, function(ph) {
    ph["formula"] <- list(NULL)
    ph
  })
  f <- hazard(phases = bare, ...)
  for (nm in names(phases)) {
    f$spec$phases[[nm]]["formula"] <- list(phases[[nm]]$formula)
    if (!is.null(f$fit$phases)) {
      f$fit$phases[[nm]]["formula"] <- list(phases[[nm]]$formula)
    }
  }
  f
}
