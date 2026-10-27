* **A Conservation of Events fit could report wrong standard errors when its
  full-information variance could not be computed (#652).** The variance from
  the search holds the conserved phase's `log_mu` fixed, so it treats that
  parameter as known, and it was taken at the starting value of `log_mu`
  rather than the solved one. When the full-information Hessian could not be
  computed or inverted, the fit fell back to that variance and reported a
  standard error for every other parameter, while the warning named only the
  conserved `log_mu`. On `avc` with two phases, one standard error was 0.17
  of its full-information value and another 1.3 times it. Such a fit now
  reports no standard errors, warns that none are reported, and records why
  in `degraded_causes`; `vcov()` returns `NA`. The estimates are unchanged.
  The conserved `log_mu` is also now counted as estimated, as it is when the
  recompute succeeds, so the fit's AIC no longer comes out 2 lower than the
  same fit's with a successful recompute.
