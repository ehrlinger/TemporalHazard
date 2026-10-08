* **A stepwise screen on a Weibull fit whose scale cannot be represented
  failed every candidate refit, and a selection-mode `hzr_bootstrap()`
  pooled such replicates as successes (#566).** The refit's warm start
  copied `mu`, which is 0, `Inf` or subnormal on such a fit, and `hazard()`
  refuses that as a starting value. `hzr_stepwise()` said so, but
  `hzr_bootstrap()` did not read it: on 400 rows with a covariate near 1000,
  20 replicates all counted as successes and `z`, a true effect, was
  selected in 75% of them, against 100% for the same model with the
  covariate centered. The warm start now takes `log(mu)` from the fit,
  moved into the range a number can hold, and carries the difference
  through the coefficients. The same 20 replicates now select what the
  centered fit selects, replicate for replicate. Separately,
  `hzr_bootstrap()` now counts the selection-mode replicates in which a
  candidate refit failed, for any reason, and warns once: such a candidate
  was never tested. Those replicates are still pooled. The warning about a
  `mu` that cannot be represented no longer says the other parameters are
  unaffected without qualification.

* **`hzr_bootstrap()` pooled replicates that stopped on `nlm()` code 4 or 5
  and failed the relative-gradient test, with no count and no warning.**
  `hazard()` warns about such a fit, but replicates run with their warnings
  suppressed, and the count for #531 leaves those codes out because
  `hazard()` warns for them. `hzr_bootstrap()` now reads the base refit and
  the final fit of each replicate and warns once with the number of
  replicates affected. Their estimates are still pooled. A stepwise refit
  that was not kept is not read, since the screen does not return it.
