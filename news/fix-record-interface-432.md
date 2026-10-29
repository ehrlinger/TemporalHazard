* **A vector-interface fit made through a wrapper that forwards
  `formula = fml` could be refit as a formula fit (#432, #406).** `hazard()`
  takes the formula route only when `formula` evaluates to a formula, but the
  stored call records `fml` whether it held one or `NULL`, so the refitting
  tools guessed the route from the call. Where `fml` was later bound to a
  formula, `hzr_bootstrap()` resampled a fit given `x =` that the plain
  vector fit refuses, and reported 3 of 3 replicates with a row labelled
  `age` carrying `mal`'s coefficients, with no warning. Where `fml` stayed
  `NULL`, `hzr_stepwise()` stopped with a message about resolving `fml`
  instead of the check a vector fit gets, and `hzr_bootstrap()` refused a fit
  the plain vector fit resamples; it now resamples it. `hazard()` now records the route
  it took in `spec$interface`, and `hzr_bootstrap()`, `hzr_stepwise()` and
  `predict(newdata = )` read it. A fit saved before this version is read as
  before.
