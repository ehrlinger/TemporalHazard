- **`hzr_stepwise()` refits under the base fit's `control` (#688).** Every
  candidate refit, and the selected model, used `hazard()`'s defaults, so a
  base fit run with `conserve = FALSE`, a raised `maxit` or several
  `n_starts` was stepped, and its selected model fitted, under settings it
  never had, with no warning. They now run under the base fit's `control`,
  as every step of a SAS `SELECTION` runs under one `PROC HAZARD`'s options.
  A `control` passed to `hzr_stepwise()` overrides it entry by entry. This
  changes the result of a screen whose base fit set a non-default `control`.
- **`hzr_bootstrap()` refits under the `control` the fit ran under.** Each
  replicate rebuilt its fit from the stored call, where `control` is an
  expression, and a fit made at top level kept no copy of the variable it
  named. A `control = ctrl` fit whose `ctrl` was changed afterwards
  bootstrapped under the new value, with no warning, and one whose `ctrl`
  was removed failed every replicate. Replicates now use the fit's stored
  `control`.
