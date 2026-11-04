- **`hzr_stepwise()` refits under the base fit's `control` (#688).** Every
  candidate refit, and the selected model, used `hazard()`'s defaults, so a
  base fit run with `conserve = FALSE`, a raised `maxit` or several
  `n_starts` was stepped, and its selected model fitted, under settings it
  never had, with no warning. They now run under the base fit's `control`,
  as every step of a SAS `SELECTION` runs under one `PROC HAZARD`'s options.
  A `control` passed to `hzr_stepwise()` overrides it entry by entry. This
  changes the result of a screen whose base fit set a non-default `control`.
