* **A `.` inside a function call in a formula was dropped without a word.**
  `terms()` writes out `.` only where it is a term of its own, so
  `Surv(time, status) ~ log(abs(.))` reached the design as `~ 1` and the
  model was fitted with no covariates. The same held for a phase's
  `hzr_phase(formula = )`. Beside a bare `.`, as in `~ . + log(.)`, `terms()`
  substituted the sum of the columns, giving one covariate `log(x1 + x2)`.
  `hazard()` now stops on `.` anywhere but as a term of the formula
  operators, naming the term, as `lm()` stops on these formulas. `~ .`,
  `~ . - age`, `~ (.)^2` and `~ .:age` are unchanged.
