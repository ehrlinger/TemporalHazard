* **A covariate named like one of its own phase's parameters is refused.**
  A design column called `nu`, `m`, `log_t_half` or `log_mu` in a cdf or
  hazard phase (`log_tau`, `gamma`, `alpha`, `eta` or `log_mu` in a g3
  phase, `log_mu` in a constant one) gave two parameters the same name,
  such as `early.nu`. The estimates were right, but `summary()` failed and a
  lookup by name, `coef(fit)["early.nu"]` or a stepwise Wald test among
  them, read the shape instead of the covariate. `hazard()` now stops,
  naming the phase and the column, so it can be renamed. A factor level
  that expands to such a name is caught too, and so is the rarer case of
  two phases making one name (a phase `a.b` beside a phase `a` with a
  covariate `b.log_mu`). A fit saved before this release with such a name
  still prints and predicts, but `hzr_evaluate()` and `hzr_bootstrap()`
  now refuse it. `hzr_stepwise()` refuses such a candidate, or a
  `data` column that the default scope would offer, before fitting any
  candidate, under every criterion.
