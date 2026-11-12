* **A covariate named like one of its own phase's parameters is refused.**
  A design column called `nu`, `m`, `log_t_half` or `log_mu` in a cdf or
  hazard phase (`log_tau`, `gamma`, `alpha`, `eta` or `log_mu` in a g3
  phase, `log_mu` in a constant one) gave two parameters the same name,
  such as `early.nu`. The estimates were right, but `summary()` failed and a
  lookup by name, `coef(fit)["early.nu"]` or a stepwise Wald test among
  them, read the shape instead of the covariate. `hazard()` now stops,
  naming the phase and the column, so it can be renamed. A factor level
  that expands to such a name is caught too.
