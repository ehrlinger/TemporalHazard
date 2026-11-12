* **A multiphase fit no longer returns a shape that has run to a limit
  (#683).** When a cdf or hazard phase's `nu` runs to 0 (the phase becomes
  a step), its `t_half` falls below the first observed time, or its
  `t_half` rises past the data, `hazard()` used to return the drifted last
  iterate, with a warning for the first two and nothing for the third. It
  now refits from a fixed grid of starts. A start that reaches the best
  objective without running to a limit the check detects is returned, with
  a warning of class `hzr_limit_restarted`. When none does, `hazard()`
  stops with an error of class `hzr_limit_refusal` that names the phase and
  parameter, so you can fix it, as `FIXNU` or `FIXTHALF` would in SAS. The
  error carries the restarts it tried as `$grid`, and the drifted fit as
  `$fit` for inspection; `print()` marks that fit as refused. A
  weakly identified shape, or a fit whose Hessian cannot be used so it has
  no standard errors (#734), tries the same grid but is never refused: the
  original is kept unless a start does strictly better. `print()` and
  `summary()` always show a `Limit check:` line, and `fit$fit$limit`
  records the trigger and the grid; a fit saved before this release reads
  "not run". `control = list(limit_restart = FALSE)` skips the grid: a hard
  limit then refuses at once, so it does not bring back the drifted fit,
  and a weak flag is recorded with the restart "not tried".
  `hzr_stepwise()` and `hzr_bootstrap()` record a refused refit under the
  reason `"limit: ..."`, and the bootstrap counts restarted replicates. The
  model is never switched for you.
* **A blind start now reaches the SAS optimum on `hm.death.patient`
  (#729).** It used to stop in another basin, where the early phase's
  `t_half` was not identified; the limit check restarts it from there.
* **Code inside a formula, a phase formula included, runs as often as it
  would without the limit check**, however many fits the check makes, and
  `hazard()`'s arguments are evaluated in the same order as without it, so
  a random term or a `runif()` in a phase draws the same numbers.
