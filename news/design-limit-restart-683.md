* **A multiphase fit no longer returns a shape that has run to a limit
  (#683).** When an early or late phase's `nu` runs to 0 (the phase becomes
  a step), its `t_half` falls below the first observed time, or its
  `t_half` rises past the data, `hazard()` used to return the drifted last
  iterate with a warning and a `$boundary` record. It now refits from a
  fixed grid of starts. A start that stays clear of every limit and reaches
  the best objective is returned, with a warning of class
  `hzr_limit_restarted`. When none does, `hazard()` stops with an error of
  class `hzr_limit_refusal` that names the phase and parameter, so you can
  fix it, as `FIXNU` or `FIXTHALF` would in SAS. A weakly identified shape,
  or a fit whose Hessian cannot be used so it has no standard errors
  (#734), tries the same grid but is never refused: the original is kept
  unless a start does strictly better. `print()` and `summary()` always
  show a `Limit check:` line, and `fit$fit$limit` records the trigger and
  the grid. `control = list(limit_restart = FALSE)` skips the grid.
  `hzr_stepwise()` and `hzr_bootstrap()` record a refused refit under the
  reason `"limit: ..."`, and the bootstrap counts restarted replicates. The
  model is never switched for you.
* **A blind start now reaches the SAS optimum on `hm.death.patient`
  (#729).** It used to stop in another basin, where the early phase's
  `t_half` was not identified; the limit check restarts it from there.
* **Code inside a formula runs once per `hazard()` call**, including a phase
  formula, however many fits the limit check makes, so a random term draws
  once.
