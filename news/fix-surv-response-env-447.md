* **The `Surv()` response of a formula now sees the formula's environment
  (#447).** It was evaluated in `data` and then in a frame of the package's
  own, never the environment the formula was written in. A fit inside a
  function whose response used a local variable, such as
  `Surv(pmin(time, cutoff), status)`, silently used a global variable of the
  same name instead, if one existed, with no warning. Without one it failed
  with "object not found". The right-hand side and `stats::lm()` already
  looked there. A column of `data` still takes precedence. `hzr_bootstrap()`
  now checks response variables in the same place, so a per-row local in the
  response is refused instead of being left out of the resampling.
