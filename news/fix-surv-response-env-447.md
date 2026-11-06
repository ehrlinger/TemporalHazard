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
  One case is still not caught: a helper function, defined in the calling
  environment, that itself uses a per-row vector, such as
  `cap <- function(x) pmin(x, caps)`. The bootstrap check looks at the
  names in the formula, not inside the functions they call, so `caps` is
  bootstrapped without being resampled. The right-hand side has had the same
  gap since before this change.
