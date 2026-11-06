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
  Two cases are still not caught (#721). A per-row vector reached through
  another object, such as `cfg$caps` with `cfg <- list(caps = ...)`, or
  through a helper function that uses it, such as
  `cap <- function(x) pmin(x, caps)`, is bootstrapped without being
  resampled. The check measures the name in the formula (`cfg`, `cap`), not
  the values it leads to. The right-hand side, and a global in the response,
  have had the same gap since before this change.
