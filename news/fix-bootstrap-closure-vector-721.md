* **`hzr_bootstrap()` now refuses a term that reads a per-row vector from
  outside `data` through a helper function or an object member (#721).**
  The refusal from #278 read the names in each formula and skipped any name
  that held a function or did not hold one value per row. So
  `z <- function(x) x + zvec` in `~ z(age)`, `cap <- function(x)
  pmin(x, caps)` in `Surv(cap(time), status)`, and `cfg$z` with
  `cfg <- list(z = ...)` all got through. Every replicate reported success
  while the vector stayed in its original order, and the interval was wrong.
  Each term of the response, the right-hand side, the phase formulas and a
  select-mode `scope` is now also evaluated on the rows of `data` in a
  second order, in each environment the refit reads it from, and a term
  that does not give each row its own value is refused by name. A term that
  depends on the order of the rows or is random, such as `cumsum(x)` or
  `jitter(age)`, is refused too, since each replicate would recompute it. A helper that only transforms its argument, such as
  `function(x) log(x)` defined next to the data, or a term computed from the
  rows as a set, such as `scale()` or `poly()`, still bootstraps. A
  bootstrap this now refuses would have returned a wrong interval.
