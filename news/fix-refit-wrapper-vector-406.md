* **A vector fit made through a wrapper that passes its formula argument by
  variable is recognised as one (#406).** `wrap <- function(dat) { fml <-
  NULL; hazard(formula = fml, time = dat$t, status = dat$s, ...) }` stores
  `call$formula = quote(fml)`: a symbol, not `NULL`. Checks that asked
  `is.null(call$formula)` therefore took the fit for a formula fit and tried
  to resolve `fml`, so `hzr_stepwise()` stopped with "stored model formula
  (`fml`) did not resolve to a formula" and `hzr_bootstrap()` with "cannot
  count the rows to resample" -- both loud, both about a name rather than
  about the interface. A fit of that shape -- `time =` in the call, the name
  still bound to `NULL` where the fit was made, and no stored data frame
  beside any design, which a formula fit must have because `hazard()`
  requires `data` with a formula -- now behaves as the same fit made
  without the wrapper: its bootstrap replicates are identical to the plain
  fit's, a multiphase one screens instead of stopping, a single-distribution
  one is refused as the vector fit it is, and one holding a design passed as
  `x` is refused by the same rule as the plain fit, rather than resampling
  rows beside a design that cannot follow them. The stored formula argument
  is read only when it is a name, and then only for its value, so an
  expression written there is never evaluated and a bootstrap neither
  repeats its side effects nor consumes the random number stream. Fits of
  other shapes are unchanged, including a fit that stores both a design and
  a frame, which could be either interface, and the messages above: which
  interface `hazard()` used cannot be recovered from a stored call
  afterwards, and the remedy is for `hazard()` to record it, which is
  tracked separately (#432). The `predict()` half of #406 is a separate
  change (#409), which carries its own rule for the same question.
