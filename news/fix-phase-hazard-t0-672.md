* **The hazard at exactly time 0 was `NaN` for a phase with `m = 0` and
  `nu > 0` (#672).** That phase's density was formed as a product of a
  survival that underflows to 0 and a power that overflows, so at time 0 it
  read `0 * Inf`. `hzr_phase_hazard()` and `hzr_decompos()$g` returned `NaN`
  where the limit is 0, and in a multiphase fit the `NaN` took the whole
  hazard with it: `predict(type = "hazard")` at time 0 returned `NaN`, with
  no message, even for the constant phase's share. Scoring start/stop data
  reaches this, because every first segment starts at 0. The density is now
  formed on the log scale, giving 0 at time 0 and the same values, to
  rounding, elsewhere. `hazard()` drops rows at time 0, so the likelihood
  never read the `NaN`.
* `predict(se.fit = TRUE)` on a Weibull fit at time 0 gave a standard error
  of `NaN` from the same `0 * Inf`, in the shape term; it is now 0 (#619 item
  22). The confidence limits there stay `NA`, as for any prediction with a
  cumulative hazard of 0.
* `?hazard` now says that `fit$fit$objective`, and the `objective` of each
  start, is the maximized log-likelihood, larger is better (or PROC HAZARD's
  objective under `objective = "sas"`). The `reltol` entry, which
  describes the negative the optimizer minimizes, says the stored value has
  the opposite sign (#674).
