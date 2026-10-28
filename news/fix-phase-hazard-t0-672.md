* **The hazard at exactly time 0 was `NaN` for a phase with `m = 0` and
  `nu > 0` (#672).** That phase's density was formed as a product of a
  survival that underflows to 0 and a power that overflows, so at time 0 it
  read `0 * Inf`. `hzr_phase_hazard()` and `hzr_decompos()$g` returned `NaN`
  where the limit is 0, and in a multiphase fit the `NaN` took the whole
  hazard with it: `predict(type = "hazard")` at time 0 returned `NaN`, with
  no message, even for the constant phase's share. Scoring start/stop data
  reaches this, because every first segment starts at 0. The same `0 * Inf`
  struck at times above 0 when `nu` is small or `t_half` is large (1e20 is
  enough near t = 0, where `h` was `NaN` as well), and there the likelihood
  did read it. Such a point was treated as
  one where the likelihood is not defined, so a fit could be steered away
  from it or fail outright: one with an early phase at `nu = 0.01` and three
  events before t = 0.002 stopped with "no usable fit", and now converges.
  The density is now formed on the log scale, giving the limit, 0, where the
  product was `NaN`, and the same values to rounding elsewhere. A fit that
  reached such a point can therefore change.
* `predict(se.fit = TRUE)` on a Weibull fit at time 0 gave a standard error
  of `NaN` from the same `0 * Inf`, in the shape term; it is now 0 (#619 item
  22). The confidence limits there stay `NA`, as for any prediction with a
  cumulative hazard of 0.
* `?hazard` now says that `fit$fit$objective`, and the `objective` of each
  start, is the maximized log-likelihood, larger is better (or PROC HAZARD's
  objective under `objective = "sas"`). The `reltol` entry, which
  describes the negative the optimizer minimizes, says the stored value has
  the opposite sign (#674).
