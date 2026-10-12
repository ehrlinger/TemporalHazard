* **A multiphase fit no longer selects a start that ran away to an
  underflow (#600).** For an early phase with `m = 0` and `nu < 0`,
  `hzr_decompos()` computed the CDF as `1 - exp(-x)`, which rounds to
  exactly 0 once `x` is below about `1e-16`, while the density, formed
  directly, stayed positive. Where that happened the cumulative hazard
  vanished and the hazard did not, so the likelihood had a false peak that
  rose without bound in `mu`. On the `avc` data, the `m`-fixed model started
  from `t_half = 0.5`, `nu = 1`, `mu = 1` with the default five starts
  climbed it to `log(mu)` 130 and a log-likelihood of +6380.86, and returned
  that start as the fit with `converged = TRUE`; three of the other four
  reach -214.875. The CDF is now `-expm1(-x)`, as the `m > 0, nu < 0`
  branch already had it, and the same call returns -214.875.
