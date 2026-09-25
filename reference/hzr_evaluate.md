# Evaluate a hazard model at parameters you supply

Computes a model's log-likelihood, and optionally its hazard and
cumulative hazard, at parameters **you** supply rather than at
parameters fitted from the data. This is what a parity check needs: the
likelihood at another program's converged estimates, evaluated by this
package's own likelihood.

## Usage

``` r
hzr_evaluate(object, theta, times = NULL)
```

## Arguments

- object:

  A `hazard` object, fitted or built with `fit = FALSE`. Its data,
  distribution and phase specification are used; its own `theta` is not.

- theta:

  Numeric vector of parameters to evaluate at, on the internal scale the
  model uses (see
  [`hzr_theta_names()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_theta_names.md)).
  Its length must match the model's parameter count. Names, when
  supplied, must match too.

- times:

  Optional numeric vector of times, for multiphase models only. When
  given, the hazard and cumulative hazard at those times are returned
  for a covariate-free ("baseline") subject: every covariate at 0, so
  the curve is the shape the supplied parameters describe, not a
  prediction for any row of the data. The other families have no
  internal shape function that takes parameters directly, and writing
  their hazards out here would duplicate
  [`predict()`](https://rdrr.io/r/stats/predict.html); `times` is
  refused for them rather than mirrored.

## Value

An object of class `hzr_evaluation`: a list with `theta`, the parameters
the likelihood was evaluated at – the vector supplied, with any
constrained entry replaced by the value its phase derives, which is
warned about as
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
warns; `logLik`, the log-likelihood there; `dist`; `n_obs` and
`n_events`, the rows the likelihood scored and the exact events among
them (a left- or interval-censored row counts in `n_obs`, not in
`n_events`); and, when `times` was given, `curve`, a data frame of
`time`, `hazard` and `cumulative_hazard`.

## Details

The result is not a fit and does not pretend to be one. It carries no
standard errors, no convergence status and no covariance matrix, because
nothing was estimated: [`print()`](https://rdrr.io/r/base/print.html)
says so on its first line, and the parameters are labelled as supplied.
Where a fitted model is what you want, use `hazard(..., fit = TRUE)`.

A phase built with `hzr_phase(constraint = )` derives one of its shapes
from the others, and that rule is applied to `theta` here as the fit
applies it, so the derived entry you pass is replaced rather than used
as given. At a fitted model's own estimates this returns that fit's
objective, under Conservation of Events too: since \#362 the fit's
objective is recomputed at the estimates it returns, except where the
fit warns that it could not, and then the two can differ.

A `theta` that passes the input checks but that the likelihood cannot
evaluate – an overflowing rate, a shape outside the family, or a value
past a guard that stops short of where the log-likelihood itself
overflows – gives `-Inf`, with a warning of class
`"hzr_evaluate_not_finite"`, for every distribution. A Weibull `mu` or
`nu` at or below 0 is still refused outright by those input checks.

## See also

[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
to fit a model,
[`hzr_theta_names()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_theta_names.md)
for the parameter order.

## Examples

``` r
data(avc, package = "TemporalHazard")
spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
               dist = "weibull", theta = c(1, 1), fit = FALSE)
# The likelihood at parameters from somewhere else, e.g. a SAS run.
hzr_evaluate(spec, theta = c(0.05, 0.9))
#> Model evaluated at SUPPLIED parameters -- not a fit
#>   nothing was estimated here: no standard errors, no convergence,
#>   no covariance. Use hazard(fit = TRUE) to fit.
#> 
#>   distribution: weibull
#>   observations: 310 (70 events)
#>   logLik at the supplied parameters: -748.63111
#> 
#>   parameters supplied:
#> [1] 0.05 0.90
```
