# Predictions from a fit loaded out of a SAS `OUTHAZ=` dataset

Rebuilds the multiphase model the `OUTHAZ=` dataset describes – which
phases are in it, their shapes, and the fitted parameter vector – and
then predicts exactly as
[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
does.

## Usage

``` r
# S3 method for class 'hzr_outhaz'
predict(
  object,
  newdata,
  type = c("hazard", "survival", "cumulative_hazard"),
  decompose = FALSE,
  se.fit = FALSE,
  level = 0.95,
  conf.type = c("log-log", "logit"),
  ...
)
```

## Arguments

- object:

  An `hzr_outhaz` object from
  [`hzr_read_outhaz()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_read_outhaz.md).

- newdata:

  Data frame with a `time` column. Required.

- type:

  One of `"hazard"`, `"survival"` or `"cumulative_hazard"`; default
  `"hazard"`, as in
  [`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md).

- decompose:

  Accepted only to keep this method's argument list identical to
  [`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)'s,
  so that a positional call means the same thing for both. `TRUE` is an
  error: an `OUTHAZ=` dataset names its phases but carries no per-phase
  decomposition, and silently returning the total prediction instead
  would answer a question that was not asked.

- se.fit:

  Logical; add delta-method standard errors and confidence limits, as
  [`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
  does.

- level:

  Numeric confidence level in `(0, 1)`; default `0.95`. Only used when
  `se.fit = TRUE`.

- conf.type:

  Transform for `type = "survival"` confidence limits when
  `se.fit = TRUE`: `"log-log"` (default) or `"logit"`, which reproduces
  SAS `PROC HAZPRED`'s survival limits. A real argument rather than part
  of `...` because a misspelling would otherwise be swallowed silently
  and return log-log limits that disagree with the SAS job being
  reproduced. As in
  [`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md),
  the value is checked only where it is used, so an ignored one does not
  make a point or hazard prediction fail.

- ...:

  Must be empty. Anything landing here is an error rather than a
  silently ignored argument.

## Value

What
[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
returns for a fit with no covariates: a numeric vector of predictions
when `se.fit = FALSE`, and when `se.fit = TRUE` a data frame with
columns `fit`, `se.fit`, `lower` and `upper` (the confidence limits at
`level`). The decomposed long-format return of
[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
is not reachable here, because `decompose = TRUE` is refused.

## Details

`newdata` is required. An `OUTHAZ=` dataset holds a converged model and
no data, so there is no fitted time vector to fall back on; without
`newdata` the prediction would be over nothing.

## What is reconstructed, and what is refused

A phase is in the model when its intercept row (`E0`, `C0`, `L0`)
carries `_STATUS_ = 1`, which is the test `PROC HAZPRED` itself applies.
The early phase becomes `hzr_phase("cdf")` and the late phase
`hzr_phase("g3")`, with the shape estimates read off the
`DELTA`/`THALF`/`NU`/`M` and `TAU`/`GAMMA`/`ALPHA`/`ETA` rows.

These cases error rather than return a number that looks like a
prediction: a dataset carrying covariates (their coefficients cannot be
matched to `newdata` columns from the file alone); a non-zero `DELTA`
(the early-phase time transformation is not implemented); and a `G1FLAG`
that disagrees with the signs of the `M` and `NU` estimates.

With `se.fit = TRUE` there are three more, because the stored covariance
is on SAS's estimation scale and has to be mapped onto this package's:

- a fit constrained by `FIXMNU1`, which ties `M` to `1/NU`;

- a fit estimating `GAMMA`, `ALPHA` or `ETA` on one of PROC HAZARD's
  *composite* late-phase scales – `log(GAMMA*ETA - 2)` or
  `log(GAMMA*ETA/ALPHA - 2)` rather than
  [`log()`](https://rdrr.io/r/base/Log.html) of the parameter. This is
  the ordinary unconstrained late phase, not an exotic case: with
  `G3FLAG = 1` and no `FIXGE2`/`FIXGAE2`, `ETA` and `ALPHA` are always
  on a composite scale and `GAMMA` is whenever `ETA` is fixed. So an
  early + constant + late fit with any free late shape gets point
  predictions but no standard errors;

- a fit under `FIXGE2` or `FIXGAE2` where a late parameter is *derived*
  from an estimated one (`ETA = 2/GAMMA`, `GAMMA = 2/ETA`,
  `ALPHA = GAMMA*ETA/2`). The derived parameter has no covariance row,
  so its contribution to the variance would be dropped without trace.

## See also

[`hzr_read_outhaz()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_read_outhaz.md),
[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md).

## Examples

``` r
f <- system.file("extdata", "outhaz-fixture.rds", package = "TemporalHazard")
if (nzchar(f)) {
  fit <- hzr_read_outhaz(f)
  predict(fit, newdata = data.frame(time = c(1, 6, 12)), type = "survival")
}
#> [1] 0.9064419 0.7547405 0.6075586
```
