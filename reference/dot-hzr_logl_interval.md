# Interval-censored log-likelihood contribution

The single place the interval-censored contribution is written. Both
[`.hzr_logl_multiphase()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_logl_multiphase.md)
and the finite-difference closure inside
[`.hzr_gradient_multiphase()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_gradient_multiphase.md)
delegate here, so the optimizer cannot step by the gradient of a
different objective than the one it evaluates.

## Usage

``` r
.hzr_logl_interval(
  cumhaz_lower,
  cumhaz_upper,
  lower,
  upper,
  weights,
  objective = c("likelihood", "sas")
)
```

## Arguments

- cumhaz_lower:

  Cumulative hazard at the interval lower bounds, \\\Lambda(l)\\.

- cumhaz_upper:

  Cumulative hazard at the interval upper bounds, \\\Lambda(u)\\.

- lower:

  Interval lower bounds. Used only by `objective = "sas"`.

- upper:

  Interval upper bounds. Used only by `objective = "sas"`.

- weights:

  Case weights. In a SAS parity run these are the ICENSOR variable,
  which is a weight (a death count in an aggregated study), not merely
  an indicator.

- objective:

  `"likelihood"` for the interval probability – the default, and the
  only statistically consistent form – or `"sas"` for the
  interval-mean-hazard density term `PROC HAZARD` accumulates. See
  `inst/dev/SAS-INTERVAL-OBJECTIVE-DESIGN.md`.

## Value

Scalar summed contribution; `-Inf` for infeasible parameters.

## Details

Callers pass **only the interval rows** – already subset by
`status == 2` – so this helper never sees the status mask and cannot
disagree with a caller about which rows are intervals.
