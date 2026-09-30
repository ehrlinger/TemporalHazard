# Print method for hzr_gof

Print method for hzr_gof

## Usage

``` r
# S3 method for class 'hzr_gof'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An `hzr_gof` object.

- digits:

  Number of decimal places for formatting.

- ...:

  Additional arguments (ignored).

## Value

The object `x` of class `c("hzr_gof", "data.frame")`, invisibly. The
data frame has one row per time point and columns: `time`, `n_risk`,
`n_event`, `n_censor`, `km_surv` (Kaplan-Meier survival), `km_cumhaz`,
`par_surv` (parametric survival), `par_cumhaz`, `cum_observed`
(cumulative observed events), `cum_expected` (cumulative expected events
from model), `residual` (cum_expected - cum_observed). Multiphase models
additionally include `par_cumhaz_<phase>` columns for per-phase
cumulative hazard contributions. A `"summary"` attribute contains scalar
diagnostics: `total_observed`, `total_expected`, `final_residual`,
`dist`, `n`, `n_tallied`. When `n_tallied` is below `n`, the printed
totals cover only those subjects, and a note says so.
