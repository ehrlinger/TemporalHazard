# Wayne Nelson cumulative hazard estimator with lognormal confidence limits

Compute the Nelson-Aalen cumulative hazard estimate with lognormal
confidence limits. Supports weighted events for severity-adjusted
analyses of repeated/recurrent events. This is the R equivalent of the
SAS `nelsonl.sas` macro.

## Usage

``` r
hzr_nelson(time, event, weight = NULL, conf_level = 0.95)

# S3 method for class 'hzr_nelson'
print(x, digits = 4, ...)
```

## Arguments

- time:

  Numeric vector of follow-up times.

- event:

  Numeric event indicator (1 = event, 0 = censored).

- weight:

  Optional numeric vector of event weights (default 1). Weights are
  applied only to events (censored observations contribute zero weight).
  Use for severity-weighted repeated events.

- conf_level:

  Confidence level for the interval (default 0.95). SAS `%NELSONT`
  defaults to `CLEVEL = 0.68268948`, a 1-SD (68.3%) interval, so a band
  reproduced at this function's default is about 1.96 times wider than
  the SAS one. Pass `conf_level = 2 * stats::pnorm(1) - 1` to match.
  Same convention as
  [`hzr_kaplan()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_kaplan.md).

- x:

  An `hzr_nelson` object.

- digits:

  Number of decimal places for formatting.

- ...:

  Additional arguments (ignored).

## Value

A data frame with one row per unique event time and columns:

- time:

  Event time.

- n_risk:

  Number at risk.

- n_event:

  Number of events at this time.

- weight_sum:

  Sum of event weights at this time.

- cumhaz:

  Nelson cumulative hazard estimate.

- std_err:

  Standard error.

- cl_lower:

  Lower lognormal confidence limit.

- cl_upper:

  Upper lognormal confidence limit.

- hazard:

  Interval hazard rate.

- cum_events:

  Cumulative (weighted) event count.

## Details

Unlike
[`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)
which uses the Breslow estimator with Greenwood variance, this function
uses the Wayne Nelson estimator with lognormal confidence limits that
are always non-negative.

## References

Nelson W (1972). Theory and applications of hazard plotting for censored
failure data. *Technometrics* 14(4):945–966.
[doi:10.1080/00401706.1972.10488991](https://doi.org/10.1080/00401706.1972.10488991)

Aalen O (1978). Nonparametric inference for a family of counting
processes. *Ann Statist* 6(4):701–726.
[doi:10.1214/aos/1176344247](https://doi.org/10.1214/aos/1176344247)

## See also

[`hzr_kaplan()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_kaplan.md)
for survival estimation.

## Examples

``` r
data(cabgkul)
nel <- hzr_nelson(cabgkul$int_dead, cabgkul$dead)
head(nel)
#> Nelson cumulative hazard estimate with lognormal CL
#> Events: 69  | Time points: 6 
#> 
#>     time n_risk n_event weight_sum cumhaz std_err cl_lower cl_upper hazard
#>  0.03285   5880      39         39 0.0066   2e-04   0.0063   0.0070 0.2019
#>  0.06571   5841       9          9 0.0082   2e-04   0.0077   0.0087 0.0469
#>  0.09856   5832       3          3 0.0087   3e-04   0.0081   0.0093 0.0157
#>  0.13142   5829       7          7 0.0099   3e-04   0.0092   0.0106 0.0365
#>  0.16427   5822       9          9 0.0114   4e-04   0.0107   0.0122 0.0471
#>  0.19713   5813       2          2 0.0118   4e-04   0.0110   0.0126 0.0105
#>  cum_events
#>          39
#>          48
#>          51
#>          58
#>          67
#>          69
```
