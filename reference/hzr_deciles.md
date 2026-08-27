# Decile-of-risk calibration

Partition observations into groups (default 10) by predicted risk and
compare observed vs. expected event counts in each group. Good
calibration means the two track each other across the risk spectrum.

## Usage

``` r
hzr_deciles(object, time, groups = 10L, status = NULL, event_time = NULL)
```

## Arguments

- object:

  A fitted `hazard` object (with `fit = TRUE`).

- time:

  Numeric scalar: the horizon at which predicted survival is used to
  **rank subjects into risk groups** (e.g. `time = 12` ranks by 12-month
  predicted survival). It does not restrict the event/expected counts,
  which are accumulated over each subject's full follow-up.

- groups:

  Integer: number of risk groups (default 10 for deciles).

- status:

  Optional numeric vector of event indicators (1 = event, 0 = censored).
  If `NULL` (default), extracted from the fitted object's stored data.

- event_time:

  Optional numeric vector of observed event/censoring times. If `NULL`
  (default), extracted from the fitted object.

## Value

A data frame with one row per risk group and columns:

- group:

  Integer group label (1 = lowest risk, ranked by predicted survival at
  `time`).

- n:

  Number of observations in the group.

- events:

  Observed event count in the group (all events over follow-up).

- expected:

  Expected event count: the sum of each subject's predicted cumulative
  hazard at its own follow-up time.

- observed_rate:

  Observed event rate (events / n).

- expected_rate:

  Expected event rate (expected / n).

- chi_sq:

  Chi-square contribution: (events - expected)^2 / expected.

- p_value:

  Upper-tail p-value from the chi-square test for this group (1 df).

- mean_survival:

  Mean predicted survival probability at the horizon in the group.

- mean_cumhaz:

  Mean predicted cumulative hazard at follow-up in the group.

An attribute `"overall"` is attached with the overall chi-square
statistic, degrees of freedom, and p-value.

## Details

This reproduces the SAS `deciles.hazard.sas` macro. **All** subjects are
ranked by predicted survival at the horizon `time` and split into
equal-sized risk groups. Within each group the **expected** event count
is the sum of each subject's predicted cumulative hazard at its *own*
follow-up time, and the **observed** count is its number of events;
under conservation of events the group totals sum to the total observed
events. The horizon therefore only stratifies subjects into risk groups
– it does not restrict or exclude any subject, and the expected/observed
totals are independent of it.

## See also

[`predict.hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/predict.hazard.md)
for the prediction types used internally.

## Examples

``` r
# \donttest{
data(avc)
avc <- na.omit(avc)
fit <- hazard(
  survival::Surv(int_dead, dead) ~ age + mal,
  data  = avc,
  dist  = "weibull",
  theta = c(mu = 0.01, nu = 0.5, beta_age = 0, beta_mal = 0),
  fit   = TRUE
)
cal <- hzr_deciles(fit, time = 120)
print(cal)
#> Decile-of-risk calibration (risk grouped at time = 120 )
#> 305 subjects, all included.
#> 10 groups, 68 observed events, 68.1 expected
#> 
#>  group  n events expected observed_rate expected_rate  chi_sq p_value
#>      1 31      2    0.949        0.0645        0.0306 1.16000  0.2810
#>      2 30      3    3.080        0.1000        0.1030 0.00187  0.9660
#>      3 31      5    4.900        0.1610        0.1580 0.00201  0.9640
#>      4 30      1    6.760        0.0333        0.2250 4.91000  0.0267
#>      5 31      4    7.330        0.1290        0.2360 1.51000  0.2190
#>      6 30      7    7.160        0.2330        0.2390 0.00353  0.9530
#>      7 31     12    6.600        0.3870        0.2130 4.41000  0.0357
#>      8 30     10    6.860        0.3330        0.2290 1.44000  0.2310
#>      9 31     11   12.300        0.3550        0.3970 0.14100  0.7070
#>     10 30     13   12.200        0.4330        0.4060 0.05720  0.8110
#>  mean_survival mean_cumhaz
#>          0.963      0.0306
#>          0.881      0.1030
#>          0.814      0.1580
#>          0.757      0.2250
#>          0.717      0.2360
#>          0.681      0.2390
#>          0.668      0.2130
#>          0.655      0.2290
#>          0.520      0.3970
#>          0.396      0.4060
#> 
#> Overall: chi-sq = 13.6 on 9 df, p = 0.136 
# }
```
