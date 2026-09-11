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
#> 10 groups, 68 observed events, 68 expected
#> 
#>  group  n events expected observed_rate expected_rate   chi_sq p_value
#>      1 31      2    0.933        0.0645        0.0301 1.220000  0.2690
#>      2 30      3    3.040        0.1000        0.1010 0.000498  0.9820
#>      3 31      5    4.720        0.1610        0.1520 0.017200  0.8960
#>      4 30      1    6.860        0.0333        0.2290 5.000000  0.0253
#>      5 31      4    7.280        0.1290        0.2350 1.480000  0.2240
#>      6 30      7    7.120        0.2330        0.2370 0.002150  0.9630
#>      7 31     12    6.570        0.3870        0.2120 4.480000  0.0343
#>      8 30     10    6.830        0.3330        0.2280 1.470000  0.2250
#>      9 31     11   12.400        0.3550        0.3990 0.154000  0.6940
#>     10 30     13   12.300        0.4330        0.4090 0.044800  0.8320
#>  mean_survival mean_cumhaz
#>          0.963      0.0301
#>          0.882      0.1010
#>          0.815      0.1520
#>          0.759      0.2290
#>          0.719      0.2350
#>          0.683      0.2370
#>          0.669      0.2120
#>          0.657      0.2280
#>          0.518      0.3990
#>          0.393      0.4090
#> 
#> Overall: chi-sq = 13.9 on 9 df, p = 0.127 
# }
```
