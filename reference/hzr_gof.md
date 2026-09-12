# Goodness-of-fit: observed vs. predicted events

Compare a fitted hazard model against the nonparametric Kaplan-Meier
estimate by computing observed and expected (parametric) event counts at
each distinct event time. This is the R equivalent of the SAS
`hazplot.sas` macro and implements the conservation-of-events
diagnostic.

## Usage

``` r
hzr_gof(object, time_grid = NULL)
```

## Arguments

- object:

  A fitted `hazard` object (with `fit = TRUE`).

- time_grid:

  Optional numeric vector of time points at which to evaluate the
  parametric model. If `NULL` (default), uses the sorted unique event
  times from the fitted data.

## Value

A data frame with one row per time point and columns:

- time:

  Evaluation time.

- n_risk:

  Number at risk (Kaplan-Meier).

- n_event:

  Number of events at this time.

- n_censor:

  Number censored at this time.

- km_surv:

  Kaplan-Meier survival estimate.

- km_cumhaz:

  Kaplan-Meier cumulative hazard (\\-\log(\text{km\\surv})\\).

- par_surv:

  Parametric survival from the fitted model, at the covariate means for
  a model with covariates.

- par_cumhaz:

  Parametric cumulative hazard, at the covariate means for a model with
  covariates.

- cum_observed:

  Cumulative observed events to this time.

- cum_expected:

  Cumulative expected events: `par_cumhaz` times the number of
  observations exiting the risk set, summed to this time.

- residual:

  Expected minus observed (`cum_expected - cum_observed`).

For multiphase models, additional columns are appended for each phase:
`par_cumhaz_<phase>`.

An attribute `"summary"` is attached with scalar diagnostics: total
observed events, total expected events, and the final residual.

## Details

At each observed event time the function computes:

- The Kaplan-Meier survival and cumulative hazard.

- The parametric survival and cumulative hazard from the fitted model
  (and per-phase components for multiphase models).

- Cumulative observed events vs. cumulative expected events (the
  parametric cumulative hazard at each time, times the number of
  observations leaving the risk set then).

- The running residual (expected minus observed).

For an intercept-only model every patient shares one curve, so the
expected count is the sum of each patient's cumulative hazard at their
own exit time. At the maximum likelihood estimate that sum equals the
number of observed events (the conservation-of-events identity): the
final residual is zero and the printed "Conservation ratio (E/O)" is 1.

A model with covariates is different. The parametric curve, and so the
expected count, is evaluated at the covariate means, one "mean patient"
standing in for everyone. The mean patient's cumulative hazard is not
the average of the patients' cumulative hazards, so the printed E/O is
not the conservation-of-events identity and need not be near 1 at a
correct fit. Read it as a mean-patient check: how closely the curve for
a patient with average covariates follows the whole cohort.

To check conservation of events for a covariate model, sum each
patient's own cumulative hazard at their follow-up time and compare the
total with the event count. Pass the covariate columns plus a `time`
column to
[`predict.hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/predict.hazard.md)
with `type = "cumulative_hazard"`; the Examples show how.

## See also

[`hzr_deciles()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_deciles.md)
for decile-of-risk calibration,
[`predict.hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/predict.hazard.md)
for prediction types.

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
gof <- hzr_gof(fit)
print(gof)
#> Goodness-of-fit: observed vs. expected events
#> Distribution: weibull  | n = 305 
#> 
#> Total observed events: 68 
#> Total expected events: 55.583 
#> Final residual (E - O): -12.417 
#> Conservation ratio (E/O): 0.817 
#> 
#> Use plot columns: time, km_surv, par_surv, cum_observed, cum_expected, residual

# With covariates, hzr_gof() uses the covariate means, so its E/O is a
# mean-patient check.  The conservation-of-events check sums each
# patient's own cumulative hazard at their follow-up time:
nd <- avc[, c("age", "mal")]
nd$time <- avc$int_dead
c(expected = sum(predict(fit, newdata = nd, type = "cumulative_hazard")),
  observed = sum(avc$dead))
#> expected observed 
#>  67.9993  68.0000 

# Plot observed vs expected events
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  ggplot(gof, aes(x = time)) +
    geom_line(aes(y = cum_observed), colour = "#D55E00") +
    geom_line(aes(y = cum_expected), colour = "#0072B2") +
    labs(x = "Time", y = "Cumulative events") +
    theme_minimal()
}

# }
```
