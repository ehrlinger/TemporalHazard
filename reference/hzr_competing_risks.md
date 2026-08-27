# Competing risks cumulative incidence

Compute cumulative incidence functions for multiple competing event
types using the Aalen-Johansen estimator with Greenwood variance. This
is the R equivalent of the SAS `markov.sas` macro.

## Usage

``` r
hzr_competing_risks(time, event)

# S3 method for class 'hzr_competing_risks'
print(x, digits = 4, ...)
```

## Arguments

- time:

  Numeric vector of follow-up times.

- event:

  Integer vector of event type indicators: 0 = censored, 1 = event type
  1, 2 = event type 2, etc.

- x:

  An `hzr_competing_risks` object.

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

- n_event_1, n_event_2, ...:

  Events of each type at this time.

- n_censor:

  Number censored at this time.

- surv:

  Overall event-free survival (freedom from all events).

- incid_1, incid_2, ...:

  Cumulative incidence for each event type.

- se_surv:

  Standard error of overall survival.

- se_1, se_2, ...:

  Standard error of each cumulative incidence.

## Details

Unlike the naive 1 - KM estimator (which overestimates incidence when
competing risks exist), this provides the correct marginal cumulative
incidence for each event type.

## Note

This estimator is **unweighted**: every observation contributes a unit
count to the at-risk and event tallies. There is no `weights` argument,
so case or inverse-probability weights are not yet supported for
competing-risks incidence (unlike
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md),
which accepts `weights`). Pre-expand weighted rows to individual records
if an approximate weighted estimate is needed.

## References

Aalen O, Johansen S (1978). An empirical transition matrix for
non-homogeneous Markov chains based on censored observations. *Scand J
Statist* 5(3):141–150.

Kalbfleisch JD, Prentice RL (1980). *The Statistical Analysis of Failure
Time Data.* Wiley, New York.

## See also

[`hzr_kaplan()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_kaplan.md)
for single-event survival estimation.

## Examples

``` r
data(valves)
valves_cc <- na.omit(valves)
# Combine death and PVE into a competing risks event variable
# 0 = censored, 1 = death, 2 = PVE
event_cr <- ifelse(valves_cc$dead == 1, 1L,
                   ifelse(valves_cc$pve == 1, 2L, 0L))
time_cr <- pmin(valves_cc$int_dead, valves_cc$int_pve)
cr <- hzr_competing_risks(time_cr, event_cr)
head(cr)
#> Competing risks cumulative incidence
#> Event types: 2  | Time points: 6 
#> Final survival: 0.9888 
#> Final incid_1 : 0.0112 
#> Final incid_2 : 0 
#> 
#>     time n_risk n_event_1 n_event_2 n_censor   surv incid_1 incid_2 se_surv
#>  0.00068   1523         1         0        0 0.9993  0.0007       0  0.0007
#>  0.00137   1522         4         0        0 0.9967  0.0033       0  0.0015
#>  0.00171   1518         1         0        0 0.9961  0.0039       0  0.0016
#>  0.00205   1517         3         0        0 0.9941  0.0059       0  0.0020
#>  0.00274   1514         4         0        0 0.9915  0.0085       0  0.0024
#>  0.00411   1510         4         0        0 0.9888  0.0112       0  0.0027
#>    se_1 se_2
#>  0.0007    0
#>  0.0015    0
#>  0.0016    0
#>  0.0020    0
#>  0.0024    0
#>  0.0027    0
```
