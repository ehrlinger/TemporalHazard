# Parse Surv() formula for hazard modeling

Extracts time, status, time_lower, time_upper, and predictors from a
formula of the form `Surv(time, status) ~ x1 + x2 + ...`. Supports
right-censored, left-censored, interval-censored, and counting-process
(start-stop) data.

## Usage

``` r
.hzr_parse_formula(formula, data)
```

## Arguments

- formula:

  A formula object with Surv() on the LHS.

- data:

  A data frame containing variables referenced in the formula.

## Value

A list with elements: time, status, time_lower, time_upper, x

## Details

For counting-process (start-stop) data, use `Surv(start, stop, event)`.
The start times are returned as `time_lower` and stop times as `time`,
enabling the likelihood to compute `H(stop) - H(start)` per epoch.

`Surv()` and this package code censoring status differently, so the
returned `status` is translated, not passed through:

|          |                |               |                   |
|----------|----------------|---------------|-------------------|
| Meaning  | TemporalHazard | `Surv` "left" | `Surv` "interval" |
| left     | `-1`           | `0`           | `2`               |
| right    | `0`            | –             | `0`               |
| event    | `1`            | `1`           | `1`               |
| interval | `2`            | –             | `3`               |

Under `type = "interval"`, `Surv()` reuses the `time2` column to hold
the status of any non-interval row, so an upper bound is read only where
the row is genuinely interval-censored.
