# Print method for fitted hazard models

Compact one-block summary of a fitted `hazard` object: sample size,
number of predictors, distribution, theta vector, and log-likelihood. S3
dispatch only – users call `print(fit)` rather than invoking this
directly.

## Usage

``` r
# S3 method for class 'hazard'
print(x, ...)
```

## Arguments

- x:

  A `hazard` object returned by
  [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md).

- ...:

  Additional arguments (ignored).

## Value

`x`, invisibly.
