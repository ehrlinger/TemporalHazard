# Print method for hazard summary objects

Formatted console display of
[`summary.hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/summary.hazard.md)
output: distribution, phase list (for multiphase), coefficient table
with standard errors, and log-likelihood. When the post-fit Hessian is
ill-conditioned or not positive-definite, a note warns that the standard
errors may be unreliable; when the Hessian could not be inverted at all,
a note reports that standard errors are unavailable. S3 dispatch only –
users call `print(summary(fit))` rather than invoking this directly.

## Usage

``` r
# S3 method for class 'summary.hazard'
print(x, ...)
```

## Arguments

- x:

  A `summary.hazard` object returned by
  [`summary.hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/summary.hazard.md).

- ...:

  Additional arguments (ignored).

## Value

`x`, invisibly.
