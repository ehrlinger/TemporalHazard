# Extract coefficients from hazard model

Extract coefficients from hazard model

## Usage

``` r
# S3 method for class 'hazard'
coef(object, ...)
```

## Arguments

- object:

  A `hazard` object.

- ...:

  Unused; for S3 compatibility.

## Value

A named numeric vector of fitted parameter estimates, or `NULL` if the
model has not been fitted (`fit = FALSE`).

## Examples

``` r
fit <- hazard(time = rexp(30, 0.5), status = rep(1L, 30),
              theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
coef(fit)
#> [1] 0.4647163 0.9724878
```
