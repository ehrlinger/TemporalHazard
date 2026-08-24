# Extract the captured console trace from an `hzr_stepwise` fit

Every run of
[`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
records the header, per-step lines, and final summary regardless of the
`trace` flag. This accessor returns the full character vector for
display or logging.

## Usage

``` r
stepwise_trace(fit)
```

## Arguments

- fit:

  An `hzr_stepwise` object.

## Value

Character vector, one element per console line.

## See also

[`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md),
which produces the object this accessor reads.

## Examples

``` r
data(avc)
avc <- na.omit(avc)
base <- hazard(survival::Surv(int_dead, dead) ~ age,
               data = avc, dist = "weibull", fit = TRUE,
               theta = c(mu = 0.01, nu = 0.5, 0))
# \donttest{
sw <- hzr_stepwise(base, scope = ~ age + mal,
                   data = avc, direction = "forward",
                   control = list(n_starts = 1))
#> Stepwise selection (direction = forward, criterion = score, slentry = 0.30, slstay = 0.20)
#> 
#> Step 1: ENTER  mal   (p = 0.001)
#> (no further action after 1 step)
#> 
#> Final model: 2 covariates, logLik = -223.55, AIC = 455.09
cat(stepwise_trace(sw), sep = "\n")
#> Stepwise selection (direction = forward, criterion = score, slentry = 0.30, slstay = 0.20)
#> 
#> Step 1: ENTER  mal   (p = 0.001)
#> (no further action after 1 step)
#> 
#> Final model: 2 covariates, logLik = -223.55, AIC = 455.09
# }
```
