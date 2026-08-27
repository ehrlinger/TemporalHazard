# Summarize a hazard model

Returns a compact summary of a `hazard` object, including model
metadata, fit diagnostics, and coefficient-level statistics when
available.

## Usage

``` r
# S3 method for class 'hazard'
summary(object, ...)
```

## Arguments

- object:

  A `hazard` object.

- ...:

  Unused; for S3 compatibility.

## Value

An object of class `summary.hazard`.

## See also

[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
for model fitting,
[`predict.hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/predict.hazard.md)
for predictions.

[`vignette("fitting-hazard-models")`](https://ehrlinger.github.io/TemporalHazard/articles/fitting-hazard-models.md)
for fitting workflows,
[`vignette("inference-diagnostics")`](https://ehrlinger.github.io/TemporalHazard/articles/inference-diagnostics.md)
for bootstrap CIs and diagnostics.

## Examples

``` r
# -- Single-phase Weibull summary ------------------------------------
fit <- hazard(time = rexp(30, 0.5), status = rep(1L, 30),
              theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
summary(fit)
#> hazard model summary
#>   observations: 30 
#>   predictors:   0 
#>   dist:         weibull 
#>   engine:       native-r-m2 
#>   converged:    TRUE 
#>   log-lik:      -49.8926 
#>   evaluations: fn=15, gr=6
#> 
#> Coefficients:
#>     estimate std_error   z_stat      p_value
#> mu 0.5284203 0.1084224 4.873719 1.095168e-06
#> nu 0.9378551 0.1375396 6.818802 9.180267e-12

# \donttest{
# -- Multiphase model summary ----------------------------------------
set.seed(42)
n   <- 200
dat <- data.frame(
  time   = rexp(n, rate = 0.25) + 0.01,
  status = rbinom(n, size = 1, prob = 0.65)
)
fit_mp <- hazard(
  survival::Surv(time, status) ~ 1,
  data   = dat,
  dist   = "multiphase",
  phases = list(
    early = hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
                       fixed = "shapes"),
    late  = hzr_phase("cdf", t_half = 5,   nu = 1, m = 0,
                       fixed = "shapes")
  ),
  fit     = TRUE,
  control = list(n_starts = 5, maxit = 1000)
)
summary(fit_mp)
#> Multiphase hazard model (2 phases)
#>   observations: 200 
#>   predictors:   0 
#>   dist:         multiphase 
#>   phase 1:      early - cdf (early risk)
#>   phase 2:      late - cdf (late risk)
#>   engine:       native-r-m2 
#>   converged:    TRUE 
#>   log-lik:      -428.716 
#>   evaluations: fn=30, gr=6
#> 
#> Coefficients (internal scale):
#> 
#>   Phase: early (cdf)
#>                estimate std_error    z_stat      p_value
#>   log_mu     -2.1152885 0.2910719 -7.267238 3.669127e-13
#>   log_t_half -0.6931472        NA        NA           NA
#>   nu          2.0000000        NA        NA           NA
#>   m           0.0000000        NA        NA           NA
#> 
#>   Phase: late (cdf)
#>               estimate  std_error   z_stat      p_value
#>   log_mu     0.5511614 0.09532828 5.781719 7.394105e-09
#>   log_t_half 1.6094379         NA       NA           NA
#>   nu         1.0000000         NA       NA           NA
#>   m          0.0000000         NA       NA           NA
# }
```
