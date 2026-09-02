# Parameter names for a phase specification, in `theta` order

Returns the names of the `theta` vector a multiphase
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
fit builds from `phases`, in the order `theta` requires. Use it to check
a hand-written starting vector against the specification it is meant to
go with, before any fit runs.

## Usage

``` r
hzr_theta_names(phases, covariates = NULL)
```

## Arguments

- phases:

  A non-empty list of
  [`hzr_phase()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase.md)
  objects, as passed to
  [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md).

- covariates:

  Optional named list mapping phase name to that phase's covariate
  column names, e.g. `list(early = c("age", "sex"))`. Phases absent from
  the list are treated as having no covariates. Names are used verbatim,
  so they must match the columns the fit will see.

## Value

A character vector whose **order is the required `theta` order**. Its
length is the number of parameters the specification implies, so
`length(hzr_theta_names(phases))` is the length `theta` must have.

## Details

`theta` is positional and its entries are not on a common scale. For an
`early` + `late` pair the layout is

    early.log_mu, early.log_t_half, early.nu, early.m,
    late.log_mu,  late.log_tau,     late.gamma, late.alpha, late.eta

so the late phase logs `mu` and `tau` while carrying `gamma`, `alpha`
and `eta` on the natural scale. **Wrapping the wrong element in
[`log()`](https://rdrr.io/r/base/Log.html) produces a fit, not an
error**, which is why a comment describing the order is not enough and
this returns the real thing. The order is a property of the
specification, so it changes the moment a phase is added, removed or
retyped.

Phase names come from `names(phases)`; unnamed phases are labelled
`phase_1`, `phase_2`, ... by the same validation
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
applies, so the labels here are the labels a fit will use.

## See also

[`hzr_phase()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase.md)
for building a specification,
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
for the fit whose `theta` this names.

## Examples

``` r
phases <- list(early = hzr_phase("cdf"), late = hzr_phase("g3"))
hzr_theta_names(phases)
#> [1] "early.log_mu"     "early.log_t_half" "early.nu"         "early.m"         
#> [5] "late.log_mu"      "late.log_tau"     "late.gamma"       "late.alpha"      
#> [9] "late.eta"        

# Check a hand-written starting vector before fitting.
theta0 <- c(log(0.05), log(0.2), 0, -0.4, log(0.03), log(1), 1, 1, 1)
stopifnot(length(theta0) == length(hzr_theta_names(phases)))
setNames(theta0, hzr_theta_names(phases))
#>     early.log_mu early.log_t_half         early.nu          early.m 
#>        -2.995732        -1.609438         0.000000        -0.400000 
#>      late.log_mu     late.log_tau       late.gamma       late.alpha 
#>        -3.506558         0.000000         1.000000         1.000000 
#>         late.eta 
#>         1.000000 

# With covariates on one phase.
hzr_theta_names(phases, covariates = list(early = c("age", "sex")))
#>  [1] "early.log_mu"     "early.log_t_half" "early.nu"         "early.m"         
#>  [5] "early.age"        "early.sex"        "late.log_mu"      "late.log_tau"    
#>  [9] "late.gamma"       "late.alpha"       "late.eta"        
```
