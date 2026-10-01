# AVC: Atrioventricular Canal Repair

Survival data for 310 patients who underwent repair of atrioventricular
septal defects (congenital heart disease) at the Cleveland Clinic
between 1977 and 1993. Exhibits two identifiable hazard phases: an early
post-operative risk and a constant late phase.

## Usage

``` r
avc
```

## Format

A data frame with 310 rows and 11 variables:

- study:

  Patient identifier

- status:

  NYHA functional class (1–4). A covariate, not a censoring status: the
  event indicator is `dead`.

- inc_surg:

  Surgical grade of AV valve incompetence

- opmos:

  Date of operation (months since January 1967)

- age:

  Age at repair (months)

- mal:

  Malalignment indicator (0/1)

- com_iv:

  Interventricular communication indicator (0/1)

- orifice:

  Associated cardiac anomaly indicator (0/1)

- dead:

  Death indicator (1 = dead, 0 = censored)

- int_dead:

  Follow-up interval to death or last contact (months)

- op_age:

  Interaction term: opmos x age

## Source

Blackstone, Naftel, and Turner (1986)
[doi:10.1080/01621459.1986.10478314](https://doi.org/10.1080/01621459.1986.10478314)
. Cleveland Clinic Foundation.

## See also

[`vignette("fitting-hazard-models")`](https://ehrlinger.github.io/TemporalHazard/articles/fitting-hazard-models.md),
[`vignette("prediction-visualization")`](https://ehrlinger.github.io/TemporalHazard/articles/prediction-visualization.md)

Other datasets:
[`cabgkul`](https://ehrlinger.github.io/TemporalHazard/reference/cabgkul.md),
[`omc`](https://ehrlinger.github.io/TemporalHazard/reference/omc.md),
[`tga`](https://ehrlinger.github.io/TemporalHazard/reference/tga.md),
[`uslife2023`](https://ehrlinger.github.io/TemporalHazard/reference/uslife2023.md),
[`valves`](https://ehrlinger.github.io/TemporalHazard/reference/valves.md)

## Examples

``` r
data(avc)
avc <- na.omit(avc)

# Kaplan-Meier survival
km <- survival::survfit(survival::Surv(int_dead, dead) ~ 1, data = avc)
plot(km, xlab = "Months after AVC repair", ylab = "Survival",
     main = "AVC: Kaplan-Meier survival estimate")


# \donttest{
# Two-phase hazard fit (early CDF + constant -- what AVC supports)
fit <- hazard(
  survival::Surv(int_dead, dead) ~ 1, data = avc,
  dist = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1),
    constant = hzr_phase("constant")
  ),
  fit = TRUE, control = list(n_starts = 5, maxit = 1000)
)
#> Warning: Hessian is not positive-definite at the optimum; standard errors may be unreliable
#> Warning: Non-positive variance estimates; the optimum may not be a proper maximum
#> Warning: Hessian is not positive-definite at the optimum; standard errors may be unreliable
#> Warning: Non-positive variance estimates; the optimum may not be a proper maximum
#> Warning: Hessian is not positive-definite at the optimum; standard errors may be unreliable
#> Warning: Non-positive variance estimates; the optimum may not be a proper maximum
#> Warning: Hessian is not positive-definite at the optimum; standard errors may be unreliable
#> Warning: Non-positive variance estimates; the optimum may not be a proper maximum
#> Warning: Hessian is not positive-definite at the optimum; standard errors may be unreliable
#> Warning: Non-positive variance estimates; the optimum may not be a proper maximum
summary(fit)
#> Multiphase hazard model (2 phases)
#>   observations: 305 
#>   predictors:   0 
#>   dist:         multiphase 
#>   phase 1:      early - cdf (early risk)
#>   phase 2:      constant - constant (flat rate)
#>   engine:       native-r-m2 
#>   converged:    TRUE 
#>   gradient:     relative 2.69e-06 (SAS/C requires <= 6.06e-06; met)
#>   log-lik:      -211.468 
#>   Not done in this run: none
#>   evaluations: fn=5, gr=1
#> 
#> Coefficients (internal scale):
#> 
#>   Phase: early (cdf)
#>                estimate std_error     z_stat      p_value
#>   log_mu     -1.4575071 0.1412491 -10.318701 5.800266e-25
#>   log_t_half -1.7947307 0.3641622  -4.928383 8.291296e-07
#>   nu          1.4542057 0.5504477   2.641860 8.245217e-03
#>   m           0.9265408 0.7085363   1.307683 1.909809e-01
#> 
#>   Phase: constant (constant)
#>           estimate std_error    z_stat      p_value
#>   log_mu -7.523149 0.4714422 -15.95774 2.516989e-57
# }
```
