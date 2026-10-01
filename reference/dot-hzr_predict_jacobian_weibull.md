# Jacobian of Weibull predictions with respect to theta

Jacobian of Weibull predictions with respect to theta

## Usage

``` r
.hzr_predict_jacobian_weibull(type, theta, time, x, p)
```

## Arguments

- type:

  Prediction type (see `.hzr_predict_with_se`).

- theta:

  MLE parameter vector `c(log(mu), nu, beta_1, ...)` (#566).

- time:

  Prediction times (may be NULL for `hazard` / `linear_predictor`).

- x:

  Design matrix (n x p_cov) or NULL.

- p:

  Length of theta.

## Value

Numeric n x p Jacobian, with respect to `log(mu)`, `nu` and beta.
