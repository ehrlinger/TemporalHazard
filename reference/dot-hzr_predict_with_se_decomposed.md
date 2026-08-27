# Per-phase + total cumulative-hazard predictions with delta-method CLs

Long-format companion to
[`.hzr_predict_with_se()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_predict_with_se.md)
for multiphase models under `decompose = TRUE`. Computes log-scale CLs
for the total cumulative hazard and for each phase's additive
contribution `H_j(t) = mu_j(x) Phi_j(t)`. Per-phase CLs use only that
phase's parameter block, so they do NOT sum to the total CL (cross-phase
covariance contributes only to the total).

## Usage

``` r
.hzr_predict_with_se_decomposed(
  object,
  time,
  x_list,
  cov_counts,
  phases,
  level = 0.95
)
```

## Arguments

- object:

  Fitted multiphase `hazard` object.

- time:

  Numeric prediction times (length n).

- x_list:

  Named list of per-phase design matrices.

- cov_counts:

  Named integer covariate counts.

- phases:

  Named list of `hzr_phase` objects.

- level:

  Confidence level.

## Value

Long `data.frame(time, component, fit, se.fit, lower, upper)`;
`component` is an ordered factor with levels
`c("total", names(phases))`.
