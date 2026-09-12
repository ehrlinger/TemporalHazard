# Gradient of the multiphase log-likelihood

Computes the score vector \\d\ell / d\theta\\ using analytic chain-rule
formulas for `log_mu` and `beta` parameters, and finite-difference
derivatives (via
[`.hzr_phase_derivatives()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_phase_derivatives.md))
for shape parameters (`log_t_half`, `nu`, `m`): central, except that the
`m` stencil stays on the side of `m` near `m = 0`, where the phase
family has a cusp.

## Usage

``` r
.hzr_gradient_multiphase(
  theta,
  time,
  status,
  time_lower = NULL,
  time_upper = NULL,
  x = NULL,
  weights = NULL,
  phases,
  covariate_counts,
  x_list,
  objective = c("likelihood", "sas"),
  sanitize = TRUE,
  ...
)
```

## Arguments

- theta:

  Full parameter vector (internal scale).

- time:

  Numeric vector of follow-up times (n).

- status:

  Numeric event indicator: 1 = event, 0 = right-censored, -1 =
  left-censored, 2 = interval-censored.

- time_lower:

  Optional lower bounds for interval censoring.

- time_upper:

  Optional upper bounds for left/interval censoring.

- x:

  Design matrix (unused directly; kept for interface compatibility).

- phases:

  Named list of validated `hzr_phase` objects.

- covariate_counts:

  Named integer vector of per-phase covariate counts.

- x_list:

  Named list of per-phase design matrices.

- objective:

  Which interval-censored contribution to accumulate; see
  [`.hzr_logl_interval()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_logl_interval.md).
  Exact-event, right-censored and left-censored rows are unaffected.

- sanitize:

  Logical; `FALSE` returns `NA` where the default returns 0. Used by
  SAS/C's acceptance test, which must not read a zero it could not
  compute as a small gradient.

- ...:

  Ignored.

## Value

Numeric vector of length `length(theta)`: the gradient. With
`sanitize = TRUE` (the default) a component that cannot be evaluated is
0, and so is the whole vector at an infeasible point (guards the
optimizer); with `sanitize = FALSE` those are `NA`.
