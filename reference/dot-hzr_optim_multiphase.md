# Fit a multiphase additive hazard model via maximum likelihood

Assembles starting values from phase specifications, resolves per-phase
design matrices, and delegates to `.hzr_optim_generic()`.

## Usage

``` r
.hzr_optim_multiphase(
  time,
  status,
  time_lower = NULL,
  time_upper = NULL,
  x = NULL,
  theta_start = NULL,
  weights = NULL,
  control = list(),
  phases,
  objective = c("likelihood", "sas"),
  formula_global = NULL,
  data = NULL
)
```

## Arguments

- time:

  Numeric vector of follow-up times.

- status:

  Numeric event indicator vector.

- time_lower:

  Optional lower bounds for interval censoring.

- time_upper:

  Optional upper bounds for left/interval censoring.

- x:

  Global design matrix (n x p) or NULL.

- theta_start:

  Starting parameter vector (full internal scale). If NULL, assembled
  automatically from phase specs.

- control:

  Named list of control options.

- phases:

  Named list of validated `hzr_phase` objects.

- objective:

  Which interval-censored contribution to accumulate; see
  [`.hzr_logl_interval()`](https://ehrlinger.github.io/temporal_hazard/reference/dot-hzr_logl_interval.md).

- formula_global:

  The global formula (used when phases have no phase-specific formula).

- data:

  Data frame containing covariates (needed for phase-specific formula
  evaluation).

## Value

List with par (internal scale), value, convergence, vcov, etc.
