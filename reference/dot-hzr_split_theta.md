# Split the full theta vector into per-phase sub-vectors

Split the full theta vector into per-phase sub-vectors

## Usage

``` r
.hzr_split_theta(theta, phases, covariate_counts)
```

## Arguments

- theta:

  Numeric vector: full parameter vector (internal scale).

- phases:

  Named list of validated `hzr_phase` objects.

- covariate_counts:

  Named integer vector: number of covariates per phase.

## Value

Named list of numeric vectors, one per phase.
