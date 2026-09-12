# Extract starting values from a phase specification

Returns initial theta sub-vector on the estimation (internal) scale:
log(mu), then log(t_half), nu, m for `"cdf"`/`"hazard"` phases or
log(tau), gamma, alpha, eta for `"g3"` (nothing for `"constant"`),
followed by zeros for covariate coefficients.

## Usage

``` r
.hzr_phase_start(phase, n_covariates = 0L, mu_start = 0.1)
```

## Arguments

- phase:

  An `hzr_phase` object.

- n_covariates:

  Integer; number of covariate columns.

- mu_start:

  Numeric scalar; initial scale parameter (default 0.1).

## Value

Unnamed numeric vector of starting values.
