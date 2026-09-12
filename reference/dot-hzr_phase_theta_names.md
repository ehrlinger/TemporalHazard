# Generate parameter names for a phase's sub-vector

Builds the named labels for the parameter block that belongs to a single
phase in the full theta vector. The block layout is:

## Usage

``` r
.hzr_phase_theta_names(phase, phase_name, covariate_names = character(0))
```

## Arguments

- phase:

  An `hzr_phase` object.

- phase_name:

  Character label for the phase (e.g. `"early"`, `"phase_1"`).

- covariate_names:

  Character vector of covariate column names that this phase uses. Can
  be length 0 if no covariates.

## Value

Character vector of parameter names.

## Details

- For `"cdf"`/`"hazard"`:
  `[log_mu, log_t_half, nu, m, beta_1, ..., beta_p]`

- For `"g3"`:
  `[log_mu, log_tau, gamma, alpha, eta, beta_1, ..., beta_p]`

- For `"constant"`: `[log_mu, beta_1, ..., beta_p]`
