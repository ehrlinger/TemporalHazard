# Bootstrap resampling for hazard model coefficients

Resample data with replacement, refit the hazard model on each
replicate, and accumulate coefficient distributions. Returns a tidy data
frame of per-replicate estimates with summary statistics. This is the R
equivalent of the SAS `bootstrap.hazard.sas` macro.

## Usage

``` r
hzr_bootstrap(
  object,
  n_boot = 200L,
  fraction = 1,
  seed = NULL,
  verbose = FALSE,
  scope = NULL,
  direction = c("both", "forward", "backward"),
  criterion = c("wald", "aic"),
  slentry = 0.3,
  slstay = 0.2,
  max_steps = 50L,
  max_move = 4L,
  force_in = character(),
  force_out = character(),
  ...
)

# S3 method for class 'hzr_bootstrap'
print(x, digits = 4, ...)
```

## Arguments

- object:

  A fitted `hazard` object (with `fit = TRUE`).

- n_boot:

  Integer: number of bootstrap replicates (default 200).

- fraction:

  Numeric in (0, 1\]: fraction of data to sample per replicate (default
  1.0 for full bootstrap; \< 1 for bagging).

- seed:

  Optional integer random seed for reproducibility. When supplied,
  `set.seed(seed)` is called at function entry, jumping the global RNG
  to the seeded state; it is not restored on exit. Pass `NULL` (the
  default) to skip the
  [`set.seed()`](https://rdrr.io/r/base/Random.html) call and start from
  the caller's current RNG state. Note that the bootstrap consumes
  random numbers either way, so the global RNG state will advance during
  the call – `seed = NULL` avoids the *reset* at entry, not the advance
  during resampling.

- verbose:

  Logical; if `TRUE`, print progress every 50 replicates.

- scope:

  Candidate variable scope for embedded stepwise selection during each
  bootstrap replicate. `NULL` (default) preserves the original
  fixed-formula bootstrap: every replicate refits `object`'s exact
  model, and `summary$pct` is always ~100. When supplied (a one-sided
  formula, character vector, or – for multiphase fits – a named list of
  one-sided formulas keyed by phase, matching
  [`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)'s
  `scope`), each replicate runs a fresh
  [`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
  selection instead; see Details.

- direction, criterion, slentry, slstay, max_steps, max_move, force_in,
  force_out:

  Passed through to
  [`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
  on each replicate when `scope` is supplied; ignored when
  `scope = NULL`. See
  [`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
  for definitions and defaults.

- ...:

  Additional arguments forwarded to
  [`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
  (e.g. `control = list(n_starts = 1)`) when `scope` is supplied;
  ignored otherwise.

- x:

  An `hzr_bootstrap` object.

- digits:

  Number of decimal places for formatting.

## Value

A list with class `"hzr_bootstrap"` containing:

- replicates:

  Data frame with columns `replicate`, `parameter`, and `estimate` – one
  row per parameter per successful replicate.

- summary:

  Data frame with columns `parameter`, `n`, `pct`, `mean`, `sd`, `min`,
  `max`, `ci_lower`, `ci_upper` – one row per parameter. In
  `mode = "select"`, `pct` is the selection frequency and the other
  statistics are conditional on selection.

- n_success:

  Number of successfully converged replicates.

- n_failed:

  Number of replicates that failed to converge.

- mode:

  `"refit"` (fixed-formula bootstrap) or `"select"` (embedded stepwise
  selection).

- scope:

  Only present when `mode == "select"`: the candidate scope used.

## Details

When `scope` is supplied, each replicate instead runs a fresh
[`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
selection on the resampled data (starting from a fixed-shape refit of
`object`) instead of refitting `object`'s exact formula. This is the R
equivalent of the SAS `%HAZBOOT` macro: fit the hazard shape with no
covariates (fixing it via `hzr_phase(..., fixed = "shapes")`), then
bootstrap-screen candidate covariates for how often they enter the
model. `summary$pct` then reports the selection frequency across
replicates, and `summary$mean`/`sd`/`ci_*` describe the coefficient
distribution conditional on selection.

## See also

[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
for model fitting,
[`vcov.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/vcov.hazard.md)
for Hessian-based standard errors,
[`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
for the selection procedure used when `scope` is supplied.

## Examples

``` r
# \donttest{
data(avc)
avc <- na.omit(avc)
fit <- hazard(
  survival::Surv(int_dead, dead) ~ age + mal,
  data  = avc,
  dist  = "weibull",
  theta = c(mu = 0.01, nu = 0.5, 0, 0),
  fit   = TRUE
)
bs <- hzr_bootstrap(fit, n_boot = 50, seed = 123)
print(bs)
#> Bootstrap inference for hazard model
#> Mode: fixed refit 
#> Replicates: 50 successful, 0 failed
#> 
#>  parameter  n pct    mean     sd     min     max ci_lower ci_upper
#>         mu 50 100  0.0004 0.0004  0.0000  0.0023   0.0000   0.0013
#>         nu 50 100  0.2219 0.0144  0.1963  0.2618   0.1986   0.2462
#>        age 50 100 -0.0062 0.0028 -0.0129 -0.0022  -0.0125  -0.0029
#>        mal 50 100  0.8464 0.2701  0.2673  1.3909   0.4337   1.3471

# Embedded stepwise selection: screen candidate covariates for how
# often they enter the model across resamples (R equivalent of SAS
# %HAZBOOT).
base <- hazard(
  survival::Surv(int_dead, dead) ~ 1,
  data  = avc,
  dist  = "weibull",
  theta = c(mu = 0.01, nu = 0.5),
  fit   = TRUE
)
bs_sel <- hzr_bootstrap(base, n_boot = 20, seed = 123,
                         scope = ~ age + mal,
                         slentry = 0.3, slstay = 0.2)
print(bs_sel)
#> Bootstrap inference for hazard model
#> Mode: embedded stepwise selection 
#> Replicates: 20 successful, 0 failed
#> 
#>  parameter  n pct    mean     sd     min     max ci_lower ci_upper
#>        age 20 100 -0.0068 0.0031 -0.0125 -0.0022  -0.0116  -0.0027
#>         mu 20 100  0.0007 0.0009  0.0001  0.0037   0.0001   0.0028
#>         nu 20 100  0.2249 0.0143  0.2040  0.2620   0.2061   0.2543
#>        mal 14  70  0.9749 0.2623  0.5514  1.3572   0.5734   1.3368
# }
```
