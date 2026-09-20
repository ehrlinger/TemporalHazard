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
  criterion = c("score", "wald", "aic"),
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
  the caller's current RNG state. The bootstrap consumes random numbers
  either way, so the global RNG state will advance during the call;
  `seed = NULL` avoids the *reset* at entry, not the advance during
  resampling.

- verbose:

  Logical; if `TRUE`, display a text progress bar over the `n_boot`
  replicates (via
  [`utils::txtProgressBar()`](https://rdrr.io/r/utils/txtProgressBar.html)).

- scope:

  **Experimental.** Candidate variable scope for embedded stepwise
  selection during each bootstrap replicate. The argument and the shape
  of the object it returns may change in a future release; see the
  "Selection mode is experimental" section below. `NULL` (default)
  preserves the original fixed-formula bootstrap: every replicate refits
  `object`'s exact model, and `summary$pct` is always ~100. When
  supplied (a one-sided formula, character vector, or, for multiphase
  fits, a named list of one-sided formulas keyed by phase, matching
  [`hzr_stepwise()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)'s
  `scope`; a two-sided formula is an error), each replicate runs a fresh
  [`hzr_stepwise()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
  selection instead; see Details.

- direction, slentry, slstay, max_steps, max_move, force_in, force_out:

  Passed through to
  [`hzr_stepwise()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
  on each replicate when `scope` is supplied. With `scope = NULL`
  nothing reads them, so a value other than the default is an error; the
  default's own value is accepted, whether or not it was passed, so that
  a wrapper forwarding its defaults still works.
  `direction = "backward"` with a non-empty `scope` is an error too: a
  backward screen does not read `scope`. For a backward screen on each
  replicate, pass an empty scope such as `~ 1`. See
  [`hzr_stepwise()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
  for definitions and defaults.

- criterion:

  Entry / retention rule passed through to
  [`hzr_stepwise()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
  on each replicate when `scope` is supplied; with `scope = NULL`, a
  value other than the default is an error. One of `"score"` (default),
  `"wald"`, or `"aic"`. `"score"` reproduces C/SAS HAZARD's `SELECTION`
  statistic and needs no per-candidate refit, which is what makes a
  bootstrap screen over many candidates tractable. Following SAS, the
  variance used during *selection* is approximate (shaping-parameter
  covariances are ignored); final-model standard errors are unaffected.
  For single-distribution fits, `"score"` computes the observed
  information numerically via the suggested numDeriv package and errors
  if it is not installed; a multiphase fit uses the analytic Hessian
  instead and does not need it. See
  [`hzr_stepwise()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md).

- ...:

  Additional arguments forwarded to
  [`hzr_stepwise()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
  (e.g. `control = list(maxit = 500)`, which every `dist` reads) when
  `scope` is supplied. Only `control` and an `objective` equal to the
  fit's are forwarded. `trace` is accepted and ignored, since each
  replicate's screen runs quietly; use `verbose` for progress. Any other
  name is an error, and so is any `...` argument without `scope`.

- x:

  An `hzr_bootstrap` object.

- digits:

  Number of decimal places for formatting.

## Value

A list with class `"hzr_bootstrap"` containing:

- replicates:

  Data frame with columns `replicate`, `parameter`, and `estimate`, one
  row per parameter per successful replicate.

- summary:

  Data frame with columns `parameter`, `n`, `pct`, `mean`, `sd`, `min`,
  `max`, `ci_lower`, `ci_upper`, one row per parameter. In
  `mode = "select"`, `pct` is the selection frequency and the other
  statistics are conditional on selection. A free parameter whose `sd`
  is 0, to within rounding, across two or more replicates draws a
  warning naming it: the replicates did not re-estimate it. Parameters
  the fit holds fixed are exempt: those held by `hzr_phase(fixed = )`,
  shapes a constraint derives, and a conserved `log_mu`.

- n_success:

  Number of successfully converged replicates.

- n_failed:

  Number of replicates that failed: the refit stopped with an error,
  returned something other than a fit, returned a non-finite objective
  or one at the optimizer's -1e10 sentinel (which stands in for a
  likelihood that could not be evaluated), or returned a finite
  objective but no parameter estimates.

- failure_reasons:

  Named integer vector counting why replicates failed, most common
  first: the refit's error message (or `"error with an empty message"`),
  `"refit returned a <class>, not a fit object"` (`"an <class>"` when
  the class begins with a vowel; or
  `"... with no \code{fit}, not a fit object"`),
  `"refit returned no parameter estimates"`, or
  `"non-finite objective (did not converge)"`, or
  `"objective at the optimizer's -1e10 sentinel (no log-likelihood)"`.
  It sums to `n_failed`, and is an empty named integer vector, never
  `NULL`, when none failed. When every replicate fails,
  `hzr_bootstrap()` also warns, naming the most common reason.

- n_uncomputable_replicates:

  Select mode only: number of otherwise successful replicates whose
  screen stopped because no remaining candidate could be tested (its
  score statistic, or for a removal its Wald statistic, could not be
  computed), rather than because no candidate met `slentry` or `slstay`.
  A non-zero count means every reported selection frequency is biased: a
  candidate such a replicate could not test for entry counts as not
  selected, and one it could not test for removal as selected. A
  replicate that went on after deciding a variable without a Wald test
  (`wald_no_variance`) is not counted here unless it also stopped, and
  gets a warning of its own either way. Always `0` in refit mode.

- uncomputable_reasons:

  Select mode only: named integer vector counting *why* candidate scores
  were unavailable, summed over every replicate.
  `information_indefinite` is the one to read first: it marks candidates
  whose effect is too large for the score test's approximation at zero,
  typically strong variables. Those are refit and Wald-tested
  automatically, so a candidate reaching this count is one whose refit
  also failed and which therefore went untested, understating its
  selection frequency. Empty in refit mode.

- n_nonmonotone_replicates:

  Select mode only: number of otherwise successful replicates in which a
  forward step *lowered* the log-likelihood. Entered models are nested,
  so this cannot occur at the optimum; such a replicate continued from a
  refit that did not converge, and its later selections are pooled on
  the same footing as any other. Always `0` in refit mode.

- n_wald_fallback_replicates, n_wald_fallbacks:

  Select mode only: the number of otherwise successful replicates that
  entered at least one variable on a Wald test instead of the score
  statistic, and the total number of such entries across all replicates.
  The score criterion declines a candidate whose observed information is
  indefinite at `beta = 0` (which happens when the effect is *large*),
  so those candidates are refit and Wald-tested rather than dropped. A
  high count means much of the selection was decided by a different
  criterion from the one requested, which matters most here: these
  entries drive the pooled selection frequencies on the same footing as
  every other. Always `0` in refit mode.

- mode:

  `"refit"` (fixed-formula bootstrap) or `"select"` (embedded stepwise
  selection).

- scope:

  Only present when `mode == "select"`: the candidate scope used.

## Details

When `scope` is supplied, each replicate instead runs a fresh
[`hzr_stepwise()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
selection on the resampled data (starting from a fixed-shape refit of
`object`) instead of refitting `object`'s exact formula. This is the R
equivalent of the SAS `%HAZBOOT` macro: fit the hazard shape with no
covariates (fixing it via `hzr_phase(..., fixed = "shapes")`), then
bootstrap-screen candidate covariates for how often they enter the
model. `summary$pct` then reports the selection frequency across
replicates, and `summary$mean`/`sd`/`ci_*` describe the coefficient
distribution conditional on selection.

## Selection mode is experimental

Everything reached through `scope` (the selection arguments, and the
`summary$pct` selection frequencies they produce) is new and should be
treated as unstable. The fixed-formula bootstrap (`scope = NULL`) is not
affected and has been stable since 0.9.3.

Two reasons to expect change. The first is that the design is still
being read off real runs rather than settled in advance: production
screens have already moved the defaults once and turned up several ways
a screen could report success while selecting nothing.

The second is scale, and it is the one to plan around. A screen over a
large candidate pool runs for hours, and this function writes nothing
until its final replicate, so a run that dies late loses everything.
There is no built-in way to split one screen across processes and
combine the parts. If you are running at that scale, drive
`hzr_bootstrap()` in chunks from your own script and pool the replicates
yourself: deriving each chunk's seed from its chunk number, offsetting
replicate ids so a variable selected in two chunks is not counted once,
and recomputing frequencies from the pooled replicates rather than
averaging across chunks. Whatever eventually covers that inside the
package may well change this function's interface.

## See also

[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
for model fitting,
[`vcov.hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/vcov.hazard.md)
for Hessian-based standard errors,
[`hzr_stepwise()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
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
#>         mu 50 100  0.0004 0.0004  0.0000  0.0023   0.0000   0.0012
#>         nu 50 100  0.2219 0.0144  0.1957  0.2620   0.1984   0.2457
#>        age 50 100 -0.0062 0.0028 -0.0129 -0.0022  -0.0125  -0.0029
#>        mal 50 100  0.8481 0.2710  0.2720  1.3901   0.4319   1.3463

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
#> Warning: Stepwise selection stopped because no remaining candidate could be tested for entry: its score statistic, or its Wald statistic for want of a variance, could not be computed (1 candidate score(s) were NA across the run). This is not the same as no candidate meeting `slentry` or `slstay`: the screen stopped without being able to test them. Causes: 1 x the current model's information matrix could not be inverted, so no candidate could be scored at that step.
#> Warning: 9 of 20 successful replicates stopped because no remaining candidate could be tested -- the score statistic, or for a removal the Wald statistic, could not be computed -- rather than because no candidate met `slentry` or `slstay`. Every reported selection frequency is biased by them: a candidate they could not test for entry counts as not selected, and one they could not test for removal as selected. Causes: 10 x the current model's information matrix could not be inverted, so no candidate could be scored at that step.
print(bs_sel)
#> Bootstrap inference for hazard model
#> Mode: embedded stepwise selection 
#> Replicates: 20 successful, 0 failed
#> 
#>  parameter  n pct    mean     sd     min     max ci_lower ci_upper
#>         mu 20 100  0.0004 0.0006  0.0000  0.0023   0.0000   0.0021
#>         nu 20 100  0.2202 0.0148  0.1917  0.2458   0.1963   0.2432
#>        mal 18  90  0.9822 0.3429  0.5118  1.6183   0.5286   1.5530
#>        age 11  55 -0.0067 0.0030 -0.0106 -0.0022  -0.0104  -0.0025
# }
```
