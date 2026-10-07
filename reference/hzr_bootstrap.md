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
  likelihood that could not be evaluated), did not converge, or returned
  a finite objective but no parameter estimates. A single-distribution
  refit that ends where the likelihood is not defined reports no
  objective, so it fails as a non-finite objective, not at the sentinel
  (#486).

- failure_reasons:

  Named integer vector counting why replicates failed, most common
  first: the refit's error message (or `"error with an empty message"`),
  `"refit returned a <class>, not a fit object"` (`"an <class>"` when
  the class begins with a vowel; or
  `"... with no \code{fit}, not a fit object"`),
  `"refit returned no parameter estimates"`, or
  `"non-finite objective (did not converge)"` (the reason for a
  single-distribution refit that ends where the likelihood is not
  defined, \#486), or
  `"objective at the optimizer's -1e10 sentinel (no log-likelihood)"` (a
  refit whose reported objective is still the sentinel), or
  `"refit did not converge (converged = FALSE)"` (a refit that reports
  `converged = FALSE`, including one stopped on a score that is not
  finite, \#518). It sums to `n_failed`, and is an empty named integer
  vector, never `NULL`, when none failed. When every replicate fails,
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
  selection frequency. `rows_differ` (under `criterion = "aic"`) and
  `loglik_below_base` (under any criterion) mark entries a replicate
  declined without comparing them. Such an entry counts as not selected
  unless a later step of the same replicate tested and entered it, so
  these can understate a selection frequency; a warning gives how many
  replicates completed after declining one. The tally counts attempts,
  not distinct candidates. Empty in refit mode.

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

- unresolved:

  Only present when `mode == "select"`: the `$scope$unresolved` record
  of the up-front
  [`hzr_stepwise()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
  screen on the original data, the names in `force_in`, `force_out` and
  a character `scope` that matched nothing and were ignored.
  [`print()`](https://rdrr.io/r/base/print.html) shows a line when any
  is non-empty.

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

A replicate in which any fit (the base refit, a stepwise refit or the
final fit) may not be a maximum, by the rule
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
warns on with class `"hzr_possible_false_maximum"`, is counted, and
`hzr_bootstrap()` warns once with that count. Such replicates are kept
in the pooled results.

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
#> Warning: 9 of 20 successful replicates entered at least one variable on a Wald test rather than on the score statistic (9 candidate(s) in total), because the score could not test them. Those entries were decided by a different criterion from the rest of the run, and they are in the pooled frequencies on the same footing as everything else. See `$n_wald_fallbacks`.
print(bs_sel)
#> Bootstrap inference for hazard model
#> Mode: embedded stepwise selection 
#> Replicates: 20 successful, 0 failed
#> 
#>  parameter  n pct    mean     sd     min     max ci_lower ci_upper
#>        age 20 100 -0.0067 0.0030 -0.0120 -0.0022  -0.0114  -0.0027
#>         mu 20 100  0.0005 0.0006  0.0001  0.0023   0.0001   0.0021
#>         nu 20 100  0.2255 0.0142  0.2041  0.2620   0.2061   0.2543
#>        mal 19  95  0.8695 0.2897  0.4850  1.3571   0.4970   1.3289
# }
```
