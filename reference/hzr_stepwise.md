# Stepwise covariate selection for a parametric hazard model

Run forward, backward, or two-way stepwise selection on an existing
`hazard` fit using score (Q) statistics, Wald p-values, or AIC deltas as
the entry / retention criterion. Phase-specific entry is supported for
multiphase models: a covariate can enter one phase and not another.

## Usage

``` r
hzr_stepwise(
  fit,
  scope = NULL,
  data,
  direction = c("both", "forward", "backward"),
  criterion = c("score", "wald", "aic"),
  slentry = 0.3,
  slstay = 0.2,
  max_steps = 50L,
  max_move = 4L,
  force_in = character(),
  force_out = character(),
  trace = TRUE,
  ...
)

# S3 method for class 'hzr_stepwise'
print(x, ...)

# S3 method for class 'hzr_stepwise'
summary(object, ...)

# S3 method for class 'summary.hzr_stepwise'
print(x, ...)

# S3 method for class 'hzr_stepwise'
as.data.frame(x, ...)
```

## Arguments

- fit:

  A fitted `hazard` object built via the
  `formula = Surv(...) ~ predictors, data = df` interface.

- scope:

  Candidate set. `NULL` (default) uses every data-frame column not
  already in the model for every phase. For single-distribution fits,
  pass a one-sided formula (`~ age + nyha`) or a character vector of
  names. For multiphase fits, pass a named list of one-sided formulas
  keyed by phase.

- data:

  Data frame the base fit was built on. Required for refits.

- direction:

  Search strategy — one of `"both"` (default), `"forward"`, or
  `"backward"`. Controls whether variables may only enter, only leave,
  or both. See the **Selection direction and criterion** section.

- criterion:

  Entry / retention rule — one of `"score"` (default), `"wald"`, or
  `"aic"`. `"score"` and `"wald"` both apply SAS-style p-value
  thresholds (`slentry` / `slstay`) but score entry candidates
  differently, and can therefore select different variable sets;
  `"score"` reproduces C/SAS HAZARD and needs no per-candidate refit.
  `"aic"` adds or drops whenever it lowers the AIC. See the **Selection
  direction and criterion** section.

- slentry:

  Entry p-value threshold for the score / Wald criteria. Default `0.30`
  matches SAS `SLENTRY`.

- slstay:

  Retention p-value threshold for the score / Wald criteria. Default
  `0.20` matches SAS `SLSTAY`.

- max_steps:

  Hard cap on total accepted actions. Emits a
  [`warning()`](https://rdrr.io/r/base/warning.html) if hit. Default
  `50`.

- max_move:

  Per-variable oscillation cap. When a variable has entered + exited
  more than `max_move` times it is frozen for the remainder of the run.
  Default `4`.

- force_in:

  Character vector of variables that must remain in the model. Such
  variables are still scored and reported in the selection trace, but
  are never dropped.

- force_out:

  Character vector of variables that may never be considered as
  candidates.

- trace:

  Logical; print step-by-step progress to the console. Default `TRUE`.

- ...:

  Unused.

- x:

  An `hzr_stepwise` object.

- object:

  An `hzr_stepwise` object.

## Value

An object of class `c("hzr_stepwise", "hazard")` – the final fit
augmented with:

- `steps`:

  Data frame with one row per accepted / frozen action; see Details.

- `scope`:

  Record of the candidate scope, plus `force_in`, `force_out`, and the
  frozen set.

- `criteria`:

  Named list of the threshold / direction settings actually applied,
  plus, under `criterion = "score"`, `n_uncomputable_scores` (how many
  candidate scores were `NA`), `uncomputable_reasons` (a named integer
  vector of *why*) and `stopped_uncomputable`. Read
  `uncomputable_reasons` before treating an unscored candidate as a bad
  one: `information_indefinite` marks candidates whose effect is too
  large for the score test's approximation at zero, which are typically
  the strongest variables on offer rather than degenerate ones. Those
  are now refit and tested by Wald automatically, counted in
  `n_wald_fallbacks`; a candidate still reaches `uncomputable_reasons`
  only when that refit itself fails, or when the cause is one no refit
  can rescue — which is every cause except `information_indefinite` and
  `coefficient_diverging`, the two a refit exists to rescue. Read
  `uncomputable_reasons` for which one it was in any given run. For
  every criterion it also carries `refit_failures` (the `"var"` /
  `"var@phase"` tokens of candidate moves whose refit errored or failed
  to converge), `n_refit_failures`, and `stopped_refit_failed` — `TRUE`
  when the run ended on an iteration in which refits failed, which is a
  screen that could not test its candidates rather than one that tested
  them and liked none. Check it before reading a zero-row `steps` as an
  honest null result.

- `trace_msg`:

  Character vector of the trace lines, captured regardless of the
  `trace` flag.

- `elapsed`:

  `difftime` from start to finish.

- `final_call`:

  The call that produced this result.

`print.hzr_stepwise` returns `x` invisibly.

`summary.hzr_stepwise` returns a `summary.hzr_stepwise` object (extends
`summary.hazard`) with `$stepwise_steps` and `$stepwise_trace` appended.

`print.summary.hzr_stepwise` returns `x` invisibly.

`as.data.frame.hzr_stepwise` returns the `$steps` data frame.

## Details

The `steps` data frame has columns:

- `step_num`:

  Integer sequence starting at 1.

- `action`:

  `"enter"`, `"drop"`, or `"frozen"`.

- `variable`:

  Variable affected.

- `phase`:

  Phase name (multiphase) or `NA_character_`.

- `criterion`:

  The criterion actually applied to this step — `"score"`, `"wald"`, or
  `"aic"`. Under `criterion = "score"` the drop rows read `"wald"`,
  because score is entry-only.

- `score`:

  Winning score used for the decision.

- `stat`, `df`:

  Test statistic and degrees of freedom.

- `stat_type`:

  What `stat` is on this row, and so which reference distribution
  recomputes its p-value: `"score_q"` (chi-square on `df`), `"wald_z"`
  (standard normal) or `"wald_chisq"` (chi-square on `df`). `df` alone
  does not distinguish them — a scalar Wald is reported as a *z*, not as
  its square, so it and a score Q are both recorded at `df = 1` while
  calling for different distributions. It also identifies the rows the
  Wald fallback rescued, but only among *entry* rows: under
  `criterion = "score"` those are the rows with
  `action == "enter" & stat_type == "wald_z"`. Drop rows are always
  Wald-tested under that criterion — removal follows SAS and is tested
  on the current model's Wald p-value — so they read `"wald_z"` whether
  or not the fallback ever fired. See `$criteria$n_wald_fallbacks` for
  the count.

- `p_value`, `delta_aic`:

  Always populated when computable, regardless of the active criterion.

- `logLik`, `aic`, `n_coef`:

  Goodness-of-fit diagnostics of the model *after* this step.

## Selection direction and criterion

Two arguments shape the search. `direction` decides which moves are
allowed at each step; `criterion` decides how a candidate move is scored
and whether it is accepted.

- `direction = "forward"`:

  Start from the base model and only *add* variables — the best eligible
  candidate enters each step until none clears the entry rule. Variables
  never leave once in.

- `direction = "backward"`:

  Start from the full candidate model and only *drop* variables — the
  weakest term leaves each step until all survivors clear the retention
  rule.

- `direction = "both"` (default):

  Two-way stepwise: after each entry, already-selected variables are
  re-tested and may be dropped. This is the SAS `SELECTION = STEPWISE`
  strategy. `max_move` caps how often a single variable may oscillate
  before it is frozen.

&nbsp;

- `criterion = "score"` (default):

  Accept moves on SAS-style significance thresholds, using the score (Q)
  statistic of the candidate coefficient — this reproduces C/SAS
  HAZARD's `SELECTION` statistic. Q is evaluated at the *current*
  model's MLE with the candidate's coefficient pinned at zero, so **no
  candidate refit is needed**: the reduced-model information is inverted
  once per step and reused across every candidate. Only the winner is
  refit. A candidate enters if its p-value is below `slentry`.

  Score is an *entry* criterion; the drop path never refit per candidate
  in the first place, so removals are tested on the current model's Wald
  p-value against `slstay`, as SAS does.

  For single-distribution fits, the score criterion computes the
  observed information numerically via the suggested numDeriv package
  and errors with a clear message if it is not installed; a multiphase
  fit uses the analytic Hessian instead and does not need it.

  Following SAS, the variance used during *selection* is approximate:
  shaping-parameter covariances are ignored. This affects selection only
  — final-model standard errors are unchanged and still come from the
  full Hessian. Candidates must be single-column numeric main-effect
  terms; a factor is rejected with an error rather than skipped.

- `criterion = "wald"`:

  Accept moves on SAS-style significance thresholds, using the Wald
  \\\chi^2\\ of the affected coefficient(s): a candidate enters if its
  p-value is below `slentry`, and a term is dropped if its p-value rises
  above `slstay`. Entry candidates are scored from a refit that adds the
  candidate (so its new coefficient can be tested); drop candidates are
  scored from the *current* model's Wald p-values without a
  per-candidate refit, and a single refit is run only after a drop is
  chosen. This was the default before version 1.2.0. It differs
  algorithmically from C/SAS HAZARD, so the two criteria can take
  different step paths — and select different variable sets — even when
  they converge to a similar final model.

- `criterion = "aic"`:

  Accept any move with \\\Delta\mathrm{AIC} \< 0\\ (a strictly better
  penalised fit), ignoring `slentry` / `slstay`. Entry candidates use
  the actual \\\Delta\mathrm{AIC}\\ from the candidate refit; drop
  candidates use a Wald-to-likelihood-ratio approximation,
  \\\Delta\mathrm{AIC} \approx W - 2\\\mathrm{df}\\, computed from the
  current model without a per-candidate refit (the chosen drop is refit
  afterwards). Use this for a non-significance-based,
  information-criterion search.

## See also

[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
for the base model and
[`hzr_phase()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase.md)
for multiphase scopes;
[`stepwise_trace()`](https://ehrlinger.github.io/TemporalHazard/reference/stepwise_trace.md)
to retrieve the captured selection log.

## Examples

``` r
data(avc)
avc <- na.omit(avc)
base <- hazard(survival::Surv(int_dead, dead) ~ age,
               data = avc, dist = "weibull", fit = TRUE,
               theta = c(mu = 0.01, nu = 0.5, 0))
# \donttest{
sw <- hzr_stepwise(base, scope = ~ age + mal,
                   data = avc, direction = "forward",
                   control = list(n_starts = 1))
#> Stepwise selection (direction = forward, criterion = score, slentry = 0.30, slstay = 0.20)
#> 
#> Step 1: ENTER  mal   (p = 0.001)
#> (no further action after 1 step)
#> 
#> Final model: 2 covariates, logLik = -223.55, AIC = 455.09
print(sw)
#> Stepwise selection (direction = forward, criterion = score, slentry = 0.30, slstay = 0.20)
#> 
#> Step 1: ENTER  mal   (p = 0.001)
#> (no further action after 1 step)
#> 
#> Final model: 2 covariates, logLik = -223.55, AIC = 455.09
# }
```
