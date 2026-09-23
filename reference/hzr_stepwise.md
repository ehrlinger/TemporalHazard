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
  already in the model for every phase. A candidate enters by writing
  its name into the formula text, so under this default too a column
  whose name reads as a different term enters as that term, with no
  warning: a column literally named `age:mal` enters as the interaction,
  and a column `"age "` beside `age` enters as `age` (#449). For
  single-distribution fits, pass a one-sided formula (`~ age + nyha`) or
  a character vector of names. Each name in a character `scope` is
  looked up, not parsed: a name that is exactly a column of `data` is
  that column, otherwise a name that is exactly a term label as
  [`terms()`](https://rdrr.io/r/stats/terms.html) writes it
  (`` "`_X1`" ``, `"log(age)"`, `"age:mal"`) is that term, and any other
  name is ignored with a warning naming it. The column is looked up
  first, so when `data` has a column literally named `age:mal`,
  `"age:mal"` puts that column in the scope rather than the interaction.
  Resolution decides which variables are in the scope; it does not
  change how a candidate is entered. The refit writes the name as you
  spelled it into the formula text. When that text reads as a different
  term, as `"age:mal"` reads as the interaction and `"age "` as `age`,
  the screen enters that other term with no warning (#449); when it does
  not parse, as a bare `"_X1"` does not, the candidate cannot enter
  (#441, \#438). For multiphase fits, pass a named list of one-sided
  formulas keyed by phase, naming each phase once. `scope` lists what
  may enter; a drop considers every term in the model except `force_in`
  and terms frozen by `max_move` before the iteration began (see the
  **Known limitation (the frozen set)** section). A two-sided formula is
  an error, since its left-hand side would never be a candidate, and so
  is a non-empty `scope` under `direction = "backward"`, which does not
  read it. An empty scope (`~ 1`,
  [`character()`](https://rdrr.io/r/base/character.html), or a list of
  `NULL`s and `~ 1`s) is accepted there.

- data:

  Data frame the base fit was built on. Required for refits.

- direction:

  Search strategy: one of `"both"` (default), `"forward"`, or
  `"backward"`. Controls whether variables may only enter, only leave,
  or both. See the **Selection direction and criterion** section.

- criterion:

  Entry / retention rule: one of `"score"` (default), `"wald"`, or
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
  Default `4`. In a two-way screen (`direction = "both"`), a variable
  frozen on entry can still be dropped in the same iteration; see the
  **Known limitation (the frozen set)** section.

- force_in:

  Character vector of variables that must remain in the model. Such
  variables are still scored and reported in the selection trace, but
  are never dropped. Each name is looked up, not parsed, once, when the
  screen starts: a name that is exactly a column of `data` is that
  column, so the bare `"_X1"` pins the column `_X1` although
  [`terms()`](https://rdrr.io/r/stats/terms.html) labels it `` `_X1` ``,
  and `"TRUE"` pins a column named `TRUE`. Otherwise a name that is
  exactly a term label of the model or `scope` is that term, so
  `` "`_X1`" `` and `"age:mal"` work too. Any other name, `"age "` with
  a trailing space when there is no such column, say, matches nothing
  and is ignored with a warning naming it. The column is looked up
  first: when `data` has a column literally named `age:mal`, `"age:mal"`
  resolves to that column and not to the interaction.

- force_out:

  Character vector of variables that may never be considered as
  candidates. Names are looked up as for `force_in`: a column of `data`
  first, then a term label of the model or `scope`, and a warning for a
  name that is neither.

- trace:

  Logical; print step-by-step progress to the console. Default `TRUE`.

- ...:

  Passed to every candidate refit. Only `control` (e.g.
  `control = list(maxit = 500)`) and an `objective` equal to the base
  fit's are accepted. Any other name is an error: a misspelling such as
  `slentyr` would be stored by
  [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
  without being read, and every other
  [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
  argument (the response, `weights`, `time_windows`, `dist`, `theta`,
  `phases`, `fit` and so on) is set by the refit itself, from the base
  model and `data`, so that each candidate is compared with the model it
  extends. The [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  ignore `...`.

- x:

  An `hzr_stepwise` object.

- object:

  An `hzr_stepwise` object.

## Value

An object of class `c("hzr_stepwise", "hazard")`, the final fit
augmented with:

- `steps`:

  Data frame with one row per accepted / frozen action; see Details.

- `scope`:

  Record of the candidate scope, plus `force_in`, `force_out`, and the
  frozen set. In a two-way screen, `frozen` can name a variable the
  final model does not contain; see the **Known limitation (the frozen
  set)** section. `unresolved` is a list with elements `force_in`,
  `force_out` and `scope`, each the names that matched neither a column
  of `data` nor a term label and were therefore ignored
  ([`character()`](https://rdrr.io/r/base/character.html) when none
  were). The trace, and so
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html), carries a line for
  each non-empty one, and a screen whose character `scope` was emptied
  this way says so where it stops. `candidates`, `force_in` and
  `force_out` here are the arguments as given, so a character
  `candidates` still lists the names that were ignored, as `force_in`
  and `force_out` do; read `unresolved` for which those were (#451).

- `criteria`:

  Named list of the threshold / direction settings actually applied,
  plus `n_uncomputable_scores` (how many candidate scores were `NA`,
  counted once per step: an entry the score test could not score under
  `criterion = "score"`; an entry under `criterion = "wald"`, or a
  removal under any criterion, whose Wald statistic could not be
  computed for want of a variance, `wald_no_variance`; or an entry under
  `criterion = "aic"` whose fit had no finite objective, `nonfinite`),
  `uncomputable_reasons` (a named integer vector of *why*),
  `wald_untested_removals` and `wald_untested_entries` (the `"var"` /
  `"var@phase"` tokens of variables kept in, or left out, on a step
  whose Wald test for them could not be computed; a variable tested at a
  later step is not listed. Entries are listed under
  `criterion = "wald"` only: under `"score"` an entry no test could
  reach is reported by its reason, such as `fallback_no_variance`) and
  `stopped_uncomputable` (`TRUE` when the last iteration had candidates
  for entry or for removal and could test none of them). Read
  `uncomputable_reasons` before treating an unscored candidate as a bad
  one: `information_indefinite` marks candidates whose effect is too
  large for the score test's approximation at zero, which are typically
  the strongest variables on offer rather than degenerate ones.
  Candidates with that cause, or with `coefficient_diverging`, are refit
  and tested by Wald automatically, counted in `n_wald_fallbacks`. A
  candidate still reaches `uncomputable_reasons` when that refit fails,
  or when its cause is any other, which no refit can rescue. Read
  `uncomputable_reasons` for which one it was in any given run. For
  every criterion it also carries `refit_failures` (the `"var"` /
  `"var@phase"` tokens of candidate moves whose refit errored, failed to
  converge, or was refused because the move would not change the model –
  a drop that removes no design column, \#320), `refit_failure_reasons`
  (why each one failed or was refused: the refit's error message, that
  it did not converge, or that the move changes nothing; named by the
  same tokens), `n_refit_failures`, and `stopped_refit_failed` (`TRUE`
  when the run ended on an iteration in which a refit failed or a move
  was refused. A refit failure is a screen that could not test its
  candidates, rather than one that tested them and liked none; a refusal
  is determinate – the move was tested and would have left the model as
  it was. Read `refit_failure_reasons` for which it was). Check it
  before reading a zero-row `steps` as an honest null result.

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

  The criterion actually applied to this step: `"score"`, `"wald"`, or
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
  does not distinguish them; a scalar Wald is reported as a *z*, not as
  its square, so it and a score Q are both recorded at `df = 1` while
  calling for different distributions. It also identifies the rows the
  Wald fallback rescued, but only among *entry* rows: under
  `criterion = "score"` those are the rows with
  `action == "enter" & stat_type == "wald_z"`. Drop rows are always
  Wald-tested under that criterion (removal follows SAS and is tested on
  the current model's Wald p-value), so they read `"wald_z"` whether or
  not the fallback ever fired. See `$criteria$n_wald_fallbacks` for the
  count.

- `p_value`, `delta_aic`:

  Always populated when computable, regardless of the active criterion.

- `logLik`, `aic`, `n_coef`:

  Goodness-of-fit diagnostics of the model *after* this step.

## Selection direction and criterion

Two arguments shape the search. `direction` decides which moves are
allowed at each step; `criterion` decides how a candidate move is scored
and whether it is accepted.

- `direction = "forward"`:

  Start from the base model and only *add* variables; the best eligible
  candidate enters each step until none clears the entry rule. Variables
  never leave once in.

- `direction = "backward"`:

  Start from the base model, which must already hold every candidate,
  and only *drop* variables; the weakest term leaves each step until all
  survivors clear the retention rule. `scope` is not read, so a
  non-empty one is an error: protect terms with `force_in`.

- `direction = "both"` (default):

  Two-way stepwise: on every iteration, whether or not a variable
  entered, every term in the model, the base model's included, is
  re-tested and may be dropped unless it is in `force_in` or was frozen
  by `max_move` before the iteration began; a variable frozen on entry
  can still be dropped in the same iteration (see the **Known limitation
  (the frozen set)** section). `scope` limits what may enter, not what
  may leave. This is the SAS `SELECTION = STEPWISE` strategy. `max_move`
  caps how often a single variable may oscillate before it is frozen.

&nbsp;

- `criterion = "score"` (default):

  Accept moves on SAS-style significance thresholds, using the score (Q)
  statistic of the candidate coefficient; this reproduces C/SAS HAZARD's
  `SELECTION` statistic. Q is evaluated at the *current* model's MLE
  with the candidate's coefficient pinned at zero, so **no candidate
  refit is needed**: the reduced-model information is inverted once per
  step and reused across every candidate. Only the winner is refit. A
  candidate enters if its p-value is below `slentry`.

  Score is an *entry* criterion; the drop path never refit per candidate
  in the first place, so removals are tested on the current model's Wald
  p-value against `slstay`, as SAS does.

  For single-distribution fits, the score criterion computes the
  observed information numerically via the suggested numDeriv package
  and errors with a clear message if it is not installed; a multiphase
  fit uses the analytic Hessian instead and does not need it.

  Following SAS, the variance used during *selection* is approximate:
  shaping-parameter covariances are ignored. This affects selection
  only; final-model standard errors are unchanged and still come from
  the full Hessian. Candidates must be single-column numeric main-effect
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
  different step paths (and select different variable sets) even when
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

## Known limitation (the frozen set)

In a two-way screen (`direction = "both"`), `$scope$frozen` can name a
variable that the final model does not contain. Each two-way iteration
makes a forward step and then a backward step, and the sets of protected
variables are fixed at the start of the iteration. A variable that the
forward step freezes can therefore still be dropped by the backward step
that follows it, so it is reported as frozen while the final model
excludes it. Nothing warns when this happens.

Forward-only and backward-only screens are not affected: a forward-only
screen never makes a backward step to drop the frozen variable, and a
backward-only screen never makes a forward step to freeze it on entry.
In those, `$scope$frozen` and the final model agree.

**When the two disagree, trust the final model and `$steps`.** The final
model is what was selected, and `$steps` records both the `"frozen"` row
and the `"drop"` that followed it. Read `$scope$frozen` only as the list
of variables that reached the `max_move` cap, not as a list of variables
held in the model.

This is a known limitation of this release. The analysis, including why
fixing it changes which variables are selected, is in
<https://github.com/ehrlinger/TemporalHazard/issues/378> and
<https://github.com/ehrlinger/TemporalHazard/issues/379>.

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
#> Warning: 'control' element(s) with no effect on this dist = "weibull" fit, ignored: control$n_starts (it applies only to dist = "multiphase").
#> Stepwise selection (direction = forward, criterion = score, slentry = 0.30, slstay = 0.20)
#> 
#> Step 1: ENTER  mal   (p = 0.001)
#> (no further action after 1 step)
#> 
#> Final model: 2 covariates, logLik = -218.60, AIC = 445.19
print(sw)
#> Stepwise selection (direction = forward, criterion = score, slentry = 0.30, slstay = 0.20)
#> 
#> Step 1: ENTER  mal   (p = 0.001)
#> (no further action after 1 step)
#> 
#> Final model: 2 covariates, logLik = -218.60, AIC = 445.19
# }
```
