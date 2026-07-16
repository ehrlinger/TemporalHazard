# `criterion = "score"` — score-test stepwise selection — design

Date: 2026-07-16
Status: approved (design), not yet implemented

Development-only document: `inst/dev/` is `.Rbuildignore`d and does not ship.

## Purpose

Add a score-test (`Q`-statistic) criterion to `hzr_stepwise()` for multiphase
fits, reproducing what SAS/C HAZARD's `SELECTION` statement actually computes,
and making SAS-scale variable screens computationally feasible.

## Why

Two independent problems, one cause.

**Speed.** Profiling `hzr_stepwise()` over a 92-variable candidate pool
(2026-07-16) put **99.77% of the time inside `.hzr_refit_with_scope()`**, 94% of
that in `stats::optim`. `R/stepwise-step.R`'s forward step performs **one full
model refit per candidate per step** — 92 variables × 2 phases = 184
optimizations per step. Measured **~25 min per bootstrap replicate**, so a
1000-resample screen matching SAS's `resampl=1000` is **~17 days**.

**Correctness of comparison.** SAS does not refit. Its listing header reads
*"Summary of Variables not in the Model (Based on Q Statistics)"* — entry
candidates are scored with a score test evaluated at the current estimates.
R and SAS therefore run **different selection algorithms**, take different step
paths, and `tests/testthat/test-stepwise-parity.R` has to settle for asserting
a "competitive log-likelihood" rather than comparing steps.

One change fixes both: score the entry candidates instead of refitting them.

## Scope

* **Multiphase only.** `dist = "multiphase"` is the case with the demonstrated
  need and the target of the SAS parity work. Other distributions keep
  `criterion = "wald"`; requesting `"score"` for them is an error with a clear
  message.
* **Entry path only.** The drop path already scores from the current model's
  Wald p-values with no per-candidate refit (`R/stepwise-step.R`), so it is
  already cheap and is left untouched.

## The statistic

For candidate covariate `v` on phase `j`, with the current model at its MLE
`θ̂` and the candidate's coefficient pinned at `β = 0`:

* **Score:** `U_β = ∂logL/∂β` evaluated at `(θ̂, β = 0)` — one call to the
  existing `.hzr_gradient_multiphase()` on the expanded parameter vector.
* **Variance:** `V_β`, the candidate's information block, adjusted for the
  other `μ`/`β` parameters **but not the shape parameters** (see below).
* **Statistic:** `Q = U_β² / V_β`, compared against `χ²(1)`. The entering
  candidate is the one with the largest `Q` whose p-value beats `slentry`.

Because `θ̂` is unchanged across candidates within a step, the nuisance block is
computed **once per step** and reused for all candidates. Each candidate then
costs one gradient evaluation plus a small rank-one update instead of a full
`optim()` run.

### SAS's approximation is deliberate, and we match it

The captured SAS listing states:

```
* During variable selection steps, variances are approximate because
  shaping parameter covariances are ignored.
```

So SAS's `Q` is **not** the textbook efficient score statistic: it ignores the
covariance between the candidate's `β` and the shape parameters. Computing the
statistically-preferable efficient score would produce a different number and
**could not match SAS**, which would leave the parity fixture unable to validate
this feature at all.

Decision (owner, 2026-07-16): **match SAS's approximation.** This package exists
to reproduce HAZARD; matching is what makes per-step `Q` checkable against the
reference. R therefore ships a knowingly approximate statistic during selection,
documented as such in `?hzr_stepwise`. Note this affects *selection* only —
final-model standard errors are unchanged and still use the full Hessian.

**Risk:** the approximation is inferred from that one line of listing text. If
the reading is wrong, R's `Q` will not match the fixture. That is the point of
making the fixture an acceptance gate rather than a smoke test — a wrong
reading fails loudly and early, before the default flips.

## Architecture

| File | Change |
|------|--------|
| `R/score-test.R` *(new)* | `.hzr_score_q_multiphase()` — computes `Q`, `df`, `p_value` for one candidate. One responsibility; unit-testable in isolation. |
| `R/candidate-score.R` | `mode = "entry"` accepts a precomputed score result instead of requiring a fitted `candidate` object. |
| `R/stepwise-step.R` | forward step dispatches to the score path when `criterion = "score"`; the existing refit loop remains for `"wald"`. |
| `R/stepwise.R` | `criterion = c("score", "wald", "aic")` for multiphase — score first, i.e. the new default. Non-multiphase keeps `c("wald", "aic")` and errors on `"score"`. |

### Interfaces

* `.hzr_score_q_multiphase(current, var, phase, data, nuisance = NULL)`
  → `list(stat, df, p_value)`.
* `.hzr_score_nuisance_multiphase(current)` → the per-step reusable block,
  computed once and passed to every candidate in that step.

## Correctness

**Acceptance gate: SAS's own `Q` values, in CI.**

`inst/fixtures/stepwise-avc-forward-wald.rds` is already committed and built on
the **bundled `avc` data** — no PHI, no secure volume, so it runs on every PR.
It carries SAS's per-step trace: 8 steps (`enter` and `drop`), `stat` / `df` /
`p_value`, phases `early`/`constant`, 5 candidates, `SLE = 0.30`, `SLS = 0.20`,
plus the final model's 8 coefficients and log-likelihood.

Its `stat` column **is** SAS's `Q`. (The `avc-forward-wald` filename and
`meta$criterion = "wald"` are misnomers — HAZARD's `SELECTION` always uses `Q`;
"wald" described the R side being compared. Correcting that metadata is part of
this work.)

Assertions:

1. **Per-step `Q` matches SAS** for the `enter` steps, to tolerance. This is the
   gate: if it fails, the statistic is wrong and the feature does not land.
   Do not widen the tolerance to pass — a mismatch is a finding.
2. `test-stepwise-parity.R` is upgraded from "competitive log-likelihood" to a
   step-by-step comparison, which is what it was always reaching for.
3. A speed assertion, so the 17-day regression cannot silently return.

## Breaking change

`criterion` defaults to `"score"` for multiphase. `hzr_stepwise()` has been on
CRAN since 0.9.8, so **re-running an existing multiphase analysis can select a
different variable set with no code change**. That is a real hazard for a
package people publish from.

Decision (owner, 2026-07-16): ship it as the default anyway. Wald-with-refit was
a deviation from the C/SAS reference this package exists to reproduce; matching
HAZARD is a correction, not a preference. `criterion = "wald"` restores the
previous behaviour exactly.

`NEWS.md` gains a **Breaking changes** entry at the top of 1.2.0 stating the
symptom, the reason, and the opt-out.

**Version digits are not touched.** Whether this forces a minor bump is the
owner's decision under the project's versioning rule.

## Out of scope

* Score criterion for the four single-distribution families. Each has an
  analytic gradient, so it is tractable, but only multiphase has a demonstrated
  need.
* The efficient (non-approximate) score as an option.
* Changing the drop path, which already avoids per-candidate refits.
