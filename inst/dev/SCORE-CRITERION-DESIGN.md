# `criterion = "score"` — score-test stepwise selection — design

Date: 2026-07-16
Status: approved (design), not yet implemented

Development-only document: `inst/dev/` is `.Rbuildignore`d and does not ship.

## Purpose

Add a score-test (`Q`-statistic) criterion to `hzr_stepwise()`, reproducing what
SAS/C HAZARD's `SELECTION` statement computes for multiphase fits, and making
SAS-scale variable screens computationally feasible.

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

* **All five distributions** — `exponential`, `weibull`, `loglogistic`,
  `lognormal` and `multiphase`. Each already has an analytic score vector
  (`.hzr_gradient_<dist>()`), so the statistic is tractable for all of them, and
  a uniform `criterion = "score"` avoids an API where the same argument is valid
  for one `dist` and an error for four.
* **Entry path only.** The drop path already scores from the current model's
  Wald p-values with no per-candidate refit (`R/stepwise-step.R:341`), so it is
  already cheap and is left untouched.

### Validation is two-tier, because a reference exists for only one family

`PROC HAZARD` **is** the Blackstone multiphase decomposition. It has no
single-distribution stepwise, so there is no SAS `Q` to reproduce for the other
four families — "match SAS's approximation" has no referent there. The one SAS
`Q` fixture (`stepwise-avc-forward-wald.rds`) is a **two-phase** model despite
its `meta$dist = "weibull"` label, which records the *shaping* family, not
`dist`.

Consequently:

| Family | Gate |
|--------|------|
| `multiphase` | **SAS's per-step `Q`** from the committed fixture. Acceptance gate; a mismatch is a finding, not a tolerance to widen. |
| the other four | **Numeric oracle:** the analytic score matches `numDeriv`'s gradient at the MLE, and `Q` agrees with the Wald chi-square from an actual refit where theory requires it (small effects, `df = 1`). |

The numeric oracle is weaker than a reference implementation, but it is the same
standard this package already holds its analytic Hessians to (see
`tests/testthat/test-*-hessian.R`, which validate against `numDeriv`).

**The shape-covariance approximation is applied uniformly.** For the four
single-distribution families the analogue of "ignore shaping parameter
covariances" is to ignore the covariance between the candidate's `β` and that
family's shape parameter (e.g. `nu` for Weibull). Applying it consistently keeps
one statistic across all five rather than two subtly different ones — but note
this means the four unreferenced families inherit an approximation chosen to
match SAS, not one independently justified for them. Documented in
`?hzr_stepwise`.

## The statistic

For candidate covariate `v` (on phase `j` for multiphase fits; `phase` is `NULL`
otherwise), with the current model at its MLE `θ̂` and the candidate's
coefficient pinned at `β = 0`:

* **Score:** `U_β = ∂logL/∂β` evaluated at `(θ̂, β = 0)` — one call to the
  family's existing analytic score, `.hzr_gradient_<dist>()`, on the expanded
  parameter vector.
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
| `R/score-test.R` *(new)* | `.hzr_score_q()` — computes `Q`, `df`, `p_value` for one candidate, dispatching on `dist` to the family's analytic score. One responsibility; unit-testable in isolation. |
| `R/candidate-score.R` | `mode = "entry"` accepts a precomputed score result instead of requiring a fitted `candidate` object. |
| `R/stepwise-step.R` | forward step dispatches to the score path when `criterion = "score"`; the existing refit loop remains for `"wald"`. |
| `R/stepwise.R` | `criterion = c("score", "wald", "aic")` for every `dist` — score first, i.e. the new default across the board. |

### Interfaces

* `.hzr_score_q(current, var, phase = NULL, data, nuisance = NULL)`
  → `list(stat, df, p_value)`. `phase` is `NULL` for single-distribution fits.
* `.hzr_score_nuisance(current)` → the per-step reusable block, computed once
  and passed to every candidate in that step.

## Correctness

### Tier 1 — multiphase: SAS's own `Q` values, in CI (acceptance gate)

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

### Tier 2 — the other four families: numeric oracle

No SAS reference exists (see Scope). For `exponential`, `weibull`,
`loglogistic` and `lognormal`:

1. **Score vs `numDeriv`:** the analytic `U_β` at `(θ̂, β = 0)` matches
   `numDeriv::grad()` of the log-likelihood, to the tolerance the existing
   Hessian tests use.
2. **Q vs refit Wald:** for a candidate with a small true effect, `Q` agrees
   with the Wald chi-square from an actual `.hzr_refit_with_scope()` — both are
   asymptotically `χ²(1)` and must converge. This uses the existing, trusted
   refit path as the oracle.
3. **Degenerate guards:** a candidate that is constant, collinear with an
   existing term, or all-`NA` yields a non-finite `Q` that is skipped rather
   than selected.

This tier cannot catch an error shared by both the analytic score and the refit
path. That limitation is the price of there being no reference implementation
for these families, and is why only multiphase carries an acceptance gate.

## Breaking change

`criterion` defaults to `"score"` for **every** `dist`. `hzr_stepwise()` has been
on CRAN since 0.9.8, so **re-running any existing stepwise analysis can select a
different variable set with no code change**. That is a real hazard for a package
people publish from, and widening the scope widened this blast radius from
multiphase users to all `hzr_stepwise()` users.

Decision (owner, 2026-07-16): ship it as the default anyway. Wald-with-refit was
a deviation from the C/SAS reference this package exists to reproduce; matching
HAZARD is a correction, not a preference. `criterion = "wald"` restores the
previous behaviour exactly.

`NEWS.md` gains a **Breaking changes** entry at the top of 1.2.0 stating the
symptom, the reason, and the opt-out.

**Version digits are not touched.** Whether this forces a minor bump is the
owner's decision under the project's versioning rule.

## Out of scope

* The efficient (non-approximate) score as an option.
* Capturing SAS fixtures for the single-distribution families — `PROC HAZARD`
  is the multiphase decomposition and has no such models to capture.
* Changing the drop path, which already avoids per-candidate refits.
