# Pre-CRAN Parity Test Inventory

**Source:** `~/Documents/GitHub/hazard/examples/` (18 paired `.sas`/`.lst` fixtures)
**Capture vintage:** v4.4.6 macOS (refreshed 2026-05-12)
**Status:** Complete. Captures finalized 2026-05-12; all 54 parity expectations passing.

## Scope summary

| Group | Count | Disposition |
|---|---|---|
| Hazard model fits (in scope for `hazard()` parity) | 8 fixtures, 10 fits | Primary parity targets |
| Utility-only fixtures (Kaplan, Nelson, deciles, predict, bootstrap, stratified) | 7 fixtures | Map to `hzr_*` helpers; secondary parity |
| Out of scope (logistic / gompit / cox) | 3 fixtures (`lg.*`) | Excluded — not parametric hazard |

## Fixture inventory

### A. Primary parity targets — PROC HAZARD fits

Numerical reference values: final log-likelihood, MLE table (Estimate / Std error / Z / Prob>|Z|), natural-scale parameter table, asymptotic vcov + correlation matrix. All anchored on the strings `Log likelihood =`, `Parameter Estimate Summary`, `Estimates for Model Parameters`, `Asymptotic Variance-Covariance Matrix`, `Asymptotic Correlation Matrix`.

| Fixture | Dataset | Time/Event | Phases | Distribution | Covariates | CoE | Fits in `.lst` | Notes |
|---|---|---|---|---|---|---|---|---|
| `hz.deadp.KUL` | CABG (n=5880, 545 ev) | INT_DEAD / DEAD | E+C+L, all fixed shapes | E: M=Nu=1, L: Weibull(γ=3) | none | NOCONSERVE | 1 | Simplest in-scope fit — 3 free params (E0, C0, L0). LL=−3740.52. **Use as parser anchor.** |
| `hz.death.AVC` | AVCS (n=310, 70 ev) | INT_DEAD / DEAD | E+C, fixed M=1, FIXTHALF FIXNU | E + C | none (covariates separate run) | CONSERVE | 1 | LL=−210.501. AVCS dataset — same as `data/avc` in temporal_hazard. |
| `hz.te123.OMC` | PRIMISOL | INT_TE / TE | E+L, Weibull, NU=−1.74 | none | CONSERVE | 2 | Two PROC HAZARD calls (initial + refined). Negative NU exercises `hzr_decompos` sign dispatch. |
| `hz.tm123.OMC` | PRIMISOL | INT_TE / TE | E+L, Weibull | none | CONSERVE | 1 | TM endpoint (vs TE in te123); paired competing-risks pattern. |
| `hm.deadp.VALVES` | VALVES | INT_DEAD / DEAD | E+C, FIXNU FIXTHALF FIXM, M=−2.156 | none | (NOCOV NOCOR) | 1 | Negative-M exercise. Known divergence vs v4.3.0 from XPORT byte-swap fix (`tests/corpus/FINDINGS.md §2a`). |
| `hm.death.patient` | VALVES | INT_DEAD / DEAD | E+C, FIXM | covariate fit (per-patient) | CONSERVE | 1 | Patient-level (vs deadp = aggregated). |
| `hm.death.AVC.deciles` | AVCS | INT_DEAD / DEAD | E+C, FIXM | EARLY: AGE, COM_IV, MAL, OPMOS, OP_AGE, STATUS | CONSERVE | 1 | Multivariable EARLY-phase fit. Followed by `%DECILES` → maps to `hzr_deciles()` for secondary parity. |
| `hm.death.AVC` | AVCS | INT_DEAD / DEAD | E+C, FIXM | varies | CONSERVE | 2 | **Fit 1 is stepwise** (`SELECTION SLE=0.2 SLS=0.1` over 6 EARLY + 6 CONSTANT candidates with `/S`, `/I`, `/E` modifiers) — route to `hzr_stepwise()` not `hazard()`. **Fit 2** uses pre-fit covariates with shaping modifiers (`/S`, `/I`, `/E`) — *not currently supported in `hazard()` R API*. See gap §1 below. |

**Free-parameter count summary** (informs vcov dim):
- 3 free: `hz.deadp.KUL`
- 2–3 free: `hz.death.AVC`, `hz.te123.OMC`, `hz.tm123.OMC`, `hm.deadp.VALVES`
- 8+ free (covariates): `hm.death.AVC.deciles`, `hm.death.patient`, `hm.death.AVC` fit 2

### B. Utility-function targets — auxiliary macros

| Fixture | Macros | Maps to | Notes |
|---|---|---|---|
| `ac.death.AVC` | `%KAPLAN`, `%NELSONT` | `hzr_kaplan()`, `hzr_nelson()` | Actuarial Kaplan-Meier + Nelson-Aalen. No hazard fit — pure non-parametric. |
| `hp.death.AVC` | `%KAPLAN` + `%HAZPRED` + `%PLOT` | `hzr_kaplan()` + `predict.hazard()` | Hazard predictions overlaid on KM. Needs upstream `OUTHAZ` fit + predict call. |
| `hp.death.AVC.hm1` | `%HAZPRED` only | `predict.hazard()` | Pure prediction from saved estimates. |
| `hp.death.AVC.hm2` | `%HAZPRED` only | `predict.hazard()` | Variant of hm1. |
| `hs.death.AVC.hm1` | `%KAPLAN` + 2x `%HAZPRED` | `hzr_kaplan()` + stratified `predict.hazard()` | Stratified predictions. |
| `bs.death.AVC` | `%HAZBOOT(RESAMPL=5, SLE=0.12, SLS=0.1)` | `hzr_bootstrap()` + `hzr_stepwise()` | Bootstrap with stepwise selection. R=5 — small for testing; output will be variance summary, not deterministic estimates. |
| `bs.death.AVC.summary` | PROC UNIVARIATE/MEANS/TRANSPOSE on bootstrap output | (none) | SAS-side post-processing only. Defer — could later validate against `hzr_bootstrap()$summary`. |

### C. Out of scope

| Fixture | Macros | Why excluded |
|---|---|---|
| `lg.dead30d.AVC` | `%LOGIT`, `%LOGITGR` (GOMPIT=1) | Logistic / gompit regression — not parametric hazard. |
| `lg.deadlate.AVC` | `%LOGIT (COX=1)` | Cox PH — different model family. |
| `lg.death.AVC` | `%LOGIT` | Plain logistic. |

## Known gaps to resolve before parity tests run

1. **Shaping modifiers `/S`, `/I`, `/E`** — used in `hm.death.AVC` fit 2 and the stepwise candidates. Not exposed in `hazard()` R API. Options: (a) extend the API; (b) hand-translate the affected covariates and pre-apply transforms; (c) drop fit 2 from primary parity and note as deferred. Recommendation: (c) for v1.0 CRAN — keep parity targets to fits that exercise the public API as-shipped.

2. **Stepwise fixtures** (`hm.death.AVC` fit 1, `bs.death.AVC` via HAZBOOT) — parity belongs in a separate harness against `hzr_stepwise()` rather than `hazard()`. Use **selected variables** and **final-model LL/MLEs** as the parity quantities, not the iteration trace (which depends on candidate ordering).

3. **CRLF + form-feed encoding** — all `.lst` files have CRLF line endings and `\f` page separators. Parser must run `tr -d '\r' | tr '\f' '\n'` (or R equivalent: `readLines(..., warn = FALSE)` then `gsub("\f", "\n")`) and `grep -a` to avoid binary detection. Verified across all 18 files.

4. **Multiple fits per `.lst`** — `hm.death.AVC` (2 fits) and `hz.te123.OMC` (2 fits) must be split by `Final Results:` anchor. Iteration-trace LLs (`Log likelihood =` appearing before `Final Results:`) belong to the run that *follows*, not precedes.

5. **v4.4.6 capture refresh** — reference values in this inventory are placeholders pending final captures. The mapping table above is version-stable; only the LL/MLE column will update.

## Proposed parity tolerances

| Quantity | Tolerance |
|---|---|
| Final log-likelihood | absolute 1e-2 (matches SAS print precision of 6 sig figs at LL ~ −10³) |
| MLE natural-scale | relative 1e-3 |
| MLE log-scale (E0, C0, L0) | absolute 1e-3 |
| Std error | relative 1e-2 (numerical Hessian differences expected) |
| vcov / correlation | relative 1e-2 |
| Z / p-value | derived, not directly compared |

Tolerances widen for `hm.deadp.VALVES` per documented v4.3.0→v4.4 divergence (FINDINGS.md §2a).

## Prototype status (2026-05-12)

**Parser shipped** at `tests/testthat/helper-sas-parity.R` (≈450 LOC). API:
`.hzr_parse_sas_lst(path) -> list(fits = list(...))`;
`.hzr_derive_primisol(raw_path) -> list(te, te_mod, tm)` (new).
**P1 gap #1 closed 2026-05-12** — cursor-based row routing + form-feed
title-line skip; all 10 in-scope fits now produce symmetric vcov +
correlation matrices.  **P1 gap #4 closed 2026-05-12** — PRIMISOL DATA
step reproduced in `.hzr_derive_primisol()`; raw file discovered at
`~/Documents/GitHub/hazard/examples/data/omc`; n_obs=382 / n_events=44
verified against SAS.

**Parity tests** at `tests/testthat/test-sas-parity.R` — **54/54 PASS**:
- `hz.deadp.KUL` (3-phase all-fixed, 3 free): LL=−3740.52, MUE/MUC/MUL
  natural-scale, SE(E0/C0/L0), 3 off-diagonal vcov entries (14 expectations).
- `hz.death.AVC` (2-phase free-shape Early): LL=−210.501,
  THALF/NU/MUE/MUC natural-scale, SE(E2) [SE(E0) skipped per P2 #11]
  (9 expectations).
- `hz.te123.OMC` fit 1 (left-truncated 2-phase): CoE-applied guard, LL within
  1e-5, THALF/NU/ETA natural-scale (5e-3), MUE (1e-4) and MUL as a ratio to
  SAS (5e-4) (8 expectations).  P1 #6 closed; tolerances re-measured
  2026-09-08.  MUL is asserted as a ratio because at 2.1e-04 a direct
  `expect_equal()` falls into waldo's absolute branch and cannot fail.
- `hz.te123.OMC` fit 2 (modulated renewal + late covariates): LL within 1e-5,
  THALF/NU/MUE/MUL natural-scale (1e-3), NOPREVTE/NOTEE coefficients,
  SE(E2/NOPREVTE/NOTEE) (14 expectations).
- `hz.tm123.OMC` fit 1 (morbidity-weighted 2-phase): CoE-applied guard, LL
  within 1e-5, THALF/ETA natural-scale (5e-3), MUE (1e-3) and MUL as a ratio
  to SAS (1e-3) (7 expectations).  P1 #6 closed; tolerances re-measured
  2026-09-08.
- `hm.deadp.VALVES` null model (Case 2L: m<0, nu=0): LL=−1864.76 exact,
  finite MUE/MUC (3 expectations).  Confirms P1 #5 closed (misdiagnosis).

### Parser validation across all 8 in-scope fixtures

| Fixture | Fits | LL | params | vcov dim | vcov symmetric |
|---|---|---|---|---|---|
| `hz.deadp.KUL` | 1 | ✅ | 3 | 3×3 | ✅ |
| `hz.death.AVC` | 1 | ✅ | 4 | 4×4 | ✅ |
| `hz.te123.OMC` | 2 | ✅ | 5, 6 | 5×5, 6×6 | ✅ |
| `hz.tm123.OMC` | 1 | ✅ | 4 | 4×4 | ✅ |
| `hm.deadp.VALVES` | 1 | ✅ | 7 | NULL (NOCOV) | n/a |
| `hm.death.patient` | 1 | ✅ | 12 | 12×12 | ✅ |
| `hm.death.AVC.deciles` | 1 | ✅ | 13 | 13×13 | ✅ |
| `hm.death.AVC` | 2 | ✅ | 11, 13 | NULL, 13×13 | n/a, ✅ |

All 10 fits across 18 fixtures parse log-likelihood, observation counts,
conserved-events counts, parameter-estimate-summary tables, natural-scale
parameter tables, and (where present) symmetric vcov + correlation
matrices. P1 gap #1 closed 2026-05-12 — form-feed page-break title lines
were being mis-parsed as row labels, advancing the cursor incorrectly;
fix skips lines whose leading identifier is not in `col_names`.

### Gaps surfaced during prototyping

The following are **logged for closure after v4.4.6 captures ship.** None
block the KUL prototype.

**P1 — must-fix before extending parity beyond simple fits:**

1. **Wide-vcov asymmetry (12+ free params).** The duplicate-name routing
   in `.hzr_extract_matrix()` doesn't fully disambiguate parameters that
   appear in both phases (e.g., `STATUS` in Early and Constant of
   `hm.death.AVC.deciles`). Symptom: matrices fill but `M != t(M)`.
   Affects: 3 fixtures (`hm.death.patient`, `hm.death.AVC.deciles`,
   `hm.death.AVC` fit 2). Fix: thread phase context into row routing so
   `Early:STATUS` and `Constant:STATUS` get distinct row indices.

2. **Shaping modifiers `/S`, `/I`, `/E`.** Used in
   `hm.death.AVC.sas` fit 2 covariate spec. Not exposed in R API.
   Options: (a) extend `hazard()` API; (b) hand-translate affected
   covariates; (c) defer fit 2 from primary parity. Recommendation (c)
   for v1.0 CRAN — primary parity uses fits that exercise the public
   API as-shipped.

3. **Stepwise fixtures need a separate harness.** `hm.death.AVC` fit 1
   uses `SELECTION SLE=0.2 SLS=0.1`. Route to `hzr_stepwise()` and
   compare on (selected vars, final LL, final MLEs) — not iteration
   trace, which depends on candidate ordering.

**P2 — nice-to-have:**

4. **`fit$fit$se` is `NA`** when shape params are fixed (KUL case).
   Diagonal SEs are recoverable from `sqrt(diag(fit$fit$vcov))[free]`,
   which the test does. Worth surfacing as a top-level convenience.

5. **`vcov(fit)` errors** with "invalid 'nrow' value" on a fit with
   fixed shapes. Direct access `fit$fit$vcov` works. The S3 method
   needs to handle the embed-with-NA-padding layout. (Probably an
   easy fix in `R/vcov-hazard.R`.)

6. **`logLik()` S3 method missing** — `logLik(fit)` errors with "no
   applicable method." Workaround: `fit$fit$objective`. Trivial to
   add a method that wraps the objective + sets `df` and `nobs`
   attributes for `AIC()` compat.

**P1 (new gaps surfaced 2026-05-12 during extension to additional fixtures):**

**P1 #4 — CLOSED 2026-05-12.** OMC PRIMISOL DATA step reproduced in
`.hzr_derive_primisol()` in `helper-sas-parity.R`.  Raw flat file
(`~/Documents/GitHub/hazard/examples/data/omc`, 339 rows × 144 chars)
parsed via fixed-width column positions from the SAS INPUT statement.
Derivation: filter `TE1=1 & TE1DJUL=NA → DELETE`; replacement censoring
via ordered RP3→RP2→RP1 overwrites; multi-row expansion per patient
(1 censored + up to 3 event rows); zero-duration row for study 173
(died same day as TE) dropped to match SAS LCENSOR auto-exclusion.
Final n_obs=382, n_events=44 verified.

**P1 #6 — CLOSED 2026-06-05 (PR #65); diagnosis above was wrong.** The
gap was never Lagrange-vs-closed-form: the C binary's `consrv` also
fixes one mu, so R already matched it.  The cause was left truncation.
Under `LCENSOR STARTTME` the CoE constraint summed `CF(T)` from zero and
omitted the entry-time term, so the conserved phase absorbed a spurious
`Sum CF(ST)` — biasing its intercept (MUE 0.01601 vs SAS 0.01740) and
offsetting the LL by ~0.14 while leaving the shapes correct.  Fixed by
threading `time_lower` into `.hzr_conserve_events()` and
`.hzr_select_fixmu_phase()`.  Since then `hz.te123.OMC` fit 1 agrees with
SAS to 1.3e-04 in LL and 1.1e-05 relative in MUE, and `hz.tm123.OMC` to
4.6e-04 and 9.5e-05.  Tolerances tightened accordingly 2026-09-08; the
old ones (LL within 0.2 / 0.5, MUE not asserted) could not fail.

**P1 #5 — CLOSED 2026-05-12 (misdiagnosis).** The alleged ~17-unit LL
gap at `nu = 0`, `m < 0` was a spurious comparison: the "SAS LL" used
was the final stepwise covariate model LL (−1804.78, 7 free params),
while the "R LL" was the null model (no covariates).  Investigation
confirmed `hzr_decompos()` Case 2L (`m < 0`, `nu = 0`) is mathematically
correct: G(t_half) = 0.5 by construction, and R's null model LL matches
SAS's Initial Summary LL exactly (both −1864.76).  The covariate model
discrepancy is a multiple-local-optima issue — R's multi-start optimizer
finds a better basin (LL ≈ −1786) than SAS's stepwise-initialized BFGS
(LL = −1804.78) on the same likelihood function.  This is expected
behavior, not a bug.  The final covariate / stepwise model parity belongs
to gap #3 (stepwise harness).

**P2 #11 — CoE-projected vcov mismatch on free `mu` SEs.** When CoE
`conserve = TRUE` is on, SAS reports asymptotic vcov projected onto the
conservation manifold (constraint absorbs one DoF, increasing variance
of remaining `mu` params); R returns the raw inverse-Hessian on the
free-parameter submatrix. Symptom: on `hz.death.AVC`, R SE(log_mu_early)
= 0.059 vs SAS 0.133 (factor of ~2.25); other shape SEs match within
5e-3. Affects all CoE-enabled fits' free-mu SE comparisons. Workaround
in test: assert log_t_half SE (E2) but skip log_mu_early SE (E0).
Closure: add delta-method projection in
`R/vcov-hazard.R` or surface SAS-style CoE-vcov via a helper.

**P3 — parser tightening (cosmetic):**

7. **Field-width parsing in param-summary table** initially missed
   p-values like `0.0032` because the regex pattern `[0-9]+` didn't
   include `.`. Fixed in prototype; flagged here as a regression-test
   target.

8. **Character class collision** `[0-9.Ee]` matches parameter names
   starting with `E` (E0, E2, E3). Fixed in prototype.

9. **Closure `<<-` for subset-assign** doesn't propagate `M[i,j] <- x`
   to parent scope; collected writes into a `pending` list and flushed
   after the loop. Fixed.

10. **Form-feed + CRLF preprocessing** mandatory before any grep on
    `.lst` files (binary detection silently zeroes `grep -c` without
    `-a`). Documented in helper header.

### Files added this session

- `tests/testthat/helper-sas-parity.R` — `.lst` parser (≈230 LOC)
- `tests/testthat/test-sas-parity.R` — KUL parity test (14 expectations, all pass)

### Files to add when extending

For each new fixture, append a `test_that()` block to `test-sas-parity.R`.
Pattern: parse reference, build matching `hazard()` call, assert LL +
free-param MLEs + free-param SE + off-diagonal vcov within tolerance.

## Next step (after v4.4.6 captures ship)

1. Re-parse all 18 `.lst` files with the existing helper — no parser changes
   needed since format is version-stable.
2. Close P1 gap #1 (wide-vcov asymmetry) so the wider fits can validate
   off-diagonals.
3. Decide on shaping-modifier scope (P1 #2): extend API or defer.
4. Add `test_that()` blocks for the remaining 7 simple-fit fixtures
   (`hz.death.AVC`, both `hz.te123.OMC` fits, `hz.tm123.OMC`,
   `hm.deadp.VALVES`, `hm.death.patient`).
5. Build the stepwise harness (P1 #3) for `hm.death.AVC` fit 1.
6. Optional polish: add `logLik.hazard` S3 method and fix `vcov.hazard` to
   handle fixed-shape padding (P2 #5–6) — these aren't blockers, the test
   uses `fit$fit$objective` / `fit$fit$vcov` directly.
