# Correctness Strategy — TemporalHazard ⇄ HAZARD

**Status:** adopted (2026-10-12, John Ehrlinger); proposed 2026-06-05
**Scope:** how the R (`TemporalHazard`) and SAS/C (`hazard`) codebases are used
together to assure correctness, and how to do it *efficiently* as both mature
into maintained peer engines.
**Related:** `hazard/docs/superpowers/specs/2026-04-28-hazard-v5-design.md`
(peer-engine architecture, Parquet contract, V9/V10 gates);
[FIXTURE-GAP-LIST.md](FIXTURE-GAP-LIST.md);
[PRE-CRAN-PARITY-INVENTORY.md](PRE-CRAN-PARITY-INVENTORY.md).

---

## 1. The reframe: parity ≠ correctness

We have been running **direct parity** — R reproduces SAS's printed `.lst`
output to print precision (~5e-6). That was the right tool for **porting and
debugging**: it found real bugs that no single-engine test would have
(multiphase `vcov()` unusable, CoE ignoring left-truncation, `hzr_deciles()`
implementing a different statistic than `%DECILES`).

But exact reproduction is the wrong shape for a **steady-state correctness
regime**, and the marginal fixture already shows it. "Does R match SAS's
printout?" and "Is the answer correct?" are different questions. We should keep
asking the second and stop paying full price for the first.

---

## 2. Why exact parity does not scale

Each new direct-parity fixture costs artisanal effort that is mostly *not* about
correctness:

- **A bespoke parser per output format** — life table, nomogram, vcov matrix,
  decile table, parameter summary each needed their own `.lst` reader.
- **Reconstructing inputs SAS never printed** — e.g. the `hp.death.AVC` nomogram
  times are `n × (12/365.2425)`, not the 3-decimal values printed; the model
  starting values come from the `.sas` `PARMS`, not the `.lst`.
- **Matching SAS *conventions* that are not "more correct"** — print rounding,
  asymmetric confidence-limit transforms, the 1-SD (`CLEVEL=0.68268948`)
  default, internal scale parameterizations.
- **External gating** — capturing reference `.lst` requires SAS/CCF access (the
  Rajes/Blackstone P0 dependency); the loop is slow and not self-service.

And when the engines genuinely **disagree** (conserved-μ standard error;
survival CL transform), exact parity only *flags* the difference — it cannot say
which engine is right. Those disagreements are the high-value signal;
mechanically reproducing the agreeing 95% by hand is not.

**Not in scope of this section: bit-identical oracles for performance work.**
For a performance refactor of code already known to be correct, an oracle
comparing the old code's outputs with the new code's, bit for bit, is the
accepted instrument. It compares R with R, so it needs no parser, no
reconstructed input and no SAS access, and none of the costs above apply. And
because such a refactor claims to change no answer, exact agreement is the
right tolerance, not a print convention. The 1.2.14 speed-up (#631) was gated
this way: an oracle over 1,344 `decompos` cells, whole fits and bootstrap
replicates, with 1-ulp mutants planted to show it could fail.

---

## 3. The reference-authority principle (the decision underneath everything)

As R becomes a maintained peer engine rather than a one-way port, the definition
of "correct" must shift:

> **Correct = satisfies the mathematics, and agrees with the other engine where
> the other engine is right.**

Not "correct = matches SAS." This reframing is what licenses us to stop chasing
SAS's print conventions, and it forces a first-principles decision on each real
disagreement instead of a reflexive "copy SAS." Current open disagreements that
must be *decided*, not just matched:

- **Conserved-phase `log_mu` variance** — SAS reports a delta-method SE for the
  CoE-conserved intercept; R returns `NA`. Pick the right behaviour on theory.
- **Survival confidence-limit transform** — R's `predict(se.fit=TRUE)` CLs do
  not reproduce HAZPRED's even at the same level (~0.02 off); different
  delta-method construction. Decide the canonical transform.

Rule of thumb: a fixture exists to **lock in a decision**, not to chase a digit.

---

## 4. A tiered correctness model

Match the assurance tool to the question. Most assurance should live in the
cheapest tier.

### Tier 1 — Property / invariant tests (R-only, no SAS, every CI run)
The workhorse. Both engines must satisfy mathematics that needs no
cross-reference. These catch the *kinds* of bugs we actually found, run in
seconds, and have zero external dependency:

- **CoE identity:** `Σ events = Σ [H(stop) − H(start)]` at the optimum.
- **Stationarity:** score (analytic gradient) ≈ 0 at the reported MLE.
- **Derivative agreement:** analytic gradient / Hessian == `numDeriv`.
- **Distributional sanity:** `S(t) ∈ [0,1]` and monotone ↓; `H(t)` monotone ↑;
  `S = exp(−H)`; hazard ≥ 0.
- **Invariance:** log-likelihood invariant to reparameterization / time scaling.
- **Self-consistency:** LL recomputed at returned coefficients == reported
  objective (caught the CoE final-adjustment bug).
- **Non-parametric anchors:** KM / Nelson-Aalen == `survival::survfit` (free via
  the existing wrapper).

We already write several of these ad hoc; they should become a **systematic,
dataset-driven suite** — the primary correctness regime.

### Tier 2 — Cross-engine differential testing on a shared structured contract
The *efficient* form of parity, and the operational form of the v5 **V9 gate**:

1. Feed **both** engines the same dataset + model spec.
2. Each writes a **structured** result (the v5 Parquet contract: estimates,
   vcov, LL, predicted S/H on a grid — not a printed listing).
3. **One generic comparator** asserts agreement within a *statistical* tolerance
   (optimizer precision / ~4 significant figures), plus the Tier-1 invariants on
   each output.

This replaces N bespoke `.lst` parsers with a single diff, and runs over **many
datasets** — synthetic across the feature space + real de-identified — rather
than the fixed set of captured `.lst`. Tolerances are defined by *statistical
meaning*, not print precision.

### Tier 3 — A small, frozen golden set (curated, regression anchors)
Keep ~7 representative SAS-captured fixtures — one per structural feature — as
regression anchors against SAS's exact numbers. We have essentially built this
already (see §5). Declare the feature matrix covered and **stop hand-porting the
long tail**.

---

## 5. Feature-coverage matrix

The golden set earns its keep by covering structural features once each, not by
exhaustively porting every example. Current coverage:

| Feature | Anchor fixture | Status |
|---|---|---|
| 3-phase, fixed shapes | `hz.deadp.KUL` | ✅ |
| 2-phase, free-shape early | `hz.death.AVC` | ✅ |
| Both-phase covariates (shared covariate) | `hm.death.patient` | ✅ |
| Covariates + GOF deciles | `hm.death.AVC.deciles` | ✅ |
| `nu=0, m<0` parametric branch | `hm.deadp.VALVES` | ✅ |
| Left-truncation + CoE (counting process) | `hz.te123.OMC` | ✅ |
| Weighted events | `hz.tm123.OMC` | ✅ |
| KM / Nelson-Aalen life tables | `ac.death.AVC` | ✅ |
| HAZPRED survival/hazard prediction | `hp.death.AVC` | ✅ (#69) |
| Interval censoring | — | ⛳ gap (R-tested, no SAS anchor) |
| 3-phase free-shape | — | ⛳ gap |
| 4-phase (R-only; no C reference — C is 3-phase) | — | ⛳ R-vs-R only |

Genuinely worth one anchor each: **interval censoring**, **3-phase free-shape**.
The remaining hand-port candidates buy little new coverage:

- **Predict variants** (`hp.*.hm1/hm2`, `hs.*`) — same prediction engine as
  `hp.death.AVC`, different covariate profiles. Tier-1/Tier-2 cover the engine.
- **Stepwise** (`hz.death.AVC` selection) — blocked on a SAS-materialized `.rds`;
  R already finds *better* optima than SAS's stepwise basin (not a bug). Belongs
  in the dedicated `hzr_stepwise()` harness, not the parity set.
- **Bootstrap** (`bs.death.AVC`) — non-deterministic; the only sensible metric is
  variable-selection frequency + mean LL, which is a Tier-1 property, not a
  point-parity.

---

## 6. What to keep, freeze, and stop

- **Keep & grow:** Tier-1 invariant suite (highest assurance per unit effort).
- **Build:** Tier-2 structured-contract comparator (rides on v5 Parquet; the
  durable cross-engine gate). One comparator, many datasets.
- **Freeze:** the Tier-3 golden anchors as-is; refresh values only on a new
  binary capture, never hand-expand.
- **Stop:** reconstructing exact inputs and matching print/CL/scale conventions
  for *new* fixtures; hand-porting predict variants / stepwise / bootstrap.
- **Decide (don't match):** conserved-μ variance; survival-CL transform.

---

## 7. Efficiency moves, concretely

1. **Define tolerances by statistical meaning**, not print precision — e.g.
   "agree to the optimizer's `reltol`" or "4 significant figures," so a CL
   convention or a rounding artifact never fails a correctness test.
2. **Generate data programmatically** — a synthetic generator spanning the
   feature space (phase counts, shapes, covariates, censoring/truncation,
   weights, sample sizes) feeds both Tier-1 and Tier-2 without SAS capture.
3. **Make the I/O structured** so comparison is generic (v5 Parquet), retiring
   the per-format `.lst` parser tax.
4. **Surface, don't bury, differences** — when the engines disagree, the test
   records the disagreement and points at the §3 decision, rather than loosening
   a tolerance until it passes.

---

## 8. Roadmap / sequencing

1. **Now (no SAS needed):** stand up the Tier-1 invariant suite as a first-class,
   dataset-driven test layer. Pull the ad-hoc invariants already scattered in
   `test-*` into it. *Highest value, zero external dependency.*
2. **Near term:** add the two worthwhile anchors (interval censoring, 3-phase
   free-shape); resolve the two §3 disagreements (conserved-μ variance,
   survival-CL transform) on first principles.
3. **With v5:** build the structured-output contract + the single generic
   comparator (Tier 2 / V9). Wire both engines to emit Parquet results; run the
   comparator over the synthetic + de-identified corpus in CI.
4. **Ongoing:** treat the golden set as frozen; route new "does R match SAS?"
   questions through Tier 2, and new "is R correct?" questions through Tier 1.

---

## 9. Open decisions (for John / Rajes / Blackstone)

- ✅ **Reference authority:** "correct = math + agree-where-right," not
  "match SAS." (§3) Ratified 2026-10-12 by adopting this document.
- **Conserved-μ variance** and **survival-CL transform:** which construction is
  canonical? (§3)
- **Coverage line:** confirm the §5 matrix is "enough" and the long-tail ports
  are explicitly out of scope.
- **Where the contract lives:** confirm the v5 Parquet result schema as the
  Tier-2 interface (vs. an interim R-side structured export to unblock Tier 2
  before v5 lands).
