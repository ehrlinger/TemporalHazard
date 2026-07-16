# bh.dead SAS parity harness — design

Date: 2026-07-16
Status: approved (design), not yet implemented

Development-only: `inst/dev/` is `.Rbuildignore`d and does not ship.

## Purpose

Compare `TemporalHazard`'s multiphase fit and `hzr_bootstrap(scope=)` embedded
stepwise screen against the SAS/C HAZARD originals (`%HAZARD`, `%HAZBOOT`) on
the `bh.dead` esophageal analysis, to debug the R implementation and to guard
against regressions.

Two SAS runs are the reference:

* `distributions/hz.dead.sas` — the shape fit (`PROC HAZARD ... noconserve`).
* `analyses/bh.dead.sas` — the bootstrap variable screen
  (`%hazboot(resampl=1000, sle=0.12, sls=0.1)`).

## Findings that motivate this work

Comparing the existing draft qmd (`analyses/bh.dead_r_bootstrap.qmd`, on the
`qhsstudies` volume) against the SAS sources and outputs surfaced three
fidelity defects. All three must be fixed for any comparison to be meaningful.

### 1. The draft screens four variables SAS never considered

`HIST_AD`, `HIST_SQ`, `GRADE`, `RESM1`, `LVI` return zero rows in
`bh.dead.lst`. They are commented out in `bh.dead.sas`'s `%hmmodel` macro via a
nested-comment trap: SAS `/* */` blocks do not nest, so

```sas
/*hist_ad, hist_sq,
grade,
/* tumloc */
```

terminates at the **first** `*/`, silently dropping `hist_ad`, `hist_sq`,
`grade` and `tumloc`. The same pattern drops `resm1` and `lvi`. The draft qmd's
scope includes `hist_ad + grade + resm1 + lvi`.

SAS's real candidate pool is **92 variables** (93 rows in each phase table,
minus the `E0`/`C0` intercept row). The draft has 19.

Selection frequency depends on the candidate pool, so a 19-variable scope
cannot be compared to a 92-variable SAS table at all.

### 2. The draft's `fixed=` matches neither SAS run

From `hz.dead.lst`'s "Estimates for Model Parameters" block (values withheld —
see Data governance below):

| Parameter | `hz.dead` (shape fit) | `bh.dead` (bootstrap base) | draft qmd |
|-----------|----------------------|----------------------------|-----------|
| THALF     | estimated            | free (`fixnu fixm`)        | fixed     |
| NU        | estimated            | fixed at its `parms` value | fixed     |
| M         | fixed at 0           | fixed at 0                 | fixed     |
| MUE       | estimated            | free                       | free      |
| MUC       | estimated            | free                       | free      |

`hz.dead` fixes **only M**. `bh.dead` fixes **nu and m**. The draft's
`fixed = "shapes"` fixes all three.

The draft also starts `t_half` at SAS's *converged* value rather than its
`parms` starting value.

### 3. The draft omits `noconserve`

`hz.dead.sas:56` reads `proc hazard data=built noconserve ...`. R defaults to
`conserve = TRUE`, so the draft runs Conservation of Events where SAS does not.
This is also why the CoE full-information-variance path was being exercised at
all.

### Non-issue (verified)

`bhblt.sas7bdat` is written **after** SAS's scaling (`bh.dead.sas:123-135`
applies `age_ln*10`, `age_e*10`, `age_2/10`, `iv_ind/50`, then writes
`library.bhblt`). `read_sas()` therefore hands R exactly the values SAS
bootstrapped; no transformation mirroring is needed. This is why the draft's mu
estimates already land within ~0.5% of SAS despite defects 2 and 3.

## Constraints

* **No PHI in the repo.** `bhblt.sas7bdat` is patient-level clinical data and
  must never enter `temporal_hazard`.
* **Data governance: nothing study-specific enters the repo either.**
  `temporal_hazard` is a **public** GitHub repository, and `.Rbuildignore` only
  excludes a file from the CRAN tarball — it does nothing about git. The SAS
  `.lst` aggregates are not PHI, but they are results from an apparently
  unpublished CCF cohort study. Publishing them, or hard-coding SAS's estimates
  and selection frequencies into committed test code, would put those results
  on the open internet. Therefore:
  * The built fixture lives **beside the data on the `qhsstudies` volume**, not
    in `inst/fixtures/`.
  * Committed code contains **no study values**. Every expectation is derived
    from the fixture at run time.
  * This design document likewise carries no estimates or frequencies.
* **Must not ship to CRAN.** The test skips unless both the data and the
  fixture are readable, so it never executes on CRAN or CI.
* **SAS's bootstrap is not reproducible.** `bh.dead.sas:338` calls
  `%hazboot(..., seed=-1, ...)`; per the macro documentation, `seed=-1` uses
  time of day. The `PCT` column is one random realization over 1000 resamples.
  Exact selection-frequency parity is impossible in principle.
* **Runtime.** A 92-variable × 1000-resample stepwise bootstrap far exceeds the
  10-minute CRAN check budget. The full run belongs in the qmd, not the suite.

## Architecture

Parsing is the single source of truth, shared by the test and the qmd, so the
two cannot drift apart on what "SAS said".

| File | Location | Public? | Purpose |
|------|----------|---------|---------|
| `parse-bhdead-lst.R` | `inst/dev/bhdead-parity/` | yes (generic code) | Parse the two `.lst` files → write `bhdead.rds` to the volume |
| `bhdead-fixture.R` | `inst/dev/bhdead-parity/` | yes (generic code) | Resolve/load/validate the fixture; no study values |
| `README.md` | `inst/dev/bhdead-parity/` | yes | How to re-capture and where the fixture lives |
| `bhdead.rds` | `qhsstudies` volume | **no** | SAS reference (stays off the repo) |
| `test-bhdead-sas-parity.R` | `tests/testthat/` | yes (generic code) | Parity test; skips unless data **and** fixture present |
| `bh.dead_r_bootstrap.qmd` | `qhsstudies` volume | no | Presentation/debug layer |

All harness code lives under `inst/dev/` (`.Rbuildignore`d, so it does not ship
to CRAN) except the test itself, which must live in `tests/testthat/` to run.
The test is written so that it is inert — a silent skip — for anyone without
the volume mounted.

### Fixture contents (`bhdead.rds`, off-repo)

A named list:

* `shape` — from `hz.dead.lst`: specified and converged THALF, NU, M, MUE, MUC.
* `early` — from `bh.dead.lst` Early Phase table: `name`, `n`, `pct`, `min`,
  `max`, `mean`, `sd` (93 rows incl. the `E0` intercept).
* `constant` — same, Constant Phase (93 rows incl. `C0`).
* `meta` — `resampl`, `sle`, `sls`, `n_obs`, capture date, source paths.

The 92-variable candidate pool is derived from `early$name` / `constant$name`
minus the intercept row — never hand-transcribed.

### Data access / skip gate

Two environment variables, each with a default pointing at the `qhsstudies`
mount:

* `TEMPORALHAZARD_BHBLT` — path to `bhblt.sas7bdat`
* `TEMPORALHAZARD_BHDEAD_FIXTURE` — path to `bhdead.rds`

The test calls `skip_if_not()` when either is unreadable, and
`skip_if_not_installed("haven")`. It therefore runs on the maintainer's machine
and skips silently on CRAN, CI, and every other machine. `haven` moves to
`Suggests`.

This mirrors `test-weights-sas-parity.R`, which skips when its fixture is
absent — extended to gate on the data as well.

## Assertions

**No SAS value is written into committed code.** Every starting value, every
expectation, and every candidate name is read from the fixture at run time.

### Shape fit — deterministic, tight

Refit `hz.dead`'s spec in R — `fixed = "m"`, `conserve = FALSE`, shape and mu
starting values taken from `fixture$shape$specified` — and assert the fitted
THALF, NU, MUE and MUC against `fixture$shape$converged` at ~1e-3 relative.

This is the primary regression guard: it is exactly where the analytic-Hessian
boundary bug and the single-free-parameter `drop=FALSE` bug surfaced.

### Bootstrap screen — stochastic, statistical

Refit `bh.dead`'s spec — `fixed = c("nu", "m")`, `conserve = FALSE`, the scope
derived from the fixture's candidate pool, `slentry = fixture$meta$sle`,
`slstay = fixture$meta$sls` — and assert:

* Spearman rank correlation of R `pct` vs the fixture's `pct` exceeds a floor.
* Of the fixture's top-10 by `pct`, R recovers at least half.
* The intercept rows (`E0` / `C0` ↔ `early.log_mu` / `constant.log_mu`) are
  selected 100% of the time.

**Thresholds are calibrated from the first real run, not invented now.** It is
not yet known whether R and SAS agree well enough for any given threshold to
pass. If they diverge substantially, that is a finding to investigate — not a
number to loosen until the test goes green.

### Staged `n_boot`

* `5` — testthat smoke check: proves the path runs end-to-end, fast.
* `50` — qmd default, for iterating on the comparison.
* `1000` — the real SAS-comparable run (matches `resampl=1000`), driven from
  the qmd with `verbose = TRUE` for the progress bar.

`n_boot` is a Quarto param so escalating is a one-line change.

## qmd rewrite

`analyses/bh.dead_r_bootstrap.qmd` mirrors the SAS step-for-step:

1. **Shape fit** (`hz.dead.sas`) — `fixed = "m"`, `conserve = FALSE`, starting
   values from the fixture. Print alongside SAS's converged values.
2. **Bootstrap base** (`bh.dead.sas` `parms ... fixnu fixm`) —
   `fixed = c("nu", "m")`, `conserve = FALSE`.
3. **Scope** — 92 variables read from the fixture, lowercased.
4. **Screen** — `n_boot` param, `slentry = 0.12`, `slstay = 0.10`,
   `verbose = TRUE`.
5. **Compare** — join R `pct` ↔ SAS `PCT`, sorted by descending absolute
   discrepancy. This table is the actual debugging deliverable.

Guards:

* **Fail loud on missing columns.** Validate all 92 scope variables exist in
  `caa` before fitting and name any that do not, rather than letting
  `hzr_stepwise()` degrade each into a per-candidate warning 92 times.
* **The qmd is presentation only.** All `.lst` parsing lives in the fixture
  script shared with the test.

## Out of scope

* An `avc`-based version of this parity check that could run in CI. Considered
  and deferred: it would require a fresh SAS run against `avc` and would not be
  literally `bh.dead`.
* Fixing `conserve`'s missing `@param control` documentation (only `n_starts`
  is currently documented). Tracked separately.
* The unexplained empty-coefficient-table print behaviour on the maintainer's
  machine. Under investigation separately; the summary object itself is
  verified correct.
