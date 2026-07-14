# `scope=` selection mode for `hzr_bootstrap()`

Date: 2026-07-14
Status: proposed

## Problem

The SAS `%HAZBOOT` macro (called from jobs like `analyses/bh.dead.sas` in
`neo_therapy`) draws bootstrap resamples and runs a **fresh stepwise variable
selection** on each one, reporting a selection *frequency* per candidate
variable (e.g. "STATUS entered 5/5 resamples, OP_AGE 3/5"). This is used to
screen out covariates that enter the model too infrequently to trust, after
first fitting a no-covariate model (`hz.dead.sas`) to establish the hazard
shape.

`hzr_bootstrap()` in this package has no equivalent: it resamples and refits
a single **fixed** formula on every replicate (covariate-stability
bootstrap). This gap is already documented in
`inst/dev/FIXTURE-GAP-LIST.md` (line 51) and in the `bs.death.AVC` test block
in `tests/testthat/test-sas-parity.R` (lines 667-709), which parses the SAS
reference frequencies but asserts nothing about R producing them.

`hzr_stepwise()` already implements SAS-style forward/backward/both Wald
selection (`slentry`/`slstay` ≈ SAS `SLE`/`SLS`). The missing piece is purely
composing it with the bootstrap resampling loop.

## Design

### API

Add one argument, `scope`, to `hzr_bootstrap()`, plus the `hzr_stepwise()`
controls needed to drive it. Signature grows from:

```r
hzr_bootstrap(object, n_boot = 200L, fraction = 1.0, seed = NULL, verbose = FALSE)
```

to:

```r
hzr_bootstrap(object, n_boot = 200L, fraction = 1.0, seed = NULL, verbose = FALSE,
              scope = NULL, direction = c("both", "forward", "backward"),
              criterion = c("wald", "aic"), slentry = 0.30, slstay = 0.20,
              max_steps = 50L, max_move = 4L,
              force_in = character(), force_out = character(), ...)
```

- `scope = NULL` (default): **unchanged behavior.** Every replicate refits
  `object`'s exact formula, as today.
- `scope` supplied (formula, character vector, or named list of formulas for
  multiphase — same shapes `hzr_stepwise()` already accepts): each replicate
  instead runs `hzr_stepwise(object, scope = scope, data = boot_data,
  direction =, criterion =, slentry =, slstay =, max_steps =, max_move =,
  force_in =, force_out =, trace = FALSE, ...)`.

`object` plays the role of `hz.dead.sas`'s output: a base fit with the
hazard shape already estimated and fixed (`hzr_phase(..., fixed =
"shapes")`). No new plumbing is needed to keep the shape fixed across
replicates — `hzr_stepwise()`'s internal refit (`.hzr_refit_with_scope()`)
already rebuilds every candidate model from `current$spec`, which carries
the `fixed` setting forward automatically.

### Accumulation logic (the actual code change)

Today, `param_names` is computed **once**, before the replicate loop, from
`object$fit$theta`, and every replicate's rows are labeled by positional
index into that fixed vector. That assumption breaks under `scope`-mode:
different replicates select different variables, in different numbers and
orders, and multiphase coefficient names are already phase-qualified
(`"early.age"` vs `"constant.age"` — confirmed via
`.hzr_candidate_coef_name()` in `stepwise-step.R`), so no cross-phase name
collision needs handling.

Change: in `scope`-mode, derive parameter names **per replicate** from that
replicate's own `boot_fit$fit$theta` (or `boot_fit$theta`, whichever
`hzr_stepwise()` returns — same field the `hazard` class already exposes),
not from a shared pre-computed vector. `fraction`, `seed`, `verbose`,
resampling of `data`/`weights`, and the success/failure accounting all stay
exactly as they are today; only the body of the per-replicate branch and the
name lookup change.

### Why summary/print need no changes

`summary$pct` is already `100 * n / n_boot` where `n` is the count of
replicates in which that parameter's row appears. In today's fixed-refit
mode every successful replicate contributes every parameter, so `pct` is
always ~100%. In `scope`-mode, a variable simply has no row in a replicate
where it wasn't selected — so `pct` becomes exactly the SAS
selection-frequency number, and `mean`/`sd`/`ci_lower`/`ci_upper` become the
coefficient distribution *conditional on selection* (matching `%sumboot`
semantics). `print.hzr_bootstrap()` needs one added line noting which mode
produced the result, so a 100% `pct` in fixed-refit mode isn't misread as
"always selected."

### Result object

Add `mode` (`"refit"` or `"select"`) to the returned list. In `"select"`
mode, also store `scope` (the candidate scope actually used) for reference,
mirroring `hzr_stepwise()`'s own `$scope` field. `class` stays
`"hzr_bootstrap"` — no new S3 class, since summary/print already generalize.

### Error handling

Per-replicate failures (non-convergence, non-finite objective) are still
caught and counted in `n_failed`, as today. But a **structurally bad**
`scope`/`slentry`/`force_in` argument (e.g. a typo'd variable name not in
`data`) should not silently manifest as "every one of 200 replicates
failed" with no diagnostic. Before entering the loop, validate `scope`-mode
arguments with one call to `hzr_stepwise()` on the *original* (unresampled)
data; let that call's own argument checks raise immediately. Only errors
inside the resampling loop itself (convergence-shaped failures) get
swallowed into `n_failed`.

### Out of scope

- Matching SAS's exact selection *path* or final variable set. That
  divergence (Wald vs. score statistic, approximate vs. full Hessian
  variance) is already documented and accepted for `hzr_stepwise()` itself
  (`hm.death.AVC` test block). This feature inherits that same documented
  gap — it is a regression-guarded capability, not a SAS-parity target.
- Any change to fixed-refit (`scope = NULL`) behavior.
- A new S3 class or print method.

## Testing (in `temporal_hazard`, using only the public `avc` fixture)

1. Regression: `scope = NULL` behavior is byte-for-byte unchanged (existing
   `bs.death.AVC` "R fixed-model bootstrap runs" test keeps passing
   untouched).
2. New test: on a small single-distribution toy fit, `scope` given as a
   one-sided formula selects a subset of variables across replicates;
   `summary$pct` values fall in `[0, 100]`; variables never entering any
   replicate are absent from `summary` (not present with `pct = 0`, since
   `split()` over zero rows can't materialize a group — document this).
3. Multiphase: extend the existing `bs.death.AVC` block
   (`test-sas-parity.R` lines 687-709) so it *exercises* R's new `scope=`
   path (small `n_boot`, e.g. 10-20, real AVC data) instead of only parsing
   the SAS reference frequencies. Assert phase-qualified names
   (`"early.*"` / `"constant.*"`) appear distinctly, `n_success > 0`, and
   `mode == "select"`. Continue to **not** assert R-vs-SAS frequency parity
   (still a documented gap) — this remains a regression guard plus a
   recorded SAS reference, upgraded from "parses SAS output" to "parses SAS
   output AND runs the equivalent R path."
4. Bad-argument test: a `scope` referencing a nonexistent column raises
   immediately (not silently reported as 200 failed replicates).

No PHI-adjacent data is used in any `temporal_hazard` test or fixture.

## Worked example (in `neo_therapy`, not committed to any git repo)

Add an R/quarto script alongside `analyses/bh.dead_r_test.qmd` that mirrors
the `hz.dead.sas` → `bh.dead.sas` hand-off:

1. Fit the shape (`hz.dead`-equivalent): `hazard(Surv(iv_dead, dead) ~ 1,
   data = caa, dist = "multiphase", phases = list(early = hzr_phase("cdf",
   ...), constant = hzr_phase("constant")), fit = TRUE)` — reusing the
   already-fitted values shown in `bh.dead_r_test.qmd` lines 379-391
   (`t_half = 1.666871, nu = -.7462, m = 0`) as the fixed shape, via `fixed
   = "shapes"`.
2. Call `hzr_bootstrap(base_fit, n_boot = 50-100, scope = <trimmed ~15-20
   variable subset drawn from bh.dead.sas's candidate list>, slentry = 0.12,
   slstay = 0.10)`, matching `bh.dead.sas`'s actual `SLE`/`SLS` call values
   (line 326).
3. Print/inspect `summary$pct` per variable as the screening output.

This uses the real `datasets/bhblt.sas7bdat` (already read by
`bh.dead_r_test.qmd` line 305) locally only. Trimmed variable list and low
`n_boot` are a deliberate first pass — scaling to the full ~90-variable list
and `n_boot = 1000` (matching SAS) is a follow-up once the mechanics are
confirmed to work, left for you to run.
