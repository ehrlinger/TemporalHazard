# What a Fit Did Not Do — Design

**Status:** Approved design (2026-09-10); pending implementation plan.
**Issue:** #242, follow-up to #197 (closed 2026-09-01).
**Version:** folds into 1.2.10 (no bump; PR #241 folds into the same patch).
**Reference implementation:** `medical-manuscript-copyedit` (Python):
`scripts/run_all.py` (`degraded`, `degraded_causes`),
`references/report-template.md` ("Not done in this run"),
`scripts/validate_report.py` (`format_degraded`, the line rebuilt from the record).

## Problem

#197 gave `fit$fit$weak` three states: a list (ridge found), `NULL` (examined,
clean), `NA` (not examined). `print()` and `summary()` write a note only for
the list and `NA`. A clean fit prints nothing, so "examined, clean" cannot be
told from "the line was never written". Nothing records *why* a step was
skipped, and nothing checks printed output against the object.

The same silence covers every other step a fit can lose. Verified in the code
on 2026-09-10 at `84931a7`:

| Step lost | Where it happens | Recorded on the object today? |
|---|---|---|
| Standard errors | `.hzr_optim_generic()` (`R/optimizer.R`) and `.hzr_safe_solve()` (`R/hessian-invert.R`), four warning sites | State only (`vcov` not a matrix). The cause lives only in the `warning()` |
| Variance of the conserved phase's `log_mu` under CoE | `.hzr_optim_multiphase()`, the `!recomputed` branch (`R/likelihood-multiphase.R`) | No. Warning only; SEs and confidence limits for that phase are understated |
| Weak-direction check | `.hzr_weak_direction()` returns `NA` from five guards | `weak = NA`. The printed note always blames "no usable Hessian", which is wrong for a SAS import |
| Conservation of Events | `.hzr_optim_multiphase()` | Yes: `spec$control$conserve_applied` + `conserve_disabled_reason`. Never printed |
| Fitting | `fit = FALSE`; `hzr_read_outhaz()` | No |

The rule behind this work: "a line that appears only on bad news cannot be
told from a line its author forgot to write" (the reference implementation's
`report-template.md`). And by the third-decay rule in the vault's
`fail-loud-engineering.md`, this is the instance that builds the check rather
than writing another note.

## Decisions (maintainer, 2026-09-10)

1. **The record sits beside `weak`; it does not replace it.** `fit$fit$weak`
   keeps its documented list / `NULL` / `NA` meaning, so `is.list(fit$fit$weak)`
   and `is.na()` keep working. The validator enforces
   `weak` is `NA` ⟺ `"weak_direction_check"` is in `degraded`.
2. **The record lives at the top level:** `fit$degraded`, `fit$degraded_causes`.
   It spans `spec` (CoE) and `fit` (SEs, weak), so neither branch owns it.
3. **Opt-outs are listed, with the cause "not requested".** One rule, no
   judgement calls: anything a default fit would have done and this one didn't
   is listed.
4. **Architecture A:** the causes are carried up as reason fields from where
   they arise, the same pattern as `conserve_disabled_reason`. They are not
   rebuilt from classed warning conditions.
5. **All five steps are in v1.**
6. **The validator `stop()`s** when the record disagrees with the object: that
   is a package bug.
7. **The block replaces two existing notes** that can give the wrong cause: the
   `weak = NA` note, and "standard errors unavailable; the Hessian could not be
   inverted".

## The record

```r
fit$degraded          # character; capability names, in canonical order; character(0) if none
fit$degraded_causes   # named character; identical(names(fit$degraded_causes), fit$degraded)
```

Canonical order (`.hzr_capabilities`):
`fitting`, `standard_errors`, `conserved_phase_variance`,
`weak_direction_check`, `conservation_of_events`.

**The object's state decides which entries appear. The reason fields supply
only the cause.** When an entry is forced by state but no reason reached the
builder, its cause is `"cause not recorded"`. That is printed, not hidden, and
a test asserts that no fixture produces it. So a lost-SE path the reason fields
miss still produces a truthful entry, and the `stop()` fires only when the
builder itself is wrong.

### Entries and causes

| Capability | Listed when (state) | Causes (exact strings) |
|---|---|---|
| `fitting` | the optimizer did not run, or the object was imported from SAS | `"not requested (fit = FALSE)"`; `"imported from SAS output; not fitted in R"` |
| `standard_errors` | `fit$fit$vcov` is not a matrix, **or** it is a matrix in which an estimated (not fixed) parameter has no finite, positive variance (Amendment 7) | `"no finite positive variance for: <names>"` (the matrix case); `"numDeriv not installed and no analytic Hessian"`; `"numDeriv::hessian() failed"`; `"Hessian has non-finite entries"`; `"Hessian not invertible"`; `"model not fitted"`; `"covariance not imported"` |
| `conserved_phase_variance` | CoE applied **and** the full-information recompute failed | `"numDeriv not installed"`; `"Hessian could not be computed"`; `"Hessian has non-finite entries"`; `"Hessian not invertible"` |
| `weak_direction_check` | `fit$fit$weak` is `NA`, exactly | `"model not fitted"`; `"imported from SAS output; no R Hessian"`; `"standard errors unavailable"`; `"Hessian condition number unavailable"`; `"covariance has non-finite entries"`; `"eigendecomposition failed"` |
| `conservation_of_events` | `dist = "multiphase"` and `isFALSE(spec$control$conserve_applied)` | from `conserve_disabled_reason`: `not_requested` → `"not requested (conserve = FALSE)"`; `unsupported_censoring` → `"status outside {0, 1} (left or interval censoring)"`; `single_phase` → `"fewer than two phases"`; `no_events` → `"no events"`; `setup_failed` → `"setup did not complete"` |

Every `NA` return in `.hzr_weak_direction()` maps to exactly one weak cause,
taken from the same code path as the return itself. The function's return
value is unchanged, so the guard cannot drift from its reason.

**Not listed, deliberately.** An ill-conditioned Hessian and a Hessian that is
not positive-definite: the step ran and gave a caveated answer. Those stay as
"Note:" lines. Bootstrap, stepwise and translator records are out of scope.

## Components

**New file `R/degraded.R`** (internal only; the API is the two documented
fields):

- `.hzr_capabilities`: the canonical order.
- `.hzr_degraded_record(...)`: a pure function of the fit state, control,
  `dist`, and whether the object was fitted or imported. Returns
  `list(degraded, degraded_causes)`.
- `.hzr_validate_degraded(object, fitted, imported = FALSE)`: `stop()`s
  unless all of these hold:
  - the record's shape: `degraded` is a character subset of
    `.hzr_capabilities`, in canonical order, with no duplicates;
    `identical(names(degraded_causes), degraded)`; every cause is a non-empty
    string;
  - agreement with state: each "listed when" condition in the table above holds
    in both directions for `fitting`, `standard_errors`, `weak_direction_check`
    and `conservation_of_events`; and `conserved_phase_variance` appears only
    when `conserve_applied` is TRUE.
- `.hzr_format_not_done(degraded, causes)`: the one formatter, used by both
  printers.
- `.hzr_check_not_done_output(lines, object)`: the output validator. It parses
  the printed block and compares it with the **raw record, never through
  `.hzr_format_not_done()`**, because a check built on the printer's formatter
  agrees with a broken printer. It returns a character vector of problems
  (`character(0)` means clean): the heading is missing; "none" was printed over
  a non-empty record; an entry is missing; a cause differs; an entry is not in
  the record.

**Reason fields carried up:**

- `.hzr_safe_solve()` gains a `reason` element on its "standard errors
  unavailable" branches.
- `.hzr_optim_generic()` returns `se_unavailable_reason` (`NA_character_` when
  SEs exist), set at each of its warning sites.
- `.hzr_optim_multiphase()` returns `conserved_variance_reason` on the
  `!recomputed` branch. The warning's current "numDeriv unavailable or the
  Hessian was not invertible" is split into its real causes.

**Stamping:** `hazard()` (`R/hazard_api.R`, where the object is assembled) and
`.hzr_outhaz_to_spec()` (`R/read-outhaz.R`, the SAS import) build the record and run
`.hzr_validate_degraded()` before returning. These are the only two places a
`hazard` object is built; `hzr_bootstrap()`, stepwise and the SAS translator
all refit through `hazard()`. `summary.hazard()` carries both fields.

## Printing

`print.hazard()` and `print.summary.hazard()` **always** print the block,
"none" included:

```
  Not done in this run: none
```
```
  Not done in this run:
    standard_errors: numDeriv not installed and no analytic Hessian
    weak_direction_check: standard errors unavailable
```

One entry per line, never `strwrap()`ped: the lines stay readable, and the
output check never has to rejoin continuation lines. The notes for a
weak-direction *list*, ill-conditioning and non-PD are unchanged.

## Testing

Every test must be able to fail.

- **Builder:** a table-driven test covering every capability × cause row
  above.
- **Validator:** hand-built rejection cases, one per rule: out of order, a
  duplicate, names not matching, an empty cause, and each state disagreement
  (for example, `weak = NA` with no entry).
- **End to end:**
  - a clean Weibull fit prints "none";
  - numDeriv mocked absent, together with a declined analytic Hessian, gives
    `standard_errors` + `weak_direction_check`;
  - `status = 2` in a multiphase fit gives CoE `unsupported_censoring`;
  - `conserve = FALSE` gives "not requested";
  - `fit = FALSE` gives `fitting` + `standard_errors` + `weak_direction_check`;
  - a SAS import fixture read with its covariance gives `fitting` +
    `weak_direction_check`; read without it, it also gives
    `standard_errors` ("covariance not imported").

  Every print/summary test runs `.hzr_check_not_done_output()` and asserts
  `character(0)`. No fixture may produce "cause not recorded".
- **Falsifiability:**
  - mock `.hzr_format_not_done()` so it drops an entry, and assert the output
    check reports it;
  - mock it so it always prints "none", and assert the same.
- **Hand mutation during development**, with the kill results reported in the
  PR: drop the `weak_direction_check` branch from the builder; skip the block
  in `print.hazard()`; swap two causes.
- `test-weak-direction.R:366` asserts the removed `NA` note; it is rewritten to
  assert the block.

## Documentation and release

- `hazard()` `@return`: document `degraded` / `degraded_causes`, and point the
  `weak` paragraph at them.
- `print.summary.hazard()` `@description`: the block replaces the two notes.
- `NEWS.md`: a bullet under 1.2.10. `DESCRIPTION` stays at 1.2.10.
- No new dependency.

## Amendments made while planning (2026-09-10)

Found while tracing the code for the implementation plan
(`inst/dev/DEGRADED-RECORD-PLAN.md`). The tables above already include them.

1. **`fit = TRUE` without `theta` skips fitting silently.** On a
   single-distribution model, `hazard()` runs the optimizer only when
   `fit && !is.null(theta)` (`R/hazard_api.R`, the dispatch). It raises no
   error and no warning. This change leaves that behaviour alone and records
   it, under its own cause rather than "fit = FALSE". Whether it should error
   is a separate decision.

   **Superseded 2026-09-10 by #243**, which makes that call `stop()`. The
   cause was dropped from #244 at the maintainer's request; `fit_ran` is now
   FALSE exactly when `fit = FALSE`.
2. **Objects without a record.** Fits saved with `saveRDS()` by an earlier
   version, and `hazard` objects that tests build by hand, have no
   `degraded`. Printing "none" for them would be the defect itself, so they
   print `Not done in this run: not recorded (object built before this record
   existed)`.
3. **The conserved-variance cause comes from `.hzr_safe_solve()`** when the
   recompute has a Hessian but cannot invert it, so it can also read
   "Hessian has non-finite entries".
4. **The validator is told whether the object was fitted or imported.**
   Nothing in a `hazard` object's state reliably says whether the optimizer
   ran, so `hazard()` tracks it (`fit_ran`) and passes it in.
5. **`hzr_read_outhaz()` needs no doc change.** It returns an `hzr_outhaz`,
   not a `hazard` object; the `hazard` object is built by the internal
   `.hzr_outhaz_to_spec()`.
6. **The weak-direction check's `eigen()` failure exit has no test.** `eigen()`
   is base R and cannot be mocked, and no input is known to make it fail on a
   finite symmetric matrix. Its cause string exists and is untested; the PR
   says so.

## Amendment 7: partial standard-error loss (2026-09-10, after review)

Found by the `r-reviewer` pass over the finished branch; the maintainer chose
to fix the whole class in this PR.

`.hzr_safe_solve()`'s non-positive-variance guard sets the offending row and
column of the covariance to `NA`, but still returns a matrix with
`reason = NA`. So an estimated parameter can finish with no standard error
while `vcov` is a matrix, and "`vcov` is not a matrix" misses it. The block
then prints "none" above a coefficient row whose SE is `NA`, which is the
defect this feature exists to prevent. The sharpest case is Conservation of
Events: a full-information recompute that "succeeds" can mask the conserved
`log_mu`'s own variance, and the recompute counts as a success, so
`conserved_phase_variance` is not listed either.

**Rule.** `standard_errors` is listed when `vcov` is not a matrix, *or* when
some **estimated** parameter has no finite, positive variance on the
diagonal. "Estimated" means not fixed: `fit$fit$fixed_mask` is `TRUE` for a
fixed parameter, and a fixed parameter's `NA` row is by design, not a loss.
Without a `fixed_mask` of the right length, every parameter counts as
estimated. The cause names the parameters:
`"no finite positive variance for: early.log_mu, late.nu"`, using the fit's
parameter names, or `par1`, `par2`, ... when none are available.

One internal helper, `.hzr_params_missing_variance(vcov, fixed_mask,
param_names)`, computes the set. The builder uses it for the entry and the
cause; the validator uses it for the state check. It is state-derived, so it
needs no carried-up reason.

**Conservation of Events follows from the mask.** After a successful
recompute the conserved `log_mu` is estimated (`fixed_mask <- !free_unc`),
so a masked variance there lists `standard_errors`. After a failed recompute
it stays fixed in the mask, and `conserved_phase_variance` records it, as
before. Nothing is counted twice.
