# Design: `hzr_repeated_events()` — an R reimplementation of the SAS `%repeat` macro

- **Date:** 2026-09-09
- **Status:** design approved; not yet implemented
- **Branch:** `feat/hzr-repeated-events`, cut from `main` at `f93d00c`

## Purpose

Repeated-events HAZARD jobs in the corpus build their fit input with the SAS macro
`%repeat`. That input was never saved, so those jobs cannot be reproduced today. This
adds a native R implementation so they can be, which is TemporalHazard's remit: SAS
reproduction.

The immediate parity target is
`/Volumes/qhsstudies/cardiac/rhythm/maze/permanent/distributions/hz.ce_cardioversion_repeated.ehb.sas`
— early CDF plus late G3, `LCENSOR IV_START`, `CONSERVE`, 9 free parameters, log
likelihood -267.885 — whose fit input `events` is built by `%repeat(...)` into SAS WORK.

## Source

`~/Documents/macro.library/repeat.sas`, 161 lines, mtime 2003-10-17:

```sas
%macro repeat(in=built, out=events, eventype=eventype,
              iv_event=iv_event, iv_end=iv_end, id=id,
              event=event, event_no=event_no, rcensor=rcensor,
              iv_start=iv_start, iv_seg=iv_seg, renewal=renewal);
```

⚠️ **Macro identity is NOT established.** See Open Items 1.

## What the macro does

`%repeat` is **long in, long out**. Its input is already one-to-many — one row per
candidate event per subject, with a text indicator naming the event of interest. It is
not a wide-to-long reshape.

The transform is gap-filling and segmentation: keep the event rows, guarantee at least
one row per subject, append a terminal censored row at end of follow-up when the last
event is earlier, then lag `iv_event` within subject to turn a set of event *times* into
a set of left-truncated *intervals*.

## Contract

```r
hzr_repeated_events(data, id, time, followup, indicator)
```

`data` is a data frame in the long shape above. The remaining arguments are
`character(1)` **column names**, matching the macro's string-valued keyword parameters
and avoiding an NSE dependency.

| argument | macro parameter | meaning |
|---|---|---|
| `id` | `id=` | subject identifier |
| `time` | `iv_event=` | interval from t=0 to each event |
| `followup` | `iv_end=` | interval from t=0 to end of follow-up |
| `indicator` | `eventype=` | 1 = event; 0 or `NA` = not an event |

**Returns** a data frame: every input column retained, with the `time` and `indicator`
columns **mutated in place** (the macro overwrites both), plus:

| column | meaning |
|---|---|
| `event` | 1 when this row is the event of interest, else 0 |
| `event_no` | running count of event repeats within subject |
| `rcensor` | right-censored indicator (see caveat below) |
| `iv_start` | interval from t=0 to the start of this segment |
| `iv_seg` | segment duration, `time - iv_start` |
| `renewal` | modulated-renewal count |
| `first` | group flag: this was the subject's first row at stage 4 |
| `last` | group flag: this was the subject's last row at stage 4 |

### Argument decisions and what they cost

**Exported, not an `inst/sas-parity/` helper.** `.hzr_derive_primisol()` is the stated
precedent, but it hardcodes *one* study's DATA step. `%repeat` is general, parameterised
and reused by every repeated-events job in the corpus, so it earns a real API with
roxygen, `\value` and examples, and a permanent commitment.

**Four input arguments, fixed output column names.** The macro's six output-name
parameters (`event=`, `event_no=`, `rcensor=`, `iv_start=`, `iv_seg=`, `renewal=`) are
not mirrored. Renaming a returned column is one line of R, so six extra arguments buy
nothing but documentation burden. If the cardioversion job passes non-default output
names, `hzr_translate_sas()` renames after the call — a translator change, not an API
change.

**`renewal` is in scope for v1.** It is two lines inside a stage we write regardless,
and the hard part — the stale `first`/`last` carry-forward, below — must be implemented
correctly whether or not the column is returned. Omitting the column while keeping the
machinery would be the worst of both. It also makes the corpus renewal variant
(`hz.te123.OMC.renewal`) a second parity target at no extra cost.

## Column-count discrepancy, asserted not hidden

SAS's macro output is 367 columns. Two of them, `lag_iv` and `number`, are pure DATA-step
loop state — a retained lag and a running counter. Returning them from a CRAN export
would be noise, so this implementation **drops those two and returns 365**.

The parity test asserts `ncol(out) == 367 - 2` with both dropped names spelled out in the
test, so the gap is an assertion that can fail rather than a silent difference a future
reader mistakes for a defect.

## Pipeline

Eight stages, the first seven mirroring the macro's seven DATA steps. Each is a small
internal function, so each is testable alone and the row ladder is observable stage by
stage rather than only at the end.

1. `rcensor <- 0`.
2. Sort by `(id, time)`; drop non-first rows that are not events.
3. Pad a subject with no first event (`time <- followup`, `rcensor <- 1`); drop leading
   non-events for multi-row subjects.
4. Record `first`/`last` as data columns; append a terminal censored row at `followup`
   for each subject whose last row has `rcensor == 0`.
5. Subset to `rcensor == 1 | indicator == 1`; derive `event`.
6. Lag `time` within subject to get `iv_start`; derive `event_no` and
   `iv_seg <- time - iv_start`; set `rcensor <- 1` where `time == followup`; derive
   `renewal` from the **stage-4** `first`/`last`.
7. Drop zero-duration non-event rows that are not a subject's first row.
8. Apply the job's post-macro `iv_start` guard (see Open Items 2).

Base R throughout — `order()`, `split()`, group flags computed from run boundaries.
**No new dependency.**

## Two traps this design exists to avoid

### The stale `first`/`last` carry-forward

The `renewal` assignment reads `first` and `last`, which are **plain data columns**
created at stage 4 — not the `first.&id` / `last.&id` automatic variables of the DATA
step it actually runs in. Between stage 4 and stage 6 the data is subset
(`if &rcensor=1 or &eventype=1`) and re-sorted.

So `first`/`last` are stale group flags carried forward. Recomputing them at stage 6
would silently give different answers for any subject whose rows were removed by that
subset. **They must be carried, never recomputed.** This produces a plausible-looking
wrong column rather than an error, which is precisely the failure mode `AGENTS.md`
names as this package's signature defect.

### Sort stability at tied times

Every stage re-runs `proc sort by &id &iv_event`. SAS's PROC SORT is **not** guaranteed
stable; R's `order()` is. Where a subject has two rows at the same `(id, time)`, R and
SAS can disagree about which row `first.&id` selects, and that propagates into
`event_no` and `renewal`.

This is not fixable in general. It is a documented parity caveat, and the parity test
asserts whether any tied `(id, time)` pairs exist in `bd_card` at all. If none do, the
trap is inert for this job and we say so **with evidence**, rather than assuming it.

## Acceptance

The reference ladder is recorded in
`hz.ce_cardioversion_repeated.ehb.log`, from `bd_card` (3246 obs, 357 vars):

```
3246 x 357  -> bd_card, the macro input
3246 x 358
 709 x 358
 963 x 360
 963 x 361
 963 x 367  -> macro output
 962 x 367  -> after the job's own post-macro iv_start guard
```

and the fit input must then match the listing's independent tallies: 962 obs, 388
events, 574 right-censored, 387 left-censored, `IV_EVENT` in
[0.0001140795, 12.99137], `IV_START` in [0.0001140795, 9.716832].

### PHI constraint

The cardioversion data is PHI. It stays on the study volume and never enters the repo.
The ladder therefore **can never run in CI**. The hardcoded counts and ranges above are
aggregates, not PHI, and may be written into a test.

### Test structure

**Unconditional synthetic tests** (`tests/testthat/test-repeated-events.R`), which run
everywhere and carry the real burden:

- a subject with no first event;
- a last event falling exactly at `followup`;
- a single-row subject;
- a zero-duration segment;
- a subject whose rows are removed at stage 5, so the stale `first`/`last` differ from
  flags recomputed at stage 6 — this is the test that fails if the carry-forward is
  reimplemented naively;
- tied `(id, time)` rows.

Assertions state exact expected vectors. Never
`expect_true(all(x %in% <every possible value>))` — see `AGENTS.md`.

**Volume-gated parity test**, `skip_if` on the presence of `/Volumes/qhsstudies`:
asserts the row ladder stage by stage, then the tallies and ranges above. Row coverage
is asserted **before** any value is compared.

Because a gated test prints green when it never ran, it is not evidence. The synthetic
set must be able to fail on its own.

## Open items

1. **Macro identity is unconfirmed.** The evidence is suggestive, not proof: this copy
   accepts all seven keyword parameters the job passes, where a rival copy at
   `/Volumes/qhsstudies/cardiac/valves/mitral/mr_severe_lvf/1999/datasets/repeat.sas`
   declares only `in=`/`out=` and would have produced five SAS errors the job's log does
   not contain; and the log numbers 162 source lines for the included file against this
   file's 161. Confirm by diffing the log's echoed source against the file. **If it is a
   different `repeat.sas`, the stage list above is wrong.**

   ⚠️ That diff settles **identity only**. It is not the acceptance test and must not be
   allowed to stand in for one: a correctly identified macro can still be reimplemented
   wrongly, and the diff would be just as clean. The log's shape ladder under
   **Acceptance** is what says the reimplementation is right. Read the source to know
   *which* macro ran; read the ladder to know whether the R code reproduces it.
2. **The guard predicate at stage 8 is unknown.** The volume was unmounted when this was
   written, so `hz.ce_cardioversion_repeated.ehb.sas` could not be read. "Post-macro
   `iv_start` guard" narrows it to roughly `iv_start < iv_end` or a zero/missing test,
   but which one determines *which* of the 963 rows is dropped. Read the job and fix the
   predicate before implementing stage 8.

3. **The ladder-to-stage mapping is not established.** The log records seven shapes and
   the macro has seven DATA steps, but they do not line up one-to-one: stages 2 and 3
   both leave 709 x 358, and the counts imply at least one step's output is not
   separately listed. The parity test is specified above as asserting the ladder "stage
   by stage", which presumes a mapping. Derive it from the log's actual step boundaries
   before writing that test; until then, only the endpoints (3246 x 357 in, 962 out) are
   known to be anchored to a specific stage.

All three open items require `/Volumes/qhsstudies` to be mounted. None blocks
implementing stages 1–7; item 3 blocks only the shape of the gated parity test.

## Gates and mechanics

The definition of done is `AGENTS.md`'s, with two corrections that apply to this change
specifically.

**Run `spelling`, which the definition of done omits.** `spelling.yaml` is a *required*
blocking check on `main`, but `AGENTS.md`'s definition-of-done list does not name it, so
following that list alone reddens CI. This change is unusually exposed: it introduces a
crowd of SAS identifiers (`iv_start`, `iv_seg`, `iv_event`, `iv_end`, `eventype`,
`rcensor`, `event_no`) into roxygen prose. Run

```r
spelling::spell_check_package(use_wordlist = TRUE)
```

as part of the gate, and add the identifiers it flags to `inst/WORDLIST`. Note that
`inst/WORDLIST` feeds only this check — it has no effect on CRAN's own incoming aspell
check, for which `devtools::check_win_devel()` is the source of truth.

**Build the check tarball under the session scratchpad, not `$TMPDIR/tree`.**
`AGENTS.md`'s recipe hardcodes a shared path that collides when two sessions run in the
same tree, and the collision surfaces as a fabricated `R CMD check` ERROR that looks like
a package defect. Use a session-unique directory.

Otherwise unchanged, and in this order: `devtools::document()`, then
`lintr::lint_package()` (0 lints), then `devtools::test()` (0 failures), then
`spelling::spell_check_package()`, then once per PR `R CMD check --as-cran` with the
manual from a clean `git archive` export of a **committed** tree.

## Out of scope

- Any change to the likelihood, the optimizer, or the existing `hzr_*` surface.
- Wiring `hzr_repeated_events()` into `hzr_translate_sas()`. The translator will need to
  recognise a `%repeat` call and emit this function, but that is a separate change with
  its own parity evidence.
- Running the cardioversion fit itself and reproducing LL -267.885. That is the payoff,
  and it is the *next* piece of work; this spec covers building the input it needs.
