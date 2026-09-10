# Design: `hzr_repeated_events()` — an R reimplementation of the SAS `%repeat` macro

- **Date:** 2026-09-09
- **Status:** implemented on this branch. SAS parity verified 2026-09-10 against the reference
  job's log and listing -- see Acceptance. The post-macro step once planned as stage 8 does not
  exist as specified -- see Open Items 2. Both of the job's fits reproduce the listing,
  verified 2026-09-10 -- see Fit parity.
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
| `rcensor` | right-censored indicator (see the caveat below, "`rcensor` and `event` overlap") |
| `iv_start` | interval from t=0 to the start of this segment |
| `iv_seg` | segment duration, `time - iv_start` |
| `renewal` | modulated-renewal count |
| `first` | group flag: this was the subject's first row at stage 4 |
| `last` | group flag: this was the subject's last row at stage 4 |

### `rcensor` and `event` overlap

`rcensor` is not a plain censoring indicator: it is 1 whenever a row's event time
equals the end of follow-up, and that includes a genuine event row, not only a
padded censoring row. A subject whose last event happens to land exactly at
`followup` gets a row with `event == 1` AND `rcensor == 1` together.

This is faithful to the macro, which sets `rcensor` in exactly that case
(`.hzr_re_stage6()`, mirroring `if &iv_event=&iv_end then &rcensor=1;`) and says
so in its own header: the macro's output-column comment for `rcensor` warns this
combination is not compatible with a `HAZARD`-style fit. The behaviour is
correct and is not changed here -- what was missing was documentation. `\value`
and the `@note` on `hzr_repeated_events()` now say this explicitly, and
`test-repeated-events.R` pins a case where both columns are 1 on the same row so
a future change that "fixes" the overlap fails loudly.

Practical consequence: do not pass `rcensor` straight through as a censoring
indicator to a `HAZARD`-style fit without first accounting for this overlap.

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

## Missing-value semantics

SAS numeric missing is an ordered value that compares EQUAL TO ITSELF, so `if . = . then` is
true. A row with both event time and end of follow-up missing therefore counts as at end of
follow-up and gets `rcensor = 1`. The internal stage (`.hzr_re_stage6()`) preserves this: it
tests `event_time == followup` with both sides possibly `NA`, and treats the both-missing
case as at end of follow-up, matching the macro rather than propagating `NA`.

The exported function nonetheless REJECTS a missing `followup` value, because carrying it
through scrambles row order downstream: the appended censored row (stage 4) sorts before the
event it terminates, and `iv_start` comes out missing. This is a deliberate deviation from the
macro, chosen because a plausible-wrong data frame is worse than a refusal. Note that SAS's
own macro header assumes end of follow-up is present.

The exported function also rejects a non-numeric indicator column, a factor `id`, a `time` or
`followup` column that is not numeric (a `Date` column, notably, used to sail through as
days-since-epoch with no error), a zero-row `data`, and a reserved output column name already
present in `data`, and warns when no indicator value equals 1. A character indicator matches
neither the event test nor the non-event test, so every row is dropped and censored rows
backfilled, producing a frame indistinguishable from a legitimately all-censored cohort --
silently, with no error. A factor `id` sorts by level order rather than value, which can give a
different subject ordering than SAS's `proc sort`.

Three further preconditions are the calling program's responsibility in the macro, and
`hzr_repeated_events()` now WARNS (rather than errors) when they are broken, naming the
affected subjects: a missing `time` value leaking `NA` into `iv_start`/`iv_seg`, an event time
greater than `followup`, and a `followup` value that is not constant within a subject. These
stay warnings, not errors, because they are faithful to the macro -- it leaves the calling
program to guarantee them and does not itself refuse a violation. The hard refusals above are
reserved for inputs where the wrong answer is silent and there is no way for a caller to tell
after the fact (a `Date` column, a name clash, an empty input); these three instead produce a
recognisably degraded but inspectable result, so a warning that says so is enough.

## Column-count discrepancy, asserted not hidden

SAS's macro output is 367 columns. Two of them, `lag_iv` and `number`, are pure DATA-step
loop state — a retained lag and a running counter. Returning them from a CRAN export
would be noise, so this implementation **drops those two and returns 365**.

The parity test will assert `ncol(out) == 367 - 2`, with both dropped names spelled out in
the test, once Task 7 (the ladder-to-stage mapping, Open Item 3) is resolved and the test can
be written -- it does not exist yet. That gap is meant to be an assertion that can fail rather
than a silent difference a future reader mistakes for a defect.

## Pipeline

Seven stages, one for each of the macro's seven DATA steps. Each is a small
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
The reference job runs one more statement after the macro, its line 65. It adjusts values and
removes no rows, it belongs to the job rather than the macro, and it is not applied here --
see Open Items 2.

Base R throughout — `order()`, group flags computed from run boundaries. **No new
dependency.**

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

Every stage re-runs `proc sort by &id &iv_event`, and R's `order()` is stable. From memory,
`PROC SORT`'s default is `EQUALS`, which keeps tied rows in their input order; if so, an
earlier draft of this document was wrong to call SAS's sort unstable. That default has not
been checked against SAS's documentation here, and nothing below rests on it. What matters
is the order rows *reach* the first sort, and in the reference job that comes out of a
`PROC SQL` join whose output order is not guaranteed.

That caveat is not inert for this job: 2845 of `bd_card`'s 3246 rows sit in 545 tied
`(ccfid, iv_event)` groups. It is nonetheless harmless, and the parity test shows it rather
than assuming it. Shuffling the input changes none of the derived columns and neither
covariate of the job's second model (`maze_prc` and `iso_pvi`, both constant within every
subject). The only columns that move are six belonging to other event types -- `ce_cva`,
`ce_embol`, `dt_cva`, `dt_embol`, `iv_cva`, `iv_embol` -- which neither fit reads. Only one
tied group mixes an event with a non-event row, and none of the 303 tied groups at a
subject's first time does; that is the one place order would change which rows survive
stage 2.

## Acceptance

Verified 2026-09-10 with the study volume mounted. Each row is a DATA step's output as the
job's `.log` records it, against R:

| `.log` line | SAS | R | step |
|---|---|---|---|
| 188 | 3246 x 357 | 3246 x 357 | `bd_card`, rebuilt by `.hzr_derive_bd_card()` |
| 196 | 3246 x 358 | 3246 x 358 | stage 1 |
| 211 | 709 x 358 | 709 x 358 | stage 2 |
| 229 | 709 x 358 | 709 x 358 | stage 3, which removes nothing in this job |
| 244 | 963 x 360 | 963 x 360 | stage 4 |
| 252 | 963 x 361 | 963 x 361 | stage 5, which removes nothing in this job |
| 267 | 963 x 367 | 963 x 365 | stage 6; R returns two fewer columns by design |
| 290 | 962 x 367 | 962 x 365 | stage 7 removes one zero-duration row |
| 315 | 962 x 367 | -- | the job's line 65, which removes nothing |

The log lines in between are `PROC SORT`s and the macro's empty `data &out; set &out;`
step, which change no shape. After the job's line 65, the fit input matches the listing's
independent tallies exactly: 962 observations, 388 events, 574 right censored, 387 left
censored, `IV_EVENT` in [0.0001140795, 12.99137] and `IV_START` in
[0.0001140795, 9.716832]. The exported function raises none of its input warnings on this
data. The fits on it are recorded under Fit parity, below.

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

**Volume-gated parity test** (`tests/testthat/test-repeated-events-parity.R`): skips
unless the maze study datasets are mounted (overridable with `HAZARD_MAZE_DATASETS`) and
`haven` is installed. It rebuilds `bd_card`, asserts every shape in the table under
Acceptance **before** comparing any value, then the tallies and ranges, then that the
result does not depend on input row order.

Because a gated test prints green when it never ran, it is not evidence. The synthetic
set must be able to fail on its own.

## Fit parity

Verified 2026-09-10 with the study volume mounted. The job's two `PROC HAZARD` statements
were put through `hzr_translate_sas()`, and the emitted status and fit calls were evaluated
over the input above. The translator leaves out only `STEEPEST`, SAS's steepest-descent
warm-up (#145). Both fits start at the job's `PARMS`.

| | SAS log likelihood | R | largest estimate gap | largest SE gap |
|---|---|---|---|---|
| first model: `CONSERVE`, 9 parameters | -267.885 | -267.8850808 | 8.3e-05 SE | 8.9e-04 relative |
| second: `NOCONSERVE`, `maze_prc` and `iso_pvi` in both phases, 13 parameters | -242.254 | -242.2541227 | 4.0e-05 SE | 5.5e-04 relative |

The comparison is made on SAS's scale. `PROC HAZARD` estimates every shape parameter on the
log scale (`E3 = Ln(|Nu|)` and so on), so R's `nu`, `m`, `gamma`, `alpha` and `eta` are
logged, and their standard errors converted by the delta method, `se(x) / |x|`. An estimate
gap is measured in the listing's own standard errors. The standard errors agree to about
three figures rather than the seven printed, which is what SAS's numeric Hessian supports.
For the first model this includes `MUL`, which the listing marks as estimated in closed form
and still reports a standard error for.

**The likelihoods agree at a point, not only at an optimum.** R's log likelihood at the first
model's `PARMS`, at SAS's final estimates and at R's own optimum is -267.8850808 each time.
At the second model's final estimates from the listing it is -242.2541227. So a gap in a fit
can only come from the optimizer, and one did.

**The second model needs `reltol` tightened.** Under default control the translated call
reports `converged = TRUE` at -242.2675, 0.013 short of the listing. Its estimates are up to
0.16 SE away, and its standard errors up to 10% off. `.hzr_optim_generic()` defaults
`reltol` to 1e-5, and BFGS stops when an iteration gains less than about `reltol * |LL|`,
here about 0.0024, which a flat ridge in the late shapes does not provide. With
`reltol = 1e-12` it reaches the listing's optimum from the same start. The parity test sets
that explicitly. Whether the default should change is a separate question, because it
touches every fit and the check-time budget.

**One start, and why.** The default five starts give the first model the same fit from
start 1. The perturbed starts end at worse local optima, LL -267.914, -268.175 and -274.368,
and their Hessians raise the not-invertible warnings, which are about those starts and not
about the reported fit. The test fits from `PARMS` alone, and agreement rests on the log
likelihood at SAS's point rather than on `n_starts`.

Each tolerance in the test fails the stopped-short fit above: half a unit in the listing's
third decimal of log likelihood, 1e-3 SE per estimate, 2e-3 relative per standard error.
The fits are in `tests/testthat/test-repeated-events-parity.R` and skip without the volume,
so CI has never run them.

## Open items, resolved 2026-09-10

1. **Macro identity.** The log's source numbering jumps from the `%inc` at line 1444 to the
   job's next statement at 1606, so the included file is 161 lines -- exactly this file's
   length. (The brief that started this work said 162; that was off by one.) SAS did not
   echo the included source, so there is nothing to diff. What settles it is functional: R
   reproduces every shape in the log and every tally in the listing.
2. **The "guard" is not a filter.** The job's line 65 is
   `if iv_start ge iv_event then iv_event=iv_event+0.0001141553;`. It moves an event time
   forward by one hour, in 365-day years, where a segment would otherwise have no length. It
   touches 15 rows and removes none: 962 in, 962 out. The 963 -> 962 drop is the macro's
   own stage 7. So stage 8 as specified -- "apply the job's guard, which removes a row" --
   does not exist, and the decision to include it in the function rested on that misreading.
   The function already returns 962 rows. **Decided 2026-09-10 by the maintainer: line 65 is
   not applied.** It is job code, not macro code, specific to HAZARD needing a positive-length
   interval under `LCENSOR`. A caller fitting that way applies it, as the job does and as the
   parity test does.
3. **The ladder-to-stage mapping** is the table under Acceptance. Stages 3 and 5 remove no
   rows in this job, which is why shapes repeat.

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
- Changing the optimizer's default `reltol`. Fit parity shows the default stops the second
  model short; the fix belongs to the optimizer, not to this function.
