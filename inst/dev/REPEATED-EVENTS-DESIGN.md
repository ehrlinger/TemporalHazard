# Design: `hzr_repeated_events()` — an R reimplementation of the SAS `%repeat` macro

- **Date:** 2026-09-09
- **Status:** implemented on this branch. SAS parity verified 2026-09-10 against the reference
  job's log and listing -- see Acceptance. The post-macro step once planned as stage 8 does not
  exist as specified -- see Open Items 2.
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
machinery would be the worst of both. (This paragraph first named the corpus job
`hz.te123.OMC.renewal` as a second parity target. It is not one: the OMC jobs never call
`%repeat`, and their renewal fit reads a `NOPREVTE` count from the raw input. The SAS
evidence for `renewal` is `ac.reintervention`; see Acceptance.)

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
data. Fitting it and reproducing LL -267.885 is out of scope; see below.

### Second job: `renewal`, against `ac.reintervention`

The cardioversion job never reads `renewal`. A first search for jobs that call `%repeat` and
then use its output (2026-09-10) covered only `general/` and `thoracic/`: it silently died
before reaching `cardiac/`. What it found:
- the consulting `tp.*` templates, which read `renewal` but were never run (none of their
  108 listings prints it);
- several thoracic jobs that print `renewal`.

A rescan of `cardiac/` found 1182 callers:
- 496 read `renewal` in code, and 404 of those have a `.lst`;
- 9 of those also save `%repeat` output to a library, and 7 of those datasets come from the
  same run as their listing;
- the maze cardioversion jobs are among the 1182, and read `iv_seg` and `event_no` only.

So `ac.reintervention` is one usable reference among several, not the only one.

The shones reoperation job (`cardiac/congenital/shones/outcomes/datasets/bd.repeated_reops.sas`)
was checked against its own output without an R rebuild. Its saved `bd_reop`, 241 rows over
121 subjects, matches its `.lst` table of `renewal` by `ev_reop` cell for cell.

The search tool itself caused two false negatives:
- **A partial result that looked complete.** `grep` here is ugrep, and the first search died
  without an error. It came back with zero `cardiac/` callers even though a known one exists.
- **Silence on SAS listings.** ugrep treats some `.lst` files as binary and then prints
  nothing, not even a zero count, unless given `-a`.

`ac.reintervention.sas` (achalasia; `thoracic/esophagus/benign/achalasia/outcomes/clinical`)
is the one verified here. It prints `ev_rein` by `renewal` directly after the macro, then saves the
macro's output unchanged as `library.bdrein`. So every returned column can be compared with
SAS's own, row by row.

Verified 2026-09-10 with the volume mounted, by `.hzr_derive_bdmult1()` and
`tests/testthat/test-repeated-events-renewal-parity.R`:
- **Rebuild:** `pndil1` 29 x 7, `reint1` 15 x 7, `cmb` 44 x 9, `bdmult` 425 x 218 and
  `bdmult1` 425 x 219, each as the `.log` records it.
- **`%repeat` stages:** 425 x 220 through stage 3, then 464 x 222 at stages 4 and 5 and
  464 x 228 at stages 6 and 7 (SAS; R's widths differ as the test states).
- **Crosstab:** the `.lst` table of `ev_rein` by `renewal`, cell for cell.
- **`bdrein`:** 464 rows over 420 subjects, in the same order, with no mismatch in
  `iv_rein`, `iv_fup`, `iv_start`, `iv_seg`, `event`, `rcensor`, `event_no`, `renewal`,
  `first` or `last`.

Two points of detail:
- **The job aliases the indicator.** It passes `event=ev_rein`, the same column as
  `eventype=`, so SAS overwrites its input indicator. R's `event` is compared with SAS's
  `ev_rein`.
- **One input line names a patient.** The job overrides `iv_fup` for one subject named by
  `ccfid`, and that subject has an event. The helper reads the identifier from the job's
  source on the volume at run time, so it never enters the repository.

What this settles, and what it does not. `bdrein`'s saved `first` and `last` are the
stage-4 flags, and they differ from recomputed flags on 75 rows. R matches SAS on all 75,
which is SAS evidence that the flags are carried forward.

But renewal gives the same answer on this data from stale flags or fresh ones:
- **The 36 stale `first` rows** are censored rows appended after a subject's only event, so
  they have `event_no = 1`. The first-row bump needs `event_no = 0`.
- **The 39 stale `last` rows** are each subject's last event, so they have `event = 1`. The
  last-row bump needs `event = 0`.

Shones shows the same from SAS's own saved columns. Stale and fresh flags differ on 122 of
its 241 rows, and the bump rule gives identical `renewal` either way; SAS matches both. That
makes two production jobs, and 705 rows, where the stale-flag branch changes nothing.

So the stale-flag branch of `renewal` is still pinned only by the synthetic fixture below.

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

**Volume-gated `renewal` parity test** (`tests/testthat/test-repeated-events-renewal-parity.R`):
it skips unless the achalasia datasets are mounted (overridable with
`HAZARD_ACHALASIA_DIR`). It rebuilds `bdmult1`, asserts every shape, then the `.lst`
crosstab, then `library.bdrein` row by row. Its failure output is counts only: no row and no
`ccfid`.

Because a gated test prints green when it never ran, it is not evidence. The synthetic
set must be able to fail on its own.

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
- Running the cardioversion fit itself and reproducing LL -267.885. That is the payoff,
  and it is the *next* piece of work; this spec covers building the input it needs.
