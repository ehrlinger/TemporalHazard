# SAS job translator: `PROC HAZARD` / `PROC HAZPRED` → R + Quarto

Design doc. Status: **approved through architecture; v1 scope set by measurement.**
Date: 2026-08-20. Branch: `feat/sas-job-translator`.

## 1. Purpose

Translate legacy SAS HAZARD jobs into `.qmd` documents that reproduce them in
`TemporalHazard`. Two audiences, deliberately coupled:

1. **A shipped analyst tool** — one export, `hzr_translate_sas()`, for the
   SAS-to-R migration.
2. **A parity harness** — every translated job whose SAS `.lst` exists becomes an
   automated SAS-versus-R fixture.

The coupling is the point. The harness renders (or `purl`s) the *emitted
document* and compares the resulting fit against SAS, so the thing under test is
the artifact that ships, not a parallel code path. Without it, the tool produces
R that looks plausible; with it, the tool produces R that is checked.

Secondary uses — bulk migration of the production corpus, and reproducible
reporting — fall out of the same machinery and are not designed for separately.

## 2. Evidence base

Nothing below is inferred from reading jobs. Two independent sources:

**The grammar** is HAZARD's own lexer, `src/hazard/hazard_l.l` and
`src/hazpred/hazpred_l.l` in the GPL-2 reference implementation: **117 keyword
rules, 77 distinct parser tokens, 30 aliases, 11 context start conditions.**

**The corpus** was measured with `tools/sas-hazard-profile.R` over three
production studies plus the public examples. Counts are `PROC HAZARD` blocks
unless noted:

| | public | preserve_root | lv_function | resilia |
|---|---|---|---|---|
| `.sas` files | 110 | 349 | 194 | 337 |
| HAZARD / HAZPRED blocks | 44 / 55 | 115 / 92 | 82 / 71 | 92 / 82 |
| `MI` (→ `MAXITER`) | 7 | 113 | 75 | 86 |
| `ICENSOR` | 0 | 43 | 74 | 42 |
| `WEIBULL` (in `PARMS`) | 19 | 87 | 76 | 88 |
| `NOCONSERVE` / `CONSERVE` | 3 / 27 | 68 / 41 | 16 / 31 | 51 / 33 |
| `INHAZ` resolved in-tree | 40% | 26% | 11% | 26% |
| distinct unresolved librefs | 7 | **1** | **1** | **2** |

Three consequences drive the design:

- **Interval censoring is the production norm** (`ICENSOR` on 74 of 76 blocks in
  `lv_function`) and is absent from the public corpus entirely. Anything designed
  against `examples/` would have missed it.
- **Conservation of Events inverts in production** — `NOCONSERVE` outnumbers
  `CONSERVE` in two of three studies, the reverse of the public corpus.
- **Unresolved `INHAZ=` collapses to one or two librefs per study.** A
  libref→path map needs *one entry* to resolve ~60 jobs.

## 3. Architecture

```
     .sas file
         │
    [ sas-lex ]          comment strip; mirrors <STMT>\*[^;]*;
         │
  [ context-aware parser ]  ← .hzr_sas_grammar, GENERATED from the lex sources
         │
   hzr_sas_job  ──────────► coverage report (translated / untranslated / unresolved)
         │
   [ renderer ]           deparse() the stored calls into chunk fences
         │
     .qmd
```

One export, `hzr_translate_sas()`. Everything else internal.

| Path | Contents |
|---|---|
| `data-raw/hazard-grammar.R` | generator; reads the `.l` files, writes `R/sysdata.rda` |
| `R/sas-lex.R` | comment stripping, normalisation, block extraction |
| `R/sas-parse-job.R` | context-aware parser → `hzr_sas_job` |
| `R/sas-job.R` | class, validator, `print` method |
| `R/sas-render-qmd.R` | renderer |
| `R/translate-sas.R` | `hzr_translate_sas()` |

### 3.1 The grammar is generated, not hand-written

`.hzr_sas_grammar` is produced by `data-raw/hazard-grammar.R` from the reference
lexer. Columns: `proc`, `context`, `keyword`, `token`, `is_alias`, `r_target`,
`implementation_status`, `note`.

A hand-written table captures only the spellings the author happened to see. The
grammar has 30 aliases that fall through to one parser token — `QUASI`/
`QUASINEWTON`, `BW`/`BACKWARD`, `SLE`/`SLENTRY`, `NOH`/`NOHAZ`, `MI`/`MAXITER`,
`P`/`PRINTIT`, `PARMS`/`PARAMETERS`, and the single-letter phase options
`E`/`I`/`M`/`O`/`S`. A corpus-inferred table fails silently on every spelling the
sample lacks.

Generating it also makes totality checkable: every one of the 117 keywords has a
row and a status, asserted in `tests/testthat/test-sas-grammar.R`. Regenerating
against a newer HAZARD becomes a drift check.

**Licence.** `hazard` is GPL-2, `TemporalHazard` is MIT. The generator lives in
`data-raw/` (already `.Rbuildignore`d), reads the `.l` files from a local
checkout path, and **only the extracted keyword table ships**. No GPL source
enters the tarball.

### 3.2 Parsing must be context-aware

Not a stylistic preference — the grammar has genuine collisions:

| spelling | context | token |
|---|---|---|
| `M` | `PARM` | `M` (early-phase shape) |
| `M` | `PHOP`, `STEP` | `MOVE` |
| `NOS` | `STEP` | `NOPRINTS` |
| `NOS` | `HZPP` | `NOSURV` |

`M` occurs 76–115 times per study as a `PARMS` shape parameter. A context-free
lookup translates it as `MOVE` and produces a wrong model with no error. The
guarantee the parser relies on — *within one proc and one context, a spelling
resolves to exactly one token* — is asserted as a test rather than assumed.

### 3.3 Block bounding follows the macro, not the PROC

`hazard_l.l` puts `)` in its whitespace class, so HAZARD itself does not treat it
as a terminator: block delimitation is the `%HAZARD(...)` **macro's** job.
Bounding a block by scanning forward for `);` therefore fails whenever the block
contains a nested paren — measured at 38 of 82 blocks in `lv_function`.

The parser instead matches **balanced parentheses from `%HAZARD(` / `%HAZPRED(`**.
A block must never run to end of file: that is how comment prose reaches the
token tables, and comment prose is where study and patient detail live.

### 3.4 `hzr_sas_job` is a thin index over R calls

```r
structure(list(
  source       = list(path =, sha256 =),
  calls        = list(fit = quote(hazard(...)), pred = quote(predict(...))),
  grid         = quote(data.frame(MONTHS = ...)),
  inhaz        = "LIB.DSNAME",   # or NULL
  outhaz       = "LIB.DSNAME",
  untranslated = data.frame(line =, construct =, reason =),
  coverage     = list(tokens_seen =, tokens_mapped =)
), class = "hzr_sas_job")
```

Model state is stored as **unevaluated calls**, not a bespoke schema. Rendering
is then `deparse()` rather than string templating; the parity harness `eval()`s
the fit call directly; the round-trip test compares call objects. Storing calls
is already house style here — `hazard(formula, data)` stores one.

The intermediate exists for exactly one reason that emitted text cannot serve:
**cross-job `INHAZ=`/`OUTHAZ=` resolution**, which needs an index of what every
*other* job wrote. Coverage accumulation and harness reuse are secondary.

## 4. Coverage policy

| Outcome | Emitted `.qmd` | `hzr_translate_sas()` |
|---|---|---|
| Fully translated | runnable | returns quietly |
| Keyword with `no_r_equivalent` | `UNTRANSLATED` callout, construct + line | warns, lists them |
| Keyword `unmapped` | `UNTRANSLATED` callout | warns, lists them |
| **Unresolved `INHAZ=`** | `stop()` at the head of the chunk | warns loudly |
| Nothing recognised | **no file written** | **errors** |

The `INHAZ` row is the one the corpus forces: a document that cannot resolve its
fit must **fail to render** rather than render a report over a missing model.

### 4.1 `INHAZ=` resolution and the `librefs` contract

Resolution order: another translated job's `OUTHAZ=` first, then the
`librefs` argument, then failure. The map is worth having because the
measurement says one entry clears ~60 jobs per study.

`librefs` is a **plain named character vector**, `c(LIBNAME = "path")`, and that
is the whole contract. It is deliberately *not* a path to a config file:

- **No new dependency.** Reading a YAML file would put `yaml` in the
  dependency graph. `Imports:` is currently the single package `survival`, and
  that cheap CRAN footprint is worth protecting.
- **No coupling to a site convention.** `_study.yml` is an institutional
  convention that a public CRAN user — the package's primary persona — does
  not have and should not need. A study that has one reads it itself and passes
  the result in; the package never learns the file exists.
- **It degrades honestly.** With no `librefs`, unresolved `INHAZ=` still fails
  loudly per the table above, rather than half-working.

Two derived guarantees, both testable:

- **Round-trip** — parse → render → re-parse yields the same `calls`.
- **Totality** — every keyword in `.hzr_sas_grammar` is exercised by at least one
  fixture, so rows cannot accumulate that nothing covers.

## 5. v1 scope, set by measurement

**In scope** (used in production, ordered by measured frequency): `DATA`,
`TIME`, `EVENT`, `PARAMETERS`, `OUTHAZ`, `INHAZ`, `OUT`, `MAXITER`, `CONDITION`,
`CONSERVE`, `NOCONSERVE`, `STEEPEST`, `QUASINEWTON`, `PRINTIT`, `NOCOR`,
`NOCOV`, `NOLOG`, `NONOTES`, `NOPRINT`, `NOHAZ`, the `PARM` parameter family
(`MUE`/`MUC`/`MUL`/`THALF`/`NU`/`M`/`TAU`/`GAMMA`/`ALPHA`/`ETA` and their `FIX*`
forms), `WEIBULL`, `EARLY`/`CONSTANT`/`LATE`, `SELECTION`, and
**`ICENSOR`/`LCENSOR`/`RCENSOR`**.

**Deferred** (present in the grammar, absent from all four corpora): `BUNDLE`,
`FIXGAE2`, `FIXGE2`, `FIXMNU1`, the `PHOP` family (`EXCLUDE`/`INCLUDE`/`MOVE`/
`ORDER`/`START`), `FAST`, `MAXSTEPS`, `MAXVARS`, `NUMERIC`, `ROBUST`,
`SEMIROBUST`, `ONEWAY`, `CLIMITS`, `NOSURV`, `NOPRINTQ`, `NOPRINTS`.

### 5.1 Highest-risk mapping: the censoring statements

`ICENSOR`/`LCENSOR`/`RCENSOR` get a dedicated parity fixture and the loudest
flag in the table. This package codes censoring `-1` left, `0` right, `1` event,
`2` interval, while `survival::Surv(type = "interval")` uses `0`/`1`/`2`/`3` for
right/event/left/interval. Carrying one coding into the other is a wrong answer
with no error, and it has shipped here before.

**`LCENSOR` is left-*truncation*, not left-censoring.** Its operand is the
counting-process entry time (`hazard_y.y`'s own comment: "start of followup"),
so it maps to `time_lower` and never touches `status`. HAZARD has no
left-censoring in this package's sense: `status = -1` is never produced by
this translator. All four corpus uses are literally `LCENSOR STARTTME`,
paired with `TIME INT_TE` (see `~/Documents/GitHub/hazard/examples/hz.te123.OMC.sas`).

**`ICENSOR`'s first operand is an event *count* (OBS column 4, `C3`), not a
0/1 flag.** A row is interval-censored where `C3 > 0`; the second operand
(`CTIME`) is the interval's *lower* bound, and `TIME` is the upper bound --
the interval runs `CTIME` -> `TIME`. The likelihood in `setlik.c` uses a
Nelson-type approximation `C3 * ln([CF(T) - CF(CT)] / (T - CT))`. `TIME`
being the upper bound is exactly what `hazard()`'s `time_upper` already
defaults to when left `NULL`, so `time_upper` is never emitted by this
translator.

| SAS | `hazard()` |
|---|---|
| `LCENSOR <var>` | `time_lower = <var>` -- entry time, `status` unchanged |
| `ICENSOR <c3> = <ctime>` | rows where `<c3> > 0`: `status = 2`, `time_lower = <ctime>` |
| `RCENSOR <var>` | `status = 0`, no branch needed -- but `<var>` is `C2`, a *count*, so it is the censored rows' `weights` (#162) |
| `EVENT <var>` | rows where `<var> > 0`: `status = 1` -- and `<var>` is `C1`, a *count*, so it is those rows' `weights` (#157) |

It also touches known open work: the 2026-08-19 `preserve_root` parity run left a
named residual where SAS maximises the interval-censored likelihood but *reports*
a log-likelihood treating interval rows as exact-event densities.

### 5.2 Conservation of Events interacts with censoring

`R/likelihood-multiphase.R:1172` auto-disables CoE when
`!all(status %in% c(0, 1))`, which `ICENSOR` guarantees. R will therefore often
disable CoE **for a different reason than SAS did**. A job carrying both
`CONSERVE` and `ICENSOR` is where the two could diverge silently, so it gets an
explicit fixture. `CONSERVE`/`NOCONSERVE` is always emitted, never defaulted.

### 5.3 `SELECTION` forks to `hzr_stepwise()`

Phase statements carry up to **320** covariates in `lv_function` — those are
stepwise *candidate pools*, not fitted covariates. A job with `SELECTION` routes
to `hzr_stepwise()`; without it, to `hazard()`. Note that HAZARD's lexer collapses
`SELECT`, `FORWARD`, `FW`, `SW` and `STEPWISE` into **one token**, while
`BACKWARD` is distinct.

**`SELECTION FORWARD` maps to `direction = "both"`, not `"forward"`.** In
`hazard_y.y`, `STEPWISE` sets option 21, `BACKWARD` 22, `ONEWAY` 34 — and every
spelling of forward/stepwise lexes to the same `STEPWISE` token, so HAZARD does
not distinguish them. Option 21 is therefore two-way. Corroborated by `slstay`,
which would be meaningless under a forward-only search. The obvious mapping
(`FORWARD` → `"forward"`) is wrong on every stepwise job and fails silently, so
this pair gets an explicit fixture.

| SAS | `setopt` | `hzr_stepwise(direction =)` |
|---|---|---|
| `FORWARD` / `FW` / `SW` / `SELECT` / `STEPWISE` | 21 | `"both"` |
| `BACKWARD` / `BW` | 22 | `"backward"` |
| `ONEWAY` / `NOSTEPWISE` / `NOSW` | 34 | no stepwise; plain `hazard()` |

`direction = "forward"` is reachable in R but **not** from a translated SAS job.

### 5.4 Prediction grids

92% of `HAZPRED` `DATA=` grids **classify** as `log_grid` or `explicit_do`,
with zero unclassifiable across all three production studies.

⚠️ **That is a classification rate, not a translation rate, and the two differ
sharply.** Measured on the public corpus after the parser landed:

| profiler class | classified | parser resolved |
|---|---|---|
| `log_grid` | 8 | 4 |
| `explicit_do` | 3 | **0** |
| `derived_set` | 3 | 0 |
| other / not found | 6 | 0 |

The classifier only asks whether the DATA step contains a `DO` loop. The parser
needs the bounds to be **literal numbers**, and real grids are written with
expressions over constants defined in the same step:

```
DTY=12/365.2425;
DO MONTHS=1*DTY,2*DTY,3*DTY,...,24 TO 180 BY 12;
```

Refusing such a grid outright is correct — a partly-read grid is a partial
`newdata`, which is the hollow-result failure this design exists to prevent.
But it means `explicit_do`, the **most common** form in `lv_function` (32
occurrences), currently resolves at 0%.

**Closing that gap means resolving the DATA step's own constant assignments**
(`NAME = <numeric expression>;`) before evaluating the `DO` list. That is a
bounded extension, not open-ended DATA-step support, and it is the single
highest-leverage change available to grid coverage. Deferred pending a decision;
until then unresolved grids emit `UNTRANSLATED` and the document says so.

`derived_set` grids emit `UNTRANSLATED` by design.

### 5.5 `%repeat` translates to `hzr_repeated_events()`

Designed 2026-09-10, branch `feat/translate-repeat`. `hzr_repeated_events()` (#241)
reimplements the macro; see `inst/dev/REPEATED-EVENTS-DESIGN.md`. This section covers only
how a job's `%repeat(...)` call reaches it.

**Evidence.** A scan of `/Volumes/qhsstudies` on 2026-09-10 found **at least 526** `.sas`
files calling `%repeat(`, 138 of them `hz.*` or `tp.hz.*` HAZARD jobs. The scan was still
running when this was written, so the count is a floor. The public corpus has none. Three
consulting templates recur, eight copies each; among the copies checked, seven of each are
identical and one differs:

| job | rewrites `OUT=` after the macro | reads a renamed output |
|---|---|---|
| `tp.hz.repeated_events` | no | `LCENSOR IV_START` |
| `tp.hz.repeated_events.modulated_renewal` | no | `EVENT EV_INFD` |
| `tp.hz.repeated.event.weighted` | yes: `KEEP` plus an `IV_EVENT` nudge | `CN_RADMT`, `EV_RADMT` |
| `hz.ce_cardioversion_repeated.ehb` | yes: the line 65 nudge | `LCENSOR IV_START` |

#### Extraction and parsing

`.hzr_sas_blocks()` gains a third block kind, `REPEAT`, matched on `%REPEAT *\(` and bounded
by its own balanced parentheses. Blocks therefore come out in file order, which is the one
property this feature cannot get wrong: the `%repeat` chunk must precede the fit that reads
its `OUT=`. Each block also records `start` and `end` offsets into the normalised text.
Those fields are additive, so existing callers are unaffected. A job with `%repeat` and no
`HAZARD`/`HAZPRED` block still errors ("no HAZARD or HAZPRED block found"); the `tp.bd.*`
dataset builders are out of scope.

`.hzr_parse_repeat()` (in `R/sas-parse-job.R`) splits the body on depth-0 commas into
`KEY=VALUE` pairs and fills the macro's twelve defaults, upper-cased. It returns the same
shape as `.hzr_parse_hazard()` (`call`, `untranslated`, `tokens_seen`, `tokens_mapped`)
plus `in` and `out`. Each argument counts as one token seen and, if accepted, one mapped.

These are refused at translate time. Each becomes a `stop()` chunk in place of the call,
plus an `$untranslated` row:

- an unknown keyword, or a positional argument (SAS itself errors on both);
- a value that is not a plain SAS name, such as `&VAR` or an empty value;
- an argument column (`ID`, `EVENTYPE`, `IV_EVENT`, `IV_END`) equal to one of the
  eight output names, for example `EVENTYPE=RCENSOR`. SAS zeroes that indicator at its
  first DATA step, so the job's own answer is garbage, and the call text alone shows it;
- two output names that are equal to each other.

#### Emitted code

For the cardioversion call, `%repeat(in=bd_card, out=events, id=ccfid, eventype=ce_card,
iv_end=iv_end, rcensor=cn_card, event=ev_card)`:

```r
# label: data -- the existing guard, once per distinct IN=
if (!exists("BD_CARD")) stop("This job read BD_CARD from a SAS DATA step, ...")

# label: repeat
EVENTS <- local({
  d <- BD_CARD
  drop <- intersect(names(d), c("EV_CARD", "EVENT_NO", "CN_CARD", "IV_START",
                                "IV_SEG", "RENEWAL", "FIRST", "LAST"))
  if (length(drop)) {
    warning(...)  # names the dropped columns
    d <- d[setdiff(names(d), drop)]
  }
  out <- hzr_repeated_events(d, id = "CCFID", time = "IV_EVENT",
                             followup = "IV_END", indicator = "CE_CARD")
  names(out)[match(c("event", "event_no", "rcensor", "iv_start", "iv_seg",
                     "renewal", "first", "last"), names(out))] <-
    c("EV_CARD", "EVENT_NO", "CN_CARD", "IV_START", "IV_SEG", "RENEWAL",
      "FIRST", "LAST")
  out
})
```

`OUT=` joins the set of datasets that an emitted chunk builds, so the fit's `DATA=EVENTS`
gets no "Assign EVENTS" guard. An `IN=` that is itself an earlier `%repeat`'s `OUT=` gets
none either.

**The rename is always emitted, not only for non-default names.** The translator
upper-cases the whole job, so the fit references `IV_START`, while the function returns
`iv_start`. R is case-sensitive, so without the rename the fit reads a column that does
not exist.

**The drop is what makes the rename safe.** Suppose `IN=` already carries a column named as
a rename target, say `EV_CARD`. Without the drop, `names<-` would produce two columns named
`EV_CARD`, and `$EV_CARD` returns the first one, which is the stale input. That would be a
wrong answer with no error. Dropping reproduces SAS exactly, because the macro assigns each
of its eight outputs unconditionally on every row before any statement reads it: `rcensor`
at the first DATA step, `first` and `last` at the fourth, `event` at the fifth, and
`event_no`, `iv_start`, `iv_seg` and `renewal` at the sixth. An argument column can never be
dropped, since that case was refused at translate time. A lower-case `rcensor` in a
mixed-case frame is not dropped either; it reaches `hzr_repeated_events()`'s own refusal,
which is loud.

#### Rewrites of `OUT=` after the macro

The maintainer decided on 2026-09-10 that a job's post-macro statements are not folded into
the function. Before each `PROC HAZARD` whose `DATA=` is an earlier `%repeat`'s `OUT=`, the
text between that macro's `end` and the fit's `start` is scanned for a `DATA` statement that
names `OUT`, or for `CREATE TABLE OUT`. Each hit emits one `stop()` chunk, directly before
the first such fit that follows it; a later fit reading the same `OUT=` gets no second copy,
since the first stop already halts the render. The chunk quotes the step: the normalised text, cut at the next `DATA`, `PROC`, `%HAZ`,
`%REPEAT` or `RUN;`. It tells the reader to replace the chunk with R code that makes the same
change. The hit is also recorded in `$untranslated`, which the translate-time warning
reports.

A silently skipped rewrite would be the house failure mode: a fit that converges over data
SAS never fitted. In the cardioversion job that means 15 zero-length segments left
un-nudged. Placing the stop at the fit, rather than right after the macro, still catches a
rewrite that follows an intervening block. `PROC SORT` is not treated as a rewrite, because
reordering rows does not change the likelihood. The quote carries no line number: every
`$untranslated` row today has `line = NA`, since the parser sees only normalised text.

#### The input dataset

`IN=` is usually built by DATA steps this translator does not translate. It gets the same
`if (!exists(...)) stop()` guard that a `PROC HAZARD` `DATA=` gets. The message asks the
reader to assign the dataset as it stood at the `%repeat` call, with the columns named as
the job spells them, in upper case. That is the existing convention for `DATA=`, kept
deliberately; the emitted code does not upper-case names for the reader. A data frame read
by `haven` keeps the stored case (`ccfid`), and `hzr_repeated_events()`'s own error for a
missing column names the one it wanted.

#### ⚠️ Open: a `LAG_IV` or `NUMBER` column in `IN=`

This was found while writing this section and is not part of the approved design. The
macro's sixth DATA step keeps two loop variables with `RETAIN lag_iv 0 number 0`. A
variable that also arrives through `SET` is re-read from the input on every iteration,
which overwrites the value retained from the previous row. So if `IN=` has a `LAG_IV` or
`NUMBER` column, SAS computes `iv_start` and `event_no` from that column rather than from its
own lag. `hzr_repeated_events()` has no such coupling, so the translated document would
disagree with the SAS listing without saying so. Here SAS is the one that is wrong: this is
a finding to document, not behaviour to reproduce. The cardioversion input has neither
column; SAS's jump from 361 to 367 columns at that step is exactly the six variables the
step adds. **Proposed:** the `repeat` chunk warns at render, naming the column and saying
that SAS's `iv_start`/`event_no` for this job are not comparable.

#### Tests

A new file, `tests/testthat/test-sas-translate-repeat.R`, evaluates the emitted chunks with
`eval()` into an environment that holds a small synthetic `BD_CARD`. It compares against
vectors derived by hand from the macro, never from the function:

1. End to end: the chunk names are exactly `c("data", "repeat", ..., "fit")`; `EV_CARD`,
   `CN_CARD`, `IV_START` and `IV_EVENT` equal hand-derived vectors; and a fit that uses
   `LCENSOR IV_START` runs (`skip_on_cran()`). This test fails if the rename is removed.
2. A non-default output name read by `EVENT EV_X`: the lower-case `event` is absent, and
   `EV_X` holds the exact vector.
3. The clash drop: `BD_CARD` carries a stale `EV_CARD` of wrong values. The test expects a
   warning naming it, exactly one `EV_CARD` column, and the macro's values. It fails if the
   drop is removed.
4. Each translate-time refusal: a `stop()` chunk whose `eval()` errors with the expected
   message, and an exact `$untranslated$construct`.
5. A post-macro rewrite: evaluation halts at the stop chunk before the fit. Negative
   controls, so the detector cannot pass by always firing: `DATA OTHER; ...` and
   `PROC SORT DATA=EVENTS` produce no stop.
6. Defaults: `%REPEAT()` guards `BUILT` and calls with `ID`, `EVENTYPE`, `IV_EVENT` and
   `IV_END`.
7. Extraction: a mixed fixture yields block kinds exactly `c("REPEAT", "HAZARD",
   "HAZPRED")`.

A volume-gated case in `test-repeated-events-parity.R` translates the real cardioversion
`.sas` file. It binds `BD_CARD` from `.hzr_derive_bd_card()`, with names upper-cased, and
asserts that evaluation halts at the line 65 stop chunk. It then checks that `EVENTS` is 962
by 365 and holds 388 `CE_CARD` events, checking the shape before comparing any value. It
does not fit the model; reproducing LL -267.885 is the next piece of work.

## 6. Testing

- `test-sas-grammar.R` — table populated, statuses valid, keyword unique per
  context, alias runs terminate, the `M`/`NOS` collisions are real.
- `test-sas-lex.R` — comment stripping against constructed fixtures, including
  the `; * ... ;` inline form and apostrophes inside comments.

  Note that comment termination is **quote-agnostic**: HAZARD's rule is
  `<STMT>\*[^;]*;`, so a comment ends at the first `;` whether or not an
  apostrophe precedes it. Treating it as quote-aware silently swallowed the
  following statement — caught in Task 1 review, in both the line-initial and
  inline paths. A quote-*aware* semicolon search is still needed for splitting
  statements, where `TITLE 'a; b';` must not split; that is what `.first_semi()`
  is for, and it must not be reused for comments.
- `test-sas-parse-job.R` — golden `hzr_sas_job` objects for representative jobs.
- `test-sas-translate-parity.R` — `skip_on_cran()`; render emitted `.qmd`,
  compare against stored SAS `.lst`.

Three further fixtures come from the gap analysis in §7 and are worth having
whether or not the translator ships:

- **G3 Weibull correspondence** (§7.1) — constrained G3 in R against
  `PARMS ... WEIBULL` in SAS. Converts the most common production configuration
  from assumed to shown.
- **`CONSERVE` + `ICENSOR`** (§7.3) — the combination where R and SAS could
  diverge despite usually agreeing.
- **Censoring status round-trip** (§5.1) — `ICENSOR`/`LCENSOR`/`RCENSOR`
  through the translator and into `hazard()`, asserting the exact status vector
  rather than membership in the set of possible codes.

Fixtures are synthesised from the public `examples/` corpus or deliberately
de-identified. **No production job content enters the repository.**

## 7. Package gaps surfaced by this analysis

The grammar plus the usage counts is, incidentally, a gap analysis of
`TemporalHazard` against what production actually does. These items are
independent of the translator and stand on their own.

| SAS construct | production use | `TemporalHazard` | verdict |
|---|---|---|---|
| `CONDITION=` | every block | `control$condition` (default 14) | present |
| `MI` / `MAXITER` | 75–113 | `control$maxit` | present |
| `NOCOR` / `NOCOV` | 6–28 | documented no-ops | correct |
| `QUASI` / `QUASINEWTON` | 21–109 | `control$method = "bfgs"` | equivalent |
| `CONSERVE` / `NOCONSERVE` | 16–68 | auto-disables on `status` outside {0,1} | right answer, different mechanism |
| `WEIBULL` late phase | 76–88 | expressible, unverified, undocumented | **verification gap** |
| `STEEPEST` | 14–109 | `method` offers only `"bfgs"` / `"nm"` | **absent** |
| `DELTA` | 0 of 3 studies | unmapped | tracked as [#143](https://github.com/ehrlinger/temporal_hazard/issues/143) |
| `RESTRICT` | 0 of 3 studies | no equivalent | defer |

### 7.1 The G3 Weibull variant is unverified, not missing

`PARMS ... WEIBULL` is `setopt(6)` → `SETG3_weibull()`, which sets
`Late.g3flag = Late.g3flag + 2`. HAZARD's late phase has four structural
variants (`g3flag` 1–4) and its likelihood branches on `1||3` versus `2||4`.
`hzr_decompos_g3()` has no variant argument.

That looks like a missing feature and is not one. The R general form

```
G3(t) = (((t/tau)^gamma + 1)^(1/alpha) - 1)^eta
```

collapses at `alpha = 1, eta = 1` to `(t/tau)^gamma` — a Weibull cumulative
hazard. The corpus corroborates the reading: `WEIBULL`, `FIXTAU`, `FIXGAMMA`
and `FIXALPHA` occur the *same* number of times in each study (87/87/87/87 in
`preserve_root`, 76/76/76/76 in `lv_function`).

So the most common production configuration is one the package can probably fit
and has **never been shown to fit**. HAZARD branches on `g3flag` for numerical
reasons — at the constraint it computes `expm1(log1p(exp(y)))` where the answer
is `exp(y)` — and R always takes the general path. Whether that costs precision
at the boundary is untested.

**Verified 2026-08-20 (R side only).** `hzr_decompos_g3(t, tau, gamma, alpha = 1,
eta = 1)` equals `(t/tau)^gamma` to 8.6e-16 at moderate times and 5.0e-13 across
1e-6 to 1e6; the hazard agrees to 1.9e-16. The worse extreme-times figure is the
`expm1(log1p(exp(y)))` precision loss predicted above — real, and harmless at
that magnitude. **The SAS leg is still outstanding:** this shows the identity
holds inside R, not that R agrees with `PARMS ... WEIBULL`.

**Action:** complete the fixture with a SAS reference run. Then document the
correspondence on `hzr_phase()`, because an undocumented equivalence is one a
user cannot find.

### 7.2 `STEEPEST` has no R equivalent

Jobs write `STEEPEST QUASI` *together* — steepest descent, then quasi-Newton.
`control$method` offers `"bfgs"` and `"nm"`; there is no steepest-descent option
and no two-stage strategy.

This is not cosmetic. The multiphase likelihood is multimodal, and `n_starts`
explores a neighbourhood rather than the space, so a different descent path can
land on a different optimum. A silent mapping of `STEEPEST QUASI` to `"bfgs"`
would surface later as an unexplained numerical discrepancy rather than as a
missing feature.

**Action:** decide explicitly — either add it to `control$method`, or document
that it is unsupported and that translated jobs may find a different optimum.
Either is defensible; silence is not. Tracked as
[#145](https://github.com/ehrlinger/temporal_hazard/issues/145).

### 7.3 Conservation of Events agrees by coincidence

`R/likelihood-multiphase.R:1172` disables CoE when
`!all(status %in% c(0, 1))`, which `ICENSOR` guarantees. Production writes
`NOCONSERVE` on 16–68 blocks per study, so R and SAS usually agree — but for
different reasons. A job carrying both `CONSERVE` and `ICENSOR` is where they
could diverge.

**Action:** a fixture for exactly that combination.

## 7.4 End-to-end parity, measured

`examples/hz.death.AVC.sas` translated from disk and fitted in R against
`inst/extdata/avc.csv` — the same 310-row dataset the SAS job reads:

| | |
|---|---|
| R, from the translated `.sas` file | **−210.500617** |
| SAS, per the job's own comments | **−210.5** |

| R parameter | fitted | `exp()` | SAS `PARMS` |
|---|---|---|---|
| `phase_1.log_mu` | −1.443200 | 0.236173 | `MUE=0.2361727` |
| `phase_1.log_t_half` | −1.889173 | 0.151210 | `THALF=0.1512095` |
| `phase_1.nu` | 1.438635 | — | `NU=1.438652` |
| `phase_1.m` | 1.000000 | — | `M=1 FIXM` |
| `phase_2.log_mu` | −7.517091 | 0.00054370 | `MUC=0.0005436977` |

Log-likelihood to six significant figures, every parameter to five or six.
This is one job on a two-phase `cdf` + `constant` model with no covariates and
no censoring beyond exact/right — it is not a claim about the interval-censored
or stepwise paths, which have no comparable reference run yet.

## 8. Open questions

1. ~~**`FORWARD` ≡ `STEPWISE` at the token level.**~~ **Resolved 2026-08-20** from
   `hazard_y.y`: they are genuinely synonymous (option 21, two-way). See §5.3.
2. **`_CLLSURV`/`_CLUSURV` confidence-limit type.** `predict.hazard()` supports
   `conf.type = "log-log"` and `"logit"`. Which does HAZPRED emit? Guessing
   silently produces wrong limits.
3. ~~**`LCENSOR` precedence against `EVENT`.**~~ **Resolved 2026-08-21** from
   `inst/dev/FIXTURE-GAP-LIST.md`'s answered Q1 and the reference grammar:
   `LCENSOR` is left-*truncation* (the counting-process entry time), not
   left-censoring, and never competes with `EVENT` or `ICENSOR` for
   `status` at all — it maps only to `time_lower`. `status = -1` is
   therefore never produced by this translator, and there is no
   `LCENSOR`/`ICENSOR` tie-break to resolve: `ICENSOR`'s `C3 > 0` and
   `LCENSOR`'s entry time apply to disjoint outputs (`status` vs.
   `time_lower`'s non-interval branch) and can both fire on the same row
   without conflict. `EVENT` still outranks `ICENSOR` for `status` — a row
   with `EVENT == 1` stays coded `1` even where `C3 > 0` — which remains a
   conservative reading chosen for symmetry, not confirmed against a SAS
   reference run. (`src/hazard/varterm.c`: `EVENT` and `ICENSOR` are
   alternatives — a job may specify `ICENSOR` alone — so the translator
   accepts that and rejects only a job with neither.) See §5.1.
4. **`CONSERVE` + `ICENSOR` in SAS.** Does SAS conserve events when interval rows
   are present, or silently disable as R does?

None blocks implementation; each blocks a specific fixture from being trusted.
