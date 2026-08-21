# Design: make the SAS translator emit documents that run

Status: approved design, not yet implemented.
Date: 2026-08-21. Target: TemporalHazard 1.2.2 GitHub release.
Issues: #151, #152, #153, #154, #155.

## The problem

`hzr_translate_sas()` emits Quarto documents that cannot render, and reports
success over them. Reproduced on a two-block job by evaluating the emitted
chunks in a clean environment holding only the input data frame:

```
coverage: 10 / 10 tokens mapped        untranslated rows: 0

  data     ok
  fit      ERROR: object 'INT_DEAD' not found
  grid     ok
  pred     ERROR: object 'fit' not found
  pred_haz ERROR: object 'fit' not found
```

Neither `R CMD check` (0/0/0) nor the 2528-test suite sees it. The tests were
written by reading the emitter, so they inherit its assumptions --
`test-sas-translate-fits.R` manufactures exactly the bindings the emitted
document lacks.

## Evidence gathered before designing

**Production corpus** (three studies under `/Volumes/qhsstudies`, 873 `.sas`
files, 93 containing `PROC HAZARD`). Read as code only; nothing from it is
copied into this repo.

| Construct | Count | Consequence |
|---|---|---|
| `SELECTION` | 36 / 93 | #152 promoted -- over a third of jobs |
| `HAZPRED` in the same file | 50 / 93 | #151's cluster hits the majority |
| `ICENSOR` | 64 / 93 | the censoring path is the common case |
| `LCENSOR` | 10 / 93 | |
| `LCENSOR` **and** `ICENSOR` | **0 / 93** | #155 is latent, not live |
| `LN_MAX = <num>` | 0 / 93 | #153's suspected regex trap does not occur here |

Production jobs are written in lowercase; `.hzr_sas_normalise()` upcases at
`R/sas-lex.R:161`, so case is handled. The public 110-file corpus is uppercase
and stylistically unrepresentative.

**Reference implementation** (`~/Documents/GitHub/hazard/src/llike/setlik.c`):

```c
c3w    = c3 * weight;                       /* C3 = COUNT OF INTERVAL CENSORED */
c1c2c3 = c1w + c2 + c3w;
llike  = -(c1c2c3) * (cumhaz - cumhst);
if (c3 > ZERO) llike += (c3w * lct);
```

C3 multiplies the contribution: #154 is a real defect. And `cumhst` = `H(STIME)`
is subtracted for **all** rows including interval ones -- SAS carries three
distinct times (`TIME`, `CTIME` interval lower, `STIME` entry) and supports the
combination natively.

**OUTHAZ datasets carry the phase structure** (`_NAME_` rows: `G1FLAG`,
`G3FLAG`, `FIXDEL0`, `FIXMNU1`, `FIXGE2`, `FIXGAE2`, then `DELTA/THALF/NU/M`,
`TAU/GAMMA/ALPHA/ETA`, `E0/C0/L0`), so a `predict()` method over a loaded
`OUTHAZ=` is reconstructible from the file alone.

## Design

### 1. Data masking in `hazard()`

`hazard()` evaluates `time`, `status`, `time_lower`, `time_upper` and `weights`
in `data`'s scope when `data` is supplied and `formula` is `NULL`:

```r
time <- eval(substitute(time), data, parent.frame())
```

The `subset()`/`transform()` idiom -- base R, no new dependency. Columns win;
non-columns fall through to the caller's environment, so `time = df$x`, a local
vector and a bare column name all keep working. With `data = NULL` nothing
changes: purely additive.

This is the root fix. `data =` is currently accepted and silently ignored on
the vector path, which is what makes the emitted call fail. Fixing it in
`hazard()` rather than in the translator also repairs the path for anyone who
passes `data =` by hand.

Rejected: emitting `AVCS$INT_DEAD` (noisy, and leaves `hazard()` misleading);
`with(AVCS, hazard(...))` (stores a call with an implicit data reference, which
`AGENTS.md` warns bites anything that rewrites or resamples a stored call); the
formula interface (requires translating the package's `-1/0/1/2` codes into
`Surv()`'s different integers, already a shipped wrong-answer bug, and `Surv`
cannot express left truncation combined with interval censoring).

### 2. Emitted document (#151)

```r
#| label: status
.hzr_status <- with(AVCS, ifelse(DEAD == 0 & C3FLAG > 0, 2, DEAD))

#| label: fit
fit <- hazard(
  data = AVCS, time = INT_DEAD, status = .hzr_status,
  time_lower = ifelse(.hzr_status == 2, ICTIME, STARTTME),
  weights = C3FLAG,
  dist = "multiphase", phases = list(...), theta = c(...),
  control = list(condition = 14), fit = TRUE
)
```

- **(a)** bind the fit: `fit <- hazard(...)`, keyed to the resolved chunk name
  (`fit`, `fit_2`, ...). Reuse the `call("<-", ...)` pattern already at
  `R/translate-sas.R:200`.
- **(b)** emit `fit = TRUE`. Without it `converged` is `NA` and `theta` holds
  the SAS *starting* values.
- **(c)** the `status` chunk gains `with(<data>, ...)`; masking covers only the
  inside of `hazard()`, not a free-standing chunk.
- **(d)** the prediction grid's column is named `time`, not the SAS `DO`
  variable, because `predict.hazard()` requires it (`R/hazard_api.R:1006`).
- **(e)** `hzr_read_outhaz()` gains an `hzr_outhaz` class and a
  `predict.hzr_outhaz()` that reconstructs the phase spec from the `G1FLAG` /
  `G3FLAG` / `FIX*` rows.

### 3. `SELECTION` (#152)

Emit `hzr_stepwise(fit = <fitname>, scope = ~ X1 + X2 + X3, ...)`. Today the
candidate pool is baked into the phase formula as forced-in terms, which
inverts the meaning of the SAS statement, and `fit`/`scope` are both absent so
the call cannot run at all.

### 4. Log grid (#153)

Use SAS's step exactly -- `INC = (5 + LN_MAX)/99.9`, not `length.out = 100`
over the closed interval, which implies `/99`. Divergence today reaches ~9.6%
at the last grid point. Also tighten the `MAX *= *([0-9.]+)` regex so it cannot
match `LN_MAX=`.

### 5. `ICENSOR` count (#154)

Emit `weights = <count variable>`. Confirmed against `setlik.c` above.

### 6. `LCENSOR` + `ICENSOR` (#155)

**Refuse loudly**: record it in `$untranslated` and emit a `stop()` in the
document. Supporting it properly needs a separate entry-time argument in
`hazard()` (SAS's `STIME`, distinct from `time_lower`'s interval role), which
is an API change to the core entry point. 0 of 93 production jobs use the
combination, so that gets its own issue rather than blocking this release.

### 7. Test harness -- build first

A helper that writes the `.qmd`, `parse()`s it, and `eval()`s its chunks in a
**clean environment holding only the input data frame**. It fails on every job
today; that is the failing test to start from.

- Unit tests use it per construct.
- `test-sas-translate-corpus.R` upgrades from "N jobs parse" to "N jobs render".

This single assertion covers 1a-1e and #152 together. Every defect in this
document was invisible to a suite that asserted the *shape* of emitted calls;
the fix is to execute them.

## Out of scope

- A separate entry-time argument in `hazard()` (the #155 root cause).
- CRAN submission. This targets a GitHub release; `cran-comments.md` still
  describes 1.2.1 and needs rewriting before any submission.
- Whether `C3 > 1` actually occurs in production data -- that needs the data,
  not the job text.

## Success criteria

1. The render harness passes on every construct-level unit test.
2. The public 110-file corpus test asserts documents **render**, not parse.
3. A `SELECTION` job emits a callable `hzr_stepwise()` with a real `scope`.
4. `document()` clean, `lintr` 0 lints, `devtools::test()` 0 failures,
   `R CMD check --as-cran` with the manual 0/0/0 under the 10-minute budget.
5. NEWS.md's experimental note is rewritten to describe what now works.
