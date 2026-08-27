# Translate a SAS HAZARD job into a Quarto document

Reads a SAS program containing `PROC HAZARD` and/or `PROC HAZPRED`
blocks and emits a Quarto document of the equivalent
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
and
[`predict.hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/predict.hazard.md)
calls.

## Usage

``` r
hzr_translate_sas(path, out_dir = NULL, librefs = NULL)
```

## Arguments

- path:

  Path to a `.sas` file.

- out_dir:

  Directory to write the `.qmd` into. `NULL` (default) parses without
  writing.

- librefs:

  Optional named character vector mapping SAS librefs to directories,
  e.g. `c(EX = "estimates")`, used to resolve `INHAZ=`. The resolved
  member is read as `<member>.sas7bdat`; give a librefs value that
  already ends in `.rds` or `.sas7bdat` (a specific file, not a
  directory) to name an already-converted fit directly, e.g.
  `c(EX = "estimates/hzdeath.rds")`.

## Value

An `hzr_sas_job` object, invisibly. `$calls` holds the emitted calls
keyed by chunk label, `$grid` the last prediction grid seen and `$inhaz`
the first unresolved `INHAZ=` (not all of each, when a job has several),
`$untranslated` the recorded gaps and `$coverage` the token counts.

## Details

**Experimental:** a job that translates does render – the emitted
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
chunk binds its fit and asks for an actual fit – but this is a
translation aid, not a turnkey reproduction, and some SAS constructs are
refused rather than translated. See the Experimental section below.

Constructs the translator does not cover are recorded on the returned
object and rendered as visible callouts, never dropped. A `PROC HAZPRED`
job whose `INHAZ=` fitted model cannot be located emits a
[`stop()`](https://rdrr.io/r/base/stop.html), so the document fails to
render rather than reporting over a model it did not load.

A job may contain more than one `PROC HAZARD` and/or `PROC HAZPRED`
block. Every block is preserved: the first of a kind keeps the bare
chunk name (`fit`, `pred`, `pred_haz`), later ones get `fit_2`, `fit_3`,
`pred_2`, and so on, so the emitted document contains every model, not
just the last one seen. `job$outhaz` is therefore a character vector,
one element per `PROC HAZARD` block that set `OUTHAZ=` (in the order the
blocks appear). Each `PROC HAZPRED` block's `INHAZ=` is resolved
independently against that vector, matching the most recently written
`OUTHAZ=` at that point in the file – mirroring SAS itself, where a
later `OUTHAZ=` write overwrites the dataset an earlier one wrote under
the same name. If a job's own `OUTHAZ=` values don't cover it, `librefs`
is tried next; distinct external `INHAZ=` values each get their own
loaded-fit chunk. When neither resolves and the job holds more than one
local fit, which fit a
[`predict()`](https://rdrr.io/r/stats/predict.html) call belongs to is
genuinely unknown; the emitted call falls back to referencing `fit` and
the ambiguity is recorded in `untranslated`, never guessed at silently.

## Experimental

The emitted document renders: the
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
chunk binds its fit to a name and passes `fit = TRUE`, so the
[`predict()`](https://rdrr.io/r/stats/predict.html) chunks have
something to predict from. Two SAS constructs are refused outright
rather than mistranslated, each emitting a
[`stop()`](https://rdrr.io/r/base/stop.html) in place of the fit: a
`SELECTION` statement requesting a stepwise screen (#152, \#160), and
`LCENSOR` combined with `ICENSOR`, which one `time_lower` argument
cannot express (#155). Prediction grids the parser cannot resolve are
refused whole, and the
[`predict()`](https://rdrr.io/r/stats/predict.html) chunks that would
have read such a grid become a
[`stop()`](https://rdrr.io/r/base/stop.html) naming it, rather than a
`predict(newdata = )` over a name no chunk builds. An unresolved
`INHAZ=` stops the render on purpose.

On a fit loaded from an external `INHAZ=` dataset, point predictions
work but `se.fit = TRUE` is refused when `PROC HAZARD` estimated a late
shape parameter on a composite scale – the generic unconstrained
three-phase case, not an exotic one. A translated `PROC HAZPRED` block
asks for confidence limits unless the SAS job says `NOCL`, so such a job
stops at its [`predict()`](https://rdrr.io/r/stats/predict.html) chunks.

`$coverage` counts tokens the parser recognised; it is not evidence that
the emitted calls execute, and a job can report full coverage with an
empty `$untranslated` while its document still errors on render. See the
1.2.2 `NEWS.md` entry. The function's API, the `hzr_sas_job` field
layout and the emitted document format are all expected to change.

## Examples

``` r
# \donttest{
job <- hzr_translate_sas(
  system.file("extdata", "hz-example.sas", package = "TemporalHazard")
)
# }
```
