- **`hzr_translate_sas()` no longer runs a model SAS never ran (#708).** A
  PROC HAZARD, PROC HAZPRED or `%repeat` inside a `%MACRO` that the job defines
  but never calls used to become a live chunk, so the document fitted a model
  SAS did not. Those chunks are now emitted with `eval: false` and a comment
  naming the macro, and listed in the job's new `$not_run`. The data guard and
  the fit a later prediction uses now come from the blocks that run. A call
  counts anywhere in the job's text, including a `*` comment and
  `CALL EXECUTE`. The check reads this file only, so a call from an
  `%INCLUDE`d file is not seen, and the `$untranslated` row says so. A macro
  the job does call keeps one flagged fit, as before.
