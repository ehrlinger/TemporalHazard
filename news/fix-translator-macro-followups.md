* **Three `hzr_translate_sas()` defects left by the macro-scope and
  `_X1` fixes are corrected.** A macro whose body sorts the `PROC HAZPRED`
  grid in place, such as `PROC SORT DATA=PRED NODUPKEY;`, was not counted as
  writing it, so calling the macro left the grid emitted with the rows the
  sort removes, and no `$untranslated` row. Such a call now refuses the grid,
  as a `DATA PRED;` step in the body does; a sort with no `DATA=` is recorded
  as a call that may rewrite the grid. A `PROC HAZPRED` inside a
  `%MACRO ... %MEND` definition read its grid where the text stands, not
  where the macro is called: a job that rebuilt the grid with `AGE=70` before
  calling the macro got predictions at the earlier `AGE=50`. Its grid and
  predictions are now refused, with an `$untranslated` row saying the grid
  depends on where the macro is called. The status chunk now assigns each
  derived column directly, `D[[".hzr_status"]] <- with(D, ...)`, instead of
  calling `transform()`, which before R 4.4.0 renamed a column such as `_X1`
  to `X_X1` whenever the chunk ran on data that already carried
  `.hzr_status`, so a second fit on the same dataset could not find the
  column.
