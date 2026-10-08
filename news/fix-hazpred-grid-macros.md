* **A `PROC HAZPRED` grid translated by `hzr_translate_sas()` now respects
  SAS macro scope.** A `DATA PRED;` step inside a `%MACRO ... %MEND`
  definition was read as if it ran where it stands. A job that built
  `PRED` with `AGE=50` and then defined a macro rebuilding it with
  `AGE=70` got predictions at 70, with no `$untranslated` row, although
  the macro was never called and SAS predicted at 50. A definition's
  steps now count only for a `PROC HAZPRED` inside the same definition.
  Calling a macro whose body writes the grid refuses the grid, as a
  macro's `OUT=` already did. An `%INCLUDE`, or a call of a macro the file
  does not define, between the grid's DATA step and the `PROC HAZPRED`
  may rewrite the grid out of sight: the grid the job shows is still
  emitted, and each such call now has an `$untranslated` row saying so.
  Jobs with no macros translate as before; the public corpus's calls and
  `$untranslated` rows are unchanged.
