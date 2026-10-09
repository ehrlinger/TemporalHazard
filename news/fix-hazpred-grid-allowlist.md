* **`hzr_translate_sas()` now emits a `PROC HAZPRED` grid only when it reads
  every statement that can write it.** It used to refuse a list of
  constructs known to rewrite a dataset and assume everything else was
  harmless, and each of three release reviews found another construct that
  produced a wrong grid with no `$untranslated` row. Four did so on the
  1.2.12 code: a macro call with no semicolon, which swallows the `DATA` or
  `PROC SORT` statement after it (`%SETUP` then `DATA PRED; AGE=70;` gave the
  earlier grid at `AGE=50`); a `WHERE` statement in a `PROC SORT` of the
  grid, which was ignored; a macro whose body replaces the grid through
  `PROC DATASETS`; and a grid step inside an open-code `%IF 0 %THEN %DO;`,
  which was read as if it ran. Each now refuses the grid, with an
  `$untranslated` row naming the statement. The rule is now an allow-list:
  global statements, `%LET` and `%PUT`, procedures whose output datasets are
  all named by an option or statement the translator reads, DATA steps, and
  `PROC SORT` with only `DATA=`, `OUT=` and an ascending `BY`. Anything else
  between a step the grid uses and the `PROC HAZPRED` refuses the grid. That
  includes an `%INCLUDE` or a call of a macro the file does not define,
  which until now was recorded with the grid still emitted. Such a call
  before the grid's first step cannot write a dataset the grid reads, but it
  may change session state the grid depends on, such as `OPTIONS OBS=`: the
  grid is emitted, and the call now has an `$untranslated` row saying so. A
  statement the translator cannot place outside any step, or `OPTIONS OBS=`
  itself, refuses the grid wherever it stands. The 15 grids the public corpus
  emits are unchanged; four of its jobs gain two such rows each, for an
  `%INCLUDE` and a macro call before the grid.

* **A `PROC HAZPRED` `TIME` statement with more than one variable now stops
  the block, as `PROC HAZPRED` does.** `TIME` takes exactly one name, and a
  second is a syntax error on which the procedure stops. `hzr_translate_sas()`
  took the first name and emitted the predictions. It now emits a `stop()`
  and an `$untranslated` row, as it does for a multi-operand `TIME` in
  `PROC HAZARD`.

* **A grid `DATA` step with `SET;` and no dataset name, or a numeric
  `LENGTH` below 8 bytes, now refuses the grid.** `SET;` reads the most
  recent dataset, and was read as no `SET` at all; a short numeric length
  stores fewer digits than the emitted grid holds.
