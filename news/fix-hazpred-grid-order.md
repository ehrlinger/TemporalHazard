* **`hzr_translate_sas()` could emit a PROC HAZPRED grid in an order SAS
  does not use, with no `$untranslated` row.** A `PROC SORT` with no `DATA=`
  sorts SAS's most recent dataset, `_LAST_`. The translator took that to be
  the dataset the job last named, but `DATA;`, an `OUTPUT` statement with no
  `OUT=` (PROC MEANS, SUMMARY, REG and others), `OPTIONS _LAST_=` and a macro
  whose body does one of these all change it unseen. A grid of 12, 6, 1
  followed by `DATA; ... PROC SORT; BY MONTHS;` was emitted as 1, 6, 12. A
  `PROC SORT` with no `DATA=` now refuses the grid when it sorts a dataset
  the grid reads or stands between a step the grid uses and the PROC
  HAZPRED, `OPTIONS _LAST_=` refuses it wherever it stands, and `DATA;` and
  an `OUTPUT` with no `OUT=` are recorded as writing an unnamed dataset.

* **A grid sorted by a character variable could come out in an order SAS
  does not use.** SAS compares character values in ASCII order and at the
  column's fixed length, which its first assignment sets: after `G='B'`,
  the values `'AB'` and `'AA'` are both stored as `'A'` and keep their input
  order, where R sorted the full strings, and in the session's locale.
  `OPTIONS SORTSEQ=` changed SAS's order unseen too. The emitted grid code
  now stops when a `PROC SORT` it reproduces is by a character variable, and
  `OPTIONS SORTSEQ=` or `NOSORTEQUALS` refuses the grid. Numeric sorts are
  unchanged; on the public corpus no emitted grid sorts by a character
  variable.

* **A DATA step between the grid and the PROC HAZPRED could rename or
  delete the grid through a function, and the grid was still emitted.** Any
  `CALL` routine, and functions such as `RENAME()`, `FDELETE()`, `SYSTEM()`,
  `OPEN()` and `FETCH()`, now refuse a grid the step stands inside.

* **A `DATA PGM=` statement between the grid and the PROC HAZPRED no longer
  passes unseen.** It runs or stores a compiled DATA step program, which may
  rewrite the grid, and now refuses a grid it stands inside.
