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

* **A grid sorted by a character variable followed the R session's locale,
  not SAS's ASCII order.** Under `en_US`, keys `B`, `_`, `A` sorted with `_`
  first, where SAS puts it last. The emitted `order()` now passes
  `method = "radix"`, which is byte order on every platform. Numeric sorts
  are unchanged.

* **A DATA step between the grid and the PROC HAZPRED could rename or
  delete the grid through a function, and the grid was still emitted.** Any
  `CALL` routine, and functions such as `RENAME()`, `FDELETE()`, `SYSTEM()`,
  `OPEN()` and `FETCH()`, now refuse a grid the step stands inside.
