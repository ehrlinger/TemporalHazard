* `hzr_translate_sas()` code failed on R before 4.4 when a phase statement
  named a variable that is not a syntactic R name, such as `_X1`: the fit
  stopped with "object '_X1' not found". The emitted status chunk adds its
  column with `transform()`, which renamed `_X1` to `X_X1` before R 4.4. The
  chunk now passes `check.names = FALSE`, so the column keeps its name on
  every R version (#609).
