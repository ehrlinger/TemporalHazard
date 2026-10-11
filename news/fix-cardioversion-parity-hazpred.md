* The repeated-events parity test of the cardioversion job now expects the
  `$untranslated` row for its `PROC HAZPRED`, which stands inside a
  `%MACRO` definition and has been refused since #625. Its two fits had
  failed since that change. The test needs the study datasets, which no CI
  runner has, so the failure was not reported. The package's behaviour is
  unchanged.
