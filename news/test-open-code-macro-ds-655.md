* A test now pins the `hzr_translate_sas()` refusal of an open-code `DATA &DS;`
  that rewrites a PROC HAZPRED grid (`%LET DS=PRED; DATA &DS; ...` between the
  grid and the `PROC HAZPRED`). The refusal itself shipped in #628. Before it,
  the macro-named write went unrecorded and the earlier grid was emitted with no
  `$untranslated` row (#655).
