* `hzr_translate_sas()` reads a macro call against the definition completed
  before it, as SAS does, not the last one in the file. A macro defined twice,
  whose first body rewrites a PROC HAZPRED grid, was read against its second
  body, and the grid was emitted with no `$untranslated` row (SAS predicts at
  the rewritten grid). The reverse order no longer refuses a grid over a body
  defined after the call, and a call before any definition is read as a macro
  the file does not define, which SAS does not resolve (#651).
