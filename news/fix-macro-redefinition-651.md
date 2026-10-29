* `hzr_translate_sas()` reads a macro call against the definition completed
  before it, as SAS does, not the last one in the file. A macro defined twice,
  whose first body rewrites a PROC HAZPRED grid, was read against its second
  body, and the grid was emitted with no `$untranslated` row (SAS predicts at
  the rewritten grid). The reverse order no longer refuses a grid over a body
  defined after the call, and a call before any definition is read as a macro
  the file does not define, which SAS does not resolve. A macro defined inside
  another macro, including one that redefines itself, is now one this cannot
  read: SAS compiles that definition only when the outer macro runs, and each
  such job had emitted the earlier grid with no row. Its calls follow the rule
  for any unread call, and when the enclosing definition holds a statement
  this cannot read, such as `OPTIONS OBS=`, the call refuses the grid wherever
  it stands. A `%MACRO` inside a macro call's argument list, such as quoted
  `%NRSTR()` text, is no longer read as a definition (#651).
