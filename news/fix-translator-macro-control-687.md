* `hzr_translate_sas()` now records SAS macro control flow that decides
  whether a fit runs, where it used to report a clean translation (#687). A
  `%HAZARD(`, `%HAZPRED(` or `%REPEAT(` call inside a `* ... ;` comment, which
  SAS's macro processor runs, was dropped with the comment; it now gets an
  untranslated row. A PROC HAZARD inside a `%MACRO` definition, or under an
  open-code `%IF` or `%DO`, was emitted as a plain fit, though SAS runs it
  once per call, or only when the condition holds; it is still emitted, now
  with an untranslated row and a warning chunk above the fit. Two jobs in the
  public HAZARD corpus show each shape.
