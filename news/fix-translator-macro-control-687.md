* `hzr_translate_sas()` now records SAS macro control flow that decides
  whether a step runs, where it used to report a clean translation (#687). A
  `%HAZARD(`, `%HAZPRED(` or `%REPEAT(` call inside a `* ... ;` comment, which
  SAS's macro processor runs, was dropped with the comment; it now gets an
  untranslated row. A PROC HAZARD, PROC HAZPRED or `%REPEAT` call inside a
  `%MACRO` definition, or under an open-code `%IF` or `%DO`, was translated as
  if it ran once, though SAS runs it once per call, or only when the
  condition holds; it is still translated, now with an untranslated row and a
  warning chunk above it. `%THEN%DO` written without a space is read as a
  `%DO`, here and in the HAZPRED grid's dataset scan. Two jobs in the public
  HAZARD corpus show the first two shapes.
