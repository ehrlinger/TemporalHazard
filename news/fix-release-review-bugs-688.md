* Errors from the 1.3.0 release review (#688). `hzr_stepwise()` could not step
  a fit whose `Surv()` response deparsed to more than one line, and
  a screen that could drop aborted on a multi-level factor named in `force_in`.
  `hzr_read_outhaz()` failed on a dataset written with `NOCOV`; it now reads,
  and `se.fit = TRUE` is what refuses. `INHAZ=OA` now finds `OUTHAZ=WORK.OA`.
  `hzr_repeated_events()` now refuses a negative `time` or `followup`, which
  used to build a frame without a warning. A multiphase fit without 'numDeriv'
  warns once, not once per start. `?hzr_stepwise` now says a refit does not
  inherit the base fit's `control`.
