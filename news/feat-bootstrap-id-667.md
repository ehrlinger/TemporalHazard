* **`hzr_bootstrap(id =)` resamples patients, not rows (#667).** Data with
  several rows per patient, such as time-varying covariates in
  counting-process form (`Surv(start, stop, event)` or `time_lower =`), were
  resampled row by row, so a replicate could hold some of a patient's
  segments and not others, or one of them twice. `id` names the patient,
  as a column of `data` or a vector with one value per row; each replicate
  then draws patients with replacement and keeps every row of each one
  drawn, and `fraction` counts patients. This is SAS `%HAZBOOT` with
  `_TVC=1, _IDVAR=`. The default, `id = NULL`, resamples rows as before, and
  the same `seed` gives the same replicates as before. A fit with entry
  times and no `id` now draws a warning that its rows may not be patients.
