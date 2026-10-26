* **On a multiphase `time_windows` fit saved by 1.2.9 or earlier,
  `predict(newdata = )` and `hzr_gof()` could swap a phase formula's own
  covariates for the global ones (#658).** Such a fit has no record of which
  phases took the global design, so it was read from the stored design's
  column names. A phase formula can use columns named like the
  window-expanded global ones (`x1_w1`, `x1_w2`) that hold other covariates,
  and that phase was read as inheriting. `predict(newdata = )` was then off
  by up to 77% relative in one example, and `hzr_gof()`'s curve moved, with
  no warning. A phase is now read as inheriting only when its stored design
  equals the global design, expanded at the stored times, in values as well
  as names. Fits saved with the record, and every fit made by a current
  version, are unchanged. `hzr_deciles()` drops the local copy of this check
  that #657 added.
