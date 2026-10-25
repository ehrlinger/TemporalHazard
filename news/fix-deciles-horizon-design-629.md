* **`hzr_deciles()` ranked a multiphase `time_windows` fit by the wrong
  window's coefficients (#629).** For the horizon ranking it read each
  subject's stored design, which `time_windows` had expanded at that
  subject's *exit* time, so a subject who exited in one window was ranked
  with that window's coefficient even when the horizon falls in another.
  On two 400-subject examples, 12 and 85 subjects landed in the wrong risk
  group, and `mean_survival` was wrong in every group, with no warning. The horizon
  design is now built at the horizon through `predict()`, as `hzr_gof()`
  and the single-distribution path already do. Expected events, taken at
  each subject's own exit time, are unchanged.
