* `hzr_stepwise()` under the default `criterion = "score"` no longer passes
  over a candidate holding `NA` in silence (#686). Such a candidate cannot be
  scored, and was declined as `non_numeric` with no warning from a completed
  run, so a weaker variable could enter in its place. It is now declined as
  `missing_values`, and the run warns that it was not tested.
  `hzr_bootstrap()` counts the replicates that declined one this way in its
  existing warning about declined entries, since the candidate's selection
  frequency otherwise tracks how often its `NA` rows were drawn.
