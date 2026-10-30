* Documentation corrected from the 1.3.0 release review (#688). The
  getting-started vignette said one `log_mu` under Conservation of Events has
  no standard error; all three have had one since #652. `hzr_stepwise()` and
  `stepwise_trace()` examples no longer pass a `control` that warns it is
  ignored, `@param fit` says which fits are accepted, and `max_steps` says it
  counts `"frozen"` rows too and is checked between iterations. The `uslife2023` source no longer names an
  internal path.
