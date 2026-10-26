* **`hzr_deciles()` miscounted expected events when a caller's `event_time`
  disagreed with a row's `time_lower` (#656).** Both cases were silent.
  An `event_time` at or below a left-truncated row's entry dropped the
  entry correction, H(entry), so the row counted H(event_time) expected
  events where none are defined; that is now an error naming the first such
  row. And a row with `time_lower` equal to its stored time, which the
  likelihood reads as having no entry, was read as truncated once
  `event_time` was later than that time, losing H(time). Truncation is now
  decided as the likelihood decides it: status 0 or 1 with
  `0 < time_lower < time` at the stored time. Calls that leave `event_time`
  at its default, and `hzr_gof()`, are unchanged.
