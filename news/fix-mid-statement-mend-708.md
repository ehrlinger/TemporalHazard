- **`hzr_translate_sas()` now closes a `%MACRO` whose `%MEND` follows text
  with no semicolon (#708).** A job that wraps its old job-control lines in
  `%macro skip; ... %mend;` has no semicolon before the `%MEND`, so the
  check read the definition as never closing and, playing safe, skipped
  nothing: every fit in a macro the job never calls stayed live. SAS reads
  `%MEND` wherever it stands, and so does the check now. A `%MEND` inside a
  quoted string is still text.
