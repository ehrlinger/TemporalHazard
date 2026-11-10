* **`predict()` gives the standard error of a cohort's average curve
  (#731).** With `se.fit = TRUE, average = TRUE`, the predictions of all rows
  that share a `time` are averaged, one row per distinct time, with columns
  `time`, `n`, `fit`, `se.fit`, `lower` and `upper`. Use it for the expected
  survival of a cohort scored through a fitted model, to set against its
  Kaplan-Meier. The standard error is the delta-method standard error of the
  mean, so it allows for the correlation between patients. Averaging each
  patient's `lower` and `upper` instead, the only way before, treats them as
  perfectly correlated and gives a band that is too wide.
* **A `se.fit = TRUE` result from `predict()` carries the derivatives behind
  its standard errors**, as the attribute `"jacobian"`: each row's derivative
  with respect to the parameters, their covariance matrix and their names.
  Any weighted sum of the predictions now has a standard error, such as the
  expected number of events for a one-sample log-rank test, without
  rebuilding the package's handling of fixed parameters or of the Weibull
  scale. A Weibull fit's first column is `log(mu)`.
