* **`predict()` gives the standard error of a cohort's average curve
  (#731, experimental).** With `se.fit = TRUE, average = TRUE`, the
  predictions of all rows that share a `time` are averaged, one row per
  distinct time, with columns `time`, `n`, `fit`, `se.fit`, `lower` and
  `upper`. Use it for the expected survival of a cohort scored through a
  fitted model, to set against its Kaplan-Meier. The standard error is the
  delta-method standard error of the mean, so it allows for the correlation
  between patients, and is never larger than the mean of the patients' own
  standard errors. The limits are the mean's own. Averaging each patient's
  `lower` and `upper`, the only way before, gives a different band that is
  not the mean's interval.
* **`hzr_jacobian()` returns the derivatives behind the standard errors of a
  `se.fit = TRUE` result from `predict()`** (experimental): each row's
  derivative with respect to the parameters, their covariance matrix and
  their names. Any weighted sum of the predictions now has a standard error,
  such as the expected number of events for a one-sample log-rank test,
  without rebuilding the package's handling of fixed parameters or of the
  Weibull scale. It matches the derivatives to the rows you hold after
  selecting, reordering or repeating rows, and refuses rows it cannot match,
  such as rows bound from another prediction. A Weibull fit's first column
  is `log(mu)`. It and `average` may change in a later release.
