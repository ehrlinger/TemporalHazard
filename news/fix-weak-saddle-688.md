* A fit stopped at a saddle with a well-conditioned Hessian had `$fit$weak`
  `NULL`, which reads as examined and well identified (#688). It is now `NA`,
  not examined, with the reason in `degraded_causes`. The fit's warnings about
  the Hessian are unchanged.
