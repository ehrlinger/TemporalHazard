* A fit stopped where the Hessian is not positive-definite could have `$fit$weak`
  `NULL`, which reads as examined and well identified (#688). It is now `NA`,
  not examined, with the reason in `degraded_causes`. The fit's warnings about
  the Hessian are unchanged.
