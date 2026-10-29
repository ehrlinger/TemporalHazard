* `fit$fit$weak` now names a multiphase phase's `log_mu`, `log_t_half` or
  `log_tau` when the data do not determine it on its own: the flattest
  direction of the covariance loads on it alone and its standard error exceeds
  10 log units. Before, `weak` was `NULL`, which means "checked, identified",
  on a fit whose Hessian was ill-conditioned and whose conserved constant phase
  had been switched off at the boundary (`log_mu` near -28, standard error
  3e4 to 1e5). The single-parameter reading (#415) named only g3 shapes. Fits
  whose `weak` was `NULL` can now carry a list and raise the "weakly identified
  fit" warning. Covariate coefficients are still never named on their own, and
  neither is a log-scale name that a covariate duplicates (#673).
