* **`hzr_bootstrap()` dropped every replicate of a fit with an offset
  covariate (#661 follow-up).** The check added for #654 sized its step
  along a combined direction on the covariates alone. With a covariate far
  from zero (say `x + 20000`), the intercept's move nearly cancels the
  covariate's, so the step shrank to nothing and a genuine maximum's small
  fall read as a flat likelihood: all 3 of 3 replicates were dropped as
  having no finite maximum, against 3 kept when the covariate is centred.
  The step is now sized on the whole linear predictor, intercept included,
  and for the Weibull the intercept `log(mu)` is weighted by the shape `nu`,
  as it enters the cumulative hazard.
