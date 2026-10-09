* A select-mode `hzr_bootstrap()` replicate whose base refit did not converge
  is now counted as failed, with the reason "base refit did not converge
  (converged = FALSE)". The rule that a replicate which does not converge
  fails applied only to the final fit, so under `criterion = "aic"` each
  candidate warm-started from the base's unfinished point, reached the real
  maximum and was credited with the base's shortfall. Under a forward Wald
  screen a replicate that entered nothing failed while one that entered
  something passed, so the replicates kept were biased towards selection. On 300 Weibull rows
  with two pure-noise columns, a base that converged on the full data but not
  on any resample gave 19 of 30 successful replicates with no warning, and
  selected the noise columns in 12 and 8 of them.
* `hzr_stepwise(criterion = "aic")` now refuses a base fit that did not
  converge, as `criterion = "score"` already did. A candidate's delta AIC is
  measured against the base's log-likelihood, so a base that stopped short
  entered pure noise at delta AIC -16.5 with no warning. `criterion = "wald"`
  refuses one too unless `direction = "forward"`: a Wald entry is tested at
  the candidate's own converged refit, but a removal is tested on the base's
  own estimates and variance.
* `hzr_deciles()` now refuses a multiphase fit that dropped the rows where a
  phase covariate was missing, with the same check and message as
  `hzr_gof()`. It predicted each subject's cumulative hazard from the shorter
  design, recycled onto the wrong subjects, and reported every subject
  included: on `avc` with half the covariate missing, 70 observed events
  against 159 expected and p < 2e-16, with no warning.
