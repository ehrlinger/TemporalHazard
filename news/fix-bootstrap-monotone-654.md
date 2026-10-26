* **`hzr_bootstrap()` pooled replicates with no finite maximum (#654).** A
  resample with separation, such as no events at one level of a factor, has
  a likelihood that keeps rising as a coefficient, the intercept or a
  phase's `log_mu` runs to infinity. The refit stops where the slope falls
  below the optimizer's tolerance and reports convergence with a small
  gradient, so it passed every check and its estimate was pooled. In one
  example 12 of 40 replicates ran out to between -7 and -15 against -1.4 on
  the full data, and the coefficient's bootstrap sd was 5.0 where the other
  28 give 0.51, with no warning. Each replicate is now tested by stepping
  further out from its estimate: along each coefficient, and along the
  direction its standard errors, or failing those its Hessian, say is
  flattest, which is the one a separation along a combination of
  parameters takes. When the log-likelihood does not fall, the replicate
  fails with the reason
  `"likelihood flat or still rising (no unique finite maximum)"`, it leaves
  the summary, and a warning gives the count. A replicate whose
  design is rank-deficient, such as a factor level absent from the resample
  or a column that is an exact combination of others, is dropped the same
  way: its estimates are not determined by the data. The fit being
  bootstrapped is tested the same way first, and warned about if it fails,
  as a multiphase phase with covariates but few events of its own can.
