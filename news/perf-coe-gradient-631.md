* **Multiphase fits under Conservation of Events are faster again, with
  every value unchanged (#631).** The gradient under CoE evaluates the score
  twice at each point (#565), and the second evaluation now reuses the
  first's per-phase quantities instead of recomputing them.
  `hzr_log1pexp()` forms each branch only where it is taken, where it used
  `ifelse()`, which formed both over every element. A bootstrap of a
  three-phase fit runs about 15 to 25% faster than in the previous release.
  Estimates, standard errors and bootstrap replicates are bit-identical.
  `hzr_log1pexp()` returns the same values, attributes and NAs as before
  on every vector, matrix and array; given a data frame, which it does not
  take, it now stops instead of returning a malformed list.
