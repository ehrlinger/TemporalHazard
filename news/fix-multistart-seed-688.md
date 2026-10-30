- A multiphase fit no longer writes `.Random.seed` to the global environment.
  The `n_starts` offsets are computed by the package and are the same numbers
  as before, so fits under R's default random number generators are unchanged.
  Under another `RNGkind()` they now match those default-generator fits, where
  they previously followed the session's generator.
