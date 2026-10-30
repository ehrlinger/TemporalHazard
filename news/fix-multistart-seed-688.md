- A multiphase fit no longer writes `.Random.seed` to the global environment.
  The offsets for `n_starts` come from a small generator of the package's own
  rather than from R's stream saved and restored around `rnorm()`, so a fit also
  no longer depends on `RNGkind()`. **The offsets are new, and no `start_seed`
  recovers the old ones, so a fit with `n_starts > 1` can stop at a different
  point or a different optimum than it did before.** The default `start_seed`
  moves from 3 to 1, because seed 3 of the new generator stopped short on a
  reference `avc` fit.
