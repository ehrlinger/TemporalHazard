* **`hzr_bootstrap(seed = )` now sets the session's random number stream, as
  `set.seed()` does, and no longer restores it on exit (#653).** Restoring the
  caller's stream (#611) meant writing `.Random.seed` into the global
  environment, which CRAN policy does not allow. The replicates for a given
  seed are unchanged. What changes is the stream after the call: it continues
  from where the bootstrap's draws left it, so a script that seeded once at
  the top and draws after a seeded bootstrap now gets different numbers from
  those draws. Save `.Random.seed` before the call and restore it after, or use
  `withr::with_preserve_seed()`, to keep your own stream. `withr` moves from
  Imports to Suggests.
