* **Multiphase fits are about twice as fast again, with every value
  unchanged (#631).** Since #578, each evaluation of a phase formed
  `log(1 - G)` and the hazard on the log scale, which only a `"hazard"`
  phase reads, and formed some of it four times over. A `"cdf"` phase now
  skips that work, and a `"hazard"` phase forms each piece once. On the
  `cabgkul` data a fixed-shape three-phase fit takes 1.0 seconds against
  1.75, and the package's tests run in about half the time. Estimates,
  standard errors, predictions and bootstrap replicates are bit-identical
  to 1.2.13's.
