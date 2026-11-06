* The test that compares the stored stepwise fit in
  `vignette("clinical-analysis-walkthrough")` with a refit now skips under
  `covr`. `covr` installs the package without its vignettes, so the test
  could not find the vignette source and failed the coverage run on every
  merge since the stored fit was added.
