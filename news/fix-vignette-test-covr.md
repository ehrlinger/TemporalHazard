* The test that compares the clinical walkthrough's stored stepwise fit with
  a refit now skips under `covr`. `covr` installs the package without its
  vignettes, so the test could not find the vignette source and failed the
  coverage run on every merge since the stored fit was added.
