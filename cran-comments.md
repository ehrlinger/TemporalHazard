# CRAN submission comments -- TemporalHazard 1.2.0

## Summary

This is an update to TemporalHazard 1.1.0 (accepted 2026-06-12).
No reviewer feedback to address. New version adds embedded stepwise
variable selection to `hzr_bootstrap()`.

## Changes since 1.1.0

* **`hzr_bootstrap()` gains a `scope` argument for embedded stepwise
  variable selection.** This is the R equivalent of SAS's `%HAZBOOT`
  procedure: each bootstrap replicate runs a fresh `hzr_stepwise()`
  selection (starting from a fixed-shape refit of the base model) instead
  of a plain refit, so `summary$pct` reports the variable's selection
  frequency across resamples and `summary$mean`/`sd`/`ci_*` describe the
  coefficient distribution conditional on selection. `scope = NULL`
  (the default) preserves the original fixed-formula bootstrap unchanged.

## Test environments

* **Local:** R 4.6.0 on macOS (aarch64-apple-darwin23).
  `R CMD check --as-cran` (with PDF manual) returns 0 errors, 0 warnings,
  0 notes.
* **GitHub Actions matrix:** ubuntu-latest (R-devel / R-release /
  R-oldrel-1), macos-latest (R-release), windows-latest (R-release).
  All checks passed (0 errors / 0 warnings).
* **Reverse-dependency check:** `tools::package_dependencies(reverse = TRUE)`
  returns 0.
* **`urlchecker::url_check()`:** all URLs correct.

## NOTE disposition

`R CMD check --as-cran` on the release tarball is clean (0/0/0).
CRAN's incoming-feasibility / aspell check may report:

* **Possibly misspelled words in DESCRIPTION:** `Naftel`, `Rajeswaran`
  (cited authors with `<doi:...>` references), `UAB` (acronym), `et`,
  `al` ("et al." citation), `multiphase` (the core domain term). None are
  misspellings. `inst/WORDLIST` lists these for the `spelling` unit test;
  the CRAN aspell check has no package-level suppression mechanism, so
  this NOTE recurs by design on every submission.

## Background

* Package contains five clinical reference datasets (~120 KB compressed)
  used in vignettes and parity tests against the original C/SAS HAZARD
  program.
* Vignettes use Quarto (`VignetteBuilder: quarto`).
* A subset of tests (SAS parity, multi-start multiphase fits, OMC
  raw-file derivation) are guarded with `skip_on_cran()` to stay within
  CRAN runtime budgets; the full suite runs on GitHub Actions.
* DOI URLs in vignettes (doi.org) may return HTTP 403 to automated
  checkers but resolve correctly in browsers.
