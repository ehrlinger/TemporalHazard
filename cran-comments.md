# CRAN submission comments -- TemporalHazard 1.2.0

## Summary

This is an update to TemporalHazard 1.1.0 (accepted 2026-06-12). There is no
reviewer feedback to address.

The major version bump reflects one breaking change: `hzr_stepwise()` now
defaults to a different selection criterion, so re-running an existing stepwise
analysis can select a different set of variables. Everything else in this
release is additive or a bug fix, and the previous behaviour remains available
exactly, through a documented argument.

## Changes since 1.1.0

### Breaking

* **`hzr_stepwise()` now defaults to `criterion = "score"`.** This reproduces
  the `SELECTION` statistic of the original C/SAS HAZARD program, which this
  package exists to reproduce. The previous default, `"wald"`, refit the model
  once per candidate variable and tested the refit's Wald chi-square -- a
  deviation from the reference implementation. Because the score and Wald paths
  take different step sequences, re-running an existing analysis can now select
  a different variable set. `criterion = "wald"` restores the old behaviour
  exactly. Both the change and the escape hatch are documented in
  `?hzr_stepwise` and in the "Breaking changes" section of `NEWS.md`.

  Removing the per-candidate refit also removed the dominant runtime cost: a
  92-variable two-phase screen fell from roughly 25 minutes per bootstrap
  replicate to seconds.

### New features

* `hzr_bootstrap()` gains a `scope` argument for embedded stepwise variable
  selection within each bootstrap replicate -- the R equivalent of SAS's
  `%HAZBOOT` procedure. `scope = NULL` (the default) preserves the original
  fixed-formula bootstrap unchanged.
* `hzr_bootstrap(verbose = TRUE)` shows a text progress bar via
  `utils::txtProgressBar()` instead of an every-50-replicates message.

### Bug fixes

* `hzr_bootstrap()` no longer fails silently (returning `n_success = 0` with no
  error and no warning) when the model was fitted inside a function. `hazard()`
  now records the environment its call was written in, so each replicate's refit
  resolves arguments passed by symbol against the caller's locals rather than
  against the package namespace.
* Multiphase models whose shape sits exactly at the `m = 0` or `nu = 0`
  limiting-case boundary no longer lose their analytic Hessian and fall back
  silently to a numerical one, or to `NA` standard errors.
* Multiphase fits with a single free parameter now use the analytic Hessian for
  standard errors, instead of falling back to a numerical one because a 1x1
  matrix had been dropped to a scalar.
* `hzr_bootstrap(scope = ..., trace = ...)` no longer errors with "formal
  argument matched by multiple actual arguments".
* `hzr_bootstrap()` no longer floods the console with per-replicate numerical
  warnings. Structural problems still surface once, up front.

`NEWS.md` carries the full detail for each item.

## Test environments

* **Local:** R 4.6.0 on macOS (aarch64-apple-darwin23).
  `R CMD check --as-cran` on the built tarball, **with the PDF manual**,
  returns 0 errors, 0 warnings, 0 notes. Overall check time is about two
  minutes, well inside the reference budget; the source tarball is 2.7 MB.
* **GitHub Actions matrix:** ubuntu-latest (R-devel / R-release / R-oldrel-1),
  macos-latest (R-release), windows-latest (R-release).
* **Reverse-dependency check:** `tools::package_dependencies(reverse = TRUE)`
  returns 0 reverse dependencies.
* **`urlchecker::url_check()`:** all 19 URLs correct.

## NOTE disposition

`R CMD check --as-cran` on the release tarball is clean (0/0/0) locally. The
CRAN incoming-feasibility aspell check is expected to report:

* **Possibly misspelled words in DESCRIPTION:** `Naftel`, `Rajeswaran`
  (cited authors, each with an accompanying `<doi:...>` reference), `UAB`
  (an acronym), `et`, `al` (from the "et al." citation), and `multiphase`
  (the core domain term, standard in the literature this package
  implements). None are misspellings.

  `inst/WORDLIST` records these for `devtools::spell_check()`. It has no
  effect on the CRAN aspell check, which offers no package-level suppression
  mechanism, so this NOTE recurs by design on every submission.

## Background

* The package contains five clinical reference datasets (~120 KB compressed)
  used in vignettes and in parity tests against the original C/SAS HAZARD
  program.
* Vignettes use Quarto (`VignetteBuilder: quarto`).
* A subset of tests -- SAS parity, multi-start multiphase fits, OMC raw-file
  derivation -- is guarded with `skip_on_cran()` to stay within CRAN runtime
  budgets. The full suite, including the SAS parity gate, runs on GitHub
  Actions.
* DOI URLs in vignettes (doi.org) may return HTTP 403 to automated checkers
  but resolve correctly in browsers.
