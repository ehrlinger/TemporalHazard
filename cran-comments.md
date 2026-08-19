# CRAN submission comments -- TemporalHazard 1.2.1

## Summary

This is an update to TemporalHazard 1.1.0 (accepted 2026-06-12). There is no
reviewer feedback to address.

This is a minor-version release that nonetheless contains one breaking change:
`hzr_stepwise()` now defaults to a different selection criterion, so re-running
an existing stepwise analysis can select a different set of variables. It is
numbered as a minor rather than a major deliberately. The previous default
deviated from the C/SAS HAZARD reference this package exists to reproduce, so
the change restores intended behaviour rather than redesigning the interface;
the signature is unchanged and the old behaviour remains available exactly,
through a documented argument. The `1.x` line is the run-up to a first
production release, and the major digit is reserved for that milestone.

Everything else in this release is additive or a bug fix. Several of the fixes
correct results that were silently wrong rather than errors that were visible;
they are grouped and described below. The breaking change
is flagged under a "Breaking changes" heading in `NEWS.md` and in
`?hzr_stepwise`, so users encounter the warning regardless of the version
number.

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
  fixed-formula bootstrap unchanged. The argument is marked experimental in
  its documentation: it ships to CRAN for the first time here, and its
  interface is still being read off production runs.
* `hzr_read_outhaz()` reads a `PROC HAZARD` `outhaz=` estimate dataset, which
  holds the converged estimates and the variance-covariance matrix at full
  double precision where the printed listing carries about seven figures.
* The SAS `.lst` parsers now ship with the installed package, under
  `sas-parity/`. They previously lived under `tests/`, which `R CMD INSTALL`
  skips unless `--install-tests` is passed, so an installed package could not
  reach them.
* `hzr_bootstrap(verbose = TRUE)` shows a text progress bar via
  `utils::txtProgressBar()` instead of an every-50-replicates message.

### Bug fixes

Most of this release is bug fixes, and they are grouped below because several
share one shape: a computation that returned a result-shaped object which meant
nothing, with no error and no warning. They were found by putting the package
through its first production analysis. `NEWS.md` carries the full detail.

**Wrong results, silently.**

* The formula interface mistranslated left- and interval-censored `Surv()`
  objects. `survival::Surv()` and this package use different integer codings
  for censoring status, and the parser passed `Surv()`'s through unchanged.
  Under `type = "left"` a left-censored row was read as right-censored -- a
  wrong answer with no outward sign. Two related faults in the same code path
  made an interval-censored formula fit return the optimizer's failure
  sentinel instead of a fit. Status codes are now translated, and a regression
  test asserts that a `Surv(type = "interval")` fit reproduces the equivalent
  vector-interface fit to 1e-8 in log-likelihood.
* `hzr_bootstrap()` did not resample fits built with the vector interface. A
  vector-interface call stores `time = d$col` as an expression, so every
  replicate re-evaluated it against the original data: the run reported
  `n_success = n_boot`, `n_failed = 0`, no warning, and `n_boot` **identical**
  replicates, with `sd` exactly 0 on every parameter. Both interfaces now
  produce identical replicates for the same model, data and seed.

**Selection that could not run, and did not say so.**

* Under the new score criterion, an interval- or left-censored multiphase fit
  could not score a single candidate. The analytic observed information is not
  defined for those rows, and the score path had no numeric fallback on that
  branch, so every candidate scored `NA` and the screen stopped having tested
  nothing. The fallback now agrees with the analytic form to 1e-4 on the
  equivalent right-censored fit.
* `hzr_stepwise()` and `hzr_bootstrap(scope = ...)` failed on a formula passed
  by symbol -- `f <- Surv(t, d) ~ 1; hazard(f, ...)` -- because the stored call
  was recovered by deparsing. Every post-entry refit threw, nothing entered,
  and a bootstrap reported full success. Four sites resolved the stored formula
  independently; they now share one helper.
* A screen that could not test its candidates is now distinguishable from one
  that tested them and found nothing. `hzr_stepwise()` reports how many
  candidate scores were unavailable and **why**, and `hzr_bootstrap()`
  aggregates the same across replicates. The distinction matters in both
  directions: a collinear candidate should be dropped, while a candidate whose
  observed information is indefinite at zero is typically a strong one that
  was passed over rather than tested.
* A select-mode `hzr_bootstrap()` that selects no covariate at all now warns
  rather than returning an empty summary that reads like a completed screen.

**Diagnostics that were silent.**

* A fit that cannot compute a Hessian now says so. The analytic Hessian
  declines for left- and interval-censored rows by design, the optimizer falls
  back to `numDeriv::hessian()`, and `numDeriv` is a `Suggests` -- so on a
  machine installed without Suggests, an interval-censored multiphase fit
  produced no standard errors and nothing named the cause. Behaviour is
  unchanged; the reason is now stated.
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

* **Local:** R 4.6.1 on macOS (aarch64-apple-darwin23).
  `R CMD check --as-cran` on the built tarball, **with the PDF manual**,
  returns `Status: OK` -- 0 errors, 0 warnings, 0 notes. Overall check time is
  about three minutes, well inside the reference budget; the source tarball is
  2.7 MB.
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
