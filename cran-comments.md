# CRAN submission comments -- TemporalHazard 1.3.0

## Summary

This is an update to TemporalHazard 1.1.0, accepted 2026-06-12.

This is a resubmission. The first 1.3.0 submission, on 2026-11-05, failed the
incoming checks on one NOTE, `Overall checktime 11 min > 10 min` on
r-devel-windows-x86_64. Its Windows and Debian pre-tests were otherwise
`Status: OK`. "Check time" below says what changed.

The versions between them, 1.2.1 to 1.2.19, were released internally and never
submitted, so everything since 1.1.0 is new to CRAN. `NEWS.md` opens with a
summary of the net change from 1.1.0, written for a user upgrading from it,
followed by the detailed record for each internal version. Several later
entries in that record reverse earlier ones; the summary nets them out.

## License change

The license is now GPL (>= 2). 1.1.0 was MIT.

The reason is one internal function, `.hzr_mt_normal()`. It reproduces the
numbers `set.seed()` and `rnorm()` give under R's default generators without
using R's random number stream, and its seeding, zero-word substitution and
normal inversion are derived from R's `src/main/RNG.c` and
`src/nmath/snorm.c`, which are GPL (>= 2). The package therefore carries that
license. `inst/COPYRIGHTS` lists the derived code and its copyright holders,
and `Authors@R` credits the R Core Team and the authors of the Mersenne
Twister, Makoto Matsumoto and Takuji Nishimura, as copyright holders. 1.1.0
had a single author and copyright holder, the maintainer, who made the change.

## Changes that can alter results

This is numbered as a minor release, but some of its changes can change the
result of a script written for 1.1.0. The "Upgrading from 1.1.0" section of
`NEWS.md` lists them. The ones a user is most likely to meet:

* **`hzr_stepwise()` now defaults to `criterion = "score"`,** the selection
  statistic of the SAS HAZARD program this package reproduces. A stepwise
  analysis re-run under 1.3.0 can select a different set of variables.
* **Several fixes move fitted values.** Among them: multiphase fits under
  Conservation of Events now reach maxima that 1.1.0 stopped short of;
  exponential, log-logistic and log-normal fits now honour left truncation;
  and left- and interval-censored `Surv()` responses on the formula interface
  were mistranslated in 1.1.0.
* **Some input 1.1.0 accepted is now refused or warned about,** such as
  `control` elements that 1.1.0 accepted and never applied.

The major digit is reserved for the first production release. These changes
correct behaviour that departed from the reference implementation, rather than
redesigning the interface, and each is documented in `NEWS.md` and on the
affected help page.

## Random number streams

* A multiphase fit with `n_starts > 1` draws its start offsets from
  `.hzr_mt_normal()` under `control$start_seed` (default 3). It neither reads
  nor writes `.Random.seed`, so fitting does not depend on or advance the
  user's stream.
* `hzr_bootstrap(seed =)` calls `set.seed(seed)` on a seed the user supplies
  and does not restore the stream afterwards, as its documentation says. It
  does not write `.Random.seed` into the global environment.

## Test environments

* **Local:** R 4.6.1 on macOS (aarch64-apple-darwin23). `R CMD check
  --as-cran` on a tarball built from a clean `git archive` export, with the
  PDF manual: `Status: OK`, 0 errors, 0 warnings, 0 notes. The source tarball
  is 4.07 MB.
* **win-builder,** on this tarball: TODO.
* **GitHub Actions, `R CMD check --as-cran` with the manual,** on the commit
  before the test skips: windows-latest
  (R-devel, R-release, R-oldrel-1), macos-latest (R-release) and
  ubuntu-latest (R-devel, R-release), all `Status: OK`, with
  `urlchecker::url_check()` clean.
* **GitHub Actions, on every merge:** ubuntu-latest (R-devel, R-release,
  R-oldrel-1 and R 4.1.0, the declared minimum), macos-latest (R-release) and
  windows-latest (R-release).

## Check time

1.1.0 took 528 to 623 s to check on CRAN's Windows flavors. The tests that run
under `R CMD check` had grown about fourfold since September 2026, and the
first 1.3.0 submission's r-devel-windows incoming check took 11 minutes:
228 s in the tests and 249 s rebuilding the vignettes.

For this resubmission:

* **Tests:** the 100 test files that took longest under the check carry a
  file-level `skip_on_cran()`, against 50 in the first submission. They still
  run on every GitHub Actions job and in the local suite. On macOS the tests
  under the check now take 23 s, against 42 s for the first submission's
  tarball.
* **Vignettes:** the most expensive vignette chunk, a stepwise selection in
  `clinical-analysis-walkthrough`, now shows its code and loads its stored
  result from `inst/extdata/vignette-fit-step.rds` (72 KB). A test that runs
  on GitHub Actions refits the vignette's code and compares it with the
  stored result.

TODO: win-builder timings for this tarball.

## Reverse dependencies

`tools::package_dependencies(reverse = TRUE, which = "all")` against CRAN
returns none.

## NOTEs

TODO after win-builder.

## Background

* The package ships six datasets, about 58 KB compressed in total. Five are
  clinical reference datasets used in the vignettes and in parity tests
  against the original C/SAS HAZARD program; the sixth, `uslife2023`, is the
  2023 United States life table from the National Center for Health
  Statistics.
* Vignettes use Quarto (`VignetteBuilder: quarto`).
* Slow tests, including the SAS parity tests and multi-start multiphase fits,
  carry `skip_on_cran()` to stay within CRAN's time budget. They run on every
  GitHub Actions job, and the SAS parity tests, which need the reference
  outputs, run on the ubuntu-latest R-release job.
* DOI URLs in the vignettes (doi.org) may return HTTP 403 to automated
  checkers but resolve in a browser.
