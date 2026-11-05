# CRAN submission comments -- TemporalHazard 1.3.0

## Summary

This is an update to TemporalHazard 1.1.0, accepted 2026-06-12. There is no
reviewer feedback to address.

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
  is 3.98 MB.
* **GitHub Actions:** ubuntu-latest (R-devel, R-release, R-oldrel-1 and R
  4.1.0, the declared minimum), macos-latest (R-release) and windows-latest
  (R-release). <!-- TODO: confirm all pass on the release commit. macOS and
  Windows check only after a merge to main, not on a pull request. -->
* <!-- TODO: win-devel (devtools::check_win_devel()) and the check-release
  workflow, with their results. -->

## Check time

<!-- TODO: replace with the idle measurement. A single local run on
2026-11-04, on a loaded machine, took 438 s overall (tests 291 s, vignette
rebuild 57 s). The 1.2.13 cut measured about 12 minutes; that slowdown (#631)
has since been fixed, and the idle runs are what confirm it. -->

## Reverse dependencies

`tools::package_dependencies(reverse = TRUE, which = "all")` against CRAN
returns none.

## URLs

<!-- TODO: urlchecker::url_check() result. The 2026-11-05 run timed out on
every CRAN.R-project.org URL at 5 seconds, which says nothing about the
URLs. -->

## Expected NOTE

The CRAN incoming aspell check is expected to flag these words in
`DESCRIPTION`, none of which is a misspelling: `Naftel` and `Rajeswaran`
(cited authors, each with a `<doi:...>` reference), `UAB` (an acronym), `et`
and `al` (from "et al."), and `multiphase` (the model's standard name in the
literature it implements). `inst/WORDLIST` records them for
`spelling::spell_check_package()`, which has no effect on CRAN's check.

<!-- TODO: confirm the list against win-devel's output. -->

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
