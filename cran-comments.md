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
* **win-builder:** R-devel (2026-11-04 r90836) and R-release (4.6.2),
  Status OK. <!-- TODO: results for the resubmitted tarball, and the
  check-release workflow. -->

## Check time

1.1.0 already took 528 to 623 s to check on CRAN's Windows flavors. A
win-builder run of the first 1.3.0 tarball, on 2026-11-05, gave Status OK on
R-devel and R-release but spent 25 minutes in the tests on each, because the
suite had grown about fourfold. 50 test files, the slowest under `R CMD
check`, now carry a file-level `skip_on_cran()`. They run on every GitHub
Actions job and in the local suite.

<!-- TODO: the win-builder R-devel and R-release times for the resubmitted
tarball, and the local idle --as-cran time. -->

## Reverse dependencies

`tools::package_dependencies(reverse = TRUE, which = "all")` against CRAN
returns none.

## URLs

<!-- TODO: urlchecker::url_check() result. The 2026-11-05 run timed out on
every CRAN.R-project.org URL at 5 seconds, which says nothing about the
URLs. -->

## NOTEs

None. The win-builder R-devel and R-release checks of 2026-11-05 gave 0
errors, 0 warnings and 0 notes, incoming feasibility included.

<!-- TODO: confirm on the resubmitted tarball. The win-builder R-oldrelease
(4.5.3) run stopped twice, with no message, inside "checking CRAN incoming
feasibility"; the same check under R 4.5.3 locally completes with no
finding. -->

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
