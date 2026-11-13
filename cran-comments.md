# CRAN submission comments -- TemporalHazard 1.3.1

## Summary

This is an update to TemporalHazard 1.3.0, accepted 2026-11-08.

It is a patch release, but some calls that ran under 1.3.0 now give a
different result. A multiphase fit whose shape parameter has run to a hard
limit (a phase becoming a step, or its half-time leaving the observed times)
is now refitted from a grid of starts with a warning, or refused with an
error that names the parameter. Under 1.3.0 such a fit returned the drifted
last iterate, often with no warning, so it was not an estimate. A weaker
flag, such as a poorly identified half-time, tries the same grid but is
never refused: when no start does better, the original fit is returned
without the restart warning, and the check's outcome is recorded on the fit. `NEWS.md`
opens the 1.3.1 section with these changes.

## Test environments

TO FILL before submitting, from the logs of this tarball's checks:

* **Local:** `R CMD check --as-cran` with the PDF manual, on a tarball built
  from a clean `git archive` export.
* **win-builder:** R-devel, R-release and R-oldrelease, with per-step timings.
* **GitHub Actions, `R CMD check --as-cran` with the manual:** the
  `check-release.yaml` run on `main` after this bump.
* **GitHub Actions, on every merge:** ubuntu-latest (R-devel, R-release,
  R-oldrel-1 and R 4.1.0, the declared minimum), macos-latest (R-release) and
  windows-latest (R-release).

## Check time

TO FILL from win-builder's `00check.log`. 1.3.0 took 5 min 3 s on win-builder
R-devel and was accepted.

## Reverse dependencies

TO FILL from `tools::package_dependencies(reverse = TRUE, which = "all")`.

## NOTEs

TO FILL from the win-builder checks.

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
