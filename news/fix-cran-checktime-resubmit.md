- The 100 test files that took longest under `R CMD check` now skip on CRAN.
  The tests that ran under the check had grown about fourfold since
  September 2026, and a win-builder check spent 25 minutes in them. The
  skipped files still run on every GitHub Actions job and in the local
  suite. No code changes.
- The stepwise example in `vignette("clinical-analysis-walkthrough")` now
  shows its code and loads the stored result of running it, instead of
  refitting during the vignette build. The refit was the most expensive
  step of the build. A test refits the vignette's code and compares, so a
  stored result the code no longer produces fails the test suite.
