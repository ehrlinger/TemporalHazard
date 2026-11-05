- The 50 test files that took longest under `R CMD check` now skip on CRAN.
  The tests that ran under the check had grown about fourfold since 1.1.0,
  and a win-builder check spent 25 minutes in them. The skipped files still
  run on every GitHub Actions job and in the local suite. No code changes.
