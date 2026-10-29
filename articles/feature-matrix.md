# Feature test matrix

This page lists each TemporalHazard feature, the tests that cover it,
and whether it is compared with an external reference. It is generated
from one saved run of the full test suite by `tools/feature-matrix.R`,
which refuses to write it unless the run is of the current tree and has
no failure.

**Measured at 016d83c4 on 2026-10-29, version 1.2.17.** That run: 2700
tests, 15800 expectations, 12 skipped. A page older than the code it
describes is still possible: regenerate it at each release.

## How to read it

- **Status** is a judgment, set by hand in
  `tools/feature-matrix-map.csv`. *Tested*: covered by the package’s own
  tests. *In progress*: a known gap, named in the notes. *Shipped* is
  reserved for a feature handed to biostatistics for critique.
- **Tests** and **expectations** count what ran in that run, across the
  feature’s test files. A test file can count towards more than one
  feature.
- **Skipped** counts tests that did not run there, usually because a
  reference file, the C binary or the study data was not available.
- **Reference** says whether the feature is compared with something
  outside the package: SAS listings, the C `hazard` binary, or an
  independent oracle. “*n* of *m* tests ran” is how many of those
  comparisons actually ran in this run, not how many exist.
- **Default start** says whether a test shows a fit started from the
  package’s defaults reaching the reference estimates. A reference test
  that starts from the estimates SAS reported shows that the likelihood
  agrees at that point, not that
  [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
  finds it on its own.

## Model

| Feature | Status | Tests | Expectations | Skipped | Reference | Default start | Notes |
|:---|:---|---:|---:|---:|:---|:---|:---|
| Exponential distribution | tested | 32 | 67 | 0 | none | none |  |
| Weibull distribution | tested | 69 | 303 | 1 | 1 file(s); 6 of 7 tests ran | none | The C-binary comparison skips where the binary is not built. |
| Log-normal distribution | tested | 35 | 77 | 0 | none | none |  |
| Log-logistic distribution | tested | 36 | 77 | 0 | none | none |  |
| Multiphase likelihood and optimizer | tested | 163 | 880 | 3 | 3 file(s); 28 of 31 tests ran | none | No criterion yet separates a genuine multiphase optimum from a stuck fit (#589). |
| Phase specification and shapes (hzr_phase) | tested | 125 | 862 | 0 | 1 file(s); 2 of 2 tests ran | none |  |
| Phase identifiability and weak directions | tested | 85 | 275 | 0 | none | none |  |
| Conservation of Events | tested | 42 | 238 | 0 | none | cited test passed | The default-start test fixes the shapes at SAS’s values, so only the coefficients start from defaults. |

## Data

| Feature | Status | Tests | Expectations | Skipped | Reference | Default start | Notes |
|:---|:---|---:|---:|---:|:---|:---|:---|
| Right censoring and data handling | tested | 111 | 507 | 0 | none | none |  |
| Interval and left censoring | tested | 112 | 369 | 0 | none | none |  |
| Entry times (left truncation) | tested | 42 | 290 | 0 | none | none |  |
| Observation weights | in progress | 57 | 173 | 2 | 1 file(s); 2 of 4 tests ran | none | SAS parity covers fractional weights; weights beside repeated events are untested (#681). |
| Repeated events | in progress | 80 | 228 | 6 | 2 file(s); 0 of 6 tests ran | none | Parity runs only with the study volume; the SAS targets are not yet set (#681). |
| Time-varying coefficients (time_windows) | tested | 23 | 107 | 0 | none | none |  |
| Time-varying covariates (start/stop epochs) | in progress | 6 | 14 | 0 | none | none | Named as under-tested at the 2026-10-29 review. The epoch tests check that splitting a row into epochs leaves the log-likelihood unchanged. Production time-varying covariate work is under way in the Hsich and COMMENCE studies; none of it is in the test suite yet. |

## Covariates

| Feature | Status | Tests | Expectations | Skipped | Reference | Default start | Notes |
|:---|:---|---:|---:|---:|:---|:---|:---|
| Formula and vector interfaces | tested | 163 | 640 | 0 | none | none |  |
| Phase-specific covariates | tested | 44 | 201 | 0 | none | none |  |

## Inference

| Feature | Status | Tests | Expectations | Skipped | Reference | Default start | Notes |
|:---|:---|---:|---:|---:|:---|:---|:---|
| Standard errors and the Hessian | tested | 91 | 267 | 0 | none | none |  |
| Convergence and degraded-fit reporting | in progress | 86 | 498 | 0 | none | none | Single-distribution false maxima warn (#531); the multiphase criterion is open (#589). |
| Prediction (predict.hazard) | tested | 253 | 1789 | 0 | none | none |  |
| Decomposition and evaluation | tested | 83 | 436 | 0 | none | none |  |
| Goodness of fit and deciles | tested | 132 | 498 | 0 | none | none |  |

## Selection

| Feature | Status | Tests | Expectations | Skipped | Reference | Default start | Notes |
|:---|:---|---:|---:|---:|:---|:---|:---|
| Stepwise selection (hzr_stepwise) | tested | 298 | 1278 | 3 | 2 file(s); 5 of 8 tests ran | none |  |
| Score and Wald tests | tested | 96 | 503 | 0 | none | none |  |
| Bootstrap (hzr_bootstrap) | tested | 82 | 614 | 0 | none | none |  |

## SAS migration

| Feature | Status | Tests | Expectations | Skipped | Reference | Default start | Notes |
|:---|:---|---:|---:|---:|:---|:---|:---|
| SAS job translator (hzr_translate_sas) | tested | 486 | 5196 | 0 | none | none |  |
| SAS output readers (hzr_read_outhaz and listings) | tested | 84 | 249 | 3 | none | none |  |

## Regenerating

From the package root, after a full run saved as `tools/release-suite.R`
saves it:

``` sh
Rscript tools/release-suite.R results.rds
Rscript tools/feature-matrix.R results.rds
```

The second command refuses if a test file belongs to no feature, if a
pattern in the map matches no file, or if the run does not match the
tests on disk.
