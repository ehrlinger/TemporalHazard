# Feature test matrix

This page lists each TemporalHazard feature, the tests that cover it,
and whether it is compared with an external reference. It is generated
from one saved run of the full test suite by `tools/feature-matrix.R`,
which refuses to write it unless the run has no failure and its recorded
commit, installed dependencies and test names all match the tree.

**Measured at 6b1588e5 on 2026-10-30, version 1.2.17.** That run: 2688
tests run, 15800 expectations, 12 tests skipped. A page older than the
code it describes is still possible: regenerate it at each release.

## How to read it

- **Status** is a judgment, set by hand in
  `tools/feature-matrix-map.csv`. *Tested*: covered by the package’s own
  tests. *In progress*: a known gap, named in the notes. *Shipped* is
  reserved for a feature handed to biostatistics for critique.
- **Tests run** and **expectations** count what ran in that run, across
  the feature’s test files and the tests cited for it from other files.
  A skipped test is not counted as run. A test file can count towards
  more than one feature.
- **Skipped** counts tests that did not run there, usually because a
  reference file, the C binary or the study data was not available.
- **Reference** counts the tests that compare the feature with something
  outside the package: SAS listings, the C `hazard` binary, or an
  independent oracle. Each is cited by name in
  `tools/feature-matrix-references.csv`. “*n* of *m* ran” means *n* of
  those *m* comparisons ran in this run. The list is selected, not
  exhaustive: “none cited” means no comparison has been cited for that
  feature, not that none exists. Tests that only check a reference file
  is well formed, or that exercise the C binary’s interface, are not
  counted.
- **Default start** says whether a test shows a fit started from the
  package’s defaults reaching the reference estimates. A reference test
  that starts from the estimates SAS reported shows that the likelihood
  agrees at that point, not that
  [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
  finds it on its own.

## Model

| Feature | Status | Tests run | Expectations | Skipped | Reference | Default start | Notes |
|:---|:---|---:|---:|---:|:---|:---|:---|
| Exponential distribution | tested | 32 | 67 | 0 | 2 of 2 ran (oracle) | none |  |
| Weibull distribution | tested | 68 | 303 | 1 | none cited | none | The C-binary file tests the binary interface (command, input and output parsing), not estimates. |
| Log-normal distribution | tested | 35 | 77 | 0 | 2 of 2 ran (oracle) | none |  |
| Log-logistic distribution | tested | 36 | 77 | 0 | 2 of 2 ran (oracle) | none |  |
| Multiphase likelihood and optimizer | tested | 148 | 775 | 1 | 8 of 9 ran (C binary, SAS) | none | No criterion yet separates a genuine multiphase optimum from a stuck fit (#589). |
| Phase specification and shapes (hzr_phase) | tested | 125 | 862 | 0 | 2 of 2 ran (oracle) | none |  |
| Phase identifiability and weak directions | tested | 85 | 275 | 0 | none cited | none |  |
| Conservation of Events | tested | 42 | 238 | 0 | none cited | cited test passed | The default-start test fixes the shapes at SAS’s values, so only the coefficients start from defaults. |

## Data

| Feature | Status | Tests run | Expectations | Skipped | Reference | Default start | Notes |
|:---|:---|---:|---:|---:|:---|:---|:---|
| Right censoring and data handling | tested | 111 | 507 | 0 | none cited | none |  |
| Interval and left censoring | tested | 112 | 369 | 0 | none cited | none |  |
| Entry times (left truncation) | tested | 43 | 302 | 0 | 7 of 7 ran (oracle, SAS) | none |  |
| Observation weights | in progress | 56 | 182 | 2 | 1 of 3 ran (SAS) | none | SAS parity covers fractional weights; weights beside repeated events are untested (#681). |
| Repeated events | in progress | 74 | 228 | 6 | 0 of 4 ran (SAS) | none | Parity runs only with the study volume; the SAS targets are not yet set (#681). |
| Time-varying coefficients (time_windows) | tested | 23 | 107 | 0 | none cited | none |  |
| Time-varying covariates (start/stop epochs) | in progress | 6 | 14 | 0 | none cited | none | Named as under-tested at the 2026-10-29 review. The epoch tests check that splitting a row into epochs leaves the log-likelihood unchanged. Production time-varying covariate work is under way in the Hsich and COMMENCE studies; none of it is in the test suite yet. |

## Covariates

| Feature | Status | Tests run | Expectations | Skipped | Reference | Default start | Notes |
|:---|:---|---:|---:|---:|:---|:---|:---|
| Formula and vector interfaces | tested | 163 | 640 | 0 | none cited | none |  |
| Phase-specific covariates | tested | 44 | 201 | 0 | none cited | none |  |

## Inference

| Feature | Status | Tests run | Expectations | Skipped | Reference | Default start | Notes |
|:---|:---|---:|---:|---:|:---|:---|:---|
| Standard errors and the Hessian | tested | 92 | 268 | 0 | 1 of 1 ran (C binary) | none |  |
| Convergence and degraded-fit reporting | in progress | 86 | 498 | 0 | none cited | none | Single-distribution false maxima warn (#531); the multiphase criterion is open (#589). |
| Prediction (predict.hazard) | tested | 256 | 1816 | 0 | 3 of 3 ran (SAS) | none |  |
| Decomposition and evaluation | tested | 83 | 436 | 0 | none cited | none |  |
| Goodness of fit and deciles | tested | 135 | 531 | 0 | 3 of 3 ran (SAS) | none |  |

## Selection

| Feature | Status | Tests run | Expectations | Skipped | Reference | Default start | Notes |
|:---|:---|---:|---:|---:|:---|:---|:---|
| Stepwise selection (hzr_stepwise) | tested | 296 | 1284 | 0 | 2 of 2 ran (SAS) | none |  |
| Score and Wald tests | tested | 96 | 503 | 0 | none cited | none |  |
| Bootstrap (hzr_bootstrap) | tested | 85 | 632 | 0 | none cited | none |  |

## SAS migration

| Feature | Status | Tests run | Expectations | Skipped | Reference | Default start | Notes |
|:---|:---|---:|---:|---:|:---|:---|:---|
| SAS job translator (hzr_translate_sas) | tested | 486 | 5196 | 0 | none cited | none |  |
| SAS output readers (hzr_read_outhaz and listings) | tested | 81 | 249 | 3 | none cited | none |  |

## Regenerating

From the package root, on a committed tree:

``` sh
Rscript tools/save-suite.R results.rds
Rscript tools/feature-matrix.R results.rds
```

`tools/release-suite.R results.rds` makes the same run with the study
data required, and its output feeds the second command too. Both attach
the run’s provenance: the commit, the installed dependency versions and
the time.

The second command refuses a run made on an uncommitted tree, a run of a
commit whose tests, code or data differ from the current one, a run
under other dependency versions, a run with a failure, and a run whose
test names differ from those on disk. It also refuses a test that no
feature counts, a pattern in the map that matches no file, and a cited
test that is not in the run. Each saving script runs the suite on a
private export of the commit, so nothing done to the working tree during
the run can reach it, and refuses to start with uncommitted changes
under the paths a test reads. It also refuses a run during which the
installed package versions changed. Files outside the repository, such
as the SAS fixtures and the study data, are not part of the record.

`Rscript tools/test-feature-matrix.R` plants each of those faults and
checks the refusal, and the `lint` workflow runs it on every pull
request.
