# Package Architecture

## 1 Overview

TemporalHazard is a pure-R reimplementation of the C/SAS HAZARD
procedure originally developed at the University of Alabama at
Birmingham (UAB) for multi-phase parametric hazard modeling (Blackstone,
Naftel, and Turner 1986). The SAS/C code and this R package are
currently developed and maintained at The Cleveland Clinic Foundation,
and the R code was wholly developed at The Cleveland Clinic Foundation.
The package provides a unified framework for fitting additive hazard
models with an arbitrary number of temporal phases, each governed by the
three-parameter `decompos(t; t_half, nu, m)` family. The generalized
temporal decomposition extends naturally to longitudinal mixed-effects
settings (Rajeswaran et al. 2018).

This vignette documents the internal architecture: how source files are
organized, how functions compose into the fitting pipeline, how golden
fixtures ensure regression-free development, and what reference datasets
ship with the package.

## 2 Source file organization

The R source lives in `R/` and follows a layered design. Lower layers
know nothing about higher layers; control flows downward from the API
through distribution-specific optimizers to shared numerical primitives.

| Layer | File | Purpose |
|:---|:---|:---|
| User API | hazard_api.R | hazard(), predict.hazard(), print/summary/coef/vcov S3 methods |
| User API | argument_mapping.R | SAS HAZARD -\> R parameter translation table |
| Multiphase engine | likelihood-multiphase.R | Additive N-phase likelihood, cumhaz/hazard evaluators, multi-start optimizer |
| Multiphase engine | phase-spec.R | hzr_phase() constructor, validators, starting-value assembly |
| Multiphase engine | decomposition.R | Unified decompos(t; t_half, nu, m) parametric family |
| Single-phase likelihoods | likelihood-weibull.R | Weibull PH likelihood, gradient, optimizer |
| Single-phase likelihoods | likelihood-exponential.R | Exponential (constant hazard) likelihood |
| Single-phase likelihoods | likelihood-loglogistic.R | Log-logistic (proportional odds) likelihood |
| Single-phase likelihoods | likelihood-lognormal.R | Log-normal (AFT) likelihood |
| Shared infrastructure | optimizer.R | Generic L-BFGS-B/BFGS optimizer with Hessian-based vcov |
| Shared infrastructure | math_primitives.R | Numerically stable log1pexp, log1mexp, clamp_prob |
| Shared infrastructure | formula-helpers.R | Surv() formula parsing for right/left/interval censoring |
| Shared infrastructure | golden_fixtures.R | Synthetic fixture generators (.rds reference outputs) |
| Shared infrastructure | parity-helpers.R | Stubs for cross-validating against C HAZARD binary |

R source files by architectural layer {.table .caption-top}

## 3 Function call graph

The primary user entry point is
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md),
which dispatches to distribution-specific optimizers. The diagram below
shows the call flow for a multiphase fit.

    Call graph for a multiphase hazard fit
    --------------------------------------

    hazard(dist='multiphase')
    +-- .hzr_optim_multiphase()
    |   +-- Resolve per-phase design matrices
    |   +-- .hzr_phase_start() per phase
    |   +-- .hzr_optim_generic()
    |       +-- stats::optim(method='BFGS')
    |           +-- .hzr_logl_multiphase()
    |               +-- .hzr_split_theta()
    |               +-- .hzr_multiphase_cumhaz()
    |               |   +-- hzr_phase_cumhaz() / hzr_decompos_g3()
    |               +-- .hzr_multiphase_hazard()
    |                   +-- hzr_phase_hazard() / hzr_decompos_g3()
    +-- predict.hazard()
        +-- .hzr_multiphase_cumhaz()
        +-- .hzr_multiphase_hazard()

For single-phase distributions (Weibull, exponential, log-logistic,
log-normal),
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
dispatches to the corresponding `.hzr_optim_<dist>()` function, which
calls `.hzr_optim_generic()` with the distribution-specific
log-likelihood and gradient.

### 3.1 Key internal functions

#### 3.1.1 Theta vector layout

The multiphase model concatenates per-phase sub-vectors into a single
flat `theta` on the **internal (estimation) scale**:

| Phase.type | Sub.vector | Count |
|:---|:---|:---|
| cdf / hazard | \[log_mu, log_t_half, nu, m, beta_1, …, beta_p\] | 4 + p |
| g3 | \[log_mu, log_tau, gamma, alpha, eta, beta_1, …, beta_p\] | 5 + p |
| constant | \[log_mu, beta_1, …, beta_p\] | 1 + p |

Internal theta layout per phase {.table .caption-top}

[`.hzr_split_theta()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_split_theta.md)
partitions the concatenated vector into a named list of per-phase
sub-vectors.
[`.hzr_unpack_phase_theta()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_unpack_phase_theta.md)
extracts named parameters (`log_mu`, `log_t_half`, `nu`, `m`, `beta`)
from each sub-vector.

#### 3.1.2 Phase specification

`hzr_phase(type, t_half, nu, m, formula)` creates a lightweight S3
object that stores the phase type, initial shape parameters, and an
optional phase-specific covariate formula. The constructor validates
inputs and returns `NA` for shape parameters of constant phases.

A list of `hzr_phase` objects is validated by
[`.hzr_validate_phases()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_validate_phases.md),
which auto-names unnamed phases (“phase_1”, “phase_2”, …) and catches
the common mistake of passing a bare `hzr_phase` instead of a list.

#### 3.1.3 Decomposition engine

`hzr_decompos(time, t_half, nu, m)` is the mathematical core. It
computes `G(t)` (CDF), `g(t)` (density), and `h(t)` (hazard) for the
three-parameter family across six valid sign combinations of `nu` and
`m`. The phase-level helpers
[`hzr_phase_cumhaz()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase_cumhaz.md)
and
[`hzr_phase_hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase_hazard.md)
wrap
[`hzr_decompos()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos.md)
and apply the phase type mapping:

| Phase type   | \\\Phi(t)\\                            | \\\phi(t)\\ |
|:-------------|:---------------------------------------|:------------|
| `"cdf"`      | \\G_1(t)\\                             | \\g_1(t)\\  |
| `"hazard"`   | \\-\log(1-G_1(t))\\                    | \\h_1(t)\\  |
| `"g3"`       | \\G_3(t; \tau, \gamma, \alpha, \eta)\\ | \\g_3(t)\\  |
| `"constant"` | \\t\\                                  | \\1\\       |

#### 3.1.4 Multi-start optimization

[`.hzr_optim_multiphase()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_optim_multiphase.md)
runs a multi-start strategy: the first start uses the assembled starting
values (from `theta` or `hzr_phase` specs); subsequent starts perturb
randomly (sd = 0.5). The best log-likelihood across all starts is kept.
This helps the optimizer escape the shallow local optima that are common
in multiphase models.

The optimizer delegates to `.hzr_optim_generic()`, which wraps
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) with method
`"BFGS"` (unconstrained; all scale parameters are log-transformed). The
Hessian is computed numerically at the converged point for
variance-covariance estimation.

## 4 Golden fixture system

Golden fixtures are pre-fitted model results saved as `.rds` files in
`inst/fixtures/`. They serve as **regression anchors**: each test run
refits the model on the same data and compares estimates to the stored
values. This catches regressions when the likelihood, gradient, or
optimizer changes.

### 4.1 Fixture format

Every fixture is a named list with a consistent structure:

``` r
list(
  description = "Human-readable description",
  data = list(time, status, [x], n, events),
  seed = 42,
  true_params = list(...),   # simulation truth (synthetic fixtures)
  fit = list(
    theta     = <named numeric>,  # converged parameter vector
    logl      = <scalar>,         # log-likelihood at convergence
    converged = <logical>,        # convergence flag
    vcov      = <matrix>          # asymptotic variance-covariance
  ),
  timestamp = <POSIXct>
)
```

For C-reference fixtures (e.g., `mp_c_reference_kul.rds`), the `fit`
element is replaced by `c_reference` containing the C binary output:
parameter estimates, standard errors, variance-covariance matrix, and
the C log-likelihood.

### 4.2 Current fixtures

| File | Distribution | Description | Source |
|:---|:---|:---|:---|
| hz_univariate.rds | Weibull | Univariable shape estimation (n=100) | Synthetic (seed=42) |
| hm_multivariate.rds | Weibull | 2 covariates (n=100) | Synthetic (seed=42) |
| hm_edge_case.rds | Weibull | Edge case: n=20, 3 covariates | Synthetic (seed=42) |
| hz_loglogistic.rds | Log-logistic | Univariable (n=80) | Synthetic (seed=42) |
| hz_lognormal.rds | Log-normal | Univariable (n=80) | Synthetic (seed=42) |
| mp_synthetic_3phase.rds | Multiphase (3) | Synthetic early CDF + constant + late hazard (n=500) | Synthetic (seed=42) |
| mp_c_reference_kul.rds | Multiphase (3) | KUL CABG dataset + C HAZARD binary reference output (n=5880) | Clinical + C binary |

Golden fixture inventory {.table .caption-top style="width:100%;"}

### 4.3 Regenerating fixtures

Fixtures must be regenerated whenever the model parameterization
changes. The generators are maintainer-only helpers kept in `data-raw/`
(outside the installed package). Run them interactively from the package
source tree after `devtools::load_all()` and sourcing the generator
script:

``` r

devtools::load_all()
source("data-raw/golden_fixtures.R")

# Each generator requires an explicit output directory (no default path
# is written). Use the project source dir to persist, or tempdir() to test.
dir <- "inst/fixtures"

# Single-phase distributions
.hzr_create_synthetic_golden_fixtures(dir, seed = 42L)   # Weibull variants
.hzr_create_loglogistic_golden_fixture(dir, seed = 42L)  # Log-logistic
.hzr_create_lognormal_golden_fixture(dir, seed = 42L)    # Log-normal

# Multiphase
.hzr_create_multiphase_golden_fixture(dir, seed = 42L)   # Synthetic 3-phase
.hzr_create_c_reference_kul_fixture(dir)                  # C binary reference
```

Pass an explicit `seed` for reproducibility; the global RNG state is
saved and restored, so the caller’s stream is unaffected. With the
default (`seed = NULL`) no seed is set. Commit the updated `.rds` files
alongside any code changes.

## 5 Testing strategy

The test suite (`tests/testthat/`) is organized into four tiers:

| Tier | Files | Purpose |
|:---|:---|:---|
| Unit tests | test-math-primitives, test-decomposition, test-phase-spec, test-argument-mapping | Verify individual functions in isolation |
| Distribution tests | test-gradient-weibull, test-exponential-dist, test-loglogistic-dist, test-lognormal-dist, test-multiphase-gradient | Likelihood, gradient, and optimizer for each distribution |
| Integration tests | test-hazard-api, test-predict-types, test-interval-censoring-*, test-time-varying-*, test-multiphase-likelihood, test-multiphase-api | End-to-end: hazard() -\> predict() -\> summary() pipeline, censoring types, multiphase wiring |
| Parity tests | test-parity-core, test-parity-edge-cases, test-parity-c-binary, test-multiphase-parity | Golden fixture round-trip, C binary cross-validation |

Test suite tiers {.table .caption-top}

#### 5.0.1 Multiphase parity tests

The multiphase parity tests (`test-multiphase-parity.R`) validate
against the C HAZARD binary output for the KUL CABG dataset:

1.  **Likelihood evaluation** —Evaluates the R log-likelihood at the C
    converged parameters and asserts it matches the C output (-3740.52).

2.  **Decomposition consistency** —Verifies phase additivity and CDF
    saturation at the C reference parameter values.

3.  **Conservation of events** —Checks that the model-implied expected
    events (\\\sum \[1 - \exp(-H(t_i))\]\\) matches the observed event
    count (545), as reported by the C binary (544.9993).

4.  **Profile standard errors** —Computes a numerical Hessian varying
    only the 3 log(mu) parameters (shapes held fixed), matching the C
    binary’s estimation strategy, and compares standard errors.

5.  **Full fit convergence** —Fits the R multiphase optimizer on the
    full dataset with informed starting values and checks that the
    log-likelihood meets or exceeds the C reference.

## 6 Dataset catalog

TemporalHazard ships five clinical reference datasets in
`inst/extdata/`, converted from the original C/SAS HAZARD test data.
These datasets are used in vignette examples and parity testing.

| File | Study | n | Events | Covariates | SAS.origin |
|:---|:---|---:|:---|:---|:---|
| avc.csv | Atrioventricular canal repair | 310 | 70 deaths | NYHA, age, anatomy, era | hz.death.AVC, hm.death.AVC |
| cabgkul.csv | Coronary artery bypass grafting (KU Leuven) | 5880 | 545 deaths | None (intercept only) | hz.deadp.KUL |
| omc.csv | Open mitral commissurotomy | 339 | thromboembolic events (repeated) | TE events (repeated measures) | hz.te123.OMC |
| tga.csv | Transposition of great arteries (arterial switch) | 470 | deaths | anatomy, coronary pattern, era | hs.dthar.TGA |
| valves.csv | Primary valve replacement | 1533 | deaths, PVE, reoperation | age, NYHA, valve position, pathology, race | hm.deadp.VALVES |

Reference datasets (lazy-loaded via
[`data()`](https://rdrr.io/r/utils/data.html)) {.table .caption-top}

### 6.1 Loading datasets

Show code

``` r

# Datasets are lazy-loaded with the package — just reference them directly.
# Raw CSVs are also available in inst/extdata/ for advanced use:
#   read.csv(system.file("extdata", "cabgkul.csv", package = "TemporalHazard"))

data(cabgkul)
str(cabgkul)
#> 'data.frame':    5880 obs. of  2 variables:
#>  $ int_dead: num  201.83 195.06 7.13 126.36 187.57 ...
#>  $ dead    : int  0 0 1 1 0 0 1 1 0 1 ...

data(avc)
str(avc)
#> 'data.frame':    310 obs. of  11 variables:
#>  $ study   : chr  "001C" "002C" "004C" "005C" ...
#>  $ status  : int  3 3 1 2 2 3 1 1 3 3 ...
#>  $ inc_surg: int  4 3 2 3 1 2 3 2 3 3 ...
#>  $ opmos   : num  9.46 34.07 51.58 55 60.65 ...
#>  $ age     : num  69.2 53.7 286.1 154.6 48.4 ...
#>  $ mal     : int  0 0 0 1 0 0 0 0 0 0 ...
#>  $ com_iv  : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ orifice : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ dead    : int  1 1 0 0 0 0 0 0 0 1 ...
#>  $ int_dead: num  0.0534 0.3778 91.5337 111.608 106.8112 ...
#>  $ op_age  : num  654 1828 14759 8505 2933 ...
```

### 6.2 Dataset details

#### 6.2.1 AVC (atrioventricular canal repair)

The AVC dataset contains 310 patients who underwent repair of
atrioventricular septal defects. It is the primary example dataset used
throughout the hazard modeling examples and has the richest set of
covariates for multivariable analysis.

| Variable | Label                                     | Type      |
|:---------|:------------------------------------------|:----------|
| study    | Study number                              | character |
| status   | NYHA functional class (I-V)               | numeric   |
| inc_surg | Surgical grade of AV valve incompetence   | numeric   |
| opmos    | Date of operation (months since Jan 1967) | numeric   |
| age      | Age (months) at repair                    | numeric   |
| mal      | Important associated cardiac anomaly      | numeric   |
| com_iv   | Interventricular communication            | numeric   |
| orifice  | Accessory left AV valve orifice           | numeric   |
| dead     | Death (event indicator)                   | numeric   |
| int_dead | Follow-up interval (months)               | numeric   |
| op_age   | Interaction: opmos x age                  | numeric   |

AVC dataset variables {.table .caption-top}

#### 6.2.2 CABG/KUL (coronary artery bypass grafting)

The KUL dataset is a large series of 5880 primary isolated CABG patients
from KU Leuven (1971–July 1987). It serves as the primary benchmark for
C binary parity testing because it has the simplest structure
(intercept-only, right-censored) combined with a large sample size that
exercises all three temporal phases.

The original SAS analysis (`hz.deadp.KUL.sas`) also includes
return-of-angina and reintervention endpoints (recorded in the original
fixed-width file with 6 columns), but only the death endpoint
(`int_dead`, `dead`) is extracted for parity testing.

#### 6.2.3 OMC (open mitral commissurotomy)

The OMC dataset contains 339 patients and is unique in the collection
because it involves **repeated thromboembolic events** (up to 3 per
patient) with **left censoring**. The SAS analysis transforms the
dataset into a repeated-events format using STARTTME and CENSORED
indicators, exercising the interval censoring likelihood.

#### 6.2.4 TGA (transposition of great arteries)

The TGA dataset contains 470 patients who underwent the arterial switch
operation. It includes derived variables (logarithmic and inverse
transforms of clinical measurements) and is primarily used for
sensitivity analysis and internal validation examples.

#### 6.2.5 Valves (primary valve replacement)

The VALVES dataset is the largest multivariable example with 1533
patients and multiple endpoints (death, prosthetic valve endocarditis,
bioprosthesis degeneration, reoperation). The SAS analysis
(`hm.deadp.VALVES.sas`) fits a multivariable 3-phase model with 10
covariates across multiple event types.

## 7 SAS/C parameter mapping

The original C/SAS HAZARD procedure uses a different parameterization
than TemporalHazard. The
[`hzr_argument_mapping()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_argument_mapping.md)
function documents the full translation.

Show code

``` r

mapping <- hzr_argument_mapping()
knitr::kable(
  mapping[, c("legacy_input", "r_parameter", "implementation_status", "notes")],
  caption = "SAS HAZARD to R parameter mapping (excerpt)"
)
```

| legacy_input | r_parameter | implementation_status | notes |
|:---|:---|:---|:---|
| TIME variable | time | implemented | Core observation time input. |
| EVENT/censor variable | status | implemented | Event indicator currently retained as numeric in object\\data\\status. |
| X covariate block | x | implemented | Future versions will support richer design encoding helpers. |
| initial parameters | theta | implemented | Used by predict.hazard as coefficient vector. |
| baseline distribution | dist | implemented | Current default is ‘weibull’; more options planned. |
| control options | control | implemented | Control list is stored and reserved for optimizer parity. |
| additional legacy options | … | implemented | Supports legacy-style pass-through options during migration. |
| t | time | implemented | Canonical SAS migration uses TIME= mapping. |
| status | status | implemented | Canonical SAS migration uses EVENT= mapping. |
| theta0 | theta | planned | SAS PARMS syntax parser not yet implemented. |
| dist | dist | implemented | SAS DIST keyword maps directly to dist. |
| phases (3-phase structure) | phases (list of hzr_phase()) | implemented | Use dist=‘multiphase’ with phases argument. N-phase generalization of legacy 3-phase model. |
| MU_1, MU_2, MU_3 | mu (via exp(log_mu) in theta) | implemented | Each phase has its own scale mu_j(x) = exp(alpha_j + x\*beta_j). Starting value via hzr_phase(). |
| THALF / RHO (early) | hzr_phase(t_half=) | implemented | Half-life: time at which G(t_half) = 0.5. Same concept as SAS RHO/THALF. |
| NU (early) | hzr_phase(nu=) | implemented | Time exponent controlling rate dynamics. Same parameter name as SAS early NU. |
| M (early) | hzr_phase(m=) | implemented | Shape exponent controlling distributional form. Same parameter name as SAS early M. |
| DELTA (early) | (absorbed by decompos) | implemented | NOT IMPLEMENTED, not absorbed. The C DELTA controls B(t) = (exp(delta\*t)-1)/delta, which enters rho, the time argument and the density Jacobian separately; R computes the delta=0 branch of each. A job with DELTA != 0 is refused or flagged, never fitted silently. |
| G2 constant phase | hzr_phase(‘constant’) | implemented | Flat background rate. No shape parameters estimated. SAS G2 equivalent. |
| TAU (late) | hzr_phase(‘g3’, tau=) | implemented | Late-phase G3 scale parameter. Maps directly to hzr_phase(‘g3’, tau=). |
| GAMMA (late) | hzr_phase(‘g3’, gamma=) | implemented | Late-phase G3 time exponent. Maps directly to hzr_phase(‘g3’, gamma=). |
| ALPHA (late) | hzr_phase(‘g3’, alpha=) | implemented | Late-phase G3 shape parameter. alpha=0 gives exponential case. Maps directly to hzr_phase(‘g3’, alpha=). |
| ETA (late) | hzr_phase(‘g3’, eta=) | implemented | Late-phase G3 outer exponent. Maps directly to hzr_phase(‘g3’, eta=). |

SAS HAZARD to R parameter mapping (excerpt) {.table .caption-top}

### 7.1 Early phase (G1) mapping

The SAS early phase uses four parameters: DELTA, RHO (or THALF), NU, M.
These collapse onto the three-parameter decompos family:

- **DELTA** —Time transformation `B(t) = (exp(delta*t) - 1)/delta`. When
  DELTA = 0 (the common case), `B(t) = t` and the transformation is
  absorbed. Non-zero DELTA is not currently supported.
- **RHO / THALF** —Scale parameter.
  `RHO = NU * THALF * ((2^M - 1)/M)^NU`. Maps directly to `t_half` in
  [`hzr_phase()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase.md).
- **NU** —Time exponent. Maps directly.
- **M** —Shape exponent. Maps directly.

### 7.2 Late phase (G3) mapping

The SAS late phase uses the G3 decomposition with four parameters, now
directly supported via `hzr_phase("g3", ...)`:

- **TAU** –\> `tau` (scale parameter)
- **GAMMA** –\> `gamma` (time exponent)
- **ALPHA** –\> `alpha` (shape parameter; `alpha = 0` gives exponential
  case)
- **ETA** –\> `eta` (outer exponent)

The G3 formula (for `alpha > 0`) is: \\G_3(t) =
\left(\left((t/\tau)^\gamma + 1\right)^{1/\alpha} - 1\right)^\eta\\

Unlike G1, G3 is unbounded —it can grow without limit, making it
suitable for late-phase rising hazards. For the KUL benchmark with
`gamma = 3, alpha = 1, eta = 1`, this simplifies to \\G_3(t) = t^3\\.

## 8 Version history

The package follows semantic versioning with a prerelease qualifier
during active development:

- **v0.1.0** —Single-phase engine: Weibull, exponential, log-logistic,
  log-normal distributions with formula interface, predict, and golden
  fixture testing.
- **v0.9.0** —Multiphase engine: N-phase additive cumulative hazard,
  [`hzr_phase()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase.md)
  specification, decomposition engine, C binary parity tests, dataset
  catalog.
- **v0.9.1** (current) —Vignette suite, roxygen multiphase examples, CI
  workflow fixes (load_pkgload), SAS missing-value handling
  (`na.strings`), print.hazard phase labels, and README refresh.

## References

Blackstone EH, Naftel DC, Turner ME Jr. The decomposition of
time-varying hazard into phases, each incorporating a separate stream of
concomitant information. *J Am Stat Assoc.* 1986;81(395):615-624. doi:
[10.1080/01621459.1986.10478314](https://doi.org/10.1080/01621459.1986.10478314)

Rajeswaran J, Blackstone EH, Ehrlinger J, Li L, Ishwaran H, Parides MK.
Probability of atrial fibrillation after ablation: Using a parametric
nonlinear temporal decomposition mixed effects model. *Stat Methods Med
Res.* 2018;27(1):126-141. doi:
[10.1177/0962280215623583](https://doi.org/10.1177/0962280215623583)
