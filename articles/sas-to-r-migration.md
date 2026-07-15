# Migrating from SAS HAZARD to TemporalHazard

``` r

library(TemporalHazard)
```

## Overview

If you’ve been running multiphase hazard analyses in SAS HAZARD — the
macro suite that wraps the C HAZARD binary written by Blackstone,
Naftel, and Turner at UAB — this vignette is the bridge to running the
same analyses in R. Every SAS `HAZARD` statement, every common macro
option, and every output field maps to something in `TemporalHazard`;
this document spells out the correspondences and flags the small handful
of features the R port doesn’t yet cover.

The mapping is faithful, not literal. SAS HAZARD is a step-driven DSL:
you write a `PROC HAZARD` block with separate `TIME`, `EVENT`, `PARMS`,
`EARLY`, `CONSTANT`, and `SELECTION` statements, and the macro assembles
them into a binary call. R is function-driven: you write a single
[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
call with named arguments. The conceptual pieces are identical; the
syntactic ergonomics differ. Use this vignette to translate an existing
SAS analysis, or as a reference when comparing the package’s output
against a SAS HAZARD reference fit for parity testing.

The full formal argument mapping table is available programmatically and
ships with the package — handy when you want to grep for a specific SAS
parameter name and find its R equivalent without scrolling through the
prose below:

``` r

knitr::kable(
  TemporalHazard::hzr_argument_mapping(),
  caption = "Formal argument map: SAS HAZARD/C → hazard()",
  col.names = c(
    "SAS Statement", "Legacy Input", "C Concept",
    "R Parameter", "Required", "Expected Type",
    "Transform Rule", "Status", "Notes"
  )
)
```

| SAS Statement | Legacy Input | C Concept | R Parameter | Required | Expected Type | Transform Rule | Status | Notes |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| HAZARD | TIME variable | obs time array | time | TRUE | numeric vector | pass through as numeric | implemented | Core observation time input. |
| HAZARD | EVENT/censor variable | event indicator array | status | TRUE | numeric/logical vector | coerce to numeric 0/1 | implemented | Event indicator currently retained as numeric in object\\data\\status. |
| HAZARD | X covariate block | design matrix | x | FALSE | numeric matrix or data.frame | data.frame -\> data.matrix | implemented | Future versions will support richer design encoding helpers. |
| HAZARD | initial parameters | parameter vector | theta | FALSE | numeric vector | length must equal ncol(x) when x is present | implemented | Used by predict.hazard as coefficient vector. |
| HAZARD | baseline distribution | phase distribution selector | dist | FALSE | character scalar | normalized lower-case label | implemented | Current default is ‘weibull’; more options planned. |
| HAZARD | control options | optimizer/control struct | control | FALSE | named list | stored in spec\$control | implemented | Control list is stored and reserved for optimizer parity. |
| HAZARD | additional legacy options | misc legacy switches | … | FALSE | named arguments | stored in legacy_args for parity | implemented | Supports legacy-style pass-through options during migration. |
| TIME | t | time vector | time | TRUE | numeric vector | pass through | implemented | Canonical SAS migration uses TIME= mapping. |
| EVENT | status | event vector | status | TRUE | numeric/logical vector | coerce to numeric | implemented | Canonical SAS migration uses EVENT= mapping. |
| PARMS | theta0 | starting coef | theta | FALSE | numeric vector | map PARMS/INITIAL to theta | planned | SAS PARMS syntax parser not yet implemented. |
| DIST | dist | dist selector | dist | FALSE | character scalar | map DIST= to dist | implemented | SAS DIST keyword maps directly to dist. |
| HAZARD | phases (3-phase structure) | 3-phase Early/Const/Late | phases (list of hzr_phase()) | FALSE | list of hzr_phase objects | list(early=hzr_phase(‘cdf’,…), constant=hzr_phase(‘constant’), late=hzr_phase(‘g3’,…)) | implemented | Use dist=‘multiphase’ with phases argument. N-phase generalization of legacy 3-phase model. |
| HAZARD | MU_1, MU_2, MU_3 | per-phase scale factors | mu (via exp(log_mu) in theta) | FALSE | numeric (per-phase) | exp(alpha_j) in internal parameterization; estimated on log scale | implemented | Each phase has its own scale mu_j(x) = exp(alpha_j + x\*beta_j). Starting value via hzr_phase(). |
| G1 | THALF / RHO (early) | early half-life | hzr_phase(t_half=) | FALSE | positive scalar | maps directly to hzr_phase(t_half=) starting value | implemented | Half-life: time at which G(t_half) = 0.5. Same concept as SAS RHO/THALF. |
| G1 | NU (early) | early time exponent | hzr_phase(nu=) | FALSE | numeric scalar | maps directly to hzr_phase(nu=) starting value | implemented | Time exponent controlling rate dynamics. Same parameter name as SAS early NU. |
| G1 | M (early) | early shape | hzr_phase(m=) | FALSE | numeric scalar | maps directly to hzr_phase(m=) starting value | implemented | Shape exponent controlling distributional form. Same parameter name as SAS early M. |
| G1 | DELTA (early) | early time transform | (absorbed by decompos) | FALSE | numeric scalar | time transform B(t) = (exp(delta\*t)-1)/delta absorbed into decompos shape | implemented | The C DELTA controlled B(t) = (exp(delta\*t)-1)/delta. This transform is absorbed by decompos(). |
| G2 | G2 constant phase | constant hazard rate phase | hzr_phase(‘constant’) | FALSE | hzr_phase(‘constant’) | hzr_phase(‘constant’) with no shape parameters | implemented | Flat background rate. No shape parameters estimated. SAS G2 equivalent. |
| G3 | TAU (late) | late G3 scale | hzr_phase(‘g3’, tau=) | FALSE | positive scalar | maps directly to hzr_phase(‘g3’, tau=) for late phase | implemented | Late-phase G3 scale parameter. Maps directly to hzr_phase(‘g3’, tau=). |
| G3 | GAMMA (late) | late G3 time exponent | hzr_phase(‘g3’, gamma=) | FALSE | numeric scalar | maps directly to hzr_phase(‘g3’, gamma=) for late phase | implemented | Late-phase G3 time exponent. Maps directly to hzr_phase(‘g3’, gamma=). |
| G3 | ALPHA (late) | late G3 shape | hzr_phase(‘g3’, alpha=) | FALSE | numeric scalar | maps directly to hzr_phase(‘g3’, alpha=) for late phase | implemented | Late-phase G3 shape parameter. alpha=0 gives exponential case. Maps directly to hzr_phase(‘g3’, alpha=). |
| G3 | ETA (late) | late G3 outer exponent | hzr_phase(‘g3’, eta=) | FALSE | numeric scalar | maps directly to hzr_phase(‘g3’, eta=) for late phase | implemented | Late-phase G3 outer exponent. Maps directly to hzr_phase(‘g3’, eta=). |

Formal argument map: SAS HAZARD/C → hazard() {.table .caption-top}

------------------------------------------------------------------------

## Statement-by-statement mapping

A SAS HAZARD analysis is a stack of statements inside a `PROC HAZARD`
block. We walk through them in roughly the order they appear in a
typical analysis script, showing the SAS form and the corresponding R
call.

### `PROC HAZARD DATA=`

The procedure-level statement that names the dataset and sets global
options. The `NOCOV` / `NOCOR` flags suppress covariance and correlation
output; `CONDITION=` is a tolerance switch for the convergence check.

``` sas
PROC HAZARD DATA=AVCS NOCOV NOCOR CONDITION=14;
```

Maps to
[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
`control` list:

``` r

fit <- hazard(
  ...,
  control = list(
    nocov      = TRUE,   # suppress covariance output
    nocor      = TRUE,   # suppress correlation output
    condition  = 14      # CONDITION= switch
  )
)
```

Additional `PROC HAZARD` options with no direct R equivalent yet are
passed through `...` as named arguments and stored in `fit$legacy_args`.

------------------------------------------------------------------------

### `TIME`

Names the follow-up time variable. In SAS HAZARD this is a separate
statement; in R it’s the first argument to
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) inside the
formula.

``` sas
TIME INT_DEAD;
```

Maps directly to `time`:

``` r

fit <- hazard(
  time   = avcs$INT_DEAD,
  ...
)
```

- Must be a non-negative numeric vector in the same time unit as the
  original analysis (months, in the AVC example).
- Missing values are not permitted.

------------------------------------------------------------------------

### `EVENT`

Names the event-indicator variable. Like `TIME`, this is a separate
statement in SAS HAZARD but enters the R formula through
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html).

``` sas
EVENT DEAD;
```

Maps to `status`:

``` r

fit <- hazard(
  status = avcs$DEAD,   # 1 = event, 0 = censored
  ...
)
```

- `TemporalHazard` coerces `status` to `numeric`; logical vectors are
  also accepted.
- The HAZARD convention uses `1` = occurred, `0` = censored/alive at
  last contact.

------------------------------------------------------------------------

### `PARMS`

Supplies starting values for the optimizer and flags which parameters to
hold fixed (`FIXM`, `FIXMU`, etc.). The starting values matter — for the
multiphase optimizer in particular, a poor starting point can park the
fit at a local minimum well away from the global MLE.

``` sas
PARMS MUE=0.3504743 THALF=0.1905077 NU=1.437416 M=1 FIXM
      MUC=4.391673E-07;
```

Maps to `theta` (coefficient/parameter vector) and `control` (fix
flags):

``` r

fit <- hazard(
  theta   = c(MUE = 0.3504743, THALF = 0.1905077, NU = 1.437416,
              M = 1,           MUC   = 4.391673e-07),
  control = list(
    fix = c("M")   # FIXM → freeze M during optimization
  ),
  ...
)
```

> **Note:** Full SAS `PARMS` syntax is not mirrored one-for-one in the
> public API. In the current package, supply the parameter vector
> directly via `theta` and record fixed-parameter intent in
> `control$fix`. The legacy parity helpers do generate SAS-style control
> text when comparing against the historical binaries.

------------------------------------------------------------------------

### `EARLY` and `CONSTANT` covariate blocks

These statements assign covariates to specific phases. In SAS HAZARD
each block lists the covariates that affect that phase and their
starting coefficients. Covariates can appear in multiple blocks with
different starting values — same column, different phase-specific
effect.

``` sas
EARLY  AGE=-0.03205774, COM_IV=1.336675, MAL=0.6872028,
       OPMOS=-0.01963377, OP_AGE=0.0002086689, STATUS=0.5169533;
CONSTANT INC_SURG=1.375285, ORIFICE=3.11765, STATUS=1.054988;
```

In HAZARD, `EARLY` and `CONSTANT` define phase-specific covariate
coefficients. In `TemporalHazard` these are unified into a single design
matrix `x` and coefficient vector `theta` during M1. Phase assignment
will be formalised in M2 when the multi-phase likelihood is implemented.

**Current convention (M1):** combine all covariates into `x` and supply
the corresponding starting coefficients in `theta`.

``` r

# Build design matrix from the AVC data set
X <- data.matrix(avcs[, c("AGE", "COM_IV", "MAL", "OPMOS", "OP_AGE",
                           "STATUS", "INC_SURG", "ORIFICE")])

# Starting values from SAS EARLY + CONSTANT blocks combined
theta_start <- c(
  AGE      = -0.03205774,
  COM_IV   =  1.336675,
  MAL      =  0.6872028,
  OPMOS    = -0.01963377,
  OP_AGE   =  0.0002086689,
  STATUS   =  0.5169533,   # EARLY phase coefficient
  INC_SURG =  1.375285,
  ORIFICE  =  3.11765
)

fit <- hazard(
  time   = avcs$INT_DEAD,
  status = avcs$DEAD,
  x      = X,
  theta  = theta_start,
  dist   = "weibull"
)
```

------------------------------------------------------------------------

### `SELECTION`

Triggers stepwise covariate selection inside the hazard fit. `SLE` is
the significance level for entry, `SLS` for staying. SAS HAZARD embedded
this in the same `PROC HAZARD` call; in R the stepwise loop is a
separate
[`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
function that wraps a fitted
[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
object.

``` sas
SELECTION SLE=0.2 SLS=0.1;
```

`TemporalHazard` now ships
[`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md),
which implements forward, backward, and two-way stepwise selection with
SAS-style SLENTRY / SLSTAY thresholds, phase-specific entry for
multiphase models, and a MOVE oscillation guard. FAST-screening
(Lawless-Singhal approximate Wald updates) is not yet implemented; every
candidate currently gets a full refit. See
[`?hzr_stepwise`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
for the full option list.

------------------------------------------------------------------------

## SAS macro equivalents

The legacy SAS distribution ships a suite of macros for nonparametric
baselines, goodness-of-fit diagnostics, bootstrap inference, and
variable calibration that sit alongside `PROC HAZARD`. Each has an R
equivalent in `TemporalHazard`:

| SAS macro | R function | Purpose |
|:---|:---|:---|
| `kaplan.sas` | [`hzr_kaplan()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_kaplan.md) | KM survival with logit-transformed exact CL |
| `nelsonl.sas` / `nelsont.sas` | [`hzr_nelson()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_nelson.md) | Nelson cumulative hazard with lognormal CL |
| `hazplot.sas` | [`hzr_gof()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_gof.md) | Parametric vs. KM overlay + CoE goodness-of-fit |
| `deciles.hazard.sas` | [`hzr_deciles()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_deciles.md) | Decile-of-risk calibration + chi-square GOF |
| `logit.sas` / `logitgr.sas` | [`hzr_calibrate()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_calibrate.md) | Quantile-group calibration (logit / Gompertz / Cox link) |
| `bootstrap.hazard.sas` | [`hzr_bootstrap()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_bootstrap.md) | Bootstrap resampling + summary |
| `markov.sas` | [`hzr_competing_risks()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_competing_risks.md) | Aalen-Johansen competing-risks incidence |

Minimal illustrations on the AVC dataset follow. The
[`vignette("inference-diagnostics")`](https://ehrlinger.github.io/temporal_hazard/articles/inference-diagnostics.md)
walkthrough builds these into a full analysis workflow.

``` r

library(survival)
data(avc)
avc <- na.omit(avc)

# Base parametric fit for downstream diagnostics
fit <- hazard(
  Surv(int_dead, dead) ~ age + status + mal + com_iv,
  data = avc, dist = "weibull",
  theta = c(mu = 0.2, nu = 1, rep(0, 4)),
  fit = TRUE, control = list(maxit = 500)
)
```

``` r

km <- hzr_kaplan(time = avc$int_dead, status = avc$dead)
head(km, 4)
#> Kaplan-Meier estimate with logit confidence limits
#> Events: 14  | Time points: 4 
#> Survival range: 0.9541 to 0.9836 
#> RMST at last event: 0.016 
#> 
#>         time n_risk n_event n_censor survival std_err cl_lower cl_upper cumhaz
#>  0.001368954    305       5        0   0.9836  0.0073   0.9612   0.9932 0.0165
#>  0.002737907    300       2        0   0.9770  0.0086   0.9527   0.9890 0.0232
#>  0.008213721    298       2        0   0.9705  0.0097   0.9443   0.9846 0.0300
#>  0.016427440    296       5        0   0.9541  0.0120   0.9240   0.9726 0.0470
#>   hazard density   life
#>  12.0744 11.9752 0.0014
#>   4.8862  4.7901 0.0027
#>   1.2298  1.1975 0.0081
#>   2.0741  1.9959 0.0160
```

``` r

nel <- hzr_nelson(time = avc$int_dead, event = avc$dead)
head(nel, 4)
#> Nelson cumulative hazard estimate with lognormal CL
#> Events: 14  | Time points: 4 
#> 
#>         time n_risk n_event weight_sum cumhaz std_err cl_lower cl_upper  hazard
#>  0.001368954    305       5          5 0.0164  0.0033   0.0109   0.0237 11.9752
#>  0.002737907    300       2          2 0.0231  0.0047   0.0153   0.0335  4.8699
#>  0.008213721    298       2          2 0.0298  0.0057   0.0201   0.0425  1.2256
#>  0.016427440    296       5          5 0.0467  0.0067   0.0350   0.0610  2.0565
#>  cum_events
#>           5
#>           7
#>           9
#>          14
```

``` r

head(hzr_gof(fit), 4)
#> Goodness-of-fit: observed vs. expected events
#> Distribution: weibull  | n = 305 
#> 
#> Total observed events: 68 
#> Total expected events: 50.433 
#> Final residual (E - O): -17.567 
#> Conservation ratio (E/O): 0.742 
#> 
#> Use plot columns: time, km_surv, par_surv, cum_observed, cum_expected, residual
```

``` r

hzr_deciles(fit, time = 120, groups = 10)
#> Decile-of-risk calibration (risk grouped at time = 120 )
#> 305 subjects, all included.
#> 10 groups, 68 observed events, 68 expected
#> 
#>  group  n events expected observed_rate expected_rate chi_sq p_value
#>      1 31      0     1.37        0.0000        0.0441 1.3700  0.2420
#>      2 30      1     1.64        0.0333        0.0546 0.2480  0.6190
#>      3 31      1     2.53        0.0323        0.0816 0.9250  0.3360
#>      4 30      2     3.49        0.0667        0.1160 0.6350  0.4260
#>      5 31      7     4.58        0.2260        0.1480 1.2800  0.2570
#>      6 30     10     5.67        0.3330        0.1890 3.3100  0.0689
#>      7 31      8     7.59        0.2580        0.2450 0.0218  0.8830
#>      8 30     14     8.71        0.4670        0.2900 3.2100  0.0734
#>      9 31     10    13.50        0.3230        0.4350 0.8960  0.3440
#>     10 30     15    18.90        0.5000        0.6320 0.8220  0.3650
#>  mean_survival mean_cumhaz
#>          0.950      0.0441
#>          0.931      0.0546
#>          0.906      0.0816
#>          0.865      0.1160
#>          0.818      0.1480
#>          0.708      0.1890
#>          0.655      0.2450
#>          0.538      0.2900
#>          0.441      0.4350
#>          0.221      0.6320
#> 
#> Overall: chi-sq = 12.7 on 9 df, p = 0.176
```

``` r

hzr_calibrate(avc$age, avc$dead, groups = 10, link = "logit")
#> Variable calibration (logit link, 10 groups)
#> 
#>  group  n events    mean     min     max  prob link_value
#>      1 30     11   3.519   1.051   5.388 0.367     -0.547
#>      2 31     11   8.665   5.421  11.532 0.355     -0.598
#>      3 30     13  15.194  11.631  18.497 0.433     -0.268
#>      4 31     11  23.077  18.990  27.828 0.355     -0.598
#>      5 30      7  43.544  28.124  57.167 0.233     -1.190
#>      6 31      3  72.066  59.730  86.408 0.097     -2.234
#>      7 30      2 101.154  86.507 117.522 0.067     -2.639
#>      8 31      3 162.739 121.169 203.733 0.097     -2.234
#>      9 30      4 247.051 205.343 297.140 0.133     -1.872
#>     10 31      3 530.623 324.573 790.981 0.097     -2.234
```

``` r

set.seed(1)
hzr_bootstrap(fit, n_boot = 20)  # small for vignette build
#> Bootstrap inference for hazard model
#> Mode: fixed refit 
#> Replicates: 20 successful, 0 failed
#> 
#>  parameter  n pct    mean     sd     min    max ci_lower ci_upper
#>         mu 20 100  0.0000 0.0000  0.0000 0.0000   0.0000   0.0000
#>         nu 20 100  0.2446 0.0170  0.2222 0.2866   0.2229   0.2798
#>        age 20 100 -0.0029 0.0016 -0.0070 0.0002  -0.0061  -0.0003
#>     status 20 100  0.7035 0.2262  0.3252 1.0922   0.3475   1.0776
#>        mal 20 100  0.4622 0.3263 -0.2049 0.9676  -0.1302   0.9326
#>     com_iv 20 100  0.7395 0.3402  0.0776 1.3490   0.1436   1.3182
```

``` r

data(valves)
# Synthesize a competing-risks indicator for illustration.
valves$ev <- with(valves, ifelse(dead == 1, 1L,
                                   ifelse(pve == 1,  2L,
                                   ifelse(reop == 1, 3L, 0L))))
head(hzr_competing_risks(valves$int_dead, valves$ev), 4)
#> Competing risks cumulative incidence
#> Event types: 3  | Time points: 4 
#> Final survival: 0.9941 
#> Final incid_1 : 0.0059 
#> Final incid_2 : 0 
#> Final incid_3 : 0 
#> 
#>     time n_risk n_event_1 n_event_2 n_event_3 n_censor   surv incid_1 incid_2
#>  0.00068   1533         1         0         0        0 0.9993  0.0007       0
#>  0.00137   1532         4         0         0        0 0.9967  0.0033       0
#>  0.00171   1528         1         0         0        0 0.9961  0.0039       0
#>  0.00205   1527         3         0         0        0 0.9941  0.0059       0
#>  incid_3 se_surv   se_1 se_2 se_3
#>        0  0.0007 0.0007    0    0
#>        0  0.0015 0.0015    0    0
#>        0  0.0016 0.0016    0    0
#>        0  0.0020 0.0020    0    0
```

------------------------------------------------------------------------

## Full worked example: AVC death after repair

Statement-by-statement mapping is fine for reference, but the full
gestalt only lands when you see a complete SAS HAZARD analysis and its R
translation side by side. The example below is the final multivariable
model from `examples/hm.death.AVC.sas` in the reference C repository —
death after atrioventricular canal repair, two-phase model with
covariates assigned to specific phases. It’s the canonical “this is what
a real SAS HAZARD analysis looks like” specimen.

### SAS (original)

The original SAS code, exactly as it appears in the reference
distribution. Notice how the statements stack inside a single
`%HAZARD()` macro call.

``` sas
%HAZARD(
PROC HAZARD DATA=AVCS P CONSERVE OUTHAZ=OUTEST CONDITION=14 QUASI;
     TIME INT_DEAD;
     EVENT DEAD;
     PARMS MUE=0.3504743 THALF=0.1905077 NU=1.437416 M=1 FIXM
           MUC=4.391673E-07;
     EARLY   AGE=-0.03205774, COM_IV=1.336675,  MAL=0.6872028,
             OPMOS=-0.01963377, OP_AGE=0.0002086689, STATUS=0.5169533;
     CONSTANT INC_SURG=1.375285, ORIFICE=3.11765, STATUS=1.054988;
);
```

### R equivalent (current runnable translation pattern)

The same model in `TemporalHazard`. Every SAS statement above has a
corresponding R argument here: `TIME` and `EVENT` collapse into
`Surv(int_dead, dead)` on the left of the formula; the global covariate
list goes on the right; the `PARMS` starting values become the `theta`
vector; phase-specific assignments use the `formula` argument inside
each
[`hzr_phase()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase.md);
`FIXM` becomes `fixed = "shapes"` (or a subset of shape names). The
line-by-line correspondence is the point — once you’ve translated one
analysis this way, the pattern carries to every other.

``` r

# Assumed: avcs is a data.frame read from the AVC flat file
# (same variables as the SAS DATA step)

avcs <- avcs |>
  transform(
    LN_AGE   = log(AGE),
    LN_OPMOS = log(OPMOS),
    LN_INC   = ifelse(is.na(INC_SURG), NA, log(INC_SURG + 1)),
    LN_NYHA  = log(STATUS)
  )

# Replace missing INC_SURG with column mean (mirrors PROC STANDARD REPLACE)
avcs$INC_SURG[is.na(avcs$INC_SURG)] <- mean(avcs$INC_SURG, na.rm = TRUE)

X <- data.matrix(avcs[, c("AGE", "COM_IV", "MAL", "OPMOS", "OP_AGE",
                           "STATUS", "INC_SURG", "ORIFICE")])

fit <- hazard(
  time    = avcs$INT_DEAD,
  status  = avcs$DEAD,
  x       = X,
  theta   = c(
    # Hazard shape parameters (from PARMS)
    MUE   = 0.3504743,
    THALF = 0.1905077,
    NU    = 1.437416,
    M     = 1,
    MUC   = 4.391673e-07,
    # Covariate coefficients (from EARLY + CONSTANT blocks)
    AGE      = -0.03205774,
    COM_IV   =  1.336675,
    MAL      =  0.6872028,
    OPMOS    = -0.01963377,
    OP_AGE   =  0.0002086689,
    STATUS   =  0.5169533,
    INC_SURG =  1.375285,
    ORIFICE  =  3.11765
  ),
  dist    = "weibull",
  control = list(
    condition = 14,
    quasi     = TRUE,
    conserve  = TRUE,
    fix       = c("M")   # FIXM
  )
)

fit
```

------------------------------------------------------------------------

## Data preparation differences

| SAS HAZARD | TemporalHazard R |
|----|----|
| `PROC STANDARD REPLACE` for missing | Replace with `mean(..., na.rm = TRUE)` manually |
| Log transforms in `DATA` step | [`transform()`](https://rdrr.io/r/base/transform.html) or `dplyr::mutate()` |
| Variable labels via `LABEL` | Column names of `x` matrix carry through |
| `FILENAME` / `INFILE` / `INPUT` | [`read.table()`](https://rdrr.io/r/utils/read.table.html), [`read.csv()`](https://rdrr.io/r/utils/read.table.html), or `haven::read_sas()` |
| SAS formats (date, currency, etc.) | Standard R numeric; apply [`as.Date()`](https://rdrr.io/r/base/as.Date.html) if needed |

------------------------------------------------------------------------

## Output object

[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
returns a list of class `hazard`:

| Slot | Contents |
|----|----|
| `$call` | Unevaluated [`match.call()`](https://rdrr.io/r/base/match.call.html) for reproducibility |
| `$spec$dist` | Baseline distribution label |
| `$spec$control` | Control options passed at construction |
| `$data$time` | Follow-up time vector |
| `$data$status` | Event indicator (numeric) |
| `$data$x` | Design matrix |
| `$fit$theta` | Coefficient vector (starting values at M1; fitted at M2+) |
| `$fit$converged` | `NA` at M1; `TRUE`/`FALSE` from M2 optimizer |
| `$fit$objective` | Log-likelihood at convergence (M2+) |
| `$legacy_args` | Named pass-through arguments for parity |

------------------------------------------------------------------------

## Prediction

In SAS HAZARD the `P` (predict / print) option on the `PROC HAZARD` line
writes predicted survival, hazard, and cumulative-hazard tables to the
output dataset. In R the same predictions come from
[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
on the fitted object, with the requested quantity chosen via the `type=`
argument.
[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
currently supports four output types — `"linear_predictor"`, `"hazard"`,
`"survival"`, and `"cumulative_hazard"` — covering every quantity the
SAS `P` option produces:

``` r

# Linear predictor (X %*% theta)
eta <- predict(fit, type = "linear_predictor")

# Hazard scale (exp(eta))
hz  <- predict(fit, type = "hazard")
```

The code chunk above shows only `"linear_predictor"` and `"hazard"`;
`"survival"` and `"cumulative_hazard"` follow the same pattern. For
multiphase fits pass `decompose = TRUE` to get per-phase cumulative
hazard contributions instead of just the total.

------------------------------------------------------------------------

## Known limitations vs. SAS HAZARD

A SAS HAZARD veteran migrating to `TemporalHazard` should be aware of
the following scope limits as of v0.9.8. Detailed status per feature is
tracked in `inst/dev/SAS-PARITY-GAP-ANALYSIS.md` and
`inst/dev/DEVELOPMENT-PLAN.md`.

### Stepwise variable selection (`SELECTION` statement)

- **Supported:**
  [`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
  implements forward / backward / two-way stepwise with SAS-style
  SLENTRY / SLSTAY thresholds, per-phase entry for multiphase, and a
  MOVE oscillation guard.
- **Not yet supported:** FAST-screening (Lawless-Singhal approximate
  Wald updates). Every candidate currently gets a full refit, which is
  slower but always correct.

### Per-phase time-varying windows

- **Supported:** `time_windows` applies one piecewise-constant cut-point
  set across all covariates.
- **Not yet supported:** distinct `EARLY`/`LATE` window sets per phase.
  Workaround: expand the design matrix manually before calling
  [`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md).

### Output datasets (`OUTEST=`, `OUTVCOV=`)

- The R equivalents are `coef(fit)` and `vcov(fit)` in memory; there is
  no automatic dataset-export mode.
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) them yourself if
  you need an on-disk artefact.

### Density / quantile prediction types

- [`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
  covers `hazard`, `cumulative_hazard`, `survival`, and
  `linear_predictor`. Density and quantile (median survival) types are
  not wired up; derive them from `cumulative_hazard` and `survival`
  predictions if needed.

### Previously listed gaps that have since been closed

For users migrating from older TemporalHazard versions or reading older
SAS parity notes:

- **`weights` on all distributions** — shipped v0.9.6. Weibull,
  exponential, log-logistic, log-normal, and multiphase all honour
  observation weights end-to-end (LL + analytic gradient + multiphase
  Conservation of Events).
- **`Surv(start, stop, event)` with `start > 0`** — shipped v0.9.7.
  Counting-process rows contribute `H(stop) - H(start)` for Weibull and
  multiphase. The previous
  [`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
  guard against non-zero starts is gone.
- **Delta-method prediction confidence limits** — shipped v0.9.8. Use
  `predict(..., se.fit = TRUE, level = 0.95)` to get a data frame with
  `fit`, `se.fit`, `lower`, and `upper` per row. Closed-form Jacobian
  for Weibull and multiphase,
  [`numDeriv::jacobian`](https://rdrr.io/pkg/numDeriv/man/jacobian.html)
  fallback for exponential / log-logistic / log-normal.

------------------------------------------------------------------------

## Rcpp acceleration note

`TemporalHazard` is currently pure R. If profiling against large real
datasets turns up bottlenecks in the likelihood kernel or the optimizer
inner loop, those specific functions will be re-implemented with Rcpp.
The public interface
([`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md),
[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md),
etc.) will not change.

## References

Blackstone EH, Naftel DC, Turner ME Jr. The decomposition of
time-varying hazard into phases, each incorporating a separate stream of
concomitant information. *J Am Stat Assoc.* 1986;81(395):615–624. doi:
[10.1080/01621459.1986.10478314](https://doi.org/10.1080/01621459.1986.10478314)
