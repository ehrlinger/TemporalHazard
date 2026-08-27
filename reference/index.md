# Package index

## Hazard model

- [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
  : Build and optionally fit a hazard model
- [`summary(`*`<hazard>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/summary.hazard.md)
  : Summarize a hazard model
- [`predict(`*`<hazard>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/predict.hazard.md)
  : Predict from a hazard model object
- [`coef(`*`<hazard>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/coef.hazard.md)
  : Extract coefficients from hazard model
- [`vcov(`*`<hazard>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/vcov.hazard.md)
  : Extract variance-covariance matrix from hazard model

## Multiphase phases

- [`hzr_phase()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase.md)
  [`print(`*`<hzr_phase>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase.md)
  : Specify a single hazard phase
- [`is_hzr_phase()`](https://ehrlinger.github.io/TemporalHazard/reference/is_hzr_phase.md)
  : Test if an object is an hzr_phase
- [`hzr_phase_cumhaz()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase_cumhaz.md)
  : Cumulative hazard contribution from a single phase
- [`hzr_phase_hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase_hazard.md)
  : Instantaneous hazard contribution from a single phase

## Variable selection

- [`hzr_stepwise()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
  [`print(`*`<hzr_stepwise>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
  [`summary(`*`<hzr_stepwise>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
  [`print(`*`<summary.hzr_stepwise>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
  [`as.data.frame(`*`<hzr_stepwise>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_stepwise.md)
  : Stepwise covariate selection for a parametric hazard model

- [`stepwise_trace()`](https://ehrlinger.github.io/TemporalHazard/reference/stepwise_trace.md)
  :

  Extract the captured console trace from an `hzr_stepwise` fit

## Decomposition

- [`hzr_decompos()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos.md)
  : Generalized temporal decomposition
- [`hzr_decompos_g3()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos_g3.md)
  : Late-phase (G3) temporal decomposition

## Datasets

- [`avc`](https://ehrlinger.github.io/TemporalHazard/reference/avc.md) :
  AVC: Atrioventricular Canal Repair
- [`cabgkul`](https://ehrlinger.github.io/TemporalHazard/reference/cabgkul.md)
  : CABGKUL: Primary Isolated Coronary Artery Bypass Grafting (KU
  Leuven)
- [`omc`](https://ehrlinger.github.io/TemporalHazard/reference/omc.md) :
  OMC: Open Mitral Commissurotomy
- [`tga`](https://ehrlinger.github.io/TemporalHazard/reference/tga.md) :
  TGA: Transposition of the Great Arteries
- [`uslife2023`](https://ehrlinger.github.io/TemporalHazard/reference/uslife2023.md)
  : US Life Table 2023: All-Interval-Censored SAS Parity Anchor
- [`valves`](https://ehrlinger.github.io/TemporalHazard/reference/valves.md)
  : Valves: Primary Heart Valve Replacement

## Migration tools

- [`hzr_argument_mapping()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_argument_mapping.md)
  : Legacy HAZARD to TemporalHazard argument mapping

- [`hzr_read_outhaz()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_read_outhaz.md)
  :

  Read a SAS `outhaz` estimate dataset

- [`predict(`*`<hzr_outhaz>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/predict.hzr_outhaz.md)
  :

  Predictions from a fit loaded out of a SAS `OUTHAZ=` dataset

- [`hzr_translate_sas()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_translate_sas.md)
  : Translate a SAS HAZARD job into a Quarto document

## Numerical Primitives

- [`hzr_log1pexp()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_log1pexp.md)
  : Numerically stable log(1 + exp(x))
- [`hzr_log1mexp()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_log1mexp.md)
  : Numerically stable log(1 - exp(-x)) for x \> 0
- [`hzr_clamp_prob()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_clamp_prob.md)
  : Clamp probabilities away from 0 and 1

## Diagnostics

- [`hzr_deciles()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_deciles.md)
  : Decile-of-risk calibration
- [`print(`*`<hzr_deciles>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/print.hzr_deciles.md)
  : Print method for hzr_deciles
- [`hzr_gof()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_gof.md)
  : Goodness-of-fit: observed vs. predicted events
- [`print(`*`<hzr_gof>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/print.hzr_gof.md)
  : Print method for hzr_gof
- [`hzr_kaplan()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_kaplan.md)
  : Kaplan-Meier survival with exact logit confidence limits
- [`print(`*`<hzr_kaplan>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/print.hzr_kaplan.md)
  : Print method for hzr_kaplan
- [`hzr_calibrate()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_calibrate.md)
  : Calibrate a continuous variable against an outcome
- [`print(`*`<hzr_calibrate>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/print.hzr_calibrate.md)
  : Print method for hzr_calibrate
- [`hzr_nelson()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_nelson.md)
  [`print(`*`<hzr_nelson>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_nelson.md)
  : Wayne Nelson cumulative hazard estimator with lognormal confidence
  limits
- [`hzr_bootstrap()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_bootstrap.md)
  [`print(`*`<hzr_bootstrap>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_bootstrap.md)
  : Bootstrap resampling for hazard model coefficients
- [`hzr_competing_risks()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_competing_risks.md)
  [`print(`*`<hzr_competing_risks>`*`)`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_competing_risks.md)
  : Competing risks cumulative incidence
