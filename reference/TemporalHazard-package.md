# TemporalHazard: Temporal Parametric Hazard Modeling

Native R implementation of the multiphase parametric hazard model of
Blackstone, Naftel, and Turner (1986), with a focus on behavioral
parity, transparent numerics, and reproducible validation against the
original ‘C’/‘SAS’ HAZARD program. The package fits time-varying hazards
as an additive sum of parametric *phases* — early, constant, and late
risk streams — each carrying its own covariate effects.

## The multiphase model

The total cumulative hazard decomposes additively across \\J\\ phases:

\$\$H(t \mid \mathbf{x}) = \sum\_{j=1}^{J} \mu_j(\mathbf{x}) \\
\Phi_j(t)\$\$

where \\\mu_j(\mathbf{x}) = \exp(\alpha_j + \mathbf{x}\_j^\top
\beta_j)\\ is the phase-specific log-linear scale (intercept plus
covariate effects) and \\\Phi_j(t)\\ is the temporal shape contributed
by phase \\j\\, with its own parameters depending on the phase `type`
(see the phase vocabulary below). The instantaneous hazard and survival
follow directly:

\$\$h(t \mid \mathbf{x}) = \sum\_{j=1}^{J} \mu_j(\mathbf{x}) \\
\varphi_j(t), \qquad S(t \mid \mathbf{x}) = \exp\\\bigl(-H(t \mid
\mathbf{x})\bigr)\$\$

with \\\varphi_j = d\Phi_j/dt\\. Each phase is specified with
[`hzr_phase()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase.md)
and the model is fit by maximum likelihood with
[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md).

## Phase vocabulary

|  |  |  |  |
|----|----|----|----|
| **Phase type** | **\\\Phi_j(t)\\** | **Domain** | **Use** |
| `"cdf"` | \\G(t)\\ | \\\[0, 1\]\\ | Early risk that resolves over time |
| `"hazard"` | \\-\log(1 - G(t))\\ | \\\[0, \infty)\\ | Late or aging risk that accumulates |
| `"g3"` | \\G_3(t)\\ | \\\[0, \infty)\\ | Unbounded late risk (original C/SAS late phase) |
| `"constant"` | \\t\\ | \\\[0, \infty)\\ | Flat background rate (no shape parameters) |

Here \\G(t)\\ is the generalized temporal decomposition CDF computed by
[`hzr_decompos()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_decompos.md),
and \\G_3(t)\\ is the unbounded late-phase intensity from
[`hzr_decompos_g3()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_decompos_g3.md).
See
[`vignette("mf-mathematical-foundations")`](https://ehrlinger.github.io/temporal_hazard/articles/mf-mathematical-foundations.md)
for the full derivation.

## SAS/C HAZARD bridge

The classic three-phase HAZARD model maps directly onto
[`hzr_phase()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase.md)
calls:

|                  |                      |                         |
|------------------|----------------------|-------------------------|
| **HAZARD phase** | **Role**             | **R equivalent**        |
| G1 (early)       | Early resolving risk | `hzr_phase("cdf", ...)` |
| G2 (constant)    | Flat background rate | `hzr_phase("constant")` |
| G3 (late)        | Rising late risk     | `hzr_phase("g3", ...)`  |

TemporalHazard generalizes the fixed three-phase structure to \\N\\
phases of any type.
[`hzr_argument_mapping()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_argument_mapping.md)
gives the full parameter translation table between the SAS/C
parameterization and the R arguments.

## Main entry points

- Model fitting:

  [`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
  — build and fit single- or multiphase models;
  [`hzr_phase()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase.md)
  — specify one phase;
  [`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
  — forward, backward, or bidirectional covariate selection.

- Prediction:

  [`predict()`](https://rdrr.io/r/stats/predict.html) on a fitted
  `hazard` object — survival, cumulative hazard, and per-phase
  decomposed hazard; [`summary()`](https://rdrr.io/r/base/summary.html)
  — coefficient tables with Wald inference.

- Parametric family:

  [`hzr_decompos()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_decompos.md)
  — the early-phase (G1) decomposition \\G(t)\\, \\g(t)\\, \\h(t)\\;
  [`hzr_decompos_g3()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_decompos_g3.md)
  — the late-phase (G3) intensity;
  [`hzr_phase_cumhaz()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase_cumhaz.md)
  and
  [`hzr_phase_hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase_hazard.md)
  — per-phase \\\Phi(t)\\ and \\\varphi(t)\\.

- Diagnostics:

  [`hzr_kaplan()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_kaplan.md),
  [`hzr_nelson()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_nelson.md)
  — nonparametric references;
  [`hzr_gof()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_gof.md)
  — goodness of fit;
  [`hzr_calibrate()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_calibrate.md),
  [`hzr_deciles()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_deciles.md)
  — calibration;
  [`hzr_bootstrap()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_bootstrap.md)
  — resampling CIs;
  [`hzr_competing_risks()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_competing_risks.md)
  — cumulative incidence.

## Vignettes

- [`vignette("getting-started")`](https://ehrlinger.github.io/temporal_hazard/articles/getting-started.md):

  First fit, end to end.

- [`vignette("mf-mathematical-foundations")`](https://ehrlinger.github.io/temporal_hazard/articles/mf-mathematical-foundations.md):

  The decomposition family and multiphase model, with derivations.

- [`vignette("fitting-hazard-models")`](https://ehrlinger.github.io/temporal_hazard/articles/fitting-hazard-models.md):

  Single-phase through multiphase fitting.

- [`vignette("prediction-visualization")`](https://ehrlinger.github.io/temporal_hazard/articles/prediction-visualization.md):

  Prediction types and decomposed-hazard plots.

- [`vignette("inference-diagnostics")`](https://ehrlinger.github.io/temporal_hazard/articles/inference-diagnostics.md):

  Bootstrap CIs and diagnostics.

- [`vignette("sas-to-r-migration")`](https://ehrlinger.github.io/temporal_hazard/articles/sas-to-r-migration.md):

  Translating SAS/C HAZARD code.

## References

Blackstone EH, Naftel DC, Turner ME Jr. The decomposition of
time-varying hazard into phases, each incorporating a separate stream of
concomitant information. *J Am Stat Assoc.* 1986;81(395):615–624.
[doi:10.1080/01621459.1986.10478314](https://doi.org/10.1080/01621459.1986.10478314)

Rajeswaran J, Blackstone EH, Ehrlinger J, Li L, Ishwaran H, Parides MK.
Probability of atrial fibrillation after ablation: Using a parametric
nonlinear temporal decomposition mixed effects model. *Stat Methods Med
Res.* 2018;27(1):126–141.
[doi:10.1177/0962280215623583](https://doi.org/10.1177/0962280215623583)

## See also

Useful links:

- <https://ehrlinger.github.io/temporal_hazard/>

- <https://github.com/ehrlinger/temporal_hazard>

- Report bugs at <https://github.com/ehrlinger/temporal_hazard/issues>

## Author

**Maintainer**: John Ehrlinger <john.ehrlinger@gmail.com> \[copyright
holder\]

Authors:

- John Ehrlinger <john.ehrlinger@gmail.com> \[copyright holder\]
