# Specify a single hazard phase

Creates an `hzr_phase` object describing one term in a multiphase
additive cumulative hazard model. Pass a list of these to the `phases`
argument of
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
when `dist = "multiphase"`.

## Usage

``` r
hzr_phase(
  type = c("cdf", "hazard", "constant", "g3"),
  t_half = 1,
  nu = 1,
  m = 0,
  tau = 1,
  gamma = 1,
  alpha = 1,
  eta = 1,
  formula = NULL,
  fixed = character(0)
)

# S3 method for class 'hzr_phase'
print(x, ...)
```

## Arguments

- type:

  Character; the phase's temporal shape — one of `"cdf"` (early
  resolving risk), `"hazard"` (accumulating G1 aging risk), `"g3"` (late
  rising risk, the original C/SAS late phase), or `"constant"` (flat
  background rate). See the **Phase types** section for what each means
  and when to use it.

- t_half:

  Positive scalar; initial half-life (time at which \\G(t\_{1/2}) =
  0.5\\). Used for `"cdf"` and `"hazard"` phases. SAS early:
  `THALF`/`RHO`.

- nu:

  Numeric scalar; initial time exponent. Used for `"cdf"` and `"hazard"`
  phases. SAS early: `NU`.

- m:

  Numeric scalar; initial shape exponent. Used for `"cdf"` and
  `"hazard"` phases. SAS early: `M`.

- tau:

  Positive scalar; scale parameter for `"g3"` phases. SAS late: `TAU`.

- gamma:

  Positive scalar; time exponent for `"g3"` phases. SAS late: `GAMMA`.

- alpha:

  Non-negative scalar; shape parameter for `"g3"` phases. When
  `alpha > 0`, the generic G3 formula is used; `alpha = 0` gives the
  exponential limiting case. SAS late: `ALPHA`.

- eta:

  Positive scalar; outer exponent for `"g3"` phases. SAS late: `ETA`.

- formula:

  Optional one-sided formula (e.g. `~ age + nyha`) for phase-specific
  covariates. When `NULL` (default), the phase inherits the global
  formula from
  [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md).

- fixed:

  Character vector naming shape parameters to hold fixed during
  optimization. Valid names for `"cdf"`/`"hazard"`: `"t_half"`, `"nu"`,
  `"m"`, or `"shapes"` (shorthand for all three). Valid names for
  `"g3"`: `"tau"`, `"gamma"`, `"alpha"`, `"eta"`, or `"shapes"`
  (shorthand for all four). Fixed parameters are held at their starting
  values; only `mu` (and covariates) are estimated. Ignored for
  `"constant"` phases. This mirrors the SAS/C HAZARD workflow where
  shapes are typically fixed and only scale parameters are estimated.

- x:

  An `hzr_phase` object (for `print.hzr_phase()`).

- ...:

  Additional arguments (ignored).

## Value

An S3 object of class `"hzr_phase"` with elements:

- type:

  Phase type string.

- t_half:

  Initial half-life (cdf/hazard phases).

- nu:

  Initial time exponent (cdf/hazard phases).

- m:

  Initial shape exponent (cdf/hazard phases).

- tau:

  Scale parameter (g3 phases).

- gamma:

  Time exponent (g3 phases).

- alpha:

  Shape parameter (g3 phases).

- eta:

  Outer exponent (g3 phases).

- formula:

  Phase-specific formula or `NULL`.

- fixed:

  Character vector of fixed parameter names (may be empty).

## Role in the multiphase model

Each phase is one term \\j\\ in the additive cumulative hazard

\$\$H(t \mid \mathbf{x}) = \sum\_{j=1}^{J} \mu_j(\mathbf{x}) \\
\Phi_j(t)\$\$

where \\\mu_j(\mathbf{x}) = \exp(\alpha_j + \mathbf{x}\_j^\top
\beta_j)\\ is the phase-specific log-linear scale and \\\Phi_j(t)\\ is
the temporal shape selected by `type` (below). The `t_half`/`nu`/`m` (or
g3 `tau`/`gamma`/`alpha`/`eta`) arguments set the starting values for
that shape; `formula` attaches the covariates \\\mathbf{x}\_j\\ that
enter \\\mu_j\\.

## Phase types

The `type` argument chooses the temporal shape \\\Phi_j(t)\\ for the
phase. Each captures a qualitatively different pattern of risk over
time; a typical clinical model combines an *early*, a *constant*, and a
*late* phase so that the total hazard can fall, level off, and rise
again.

- `"cdf"` — early, resolving risk:

  Named for the **c**umulative **d**istribution **f**unction: the phase
  contributes \\\Phi(t) = G(t)\\, the bounded CDF of the temporal
  decomposition (\\0\\ at \\t = 0\\, rising to a ceiling of \\1\\).
  Because it saturates, the *hazard* it adds, \\\mu\\g(t)\\, peaks early
  and then decays toward zero — the signature of a one-time insult that
  patients either succumb to or survive past, e.g. peri-operative
  mortality. Shape set by `t_half`, `nu`, `m`. SAS/C equivalent: the
  Early (G1) phase.

- `"hazard"` — accumulating aging risk (G1 family):

  Named because the phase contributes a **cumulative hazard** built from
  the same G1 family: \\\Phi(t) = -\log(1 - G(t))\\, which is unbounded
  and monotone increasing. Its hazard \\\mu\\h(t)\\ rises without
  leveling off, so it models risk that grows as subjects age. This is an
  alternative late-risk form derived from G1; for the original SAS/C
  late phase prefer `"g3"`. Shape set by `t_half`, `nu`, `m`.

- `"g3"` — late, rising risk (original C/SAS late phase):

  Named for the **G3** (third) decomposition family used by the original
  HAZARD program for the late phase. It contributes \\\Phi(t) = G_3(t)\\
  from
  [`hzr_decompos_g3()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos_g3.md),
  an unbounded intensity with its own four-parameter shape (`tau`,
  `gamma`, `alpha`, `eta`) that is more flexible than the G1-derived
  `"hazard"` form for capturing accelerating late mortality (e.g.
  structural valve deterioration years after surgery). Use this when
  reproducing classic three-phase HAZARD models. SAS/C equivalent: the
  Late (G3) phase.

- `"constant"` — flat background rate:

  A time-invariant hazard: \\\Phi(t) = t\\, so the added hazard \\\mu\\
  is constant (the exponential model). It represents the steady, ongoing
  risk present at all follow-up times, independent of how long ago the
  time origin was. Takes no shape parameters — only its scale \\\mu\\
  (and any covariates) is estimated. SAS/C equivalent: the Constant (G2)
  phase.

The shape derivative \\\varphi_j = d\Phi_j/dt\\ (which forms the
instantaneous hazard contribution \\\mu_j\\\varphi_j(t)\\) is \\g(t)\\
for `"cdf"`, \\h(t)\\ for `"hazard"`, \\g_3(t)\\ for `"g3"`, and \\1\\
for `"constant"`.

## See also

[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
for fitting multiphase models,
[`hzr_decompos()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos.md)
for the underlying parametric family,
[`hzr_phase_cumhaz()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase_cumhaz.md)
and
[`hzr_phase_hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase_hazard.md)
for computing \\\Phi(t)\\ and \\\phi(t)\\ from these specifications.

[`vignette("fitting-hazard-models")`](https://ehrlinger.github.io/TemporalHazard/articles/fitting-hazard-models.md)
for multiphase fitting examples,
[`vignette("mf-mathematical-foundations")`](https://ehrlinger.github.io/TemporalHazard/articles/mf-mathematical-foundations.md)
for the mathematical framework.

## Examples

``` r
# Classic 3-phase Blackstone pattern
early <- hzr_phase("cdf",      t_half = 0.5, nu = 2, m = 0)
const <- hzr_phase("constant")
late  <- hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1)

# Fix all shapes (C/SAS-style: only estimate mu)
early_fixed <- hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
                          fixed = "shapes")
late_fixed  <- hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1,
                          fixed = "shapes")

# Fix only some parameters
early_partial <- hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
                            fixed = c("nu", "m"))

# Phase with specific covariates
early_cov <- hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
                        formula = ~ age + shock)

# Use in hazard():
# hazard(Surv(time, status) ~ age, data = dat,
#        dist = "multiphase",
#        phases = list(early = early, constant = const, late = late))
```
