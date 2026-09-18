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
  fixed = character(0),
  constraint = c("none", "alpha_gamma_eta", "eta_gamma")
)

# S3 method for class 'hzr_phase'
print(x, ...)
```

## Arguments

- type:

  Character; the phase's temporal shape, one of `"cdf"` (early resolving
  risk), `"hazard"` (accumulating G1 aging risk), `"g3"` (late rising
  risk, the original C/SAS late phase), or `"constant"` (flat background
  rate). See the **Phase types** section for what each means and when to
  use it.

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

  Positive finite scalar; scale parameter for `"g3"` phases. SAS late:
  `TAU`.

- gamma:

  Positive finite scalar; time exponent for `"g3"` phases. SAS late:
  `GAMMA`.

- alpha:

  Non-negative scalar; shape parameter for `"g3"` phases. When
  `alpha > 0`, the generic G3 formula is used; `alpha = 0` gives the
  exponential limiting case. SAS late: `ALPHA`.

- eta:

  Positive finite scalar; outer exponent for `"g3"` phases. SAS late:
  `ETA`.

- formula:

  Optional one-sided formula (e.g. `~ age + nyha`) for phase-specific
  covariates. It is evaluated in the `data` given to
  [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md),
  so without `data`
  [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
  refuses it, unless it builds nothing either way: an intercept-only
  `~ 1` with no global `x`. The phase's scale parameter plays the role
  of an intercept, so the design never has one: removing it
  (`~ 0 + age`, `~ age - 1`) is ignored with a warning, and builds the
  design of `~ age`. When `NULL` (default), the phase inherits the
  global design from
  [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md):
  the global formula's covariates, or `x` on the vector interface.

- fixed:

  Character vector naming shape parameters to hold fixed during
  optimization. Valid names for `"cdf"`/`"hazard"`: `"t_half"`, `"nu"`,
  `"m"`, or `"shapes"` (shorthand for all three). Valid names for
  `"g3"`: `"tau"`, `"gamma"`, `"alpha"`, `"eta"`, or `"shapes"`
  (shorthand for all four). Fixed parameters are held at their starting
  values; only `mu` (and covariates) are estimated. Ignored for
  `"constant"` phases. This mirrors the SAS/C HAZARD workflow where
  shapes are typically fixed and only scale parameters are estimated.

- constraint:

  For `"g3"` phases, a rule that *derives* one shape from the others
  rather than estimating it:

  `"none"`

  :   (default) every shape is estimated or fixed.

  `"alpha_gamma_eta"`

  :   \\\alpha = \gamma\eta/2\\, so that \\\gamma\eta/\alpha = 2\\.
      SAS/C: `FIXGAE2`.

  `"eta_gamma"`

  :   \\\eta = 2/\gamma\\, so that \\\gamma\eta = 2\\. SAS/C: `FIXGE2`.

  The derived parameter follows the others at every step of the
  optimization, so it is not a free parameter and cannot be named in
  `fixed`; `"shapes"` leaves it out. Its starting value is computed from
  the others, and a value you supply for it, here or in the `theta`
  given to
  [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md),
  is replaced, with a warning when it differs. Under
  `hazard(fit = FALSE)` that replacement is made only when no phase
  carries covariates, since otherwise its slot in `theta` is not known
  until the design is built;
  [`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
  warns when it could not be made. Its standard error in
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) is the delta-method one,
  carried from the parameters it is derived from, so confidence limits
  from [`predict()`](https://rdrr.io/r/stats/predict.html) include its
  uncertainty. Only one constraint per phase: SAS's `FIXGE2` and
  `FIXGAE2` together force \\\alpha = 1\\ and fix `tau`, `gamma` and
  `eta`, which is better written as those fixed values.

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

- constraint:

  The shape constraint (g3 phases).

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

- `"cdf"`: early, resolving risk:

  Named for the **c**umulative **d**istribution **f**unction: the phase
  contributes \\\Phi(t) = G(t)\\, the bounded CDF of the temporal
  decomposition (\\0\\ at \\t = 0\\, rising to a ceiling of \\1\\).
  Because it saturates, the *hazard* it adds, \\\mu\\g(t)\\, peaks early
  and then decays toward zero (the signature of a one-time insult that
  patients either succumb to or survive past, e.g. peri-operative
  mortality). Shape set by `t_half`, `nu`, `m`. SAS/C equivalent: the
  Early (G1) phase.

- `"hazard"`: accumulating aging risk (G1 family):

  Named because the phase contributes a **cumulative hazard** built from
  the same G1 family: \\\Phi(t) = -\log(1 - G(t))\\, which is unbounded
  and monotone increasing. Its hazard \\\mu\\h(t)\\ rises without
  leveling off, so it models risk that grows as subjects age. This is an
  alternative late-risk form derived from G1; for the original SAS/C
  late phase prefer `"g3"`. Shape set by `t_half`, `nu`, `m`.

- `"g3"`: late, rising risk (original C/SAS late phase):

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

- `"constant"`: flat background rate:

  A time-invariant hazard: \\\Phi(t) = t\\, so the added hazard \\\mu\\
  is constant (the exponential model). It represents the steady, ongoing
  risk present at all follow-up times, independent of how long ago the
  time origin was. Takes no shape parameters; only its scale \\\mu\\
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

# Derive alpha from gamma and eta (SAS/C FIXGAE2)
late_gae2 <- hzr_phase("g3", tau = 14, gamma = 22, eta = 0.18,
                        constraint = "alpha_gamma_eta")

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
