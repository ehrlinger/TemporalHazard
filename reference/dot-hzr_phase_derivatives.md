# Derivatives of phase cumulative and instantaneous hazard w.r.t. shape params

Computes \\\Phi_j(t)\\, \\\phi_j(t)\\, and their derivatives with
respect to `log_t_half`, `nu`, and `m` using finite differences on
[`hzr_decompos()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos.md):
central in `log_t_half` and `nu`, and in `m` central except near
`m = 0`, where the stencil keeps the sign of `m`.
[`hzr_decompos()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos.md)
changes formula at `m = 0` and the two sides meet in a cusp, so for
`m >= 0` a stencil that would reach 0 becomes one-sided forward (second
order), and for `m < 0` the step is capped at 1% of `|m|` (floored at
1e-10) and becomes one-sided backward if it would still reach 0. The
`log_t_half` derivative is a central difference in `log_t_half` itself,
so its step is proportional to `t_half` at every scale. It is `NaN`
where the two points of the difference cannot both be used: `t_half` too
small or too large to step, or a point at which
[`hzr_decompos()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos.md)
fails. A `"hazard"` phase's value is `-log_surv` from
[`hzr_decompos()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos.md),
which keeps its accuracy far past saturation, so the derivative does
too.

## Usage

``` r
.hzr_phase_derivatives(
  time,
  t_half,
  nu,
  m,
  type = c("cdf", "hazard", "constant")
)
```

## Arguments

- time:

  Numeric vector of positive times.

- t_half:

  Positive scalar half-life.

- nu:

  Numeric scalar time exponent.

- m:

  Numeric scalar shape exponent.

- type:

  Phase type: `"cdf"`, `"hazard"`, or `"constant"`.

## Value

Named list:

- Phi:

  Cumulative hazard contribution \\\Phi(t)\\.

- phi:

  Instantaneous hazard contribution \\\phi(t) = d\Phi/dt\\.

- dPhi_dlog_thalf:

  \\d\Phi / d(\log t\_{1/2})\\.

- dPhi_dnu:

  \\d\Phi / d\nu\\.

- dPhi_dm:

  \\d\Phi / dm\\.

- dphi_dlog_thalf:

  \\d\phi / d(\log t\_{1/2})\\.

- dphi_dnu:

  \\d\phi / d\nu\\.

- dphi_dm:

  \\d\phi / dm\\.
