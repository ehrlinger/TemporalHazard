# Derivatives of phase cumulative and instantaneous hazard w.r.t. shape params

Computes \\\Phi_j(t)\\, \\\phi_j(t)\\, and their derivatives with
respect to `t_half`, `nu`, and `m` using central finite differences on
[`hzr_decompos()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos.md).
The `log_t_half` derivative is obtained via the chain rule:
\\d\Phi/d(\log t\_{1/2}) = t\_{1/2} \cdot d\Phi/dt\_{1/2}\\.

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
