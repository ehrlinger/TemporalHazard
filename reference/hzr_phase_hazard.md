# Instantaneous hazard contribution from a single phase

Computes \\\phi_j(t) = d\Phi_j/dt\\ for one phase – the derivative of
the cumulative hazard contribution returned by
[`hzr_phase_cumhaz()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase_cumhaz.md).

## Usage

``` r
hzr_phase_hazard(
  time,
  t_half = 1,
  nu = 1,
  m = 0,
  type = c("cdf", "hazard", "constant")
)
```

## Arguments

- time:

  Numeric vector of times (\> 0).

- t_half:

  Half-life parameter (\> 0).

- nu:

  Time exponent.

- m:

  Shape parameter.

- type:

  Phase type: `"cdf"` (early – uses \\G(t)\\), `"hazard"` (late – uses
  cumulative hazard from \\h(t)\\), or `"constant"` (flat rate – \\\Phi
  = t\\).

## Value

Numeric vector of instantaneous hazard contributions \\\phi(t)\\, same
length as `time`.

## Details

- `"cdf"`: \\\phi(t) = g(t)\\ (density).

- `"hazard"`: \\\phi(t) = h(t) = g(t)/(1-G(t))\\.

- `"constant"`: \\\phi(t) = 1\\.

## See also

[`hzr_decompos()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos.md)
for the underlying parametric family,
[`hzr_phase_cumhaz()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase_cumhaz.md)
for the cumulative version.

## Examples

``` r
t_grid <- seq(0.1, 10, by = 0.1)
phi_early <- hzr_phase_hazard(t_grid, t_half = 2, nu = 2, m = 0,
                               type = "cdf")
phi_late  <- hzr_phase_hazard(t_grid, t_half = 5, nu = 1, m = 0,
                               type = "hazard")
```
