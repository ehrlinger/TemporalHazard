# Cumulative hazard contribution from a single phase

Computes \\\Phi_j(t)\\ for one phase in the additive model \\H(t\|x) =
\sum_j \mu_j(x) \Phi_j(t)\\.

## Usage

``` r
hzr_phase_cumhaz(
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

Numeric vector of cumulative hazard contributions \\\Phi(t)\\, same
length as `time`.

## Details

- `"cdf"`: \\\Phi(t) = G(t)\\. Bounded \\\[0, 1\]\\. Models early risk
  that resolves over time.

- `"hazard"`: \\\Phi(t) = -\log(1 - G(t))\\. Monotone increasing. Models
  late or aging risk. This is the cumulative hazard derived from the
  hazard function \\h(t)\\, since \\\int_0^t h(s)\\ds = -\log(1 -
  G(t))\\.

- `"constant"`: \\\Phi(t) = t\\. Ignores `t_half`, `nu`, `m`. Equivalent
  to exponential (constant hazard rate).

## See also

[`hzr_decompos()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos.md)
for the underlying parametric family,
[`hzr_phase_hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase_hazard.md)
for the instantaneous hazard contribution.

## Examples

``` r
t_grid <- seq(0.1, 10, by = 0.1)
phi_early <- hzr_phase_cumhaz(t_grid, t_half = 2, nu = 2, m = 0,
                               type = "cdf")
phi_late  <- hzr_phase_cumhaz(t_grid, t_half = 5, nu = 1, m = 0,
                               type = "hazard")
phi_const <- hzr_phase_cumhaz(t_grid, type = "constant")
```
