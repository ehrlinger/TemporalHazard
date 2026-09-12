# Number of shape parameters for a phase

Returns 3 (t_half, nu, m) for `"cdf"` and `"hazard"` phases, 4 (tau,
gamma, alpha, eta) for `"g3"`, and 0 for `"constant"`.

## Usage

``` r
.hzr_phase_n_shape(phase)
```

## Arguments

- phase:

  An `hzr_phase` object.

## Value

Integer: 3, 4 or 0.
