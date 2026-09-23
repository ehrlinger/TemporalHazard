# Finite-difference derivatives of G3 phase Phi and phi w.r.t. shape params

Computes the G3 cumulative intensity, its time derivative, and their
partial derivatives with respect to log_tau, gamma, alpha, and eta using
central finite differences.

## Usage

``` r
.hzr_g3_phase_derivatives(time, tau, gamma, alpha, eta, h = 1e-05)
```

## Arguments

- time:

  Numeric vector of positive times.

- tau:

  Positive scalar scale parameter.

- gamma:

  Positive scalar time exponent.

- alpha:

  Non-negative scalar shape parameter.

- eta:

  Positive scalar outer exponent.

- h:

  Finite-difference step (default 1e-5). It is an ABSOLUTE step in
  `log_tau`, and so relative in `tau` only to first order; a relative
  step with a `1e-10` floor for `gamma` and `eta`; and strictly
  proportional for `alpha`. No caller passes it. Accuracy degrades where
  the function's natural scale in `log_tau`, roughly `alpha / gamma`, is
  far from 1: at `alpha = 0.05, gamma = 4` the `tau` derivative is
  accurate to about 1e-6 rather than the 1e-11 it reaches for ordinary
  shapes.

## Value

Named list with Phi, phi, and 8 derivative vectors.
