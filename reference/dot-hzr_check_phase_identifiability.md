# Warn when a phase has effectively left the model

Warns rather than stops: the fit is arithmetically fine and the other
phases' estimates are usable. It is the unidentified parameters that
must not be read as estimates – and which ones those are differs by
mode, so the message says which.

## Usage

``` r
.hzr_check_phase_identifiability(
  theta,
  time,
  phases,
  covariate_counts,
  x_list,
  tol = 1e-08
)
```

## Arguments

- theta:

  Full parameter vector (internal scale).

- time:

  Numeric vector of follow-up times (n).

- phases:

  Named list of validated `hzr_phase` objects.

- covariate_counts:

  Named integer vector of per-phase covariate counts.

- x_list:

  Named list of per-phase design matrices.

- tol:

  Threshold for both tests – the minimum share of \\\Lambda\\ a phase
  must reach somewhere, and the minimum relative variation its
  contribution must show. Default 1e-8: far above double precision, and
  orders of magnitude below any real contribution, so it fires on dead
  phases rather than merely small ones.

## Value

The share `data.frame`, invisibly; called for the warning.
