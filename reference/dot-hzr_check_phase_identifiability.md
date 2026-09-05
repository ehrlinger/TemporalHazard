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
  tol = 1e-08,
  other_times = NULL
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

  Threshold for all three tests – the minimum share of \\\Lambda\\ a
  phase must reach somewhere, and the minimum relative variation its
  contribution must show. Default 1e-8: far above double precision, and
  orders of magnitude below any real contribution, so it fires on dead
  phases rather than merely small ones.

- other_times:

  Further times the likelihood evaluates beyond `time` –
  counting-process entry times and interval bounds. The share and
  variation measures are taken over `time` alone, so when they are
  degenerate but these vary, the shapes still enter the likelihood and
  the measures are withheld rather than reported.

## Value

The share `data.frame`, invisibly; called for the warning.

## Details

Three conditions, not two. `absent` and `saturated` are per-phase. The
third is a property of the **times**: when their own relative range
falls below `tol` they separate no phase from any other, and the
per-phase messages would name a cause that did not occur (#211). It is
asked of the times rather than inferred from the phases, because a
saturated or absent phase is flat across times that are perfectly well
spread.

`absent` is measured on `share`, which does not depend on how the times
are spread, so it is always reported. The flatness verdict is withheld
when `other_times` vary, since the measures here cannot see those.
