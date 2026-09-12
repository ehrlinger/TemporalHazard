# Select the phase whose log_mu will be solved by conservation

Chooses the phase contributing the largest share of total cumulative
hazard at the current theta, matching the C HAZARD SETCOE strategy.

## Usage

``` r
.hzr_select_fixmu_phase(
  theta,
  time,
  status,
  phases,
  covariate_counts,
  x_list,
  weights = NULL,
  time_lower = NULL
)
```

## Arguments

- theta:

  Full parameter vector (internal scale).

- time:

  Numeric vector of follow-up times.

- status:

  Numeric event indicator.

- phases:

  Named list of validated `hzr_phase` objects.

- covariate_counts:

  Named integer vector.

- x_list:

  Named list of per-phase design matrices.

- weights:

  Optional numeric vector of row weights (length n). Defaults to unit
  weights. Applied when summing per-phase cumhaz so that selection
  happens on the same scale as the (weighted) observed event count.

- time_lower:

  Optional numeric vector. On a status 0/1 row with 0 \< time_lower \<
  time it is the counting-process entry (start) time, and phases are
  ranked by entry-time cumulative hazard `H(stop) - H(start)`, the scale
  on which events are conserved. Any other row has no entry
  (`H(start) = 0`); see `.hzr_multiphase_entry()`. `NULL` (the default)
  means no truncation.

## Value

Character: name of the phase to fix.

## Details

When one phase dominates by an extreme margin (\>10x the median
contribution) it is typically due to poorly-scaled shape parameters at
the starting theta rather than a genuine signal that the phase should
absorb all events. Selecting such a phase as the fixmu phase traps the
optimizer: CoE keeps rescaling that phase's log_mu upward while the
optimizer tries to drive it to near-zero. We therefore exclude extreme
outliers and select among the remaining phases. If all phases are
outliers (or only one phase exists) the largest is used as a fallback.
