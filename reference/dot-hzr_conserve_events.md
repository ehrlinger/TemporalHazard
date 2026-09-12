# Apply the Conservation of Events adjustment to one phase's log_mu

Given the current theta vector, analytically solve the fixmu phase's
log_mu so that total predicted events = total observed events.

## Usage

``` r
.hzr_conserve_events(
  theta,
  fixmu_phase,
  fixmu_pos,
  time,
  status,
  phases,
  covariate_counts,
  x_list,
  total_events,
  weights = NULL,
  time_lower = NULL
)
```

## Arguments

- theta:

  Full parameter vector (internal scale).

- fixmu_phase:

  Character: name of the phase whose log_mu is solved.

- fixmu_pos:

  Integer: position of that log_mu in theta.

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

- total_events:

  Numeric: sum of observed events (precomputed, weighted when `weights`
  is supplied).

- weights:

  Optional numeric vector of row weights (length n). Defaults to unit
  weights. Applied when summing per-phase cumhaz so Turner's adjustment
  is computed on the same scale as `total_events`.

- time_lower:

  Optional numeric vector of counting-process entry (start) times. When
  supplied, conservation is enforced on the entry-time scale
  (`Sum E = Sum [H(stop) - H(start)]`) by subtracting the entry-time
  cumulative hazard, matching the multiphase likelihood (and C HAZARD
  `setcoe` under `LCENSOR`/`STARTTME`). `NULL` (the default) means no
  truncation, i.e. `H(start) = 0`.

## Value

Updated theta vector with fixmu phase's log_mu adjusted.

## Details

This is called BEFORE each likelihood evaluation inside the optimizer,
matching the C HAZARD `CONSRV` entry point.
