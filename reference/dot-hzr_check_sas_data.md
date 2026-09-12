# Check the SAS objective's data preconditions at entry

Both conditions `objective = "sas"` imposes (no left-censored rows, and
a positive width on every interval row) are pure functions of the data,
so they hold or fail identically at every start. Evaluated inside the
objective they reach the user through
[`.hzr_optim_multiphase()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_optim_multiphase.md)'s
per-start `tryCatch`, which frames them as "produced no usable fit from
N starts" and invites raising `n_starts`, a remedy that cannot work.
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
calls this once, before any optimization, so a data defect is reported
as one.

## Usage

``` r
.hzr_check_sas_data(status, time, time_lower, time_upper, objective)
```

## Arguments

- status:

  Numeric event indicator.

- time:

  Event/censoring times.

- time_lower, time_upper:

  Optional censoring bounds; `NULL` means `time`.

- objective:

  Resolved objective, `"likelihood"` or `"sas"`.

## Value

`NULL`, invisibly; called for its side effect.

## Details

This does **not** replace the guards inside the objective and the
gradient. Those remain because the gradient is reachable without
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
(the score test calls it directly), so entry validation is not
guaranteed to have run. See
[`.hzr_check_sas_status()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_check_sas_status.md).

Bounds are normalised here exactly as
[`.hzr_logl_multiphase()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_logl_multiphase.md)
normalises them, so the check cannot disagree with the objective about
which rows are offenders. Indices are reported against the **data**, not
against the interval subset the inner guard sees.
