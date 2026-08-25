# Reject row types the SAS objective has no counterpart for

`PROC HAZARD` has no left-censoring statement, so no SAS run corresponds
to a fit containing left-censored rows. That is a data defect rather
than parameter infeasibility, so it stops rather than returning `-Inf`.

## Usage

``` r
.hzr_check_sas_status(status, objective)
```

## Arguments

- status:

  Numeric event indicator.

- objective:

  Resolved objective, `"likelihood"` or `"sas"`.

## Value

`NULL`, invisibly; called for its side effect.

## Details

Called from BOTH
[`.hzr_logl_multiphase()`](https://ehrlinger.github.io/temporal_hazard/reference/dot-hzr_logl_multiphase.md)
and
[`.hzr_gradient_multiphase()`](https://ehrlinger.github.io/temporal_hazard/reference/dot-hzr_gradient_multiphase.md).
Guarding only the objective would leave the gradient computing happily
for data the objective refuses – and the gradient is reachable on its
own, for instance from `.hzr_score_test()`, so the objective's refusal
is not guaranteed to come first.
