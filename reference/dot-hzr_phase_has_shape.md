# Does each phase have shape parameters at all?

A `constant` phase is `mu` and nothing else. The saturated
identifiability message says `mu` survives while the shape parameters go
flat, which is vacuous for a phase that has none – the wording defect in
\#211.

## Usage

``` r
.hzr_phase_has_shape(phases)
```

## Arguments

- phases:

  A list of validated `hzr_phase` objects.

## Value

Logical vector, one element per phase, in `phases` order.

## Details

This asks whether the parameters *exist*, not whether they are free. A
phase whose shapes are pinned with `hzr_phase(fixed = )` still has them,
and the message is true and useful there: it is usually why they were
pinned.
