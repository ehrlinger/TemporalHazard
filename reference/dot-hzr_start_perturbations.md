# Deterministic multi-start perturbations

Draws the offsets that shift the starting values for every multiphase
start after the first. The draw is taken from an internally seeded
stream and the ambient `.Random.seed` is restored on exit, so a fit is
reproducible without [`set.seed()`](https://rdrr.io/r/base/Random.html)
and never advances the caller's stream. Both halves matter: seeding
alone would still let an intervening draw change the fit, and restoring
alone would leave the fit dependent on whatever the stream happened to
hold when it was called.

## Usage

``` r
.hzr_start_perturbations(n_starts, n_par, seed, sd = 0.5)
```

## Arguments

- n_starts:

  Number of optimization starts.

- n_par:

  Length of the free-parameter vector being perturbed.

- seed:

  Single whole number selecting the ensemble of offsets.

- sd:

  Standard deviation of the offsets (default 0.5).

## Value

A list of `max(0, n_starts - 1)` numeric vectors of length `n_par`;
empty for a single start, which is never perturbed.
