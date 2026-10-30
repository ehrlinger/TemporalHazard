# Deterministic multi-start perturbations

Draws the offsets that shift the starting values for every multiphase
start after the first. They are the numbers `set.seed(seed); rnorm()`
gives under R's default generators, computed by
[`.hzr_mt_normal()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_mt_normal.md)
without R's stream, so a fit is reproducible without
[`set.seed()`](https://rdrr.io/r/base/Random.html), does not depend on
[`RNGkind()`](https://rdrr.io/r/base/Random.html), and neither reads nor
writes `.Random.seed`.

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
