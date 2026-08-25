# Diagnose phases that have effectively left the model

Two distinct failure modes, with different consequences, both silent:

## Usage

``` r
.hzr_phase_shares(theta, time, phases, covariate_counts, x_list)
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

## Value

`data.frame` with one row per phase: `share` (largest share of
\\\Lambda\\ at any observed time) and `variation` (relative range of the
phase's contribution across observed times, `NA` when the phase carries
covariates – `mu` then varies by row and the two sources of variation
cannot be separated from the contribution alone).

The shape is the same whatever happens: if no observed time carries a
usable total, both columns are `NA` rather than the frame being `NULL`.
`fit$phase_share` is then one type for a caller to handle rather than
two, and "could not be measured" stays distinct from "measured as zero".

## Details

- `"absent"`:

  The phase contributes essentially none of \\\Lambda\\ at any observed
  time – it has not started by the end of follow-up. Its `mu` **and**
  its shape are unidentified.

- `"saturated"`:

  The phase's \\\Phi\\ is effectively constant across the observed times
  – a `cdf` phase whose half-life is far shorter than the first
  observation has already finished. It then contributes \\\mu \cdot \Phi
  \approx \mu\\, a constant offset, so **`mu` stays well identified**
  while the shape parameters (`t_half`, `nu`, `m`) go exactly flat: the
  likelihood is unchanged whether they are pinned or fitted.

Share is taken of \\\Lambda\\, not of \\h\\, because every row type's
contribution runs through \\\Lambda(t)\\. A phase can supply almost none
of the instantaneous hazard late in follow-up and still be perfectly
well identified through the offset it already contributed – which is why
the hazard is the wrong basis for this test.

The **maximum** over times is the right summary rather than the mean: a
phase mattering only in the first week is identified, and a mean over
long follow-up would hide it.
