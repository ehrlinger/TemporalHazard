# OMC: Open Mitral Commissurotomy

Data for 339 patients who underwent open mitral commissurotomy at the
University of Alabama Birmingham. Contains repeated thromboembolic
events (up to 3 per patient) with left censoring, exercising the
interval censoring likelihood.

## Usage

``` r
omc
```

## Format

A data frame with 339 rows and 7 variables:

- study:

  Patient identifier

- te1:

  Indicator for first thromboembolic event

- te2:

  Indicator for second thromboembolic event

- te3:

  Indicator for third thromboembolic event

- int_dead:

  Follow-up interval to death or last contact (months)

- dead:

  Death indicator (1 = dead, 0 = censored)

- opdjul:

  Operation date (Julian)

## Source

University of Alabama Birmingham cardiac surgery registry.

## See also

Other datasets:
[`avc`](https://ehrlinger.github.io/temporal_hazard/reference/avc.md),
[`cabgkul`](https://ehrlinger.github.io/temporal_hazard/reference/cabgkul.md),
[`tga`](https://ehrlinger.github.io/temporal_hazard/reference/tga.md),
[`uslife2023`](https://ehrlinger.github.io/temporal_hazard/reference/uslife2023.md),
[`valves`](https://ehrlinger.github.io/temporal_hazard/reference/valves.md)
