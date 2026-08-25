# TGA: Transposition of the Great Arteries

Survival data for 470 patients who underwent the arterial switch
operation for transposition of the great arteries at Boston Children's
Hospital and Children's Hospital of Philadelphia. Used for sensitivity
analysis and internal validation examples.

## Usage

``` r
tga
```

## Format

A data frame with 470 rows and 14 variables:

- study:

  Patient identifier

- simple:

  Simple TGA indicator (0/1)

- dextroin:

  D-looped transposition indicator (0/1)

- ca_1rl2c:

  Coronary artery pattern indicator

- hyaaproc:

  Hybrid approach procedure indicator (0/1)

- no_tca:

  No total circulatory arrest indicator (0/1)

- tca_time:

  Total circulatory arrest time (minutes)

- age_days:

  Age at operation (days)

- arciopyr:

  Aortic cross-clamp time per year

- dead:

  Death indicator (1 = dead, 0 = censored)

- int_dead:

  Follow-up interval to death or last contact (months)

- source:

  Source institution (BCH or CHOP)

- ca1_2_l:

  Coronary artery configuration (1/2/L)

- opyear:

  Year of operation

## Source

Boston Children's Hospital and Children's Hospital of Philadelphia.

## See also

Other datasets:
[`avc`](https://ehrlinger.github.io/temporal_hazard/reference/avc.md),
[`cabgkul`](https://ehrlinger.github.io/temporal_hazard/reference/cabgkul.md),
[`omc`](https://ehrlinger.github.io/temporal_hazard/reference/omc.md),
[`uslife2023`](https://ehrlinger.github.io/temporal_hazard/reference/uslife2023.md),
[`valves`](https://ehrlinger.github.io/temporal_hazard/reference/valves.md)
