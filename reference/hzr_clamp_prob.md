# Clamp probabilities away from 0 and 1

Bounds a probability vector to \\\[\epsilon,\\ 1 - \epsilon\]\\ so that
downstream \\\log p\\ or \\\log(1 - p)\\ evaluations stay finite. Used
throughout the likelihood and prediction code to guard survival and CDF
values against exact 0 or 1.

## Usage

``` r
hzr_clamp_prob(p, eps = 1e-12)
```

## Arguments

- p:

  Numeric vector of probabilities.

- eps:

  Small positive tolerance in \\(0, 0.5)\\; default `1e-12`.

## Value

Numeric vector bounded to `[eps, 1 - eps]`.

## See also

[`hzr_log1pexp()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_log1pexp.md)
and
[`hzr_log1mexp()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_log1mexp.md)
for the companion stable-logarithm primitives.

## Examples

``` r
hzr_clamp_prob(c(0, 0.5, 1))
#> [1] 1e-12 5e-01 1e+00
```
