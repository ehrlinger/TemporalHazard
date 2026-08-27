# Numerically stable log(1 + exp(x))

Computes the "softplus" function \$\$\mathrm{log1pexp}(x) =
\log\\\bigl(1 + e^{x}\bigr)\$\$ without overflow. A naive
`log(1 + exp(x))` overflows to `Inf` for \\x \gtrsim 710\\ even though
the true value is \\\approx x\\. The identity \\\log(1 + e^{x}) = x +
\log(1 + e^{-x})\\ for \\x \> 0\\ keeps every intermediate finite.

## Usage

``` r
hzr_log1pexp(x)
```

## Arguments

- x:

  Numeric vector.

## Value

Numeric vector with element-wise \\\log(1 + e^{x})\\.

## References

Mächler M (2012). *Accurately Computing \\\log(1 - e^{-\|a\|})\\.* R
package Rmpfr vignette.
<https://CRAN.R-project.org/package=Rmpfr/vignettes/log1mexp-note.pdf>

## See also

[`hzr_log1mexp()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_log1mexp.md)
for the complementary \\\log(1 - e^{-x})\\,
[`hzr_clamp_prob()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_clamp_prob.md)
for boundary-safe probabilities.

## Examples

``` r
hzr_log1pexp(c(-50, 0, 50))
#> [1] 1.928750e-22 6.931472e-01 5.000000e+01
```
