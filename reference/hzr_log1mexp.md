# Numerically stable log(1 - exp(-x)) for x \> 0

Computes \$\$\mathrm{log1mexp}(x) = \log\\\bigl(1 - e^{-x}\bigr), \qquad
x \> 0\$\$ accurately across the full range. Following Mächler (2012),
the evaluation switches at \\x = \log 2\\: for small \\x\\ it uses
\\\log(-\mathrm{expm1}(-x))\\ (avoiding cancellation as \\x \to 0\\),
and for larger \\x\\ it uses \\\mathrm{log1p}(-e^{-x})\\ (avoiding loss
of precision as \\x \to \infty\\). Non-positive or non-finite inputs
return `NA`.

## Usage

``` r
hzr_log1mexp(x)
```

## Arguments

- x:

  Numeric vector with positive values.

## Value

Numeric vector with element-wise \\\log(1 - e^{-x})\\.

## References

Mächler M (2012). *Accurately Computing \\\log(1 - e^{-\|a\|})\\.* R
package Rmpfr vignette.
<https://CRAN.R-project.org/package=Rmpfr/vignettes/log1mexp-note.pdf>

## See also

[`hzr_log1pexp()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_log1pexp.md)
for the softplus \\\log(1 + e^{x})\\,
[`hzr_clamp_prob()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_clamp_prob.md)
for boundary-safe probabilities.

## Examples

``` r
hzr_log1mexp(c(0.01, 0.5, 5))
#> [1] -4.610166019 -0.932752130 -0.006760749
```
