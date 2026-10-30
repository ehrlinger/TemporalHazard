# Normal deviates from a private Mersenne Twister

Reproduces `set.seed(seed); rnorm(n)` under R's default
[`RNGkind()`](https://rdrr.io/r/base/Random.html), bit for bit, without
touching `.Random.seed`. Written out rather than called because calling
it means writing to the global environment, and changing the numbers
would change which optimum existing multiphase fits reach. The 32-bit
words are held in doubles: every product stays below 2^53, so the
arithmetic is exact.

## Usage

``` r
.hzr_mt_normal(seed, n)
```

## Arguments

- seed:

  Single whole number within integer range.

- n:

  Number of deviates.

## Value

Numeric vector of length `n`.
