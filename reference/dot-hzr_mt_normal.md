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

## Details

The generator is that of Matsumoto and Nishimura (1998). The seeding and
the inversion reproduce R's
[`set.seed()`](https://rdrr.io/r/base/Random.html) and
[`rnorm()`](https://rdrr.io/r/stats/Normal.html), and a test pins them
against those functions.

## References

Matsumoto M, Nishimura T (1998). Mersenne Twister: a 623-dimensionally
equidistributed uniform pseudo-random number generator. *ACM
Transactions on Modeling and Computer Simulation*, 8(1), 3-30.
[doi:10.1145/272991.272995](https://doi.org/10.1145/272991.272995)
