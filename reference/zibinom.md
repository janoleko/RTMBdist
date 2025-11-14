# Zero-inflated binomial distribution

Probability mass function, distribution function, and random generation
for the zero-inflated binomial distribution.

## Usage

``` r
dzibinom(x, size, prob, zeroprob = 0, log = FALSE)

pzibinom(q, size, prob, zeroprob = 0, lower.tail = TRUE, log.p = FALSE)

rzibinom(n, size, prob, zeroprob = 0)
```

## Arguments

- x, q:

  vector of quantiles

- size:

  number of trials (zero or more).

- prob:

  probability of success on each trial.

- zeroprob:

  zero-inflation probability between 0 and 1

- log, log.p:

  logical; return log-density if TRUE

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return.

## Value

`dzibinom` gives the probability mass function, `pzibinom` gives the
distribution function, and `rzibinom` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

## Examples

``` r
set.seed(123)
x <- rzibinom(1, size = 10, prob = 0.5, zeroprob = 0.5)
d <- dzibinom(x, size = 10, prob = 0.5, zeroprob = 0.5)
p <- pzibinom(x, size = 10, prob = 0.5, zeroprob = 0.5)
```
