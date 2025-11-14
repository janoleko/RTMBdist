# Zero-truncated Binomial distribution

Probability mass function, distribution function, and random generation
for the zero-truncated Binomial distribution.

## Usage

``` r
dztbinom(x, size, prob, log = FALSE)

pztbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)

rztbinom(n, size, prob)
```

## Arguments

- x, q:

  integer vector of counts

- size:

  number of trials

- prob:

  success probability in each trial

- log, log.p:

  logical; return log-density if TRUE

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return.

## Value

`dztbinom` gives the probability mass function, `pztbinom` gives the
distribution function, and `rztbinom` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

By definition, this distribution only has support on the positive
integers (1, ..., size). Any zero-truncated distribution is defined as
\$\$P(X=x \| X\>0) = P(X=x) / (1 - P(X=0)),\$\$ where \\P(X=x)\\ is the
probability mass function of the corresponding untruncated distribution.

## Examples

``` r
set.seed(123)
x <- rztbinom(1, size = 10, prob = 0.3)
d <- dztbinom(x, size = 10, prob = 0.3)
p <- pztbinom(x, size = 10, prob = 0.3)
```
