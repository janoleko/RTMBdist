
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RTMBdist

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/janoleko/RTMBdist/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/janoleko/RTMBdist/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The [`RTMB`](https://kaskr.r-universe.dev/RTMB) package enables powerful
and flexible statistical modelling with rich random effect structures
using automatic differentiation (AD). However, its built-in support for
probability distributions is limited to standard cases. `RTMBdist` fills
this gap by providing a collection of non-standard, AD-compatible
distributions, extending the range of models that can be implemented and
estimated with `RTMB`. All the distributions implemented in `RTMBdist`
allow for automatic simulation and residual calculation by `RTMB`.

The full list of distributions currently available is given in the [List
of
distributions](https://janoleko.github.io/RTMBdist/articles/distlist.html)
vignette.

## Installation

You can install the development version of `RTMBdist` from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("janoleko/RTMBdist")
```

## Example

Let’s pretend we want to do numerical maximum likelihood estimation
(MLE) for a gamma distribution that is parameterised in terms of mean
and standard deviation, which is available with the `gamma2` family:

``` r
library(RTMBdist)
#> Loading required package: RTMB
#> 
#> Attaching package: 'RTMBdist'
#> The following object is masked from 'package:RTMB':
#> 
#>     dnbinom2

# simulate data
x <- rgamma2(100, mean = 5, sd = 2)

# negative log-likelihood function
nll <- function(par) {
  x <- OBS(x) # mark x as the response
  mu <- exp(par[1]); ADREPORT(mu)
  sigma <- exp(par[2]); ADREPORT(sigma)
  -sum(dgamma2(x, mu, sigma, log = TRUE))
}

# automatically differentiable objective function object
obj <- MakeADFun(nll, c(log(5), log(2)), silent = TRUE)

# model fitting
opt <- nlminb(obj$par, obj$fn, obj$gr)

# model summary
summary(sdreport(obj))
#>        Estimate Std. Error
#> par   1.5705170 0.03480726
#> par   0.5151734 0.07757581
#> mu    4.8091338 0.16739278
#> sigma 1.6739287 0.12985637
```

Through the magic of `RTMB`, we can also immediately simulate new data
from the fitted model and calculate residuals:

``` r
# simulate new data
x_new <- obj$simulate()

# calculate residuals
osa <- oneStepPredict(obj, method = "cdf", trace = FALSE)
qqnorm(osa$res); abline(0, 1)
```

<img src="man/figures/README-sim_residuals-1.png" width="100%" />
