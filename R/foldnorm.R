#' Folded normal distribution
#'
#' Density, distribution function, and random generation for
#' the folded normal distribution.
#'
#' @details
#' This implementation of \code{dfoldnorm} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu location parameter
#' @param sigma scale parameter, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dfoldnorm} gives the density, \code{pfoldnorm} gives the distribution function, and \code{rfoldnorm} generates random deviates.
#'
#' @examples
#' x <- rfoldnorm(1, 1, 2)
#' d <- dfoldnorm(x, 1, 2)
#' p <- pfoldnorm(x, 1, 2)
#' @name foldnorm
NULL

#' @rdname foldnorm
#' @export
#' @importFrom RTMB dnorm logspace_add
dfoldnorm <- function(x, mu = 0, sigma = 1, log = FALSE) {

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dfoldnorm", x=x, mu=mu, sigma=sigma, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dfoldnorm", x=x, mu=mu, sigma=sigma, log=log))
  }

  const <- - log(sqrt(2 * pi) * sigma)
  part1 <- - (x - mu)^2 / (2 * sigma^2)
  part2 <- - (x + mu)^2 / (2 * sigma^2)

  logdens <- logspace_add(const + part1, const + part2)

  nonneg <- x >= 0

  if(log){
    logdens[!nonneg] <- -Inf
    return(logdens)
  } else{
    dens <- exp(logdens)
    dens[!nonneg] <- 0
    return(dens)
  }
}
#' @rdname foldnorm
#' @export
pfoldnorm <- function(q, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  # ensure sigma > 0
  # if (sigma <= 0) stop("sigma must be positive")

  below_zero <- q < 0

  denom <- sqrt(2) * sigma
  p <- 0.5 * (erf((q + mu) / denom) + erf((q - mu) / denom))

  p[below_zero] <- 0

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)

  return(p)
}
#' @rdname foldnorm
#' @export
#' @importFrom stats rnorm
rfoldnorm <- function(n, mu = 0, sigma = 1) {
  x <- rnorm(n, mu, sigma)
  abs(x)
}
