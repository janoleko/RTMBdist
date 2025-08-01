#' Laplace distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the Laplace distribution.
#'
#' @details
#' This implementation of \code{dlaplace} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu location parameter
#' @param b scale parameter, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dlaplace} gives the density, \code{plaplace} gives the distribution function, \code{qlaplace} gives the quantile function, and \code{rlaplace} generates random deviates.
#'
#' @examples
#' x <- rlaplace(1, 1, 1)
#' d <- dlaplace(x, 1, 1)
#' p <- plaplace(x, 1, 1)
#' q <- qlaplace(p, 1, 1)
#' @name laplace
NULL

#' @rdname laplace
#' @export
#' @import RTMB
dlaplace <- function(x, mu = 0, b = 1, log = FALSE) {
  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dlaplace", x=x, mu=mu, b=b, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dlaplace", x=x, mu=mu, b=b, log=log))
  }

  z <- abs(x - mu) / b
  logdens <- -z - log(2 * b)

  if(log) return(logdens)
  return(exp(logdens))
}
#' @rdname laplace
#' @export
plaplace <- function(q, mu = 0, b = 1, lower.tail = TRUE, log.p = FALSE) {
  # ensure b > 0
  # if (b <= 0) stop("b must be strictly positive.")

  z <- (q - mu) / b
  p <- ifelse(z < 0, 0.5 * exp(z), 1 - 0.5 * exp(-z))

  if (!lower.tail) p <- 1 - p

  if (log.p) return(log(p))
  return(p)
}
#' @rdname laplace
#' @export
qlaplace <- function(p, mu = 0, b = 1, lower.tail = TRUE, log.p = FALSE) {
  # ensure b > 0
  if (b <= 0) stop("b must be strictly positive.")

  if (log.p) p <- exp(p)

  if (!lower.tail) p <- 1 - p

  z <- ifelse(p < 0.5, log(2 * p), -log(2 * (1 - p)))
  return(mu + b * z)
}
#' @rdname laplace
#' @export
rlaplace <- function(n, mu = 0, b = 1) {
  # ensure b > 0
  if (b <= 0) stop("b must be strictly positive.")

  u <- runif(n)
  z <- ifelse(u < 0.5, log(2 * u), -log(2 * (1 - u)))
  return(mu + b * z)
}


