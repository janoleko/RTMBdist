#' Reparameterised negative binomial distribution
#'
#' Probability mass function, distribution function, quantile function and random generation for
#' the negative binomial distribution reparameterised in terms of mean and size.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mu mean parameter, must be positive.
#' @param size size parameter, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dnbinom2} gives the density, \code{pnbinom2} gives the distribution function, \code{qnbinom2} gives the quantile function, and \code{rnbinom2} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rnbinom2(1, 1, 2)
#' d <- dnbinom2(x, 1, 2)
#' p <- pnbinom2(x, 1, 2)
#' q <- qnbinom2(p, 1, 2)
#' @name nbinom2
NULL
#' @rdname nbinom2
#' @export
#' @importFrom RTMB dnbinom
dnbinom2 <- function(x, mu, size, log = FALSE) {

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dnbinom2", x=x, mu=mu, size=size, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dnbinom2", x=x, mu=mu, size=size, log=log))
  }

  # parameter transformation
  prob <- size / (size + mu)

  RTMB::dnbinom(x, size = size, prob = prob, log = log)
}
#' @rdname nbinom2
#' @export
#' @importFrom stats pnbinom
pnbinom2 <- function(q, mu, size, lower.tail = TRUE, log.p = FALSE) {
  stats::pnbinom(q, mu = mu, size = size, lower.tail = lower.tail, log.p = log.p)
}
#' @rdname nbinom2
#' @export
#' @importFrom stats qnbinom
qnbinom2 <- function(p, mu, size, lower.tail = TRUE, log.p = FALSE) {
  stats::qnbinom(p, mu = mu, size = size, lower.tail = lower.tail, log.p = log.p)
}
#' @rdname nbinom2
#' @export
#' @importFrom stats rnbinom
rnbinom2 <- function(n, mu, size) {
  stats::rnbinom(n, mu = mu, size = size)
}
