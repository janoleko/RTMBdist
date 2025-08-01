#' Zero-inflated Poisson distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the zero-inflated Poisson distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q integer vector of counts
#' @param n number of random values to return.
#' @param lambda vector of (non-negative) means
#' @param zeroprob zero-inflation probability between 0 and 1
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dzipois} gives the probability mass function, \code{pzipois} gives the distribution function, and \code{rzipois} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rzipois(1, 0.5, 1)
#' d <- dzipois(x, 0.5, 1)
#' p <- pzipois(x, 0.5, 1)
#' @name zipois
NULL
#' @rdname zipois
#' @export
#' @importFrom RTMB logspace_add dpois
dzipois <- function(x, lambda, zeroprob = 0, log = FALSE) {

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dzipois", x = x, lambda = lambda, zeroprob = zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzipois", x = x, lambda = lambda, zeroprob = zeroprob, log=log))
  }

  logdens <- numeric(length(x))
  zero_idx <- (x == 0)

  # Zero inflation part
  logdens[zero_idx] <- logspace_add(log(zeroprob), log(1-zeroprob) - lambda)
  logdens[!zero_idx] <- log(1 - zeroprob) + RTMB::dpois(x[!zero_idx], lambda, log = TRUE)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname zipois
#' @importFrom RTMB ppois
#' @export
pzipois <- function(q, lambda, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {
  # ensure lambda >= 0, zeroprob in [0,1]
  # if (any(lambda < 0)) stop("lambda must be >= 0")
  # if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  q <- floor(q)  # make sure it's integer-valued
  cdf <- numeric(length(q))

  below_zero <- q < 0
  is_zero <- q == 0
  positive <- q > 0

  cdf[below_zero] <- 0
  cdf[is_zero] <- zeroprob + (1 - zeroprob) * ppois(0, lambda)
  cdf[positive] <- zeroprob + (1 - zeroprob) * ppois(q[positive], lambda)

  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  return(cdf)
}
#' @rdname zipois
#' @importFrom stats runif rpois
#' @export
rzipois <- function(n, lambda, zeroprob = 0) {
  # ensure lambda >= 0, zeroprob in [0,1]
  if (any(lambda < 0)) stop("lambda must be >= 0")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  u <- runif(n)
  res <- ifelse(u < zeroprob, 0, rpois(n, lambda))
  return(res)
}
