#' Zero-inflated negative binomial distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the zero-inflated negative binomial distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of (non-negative integer) quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param size size parameter, must be positive.
#' @param prob mean parameter, must be positive.
#' @param zeroprob zero-inflation probability between 0 and 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dzinbinom} gives the density, \code{pzinbinom} gives the distribution function, and \code{rzinbinom} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rzinbinom(1, 2, 0.5)
#' d <- dzinbinom(x, 2, 0.5)
#' p <- pzinbinom(x, 2, 0.5)
#' @name zinbinom
NULL
#' @rdname zinbinom
#' @export
#' @importFrom RTMB dnbinom
dzinbinom <- function(x, size, prob, zeroprob = 0, log = FALSE) {

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dzinbinom", x = x, size=size, prob=prob, zeroprob=zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzinbinom", x = x, size=size, prob=prob, zeroprob=zeroprob, log=log))
  }

  logdens <- numeric(length(x))
  zero_idx <- (x == 0)

  # Zero inflation part
  logdens[zero_idx] <- logspace_add(log(zeroprob),
                                    log(1-zeroprob) + RTMB::dnbinom(0, size=size, prob=prob, log = TRUE))
  logdens[!zero_idx] <- log(1 - zeroprob) + RTMB::dnbinom(x[!zero_idx], size=size, prob=prob, log = TRUE)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname zinbinom
#' @importFrom stats runif rnbinom
#' @export
rzinbinom <- function(n, size, prob, zeroprob = 0) {
  # ensure size >= 0, prob in (0,1], zeroprob in [0,1]
  if (any(size <= 0)) stop("size must be > 0")
  if (any(prob <= 0 | prob > 1)) stop("prob must be in (0,1]")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  u <- runif(n)
  res <- ifelse(u < zeroprob, 0, rnbinom(n, size=size, prob=prob))
  return(res)
}
#' @rdname zinbinom
#' @importFrom stats pnbinom
#' @export
pzinbinom <- function(q, size, prob, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {
  # [validation code...]
  q <- floor(q)
  cdf <- numeric(length(q))

  below_zero <- q < 0
  is_zero <- q == 0
  positive <- q > 0

  cdf[below_zero] <- 0
  cdf[is_zero] <- zeroprob + (1 - zeroprob) * pnbinom(0, size=size, prob=prob)
  cdf[positive] <- zeroprob + (1 - zeroprob) * pnbinom(q[positive], size=size, prob=prob)

  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  return(cdf)
}

