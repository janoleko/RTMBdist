#' Zero-inflated gamma distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the zero-inflated gamma distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param n number of random values to return
#' @param shape positive shape parameter
#' @param scale positive scale parameter
#' @param zeroprob zero-inflation probability between 0 and 1.
#' @param log logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#'
#' @return
#' \code{dzigamma} gives the density, \code{pzigamma} gives the distribution function, \code{qzigamma} gives the quantile function, and \code{rzigamma} generates random deviates.
#'
#' @examples
#' x = rzigamma(1, 1, 1, 0.5)
#' d = dzigamma(x, 1, 1, 0.5)
#' p = pzigamma(x, 1, 1, 0.5)
#' @name zigamma
NULL

#' @rdname zigamma
#' @export
#' @importFrom RTMB dgamma logspace_add
dzigamma = function(x, shape, scale, zeroprob = 0, log = FALSE) {

  if(inherits(x, "simref")) {
    return(dGenericSim("dzigamma", x=x, shape = shape, scale = scale, zeroprob = zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzigamma", x=x, shape = shape, scale = scale, zeroprob = zeroprob, log=log))
  }

  logdens <- numeric(length(x))
  zero_idx <- (x == 0)

  # Zero inflation part
  logdens[zero_idx] <- log(zeroprob)
  logdens[!zero_idx] <- log(1 - zeroprob) +
    dgamma(x[!zero_idx], shape = shape, scale = scale, log = TRUE)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname zigamma
#' @importFrom stats runif rgamma
#' @export
rzigamma <- function(n, shape, scale, zeroprob = 0) {
  if (any(shape < 0)) stop("shape must be >= 0")
  if (any(scale <= 0)) stop("scale must be > 0")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  u <- runif(n)
  res <- ifelse(u < zeroprob, 0, rgamma(n, shape = shape, scale = scale))
  return(res)
}
#' @rdname zigamma
#' @importFrom RTMB pgamma
#' @export
pzigamma <- function(q, shape, scale, zeroprob = 0) {
  if (any(shape < 0)) stop("shape must be >= 0")
  if (any(scale <= 0)) stop("scale must be > 0")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  cdf <- numeric(length(q))

  below_zero <- q < 0
  is_zero <- q == 0
  positive <- q > 0

  cdf[below_zero] <- 0
  cdf[is_zero] <- zeroprob
  cdf[positive] <- zeroprob + (1 - zeroprob) * pgamma(q[positive], shape, scale)

  return(cdf)
}
