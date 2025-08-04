#' Zero-inflated beta distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the zero-inflated beta distribution.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param n number of random values to return.
#' @param shape1,shape2 non-negative shape parameters of the beta distribution
#' @param zeroprob zero-inflation probability between 0 and 1.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dzibeta} gives the density, \code{pzibeta} gives the distribution function, and \code{rzibeta} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rzibeta(1, 2, 2, 0.5)
#' d <- dzibeta(x, 2, 2, 0.5)
#' p <- pzibeta(x, 2, 2, 0.5)
#' @name zibeta
NULL
#' @rdname zibeta
#' @export
#' @importFrom RTMB dbeta
dzibeta <- function(x, shape1, shape2, zeroprob = 0, log = FALSE) {

  if(!ad_context()) {
    # shapes non-negative
    if (any(shape1 < 0)) stop("shape1 must be non-negative.")
    if (any(shape2 < 0)) stop("shape2 must be non-negative.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dzibeta", x=x, shape1=shape1, shape2=shape2, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dzibeta", x=x, shape1=shape1, shape2=shape2, log=log))
  }

  # not differentiable in x
  logdens <- numeric(length(x))
  zero_idx <- (x == 0)

  # Zero inflation part
  logdens[zero_idx] <- log(zeroprob)
  logdens[!zero_idx] <- log(1 - zeroprob) +
    RTMB::dbeta(x[!zero_idx], shape1 = shape1, shape2 = shape2, log = TRUE)

  if (log) return(logdens)
  return(exp(logdens))
}

#' @rdname zibeta
#' @export
#' @importFrom RTMB pbeta
pzibeta <- function(q, shape1, shape2, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # shapes non-negative
    if (any(shape1 < 0)) stop("shape1 must be non-negative.")
    if (any(shape2 < 0)) stop("shape2 must be non-negative.")
  }

  s1 <- 2 * sign(q) - 1 # gives -3 for q < 0, -1 for q == 0, and 1 for q > 0
  s2 <- sign(3 + s1) # only zero or 1

  iszero <- 0.5 * (1-s1) * s2
  ispos <- 0.5 * (1+s1) * s2

  p <- iszero * zeroprob +
    ispos * (zeroprob + (1 - zeroprob) * RTMB::pbeta(q, shape1, shape2))

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}
#' @rdname zibeta
#' @export
#' @importFrom stats rbeta
rzibeta <- function(n, shape1, shape2, zeroprob = 0) {
  if(!ad_context()) {
    # shapes non-negative
    if (any(shape1 < 0)) stop("shape1 must be non-negative.")
    if (any(shape2 < 0)) stop("shape2 must be non-negative.")
  }
  u <- runif(n)
  res <- ifelse(u < zeroprob, 0, rbeta(n, shape1, shape2))
  return(res)
}

