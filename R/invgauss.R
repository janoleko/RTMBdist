#' Inverse Gaussian distribution
#'
#' Density, distribution function, and random generation for
#' the inverse Gaussian distribution.
#'
#' @details
#' This implementation of \code{dinvgauss} allows for automatic differentiation with \code{RTMB}.
#' \code{qinvgauss} and \code{rinvgauss} are imported from the \code{statmod} package.
#'
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mean location parameter
#' @param shape shape parameter, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#' @param ... additional parameter passed to \code{statmod::qinvgauss} for numerical evaluation of the quantile function.
#'
#' @return
#' \code{dinvgauss} gives the density, \code{pinvgauss} gives the distribution function, \code{qinvgauss} gives the quantile function, and \code{rinvgauss} generates random deviates.
#'
#' @examples
#' x <- rinvgauss(1, 1, 0.5)
#' d <- dinvgauss(x, 1, 0.5)
#' p <- pinvgauss(x, 1, 0.5)
#' q <- qinvgauss(p, 1, 0.5)
#' @name invgauss
NULL

#' @rdname invgauss
#' @export
dinvgauss <- function(x, mean = 1, shape = 1, log = FALSE) {

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dinvgauss", x=x, mean=mean, shape=shape, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dinvgauss", x=x, mean=mean, shape=shape, log=log))
  }

  logdens <- 0.5 * log(shape) - 0.5 * log(2 * pi) - 1.5 * log(x) -
    (shape * (x - mean)^2) / (2 * mean^2 * x)

  if(log) return(logdens)

  return(exp(logdens))
}
#' @rdname invgauss
#' @export
#' @importFrom RTMB pnorm
pinvgauss <- function(q, mean = 1, shape = 1, lower.tail = TRUE, log.p = FALSE) {

  p <- RTMB::pnorm(sqrt(shape / q) * (q / mean - 1)) +
    exp(2 * shape / mean) * RTMB::pnorm(-sqrt(shape / q) * (q / mean + 1))

  if (!lower.tail) p <- 1 - p
  if (log.p) return(log(p))
  return(p)
}
#' @rdname invgauss
#' @export
#' @importFrom statmod rinvgauss
rinvgauss <- function(n, mean = 1, shape = 1) {
  statmod::rinvgauss(n, mean = mean, shape = shape)
}
#' @rdname invgauss
#' @export
#' @importFrom statmod qinvgauss
qinvgauss <- function(p, mean = 1, shape = 1, lower.tail=TRUE, log.p = FALSE, ...) {
  statmod::qinvgauss(p, mean = mean, shape = shape, lower.tail = lower.tail, log.p = log.p, ...)
}
