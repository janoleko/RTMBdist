#' Skewed students t distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the skew t distribution (type 2).
#'
#' @details
#' This corresponds to the skew t type 2 distribution in GAMLSS (\code{\link[gamlss.dist]{ST2}}), see pp. 411-412 of Rigby et al. (2019) and the version implemented in the \code{sn} package.
#' This implementation of \code{dskewt} allows for automatic differentiation with \code{RTMB} while the other functions are imported from the \code{sn} package.
#'
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mu location parameter
#' @param sigma scale parameter, must be positive.
#' @param skew skewness parameter, can be positive or negative.
#' @param df degrees of freedom, must be positive.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param tol a scalar value which regulates the accuracy of the result of qsn, measured on the probability scale.
#' @param method an integer value between 0 and 5 which selects the computing method; see ‘Details’ in the \code{\link[sn]{pst}} documentation below for the meaning of these values. If method=0 (default value), an automatic choice is made among the four actual computing methods, depending on the other arguments.
#'
#' @return
#' \code{dskewt} gives the density, \code{pskewt} gives the distribution function, \code{qskewt} gives the quantile function, and \code{rskewt} generates random deviates.
#'
#' @examples
#' x <- rskewt(1, 1, 2, 5, 2)
#' d <- dskewt(x, 1, 2, 5, 2)
#' p <- pskewt(x, 1, 2, 5, 2)
#' q <- qskewt(p, 1, 2, 5, 2)
#' @name skewt
NULL
#' @rdname skewt
#' @export
#' @importFrom RTMB dt
dskewt <- function(x, mu = 0, sigma = 1, skew = 0, df = 1e3, log = FALSE){
  z <- (x - mu) / sigma
  lambda <- (df + 1) / (df + z^2)
  omega <- skew * sqrt(lambda) * z

  pdf <- RTMB::dt(z, df, log = TRUE)
  cdf <- log(pt_ad(omega, df + 1))

  logdens <- log(2) - log(sigma) + pdf + cdf

  if(log) return(logdens)
  return(exp(logdens))
}

#' @import RTMB
pt_ad <- function(q, df) {
  x <- df / (df + q^2)
  q <- q + 1e-8 # avoid numerical issues with q = 0

  val <- RTMB::pbeta(x, df / 2, 0.5) / 2
  test <- 0.5 * (sign(q) + 1)  # test if q > 0
  test * (1 - val) + (1 - test) * val
}

#' @rdname skewt
#' @export
#' @usage
#' pskewt(q, mu = 0, sigma = 1, skew = 0, df = 1000,
#' method = 0, lower.tail = TRUE, log.p = FALSE)
#' @importFrom sn pst
pskewt <- function(q, mu = 0, sigma = 1, skew = 0, df = 1e3, method = 0, lower.tail = TRUE, log.p = FALSE) {
  # ensure sigma, df > 0
  if (sigma <= 0) stop("sigma must be strictly positive.")
  if (df <= 0) stop("df must be strictly positive.")

  pst(q, xi=mu, omega=sigma, alpha=skew, nu=df, method=method, lower.tail=lower.tail, log.p=log.p)
}

#' @rdname skewt
#' @export
#' @usage
#' qskewt(p, mu = 0, sigma = 1, skew = 0, df = 1000,
#' tol = 1e-8, method = 0)
#' @importFrom sn qst
qskewt <- function(p, mu = 0, sigma = 1, skew = 0, df = 1e3, tol = 1e-8, method = 0) {
  # ensure sigma, df > 0
  if (sigma <= 0) stop("sigma must be strictly positive.")
  if (df <= 0) stop("df must be strictly positive.")
  qst(p, xi=mu, omega=sigma, alpha=skew, nu=df, tol=tol, method=method)
}

#' @rdname skewt
#' @export
#' @importFrom sn rst
rskewt <- function(n, mu = 0, sigma = 1, skew = 0, df = 1e3) {
  # ensure sigma, df > 0
  if (sigma <= 0) stop("sigma must be strictly positive.")
  if (df <= 0) stop("df must be strictly positive.")

  as.numeric(rst(n, xi=mu, omega=sigma, alpha=skew, nu=df))
}

