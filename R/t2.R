#' Non-central and scaled students t distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the t distribution with non-centrality and scale parameters.
#'
#' @details
#' This implementation of \code{dt2} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mu location parameter
#' @param sigma scale parameter, must be positive.
#' @param df degrees of freedom, must be positive.
#' @param log logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#'
#' @return
#' \code{dt2} gives the density, \code{pt2} gives the distribution function, \code{qt2} gives the quantile function, and \code{rt2} generates random deviates.
#'
#' @examples
#' x <- rt2(1, 1, 2, 5)
#' d <- dt2(x, 1, 2, 5)
#' p <- pt2(x, 1, 2, 5)
#' q <- qt2(p, 1, 2, 5)
#' @name t2
NULL
#' @rdname t2
#' @export
#' @importFrom RTMB dt
dt2 = function(x, mu, sigma, df, log = FALSE){

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dt2", x=x, mu=mu, sigma=sigma, df=df, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dt2", x=x, mu=mu, sigma=sigma, df=df, log=log))
  }

  z <- (x - mu) / sigma
  logdens <- RTMB::dt(z, df, log = TRUE) - log(sigma)

  if(log) return(logdens)
  return(exp(logdens))
}
#' @rdname t2
#' @export
#' @importFrom stats pt
pt2 <- function(q, mu, sigma, df){
  z <- (q - mu) / sigma
  # stats::pt(z, df)
  pt_ad(z, df)
}
#' @rdname t2
#' @export
#' @importFrom stats rt
rt2 <- function(n, mu, sigma, df) {
  z <- stats::rt(n, df)
  return(mu + sigma * z)
}
#' @rdname t2
#' @export
#' @importFrom stats qt
qt2 <- function(p, mu, sigma, df){
  z <- stats::qt(p, df)
  return(mu + sigma * z)
}
