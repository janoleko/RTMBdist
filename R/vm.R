#' von Mises distribution
#'
#' Density, distribution function and random generation for the von Mises distribution.
#'
#' @details
#' This implementation of \code{dvm} allows for automatic differentiation with \code{RTMB}.
#' \code{rvm} and \code{pvm} are simply wrappers of the corresponding functions from \code{circular}.
#'
#' @param x,q vector of angles measured in radians at which to evaluate the density function.
#' @param mu mean direction of the distribution measured in radians.
#' @param kappa non-negative numeric value for the concentration parameter of the distribution.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#' @param n number of random values to return.
#' @param tol the precision in evaluating the distribution function
#' @param from value from which the integration for CDF starts. If \code{NULL}, is set to \code{mu - pi}.
#' @param wrap logical; if \code{TRUE}, generated angles are wrapped to the interval from -pi to pi.
#'
#' @return \code{dvm} gives the density, \code{pvm} gives the distribution function, and \code{rvm} generates random deviates.
#'
#' @examples
#' set.seed(1)
#' x <- rvm(10, 0, 1)
#' d <- dvm(x, 0, 1)
#' p <- pvm(x, 0, 1)
#' @name vm
NULL

#' @rdname vm
#' @export
#' @importFrom RTMB besselI
dvm = function(x, mu = 0, kappa = 1, log = FALSE) {

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dvm", x = x, mu = mu, kappa = kappa, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dvm", x = x, mu = mu, kappa = kappa, log=log))
  }

  logdens <- -log(2 * pi) -
    log(besselI(kappa, 0)) +
    kappa * cos(x - mu)

  if(log){
    return(logdens)
  } else{
    return(exp(logdens))
  }
}

#' @rdname vm
#' @export
#' @importFrom circular pvonmises
pvm = function(q, mu = 0, kappa = 1, from = NULL, tol = 1e-20) {
  # NA handling
  ind = which(!is.na(q))

  if(is.matrix(mu)){
    mu = mu[ind,]
  }
  if(is.matrix(kappa)){
    kappa = kappa[ind,]
  }

  probs = numeric(length(q))

  suppressWarnings(
    probs[ind] <- pvonmises(q[ind], mu, kappa, from = from, tol = tol)
  )

  probs[-ind] = NA

  as.numeric(probs)
}

#' @rdname vm
#' @export
#' @importFrom circular rvonmises
rvm = function(n, mu = 0, kappa = 1, wrap = TRUE) {
  suppressWarnings(
    angles <- as.numeric(rvonmises(n, mu, kappa))
  )

  # if generated angels should be wrapped, i.e. mapped to interval [-pi, pi], do so
  if(wrap){
    angles = (angles + pi) %% (2 * pi) - pi
  }
  angles
}
