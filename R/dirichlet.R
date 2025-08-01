#' Dirichlet distribution
#'
#' Density of the Dirichlet distribution.
#'
#' @details
#' This implementation of \code{ddirichlet} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x vector or matrix of quantiles. If \code{x} is a vector, it needs to sum to one.
#' If \code{x} is a matrix, each row should sum to one.
#' @param alpha vector or matrix of positive shape parameters
#' @param log logical; if \code{TRUE}, densities \eqn{p} are returned as \eqn{\log(p)}.
#'
#' @return
#' \code{ddirichlet} gives the density.
#'
#' @examples
#' # single alpha
#' alpha <- c(1,2,3)
#' x <- rdirichlet(1, alpha)
#' d <- ddirichlet(x, alpha)
#' # vectorised over alpha
#' alpha <- rbind(alpha, 2*alpha)
#' x <- rdirichlet(2, alpha)
#  d <- ddirichlet(x, alpha)
#' @name dirichlet
NULL

#' @rdname dirichlet
#' @export
#' @import RTMB
ddirichlet <- function(x, alpha, log = TRUE) {

  # potentially escape to RNG or produce error for CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("ddirichlet", x=x, alpha=alpha, log=log))
  }
  if(inherits(x, "osa")) {
    stop("Dirichlet does not support OSA residuals.")
  }

  # Check if x and alpha are vectors by checking if they have dimensions
  if (is.null(dim(x))) x <- matrix(x, nrow = 1)
  if (is.null(dim(alpha))) alpha <- matrix(alpha, nrow = 1)

  # Ensure x and alpha have the same number of columns
  if (ncol(x) != ncol(alpha)) {
    stop("x and alpha must have the same number of columns (categories).")
  }

  # Compute log of the multivariate beta function B(alpha) for each row
  log_B_alpha <- rowSums(lgamma(alpha)) - lgamma(rowSums(alpha))

  # Compute log density for each row
  log_density <- rowSums((alpha - 1) * log(x)) - log_B_alpha

  log_density <- as.numeric(log_density)

  # Return result based on the 'log' argument
  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}
#' @rdname dirichlet
#' @export
#' @importFrom stats rgamma
rdirichlet <- function(n, alpha) {

  if(is.vector(alpha)){
    longalpha <- rep(alpha, n)
    k <- length(alpha)
  } else if(is.matrix(alpha)){
    if(nrow(alpha) != n) stop("Number of rows in alpha must match n.")

    longalpha <- as.vector(t(alpha))
    k <- ncol(alpha)
  }

  x <- rgamma(n * k, shape = longalpha, scale = 1)
  x <- matrix(x, nrow = n, ncol = k, byrow = TRUE)
  x / rowSums(x)
}
