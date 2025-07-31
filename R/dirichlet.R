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
#' ddirichlet(c(0.2, 0.3, 0.5), c(1, 2, 3))
#' @name dirichlet
NULL

#' @rdname dirichlet
#' @export
#' @import RTMB
ddirichlet <- function(x, alpha, log = TRUE) {
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

  # Return result based on the 'log' argument
  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}
