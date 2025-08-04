# dwishart <- function(X, Sigma, nu) {
#
#   p <- nrow(X)
#   if(ncol(X) != p) stop("X must be a square matrix")
#
#   logdetX <- determinant(X, logarithm = TRUE)
#   logdetX <- logdetX$modulus * logdetX$sign
#
#   logdetSigma <- determinant(Sigma, logarithm = TRUE)
#   logdetSigma <- logdetSigma$modulus
#
#   SigmaInv <- solve(Sigma)
#
#   lognum <- (nu - p - 1) / 2 * logdetX - sum(rowSums(SigmaInv * X)) / 2
#
#   logdenom <- nu * p / 2 * log(2) + nu / 2 * logdetSigma + lgamma(nu / 2)
# }
