#' Joint density under a bivariate copula
#'
#' Computes the joint density (or log-density) of a bivariate distribution
#' constructed from two arbitrary margins combined with a specified copula.
#'
#' @param d1,d2 Marginal density values. If \code{log = TRUE}, supply the log-density.
#' If \code{log = FALSE}, supply the raw density.
#' @param p1,p2 Marginal CDF values. Need not be supplied on log scale.
#' @param copula A function of two arguments (\code{u}, \code{v}) returning
#'   the log copula density \eqn{\log c(u,v)}.
#'   You can either construct this yourself or use the copula constructors available (see details)
#' @param log Logical; if \code{TRUE}, return the log joint density. In this case,
#'   \code{d1} and \code{d2} must be on the log scale.
#'
#' @returns Joint density (or log-density) under the bivariate copula.
#'
#' @details
#' The joint density is
#' \deqn{f(x,y) = c(F_1(x), F_2(y)) \, f_1(x) f_2(y),}
#' where \eqn{F_i} are the marginal CDFs, \eqn{f_i} are the marginal densities,
#' and \eqn{c} is the copula density.
#'
#' The marginal densities \code{d1}, \code{d2} and CDFs \code{p1}, \code{p2}
#' must be differentiable for automatic differentiation (AD) to work.
#'
#' Available copula constructors are:
#' - \code{\link{cgaussian}} (Gaussian copula)
#' - \code{\link{cclayton}} (Clayton copula)
#' - \code{\link{cgumbel}} (Gumbel copula)
#' - \code{\link{cfrank}} (Frank copula)
#'
#' @export
#'
#' @examples
#' # Normal + Exponential margins with Gaussian copula
#' x <- c(0.5, 1); y <- c(1, 2)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pexp(y, 2)
#' dcopula(d1  , d2, p1, p2, copula = cgaussian(0.5), log = TRUE)
#'
#' # Normal + Beta margins with Clayton copula
#' x <- c(0.5, 1); y <- c(0.2, 0.8)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dbeta(y, 2, 1, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pbeta(y, 2, 1)
#' dcopula(d1, d2, p1, p2, copula = cclayton(2), log = TRUE)
#'
#' # Normal + Beta margins with Gumbel copula
#' x <- c(0.5, 1); y <- c(1, 2)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dbeta(y, 2, 1, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pbeta(y, 2, 1)
#' dcopula(d1, d2, p1, p2, copula = cgumbel(1.5), log = TRUE)
#'
#' # Normal + Exponential margins with Frank copula
#' x <- c(0.5, 1); y <- c(1, 2)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pexp(y, 2)
#' dcopula(d1, d2, p1, p2, copula = cfrank(2), log = TRUE)
dcopula <- function(d1, d2, p1, p2, copula, log = FALSE) {

  # ensure numeric stability of uniforms
  eps <- .Machine$double.eps
  p1 <- pmin.ad(pmax.ad(p1, eps), 1 - eps)
  p2 <- pmin.ad(pmax.ad(p2, eps), 1 - eps)

  # evaluate log copula density
  logC <- copula(p1, p2)

  # combine with marginal densities
  if (log) {
    logC + d1 + d2
  } else {
    exp(logC) * d1 * d2
  }
}

#' Gaussian copula constructor
#'
#' Returns a function computing the log density of the bivariate Gaussian copula,
#' intended to be used with \code{\link{dcopula}}.
#'
#' @param rho Correlation parameter (\eqn{-1 < rho < 1}).
#' @return Function of two arguments (u,v) returning log copula density.
#'
#' The Gaussian copula density is
#' \deqn{
#' c(u,v;\rho) = \frac{1}{\sqrt{1-\rho^2}}
#' \exp\left\{-\frac{1}{2(1-\rho^2)} (z_1^2 - 2 \rho z_1 z_2 + z_2^2) + \frac{1}{2}(z_1^2 + z_2^2) \right\},
#' }
#' where \eqn{z_1 = \Phi^{-1}(u)}, \eqn{z_2 = \Phi^{-1}(v)}, and \eqn{-1 < \rho < 1}.
#'
#' @seealso [cclayton()], [cgumbel()], [cfrank()]
#'
#' @importFrom RTMB qnorm
#'
#' @export
#'
#' @examples
#' x <- c(0.5, 1); y <- c(1, 2)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pexp(y, 2)
#' dcopula(d1  , d2, p1, p2, copula = cgaussian(0.5), log = TRUE)
cgaussian <- function(rho) {
  function(u, v) {
    # u, v assumed clipped in dcopula()
    z1 <- qnorm(u)
    z2 <- qnorm(v)

    # compute log copula density using Cholesky-style formula
    logdet <- -0.5 * log1p(-rho^2)  # log(1 - rho^2) stable
    exponent <- -0.5 / (1 - rho^2) * (z1^2 - 2 * rho * z1 * z2 + z2^2) +
      0.5 * (z1^2 + z2^2)

    logdet + exponent
  }
}

#' Clayton copula constructor
#'
#' Returns a function that computes the log density of the bivariate Clayton copula,
#' intended to be used with \code{\link{dcopula}}.
#'
#' The Clayton copula density is
#' \deqn{
#' c(u,v;\theta) = (1+\theta) (uv)^{-(1+\theta)}
#' \left( u^{-\theta} + v^{-\theta} - 1 \right)^{-(2\theta+1)/\theta}, \quad \theta > 0.
#' }
#'
#' @seealso [cgaussian()], [cgumbel()], [cfrank()]
#'
#' @param theta Positive dependence parameter (\eqn{\theta > 0}).
#'
#' @return A function of two arguments (u,v) returning log copula density.
#'
#' @export
#'
#' @examples
#' x <- c(0.5, 1); y <- c(0.2, 0.8)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dbeta(y, 2, 1, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pbeta(y, 2, 1)
#' dcopula(d1, d2, p1, p2, copula = cclayton(2), log = TRUE)
cclayton <- function(theta) {
  function(u, v) {
    log(theta + 1) -
      (theta + 1) * (log(u) + log(v)) -
      ((2 * theta + 1) / theta) * log(u^(-theta) + v^(-theta) - 1)
  }
}

#' Gumbel copula constructor
#'
#' Returns a function that computes the log density of the bivariate Gumbel copula,
#' intended to be used with \code{\link{dcopula}}.
#'
#' The Gumbel copula density
#'
#' \deqn{
#' c(u,v;\theta) = \exp\Big[-\big((-\log u)^\theta + (-\log v)^\theta\big)^{1/\theta}\Big] \cdot h(u,v;\theta),
#' }
#' where \eqn{h(u,v;\theta)} contains the derivative terms ensuring the function is a density.
#'
#' @seealso [cgaussian()], [cclayton()], [cfrank()]
#'
#' @param theta Dependence parameter (\eqn{\theta >= 1}).
#'
#' @return A function of two arguments (u,v) returning log copula density.
#'
#' @export
#'
#' @examples
#' x <- c(0.5, 1); y <- c(1, 2)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dbeta(y, 2, 1, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pbeta(y, 2, 1)
#' dcopula(d1, d2, p1, p2, copula = cgumbel(1.5), log = TRUE)
cgumbel <- function(theta) {
  function(u, v) {
    lu <- -log(u); lv <- -log(v)
    A <- (lu^theta + lv^theta)^(1/theta)

    logC <- -A + (theta - 1) * (log(lu) + log(lv)) +
      (-2 + 2/theta) * log(lu^theta + lv^theta) -
      log(u) - log(v)

    term <- 1 + (theta - 1) * (lu * lv)^theta / (lu^theta + lv^theta)^2
    logC + log(term)
  }
}

#' Frank copula constructor
#'
#' Returns a function computing the log density of the bivariate Frank copula,
#' intended to be used with \code{\link{dcopula}}.
#'
#' The Frank copula density is
#' \deqn{
#' c(u,v;\theta) = \frac{\theta (1-e^{-\theta}) e^{-\theta(u+v)}}
#' {\left[(e^{-\theta u}-1)(e^{-\theta v}-1) + (1 - e^{-\theta}) \right]^2}, \quad \theta \ne 0.
#' }
#'
#' @seealso [cgaussian()], [cclayton()], [cgumbel()]
#'
#' @param theta Dependence parameter (\eqn{\theta = 0}).
#'
#' @return Function of two arguments (u,v) returning log copula density.
#'
#' @export
#'
#' @examples
#' x <- c(0.5, 1); y <- c(1, 2)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pexp(y, 2)
#' dcopula(d1, d2, p1, p2, copula = cfrank(2), log = TRUE)
cfrank <- function(theta) {
  function(u, v) {
    # u, v assumed clipped in (0,1)
    ex <- exp(-theta * u)
    ey <- exp(-theta * v)
    e  <- exp(-theta)

    num <- theta * (1 - e) * exp(-theta * (u + v))
    den <- (1 - ex) * (1 - ey) + e - 1

    log(num / (den * den))
  }
}
