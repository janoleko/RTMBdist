# getting one step ahead resiudals and simulation functions from RTMB (not exported)
dGenericOSA <- get("dGenericOSA", envir = asNamespace("RTMB"), inherits = FALSE)
dGenericSim <- get("dGenericSim", envir = asNamespace("RTMB"), inherits = FALSE)


#' AD-compatible error function and complementary error function
#'
#' @param x vector of evaluation points
#'
#' @returns \code{erf(x)} returns the error function and \code{erfc(x)} returns the complementary error function.
#'
#' @examples
#' erf(1)
#' erfc(1)
#' @name erf
NULL
#' @rdname erf
#' @export
#' @importFrom RTMB pnorm
erf <- function(x) {
  2 * RTMB::pnorm(x * sqrt(2)) - 1
}
#' @rdname erf
#' @export
erfc <- function(x) {
  1 - erf(x)
}
