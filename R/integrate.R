#' Compute the area of a function with values \code{y} at the points \code{x}.
#'
#' Wrapper function for trapezoidal integration from pracma.
#'
#' @param x coordinates of points on the \code{x}-axis
#' @param y coordinates of function values
#'
#' @keywords internal
integrate <- function(x, y) {
  return(pracma::trapz(x, y))
}
