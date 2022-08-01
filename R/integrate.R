#' Compute the area of a function with values 'y' at the points 'x'.
#'
#' Wrapper function for trapezoidal integration from pracma.
#'
#' @param x coordinates of points on the x-axis
#' @param y coordinates of function values
#' @export
integrate <- function(x, y) {
  return(pracma::trapz(x, y))
}
