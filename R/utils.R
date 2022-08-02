voigt.superpos <- function(x, params) {
  out <- 0
  for (j in seq_along(params$pos)) {
    out <- out + params$amp[j] * RcppFaddeeva::Voigt(x, params$pos[j], params$gwidth[j], params$lwidth[j])
  }
  return(out)
}

baseline.linear <- function(x, params) {
  return(params$a * x + params$b)
}

#' Signal model for Voigt profiles
#'
#' @description
#' Implements the signal model discussed in `spectralem()`. Doesn't include any \eqn{f_i}.
#'
#' Parameters are passed in a list params with components \code{amp}, \code{pos}, \code{gwidth} and \code{lwidth} which are vectors of the same length.
#' For each the i-th entry yields the parameters of the i-th Voigt profile.
#' Furthermore named components \code{a} and \code{b} describe the slope and offset of the linear background respectively. See examples.
#' @examples
#' y <- voigt.model(
#'   x = seq(0, 1000),
#'   params = list(
#'     amp = c(10, 2, 3, 3),
#'     pos = c(332, 577, 697),
#'     gwidth = c(5, 5, 5),
#'     lwidth = c(0.1, 0.1, 0.1),
#'     a = 0.01,
#'     b = 15
#'   )
#' )
#' @param x vector of \code{x} coordinates the model is to be evaluated at
#' @param params parameter list as described in the description
#' @export
voigt.model <- function(x, params) {
  if ("a" %in% names(params)) {
    return(voigt.superpos(x, params) + baseline.linear(x, params))
  } else {
    return(voigt.superpos(x, params))
  }
}

#' A synthethic Raman signal
#'
#' Returns list of `x` and `y` and the true parameters.
#'
#' @param seed sets the random seed
#' @param K number of peaks
#' @param noise std of noise
#' @export
synthetic.signal <- function(seed, K, noise) {
  set.seed(seed)
  params <- list(
    amp = stats::rlnorm(K, 1, 1),
    lwidth = stats::rlnorm(K, 1, 1),
    gwidth = stats::rlnorm(K, 1, 1),
    pos = stats::runif(K, 100, 900),
    a = stats::rnorm(1, 0, 0.000001),
    b = stats::rnorm(1, 20)
  )
  x <- seq(0, 1000, 1.2)
  y <- stats::rnorm(length(x), voigt.model(x, params), noise)

  return(list(x = x, y = y, params = params))
}
