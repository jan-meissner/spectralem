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

voigt.model <- function(x, params) {
  return(voigt.superpos(x, params) + baseline.linear(x, params))
}

#' A synthethic Dataset.
#'
#' Return list of x and y and the true parameters.
#'
#' @param seed sets the random seed.
#' @param kp number of peaks.
#' @param noise sigma of gaussian noise.
#' @export
data.synthetic.stormyclouds <- function(seed, kp, noise) {
  set.seed(seed)
  params <- list(
    amp = stats::rlnorm(kp, 1, 1),
    lwidth = stats::rlnorm(kp, 1, 1),
    gwidth = stats::rlnorm(kp, 1, 1),
    pos = stats::runif(kp, 100, 900),
    a = stats::rnorm(1, 0, 0.000001),
    b = stats::rnorm(1, 20)
  )
  x <- seq(0, 1000, 1.2)
  y <- stats::rnorm(length(x), voigt.model(x, params), noise)

  return(list(x = x, y = y, params = params))
}
