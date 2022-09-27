calculate_responsibilities <- function(x, components) {
  resp <- sapply(seq_along(components), function(j) components[[j]]$pi_density(x))
  resp <- resp / rowSums(resp)
  resp
}

calc_and_assign_new_pi <- function(x, ry, components) {
  pi <- apply(ry, MARGIN = 2, function(col) integrate(x, col))
  if (!isTRUE(all.equal(1, sum(pi)))) {
    warning("Pi vector does not sum to 1. Proceed with caution.")
  }
  pi <- pi / sum(pi)
  mapply(function(c, j) {
    c$set_pi(pi[j])
  }, components, seq_along(components))
}

# Q is the partial loglikelihood; see "component_abstract" > function "Q"
maximize_Qk <- function(x, y, resp, components) {
  mapply(function(c, j) {
    c$fit(x, resp[, j] * y)
  }, components, seq_along(components))
}

mixture_density <- function(x, components) {
  d <- 0
  for (j in seq_along(components)) {
    d <- d + components[[j]]$pi_density(x)
  }
  d
}

total_loglikelihood <- function(x, p_y, components) {
  d <- mixture_density(x, components)
  integrate(x, p_y * log(d))
}

#' Robust, fast and fully automated peak fitting of Voigt profiles to spectral data
#'
#' @description
#' Implements the algorithm as discussed in <paperedoi>. The model fitted is defined as:
#' \deqn{f(x) = a x + b + \sum_{j=1}^K \beta_j V(x \mid \theta_j) + \sum_{i=1} \alpha_i f_i(x)}
#' Here \eqn{\beta_j > 0} are the amplitudes of the Voigt profiles \eqn{V} where \eqn{\theta_j} represents the position (\code{pos}),
#' Gaussian (\code{gwidth}) and Lorentzian (\code{lwidth}) width of the Voigt profile.
#' Additionally, it is possible to specify arbitrary positive background functions \eqn{f_i(x) > 0}, these are by default zero.
#' Note that \eqn{\alpha_i > 0} holds.
#'
#' The algorithm requires no starting parameters, but optionally some can be passed in \code{start_peaks}.
#' The most impactful hyperparameters are \code{typical_width} and \code{add_component_every_iters}.
#' The former should be a bit smaller than the typical width of Voigt profiles in the signal. If chosen too big or too small the fit will be usually suboptimal.
#' The hyperparameter \code{add_component_every_iters} defines the number of iterations to perform until a new peak is added.
#' Setting it to a bigger value comes at a greater computational cost but usually results in a better fit.
#' The hyperparameter \code{max_iter} defines the total number of iterations, impact is similar but usually weaker than \code{add_component_every_iters}.
#'
#' Implementation uses Voigt profiles as defined by 'RcppFaddeeva::Voigt()'.
#'
#' @param x the \code{x} coordinates of the signal
#' @param y the function values of the function to be fit evaluated at the \code{x} coordinates
#' @param K the number of peaks to fit
#' @param start_peaks optional initial peak parameters; list see examples
#' @param add_component_every_iters number of iterations to perform until a new peak is added
#' @param typical_width typical width of a voigt profile in the signal
#' @param background_model defines the background model
#' @param max_iter the maximal number of iterations
#' @param min_width minimal gaussian and lorentzian width
#' @param print_progress whether to print a progress bar
#' @param possible_peak_positions allowed range of positions for fitted peaks, can be increased to fit cut off peaks
#'
#' @return A list with components:
#' \itemize{
#'   \item fit_params - A list which contains named components \code{amp}, \code{pos}, \code{gwidth} and \code{lwidth} which are vectors of the same length.
#'                      For each the i-th entry yields the parameters of the i-th Voigt profile.
#'                      Further named components of the list \code{a} and \code{b} describe the slope and offset of the linear background respectively.
#'                      If \eqn{f_i(x)} are specified also contains a vector \code{background_amps} representing \eqn{\alpha_i}
#'   \item fit - A vector containing the function values of the fitted function evaluated at \code{x}.
#' }
#' @examples
#' \dontrun{
#' res <- spectralem(
#'   x, y,
#'   K = 4,
#' )
#'
#' res <- spectralem(
#'   x, y,
#'   K = 4,
#'   start_peaks = list(
#'     pos = c(332, 577, 697),
#'     gwidth = c(5, 5, 5),
#'     lwidth = c(0.1, 0.1, 0.1)
#'   )
#' )
#'
#' res <- spectralem(
#'   x, y,
#'   K = 4,
#'   background_model = list(linear = TRUE, f1, f2, f3)
#' )
#' }
#' @export
spectralem <- function(x, y, K,
                       start_peaks = list(pos = c(), gwidth = c(), lwidth = c()),
                       background_model = list(linear = TRUE),
                       typical_width = diff(range(x)) / 100 / 2,
                       add_component_every_iters = 10,
                       max_iter = add_component_every_iters * (K + 2 + ceiling(K / 2)),
                       min_width = 1e-6 * mean(diff(x)),
                       possible_peak_positions = range(x),
                       print_progress = TRUE) {

  # preprocess input to create internal data object
  # contains d$x, d$p_y, d$components which are central to the algorithm
  d <- preprocess(x, y, K,
    start_peaks = start_peaks,
    background_model = background_model,
    typical_width = typical_width,
    add_component_every_iters = add_component_every_iters,
    max_iter = max_iter,
    min_width = min_width,
    possible_peak_positions = possible_peak_positions
  )
  p_y <- d$p_y
  x <- d$x

  if (print_progress) pb <- utils::txtProgressBar(0, max_iter, style = 3)
  for (i in 1:max_iter) {
    if (print_progress) utils::setTxtProgressBar(pb, i)
    d$convergence_diagnostic[i] <- total_loglikelihood(x, p_y, d$components)

    # E-Step
    resp <- calculate_responsibilities(x, d$components)

    # M-Step, first step
    maximize_Qk(x, p_y, resp, d$components)

    # M-Step, second step
    calc_and_assign_new_pi(x, resp * p_y, d$components)

    # Place a new peak, if necessary
    if (((i %% add_component_every_iters) == (add_component_every_iters - 1))) {
      if (length(d$components) < d$final_num_components) {
        d$components <- d$placement_strategy$place_new_peak(x, p_y, d$components)
      }
    }
  }

  # create output object
  postprocess(d)
}
