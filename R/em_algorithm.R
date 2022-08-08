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

postprocess <- function(d) {
  pi_background <- 0
  for (comp in d$components) {
    if ("ComponentTrunctLinear" %in% class(comp)) {
      pi_background <- pi_background + comp$pi
    }

    if ("ComponentFixedDensity" %in% class(comp)) {
      pi_background <- pi_background + comp$pi
    }
  }
  if (pi_background > 0.95) {
    warning(
      "Strong background detected. Verify the fit with caution.
      If possible remove any background before fitting.
      Fitting a line and subtracting is often sufficient."
    )
  }

  if (!all((d$convergence_diagnostic - cummax(d$convergence_diagnostic)) <= mean(diff(d$convergence_diagnostic)) * 1e-2)) {
    warning(
      "Likelihood not monotonically increasing. Verify the fit with caution."
    )
  }
  # check if convergence is monotone?
  out <- list(components = d$components, convergence_diagnostic = d$convergence_diagnostic)
  out <- c(out, get_fit_and_fit_params(d))
  out
}

preprocess <- function(x, y, max_peaks,
                       start_peaks,
                       background_model,
                       typical_width,
                       add_component_every_iters,
                       max_iter,
                       min_width,
                       possible_peak_positions) {
  out <- list()

  # x to [0, 1] scale for better numerics
  out$true_x <- x

  out$scale <- diff(range(x))
  out$offset <- min(x)
  out$x <- (x - out$offset) / out$scale
  x <- out$x

  # rescale width parameters
  typical_width <- typical_width / out$scale
  min_width <- min_width / out$scale
  possible_peak_positions <- (possible_peak_positions - out$offset) / out$scale
  for (i in seq_along(start_peaks$pos)) {
    start_peaks$pos[i] <- (start_peaks$pos[i] - out$offset) / out$scale
    start_peaks$gwidth[i] <- start_peaks$gwidth[i] / out$scale
    start_peaks$lwidth[i] <- start_peaks$lwidth[i] / out$scale
  }

  # convert y to p_y
  p_y <- y

  out$y_min <- min(p_y)
  p_y <- p_y - out$y_min

  out$y_norm <- integrate(x, p_y)
  p_y <- p_y / out$y_norm

  out$p_y <- p_y

  # check if max_iter is consitent with max_peaks
  out$max_iter <- max_iter
  min_iter <- add_component_every_iters * (max_peaks -
    length(start_peaks$pos) -
    length(background_model) + 2)
  if (max_iter < min_iter) {
    stop("max_iter is too small! Increase max_iter or reduce max_peaks.")
  }

  # init background
  components <- list()
  if (background_model$linear) {
    components <- c(components, ComponentTrunctLinear$new())
  }

  for (ele in background_model) {
    if (mode(ele) %in% c("numeric")) {
      components <- c(components, ComponentFixedDensity$new(ele))
    }
  }
  out$final_num_components <- max_peaks + length(components)

  # init strategy
  out$placement_strategy <- StrategyMaxSmoothLossfield$new(
    min_width = min_width, typical_width = typical_width,
    burn_in_iters = min(20, max_iter / (1 + max_peaks)),
    possible_peak_positions = possible_peak_positions
  )


  # init components
  for (i in seq_along(start_peaks$pos)) {
    new_c <- ComponentTrunctVoigt$new(
      min_width = min_width,
      possible_peak_positions = possible_peak_positions
    )$set_params(
      start_peaks$pos[i],
      start_peaks$gwidth[i],
      start_peaks$lwidth[i]
    )
    components <- c(components, new_c)
  }
  out$components <- components

  if (length(out$components) == 0) {
    out$components <- out$placement_strategy$place_new_peak(x, p_y, out$components)
  }
  for (comp in out$components) {
    comp$set_pi(1 / length(out$components))
  }


  # init convergence diagnostic
  out$convergence_diagnostic <- rep(NA, max_iter)

  # warn too many x points
  if (length(x) > 9000) {
    warning("Algorithm will be slow due to large size of x samples. Consider sampling a subset of x points.")
  }

  out
}


get_fit_and_fit_params <- function(d) {
  x <- d$x
  components <- d$components
  params <- list()
  params$amp <- c()
  params$pos <- c()
  params$lwidth <- c()
  params$gwidth <- c()
  params$a <- NULL
  params$b <- NULL
  params$background_amps <- c()

  fit <- numeric(length(x))

  for (comp in components) {
    fit <- fit + comp$pi_density(x)

    if ("ComponentTrunctLinear" %in% class(comp)) {
      params$a <- comp$a * d$y_norm / d$scale
      params$b <- comp$b * d$y_norm + d$y_min - params$a * d$offset
    }

    if ("ComponentTrunctVoigt" %in% class(comp)) {
      params$pos <- c(params$pos, comp$pos * d$scale + d$offset)
      params$gwidth <- c(params$gwidth, comp$gwidth * d$scale)
      params$lwidth <- c(params$lwidth, comp$lwidth * d$scale)
      params$amp <- c(params$amp, comp$get_amp() * d$scale * d$y_norm)
    }

    if ("ComponentFixedDensity" %in% class(comp)) {
      params$background_amps <- c(params$background_amps, comp$get_amp() * d$y_norm)
    }
  }
  fit <- fit * d$y_norm + d$y_min


  list(fit = fit, fit_params = params)
}


#' Robust, fast and fully automated peak fitting of Voigt profiles to spectral data
#'
#' @description
#' Implements the algorithm as discussed in <paperedoi>. The model fitted is defined as:
#' \deqn{f(x) = a x + b + \sum_{j=1}^K \beta_j V(x \mid \theta_j) + \sum_{i=1} \alpha_i f_i(x)}
#' Here \eqn{\beta_j > 0} are the amplitudes of the Voigt profiles \eqn{V} where \eqn{\theta_j} represents the position (\code{pos}),
#' Gaussian (\code{gwidth}) and Lorentzian (\code{lwidth}) width of the Voigt profile. It requires no starting parameters, but optionally some can be passed in \code{start_peaks}.
#'
#' The most impactful hyperparameters are \code{typical_width} and \code{add_component_every_iters}.
#' The former should be a bit smaller than the typical width of Voigt profiles in the signal. If chosen too big or too small the fit will be usually suboptimal.
#' The hyperparameter \code{add_component_every_iters} defines the number of iterations to perform until a new peak is added.
#' Setting it to a bigger value comes at a greater computational cost but usually results in a better fit.
#' The hyperparameter \code{max_iter} defines the total number of iterations, impact is similar but usually weaker than \code{add_component_every_iters}.
#'
#' Additionally, it is possible to fit arbitrary positive background functions \eqn{f_i(x) > 0}, these are by default zero.
#' Note that \eqn{\alpha_i > 0} must hold. BSplines could be used here.
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

    if (((i %% add_component_every_iters) == (add_component_every_iters - 1))) {
      if (length(d$components) < d$final_num_components) {
        d$components <- d$placement_strategy$place_new_peak(x, p_y, d$components)
      }
    }
  }

  postprocess(d)
}
