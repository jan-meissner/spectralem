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

total_loglikelihood <- function(x, y, components) {
  d <- mixture_density(x, components)
  integrate(x, y * log(d))
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

  out <- list(components = d$components, convergence_diagnostic = d$convergence_diagnostic)
  out <- c(out, get_fit_and_fit_params(d))
  out
}

preprocess <- function(x, y, max_peaks, max_iter, placement_strategy,
                       start_peaks, background_model, add_component_every_iters) {
  out <- list()
  out$background_model <- background_model
  out$x <- x

  # convert y to p_y
  p_y <- y
  if (!any(p_y < 0)) {
    out$y_min <- min(p_y)
    p_y <- p_y - out$y_min
  }

  if (!isTRUE(all.equal(integrate(x, p_y), 1))) {
    out$y_norm <- integrate(x, p_y)
    p_y <- p_y / out$y_norm
  }
  out$p_y <- p_y

  # check if max_iter is consitent with max_peaks
  min_iter <- add_component_every_iters * (max_peaks -
    length(start_peaks$pos) -
    length(background_model) + 2)
  if (max_iter < min_iter) {
    stop("max_iter is too small! Increase max_iter or reduce max_peaks.")
  }

  # init components
  components <- background_model
  for (comp in background_model) {
    if (!("AbstractComponent" %in% class(comp))) {
      stop("All elements in background_model must inherit from AbstractComponent!")
    }
  }

  for (i in seq_along(start_peaks$pos)) {
    new_c <- ComponentTrunctVoigt$new()$set_params(
      start_peaks$pos[i],
      start_peaks$gwidth[i],
      start_peaks$lwidth[i]
    )
    components <- c(components, new_c)
  }
  out$components <- components

  # check strategy
  if (!("AbstractStrategy" %in% class(placement_strategy))) {
    stop("placement_strategy must inherit from AbstractStrategy!")
  }

  # convergence diagnostic
  out$convergence_diagnostic <- rep(NA, max_iter)

  out
}


get_fit_and_fit_params <- function(data) {
  x <- data$x
  components <- data$components
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
      params$a <- comp$a * data$y_norm
      params$b <- comp$b * data$y_norm + data$y_min
    }

    if ("ComponentTrunctVoigt" %in% class(comp)) {
      params$pos <- c(params$pos, comp$pos)
      params$gwidth <- c(params$gwidth, comp$gwidth)
      params$lwidth <- c(params$lwidth, comp$lwidth)
      params$amp <- c(params$amp, comp$get_amp() * data$y_norm)
    }

    if ("ComponentFixedDensity" %in% class(comp)) {
      params$background_amps <- c(params$background_amps, comp$get_amp() * data$y_norm)
    }
  }

  fit <- fit * data$y_norm + data$y_min
  list(fit = fit, fit_params = params)
}


#' Quickly fits voigt profiles to a signal.
#'
#' Requires no starting fit.
#' Most useful to get a fast and decent fit.
#' The algorithm was developed to produce a fast and decent fit that general non linear optimizers such as optim fail to produce.
#'
#' @param x the x coordinates of the Raman signal
#' @param y the function values of the function to be fit evaluated at the x coordinates
#' @param max_peaks the maximal number of components to fit
#' @param max_iter the maximal number of iterations
#' @param add_component_every_iters How many iteration to perform until a new peak is added. For all peaks to spawn max_iter >= (2+max_peaks) * add_component_every_iters needs to be satisfied.
#' @param start_peaks initial peak parameters allows user to pass know peaks. list of vectors
#' @param background_model defines the background model, support BSplines and a linear background. Instance of class which implements AbstractComponent.
#' @param placement_strategy (Optional) Defines strategy used to place/birth new peak. Instance of class which implements AbstractStrategy.
#' @param print_progress (Optional) Defines strategy used to place/birth new peak. Instance of class which implements AbstractStrategy.
#'
#' @export
spectralem <- function(x, y, max_peaks, max_iter,
                       print_progress = TRUE, add_component_every_iters = 10,
                       start_peaks = list(pos = c(), gwidth = c(), lwidth = c()),
                       background_model = list(ComponentTrunctLinear$new()),
                       placement_strategy = StrategyMaxSmoothLossfield$new()) {
  d <- preprocess(
    x, y, max_peaks, max_iter,
    placement_strategy = placement_strategy, start_peaks = start_peaks,
    background_model = background_model, add_component_every_iters = add_component_every_iters
  )
  p_y <- d$p_y

  pb <- utils::txtProgressBar(0, max_iter, style = 3)
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
      if (length(d$components) < max_peaks + length(background_model)) {
        d$components <- placement_strategy$place_new_peak(x, p_y, d$components)
      }
    }
  }

  postprocess(d)
}
