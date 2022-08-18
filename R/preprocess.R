calculate_py <- function(out){
  # convert y to p_y
  p_y <- out$y

  out$y_min <- min(p_y)
  p_y <- p_y - out$y_min

  out$y_norm <- integrate(out$x, p_y)
  p_y <- p_y / out$y_norm

  out$p_y <- p_y
  out
}

unit_scale_x <- function(out){
  # x to [0, 1] scale for better numerics
  out$true_x <- out$x
  out$scale <- diff(range(out$x))
  out$offset <- min(out$x)
  out$x <- (out$x - out$offset) / out$scale
  out
}

unit_scale_hyperparams <- function(out){
  # rescale width parameters
  out$typical_width <- out$typical_width / out$scale
  out$min_width <- out$min_width / out$scale
  out$possible_peak_positions <- (out$possible_peak_positions - out$offset) / out$scale

  for (i in seq_along(out$start_peaks$pos)) {
    out$start_peaks$pos[i] <- (out$start_peaks$pos[i] - out$offset) / out$scale
    out$start_peaks$gwidth[i] <- out$start_peaks$gwidth[i] / out$scale
    out$start_peaks$lwidth[i] <- out$start_peaks$lwidth[i] / out$scale
  }

  out
}

check_max_iters <- function(out){
  min_iter <- out$add_component_every_iters * (out$max_peaks -
    length(out$start_peaks$pos) -
    length(out$background_model) + 2)
  if (out$max_iter < min_iter) {
    stop("max_iter is too small! Increase max_iter or reduce max_peaks.")
  }
}

init_background <- function(out){
  components <- list()
  if (out$background_model$linear) {
    components <- c(components, ComponentTrunctLinear$new())
  }

  for (ele in out$background_model) {
    if (mode(ele) %in% c("numeric")) {
      components <- c(components, ComponentFixedDensity$new(ele))
    }
  }
  out$final_num_components <- out$max_peaks + length(components)
  out$components <- components
  out
}

init_strategy <- function(out){
  out$placement_strategy <- StrategyMaxSmoothLossfield$new(
    min_width = out$min_width, typical_width = out$typical_width,
    burn_in_iters = min(20, out$max_iter / (1 + out$max_peaks)),
    possible_peak_positions = out$possible_peak_positions
  )
  out
}

init_components <- function(out){
  for (i in seq_along(out$start_peaks$pos)) {
    new_c <- ComponentTrunctVoigt$new(
      min_width = out$min_width,
      possible_peak_positions = out$possible_peak_positions
    )$set_params(
      out$start_peaks$pos[i],
      out$start_peaks$gwidth[i],
      out$start_peaks$lwidth[i]
    )
    out$components <- c(out$components, new_c)
  }

  if (length(out$components) == 0) {
    out$components <- out$placement_strategy$place_new_peak(out$x, out$p_y, out$components)
  }
  for (comp in out$components) {
    comp$set_pi(1 / length(out$components))
  }
  out
}

burn_in_pis <- function(out){
  # burn-in all pi's (in case start_peaks are specified this is needed.)
  fixed_comps <- lapply(out$components, function(comp) {
    comp$to_fixed_component(out$x)
  })

  for (i in 1:20) {
    resp <- calculate_responsibilities(out$x, fixed_comps)
    maximize_Qk(out$x, out$p_y, resp, fixed_comps)
    calc_and_assign_new_pi(out$x, resp * out$p_y, fixed_comps)
  }

  mapply(function(fixed_comp, comp) {
    comp$set_pi(fixed_comp$pi)
  }, fixed_comps, out$components)
  out
}

check_numer_of_x_points <- function(out){
  if (length(out$x) > 9000) {
    warning("Algorithm will be slow due to large size of x samples. Consider sampling a subset of x points.")
  }
}

load_default_controls <- function(out){
    ## Defaults :
    con <- list(
      max_iter = out$add_component_every_iters * (out$max_peaks + 2 + ceiling(out$max_peaks / 2)),
      possible_peak_positions = range(out$x),
      print_progress = TRUE,
      min_width = 1e-6 * mean(diff(out$x)),
      y_min = min(out$y)
    )
    con[(namc <- names(control))] <- out$control
}

preprocess <- function(x, y, max_peaks,
                       start_peaks,
                       background_model,
                       typical_width,
                       add_component_every_iters,
                       max_iter,
                       min_width,
                       possible_peak_positions) {
  out <- list(
    x = x,
    y = y,
    max_peaks = max_peaks,
    start_peaks = start_peaks,
    background_model = background_model,
    typical_width = typical_width,
    add_component_every_iters = add_component_every_iters,
    max_iter = max_iter,
    min_width = min_width,
    possible_peak_positions = possible_peak_positions
  )

  # init convergence diagnostic
  out$convergence_diagnostic <- rep(NA, out$max_iter)

  out <- unit_scale_x(out)
  out <- unit_scale_hyperparams(out)
  out <- calculate_py(out)
  out <- init_background(out)
  out <- init_strategy(out)
  out <- init_components(out)
  out <- burn_in_pis(out)

  check_max_iters(out)
  check_numer_of_x_points(out)

  out
}