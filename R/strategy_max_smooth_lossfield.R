#' Max Smoothed Lossfield Placement Strategy
#'
#' @name StrategyMaxSmoothLossfield
#' @title Max Smoothed Lossfield Placement Strategy
#' @rdname StrategyMaxSmoothLossfield
#'
#' @description
#' Places a new Voigt peak where the lossfield is maximal.
#'
#' @export
StrategyMaxSmoothLossfield <- R6::R6Class(
  "StrategyMaxSmoothLossfield",
  inherit = AbstractStrategy,
  public = list(
    #' @field init_gauss_std initial gwidth of the .
    init_gauss_std = NULL,

    #' @field blur_sig radius of gaussian smoothing kernel applied to the loss.
    blur_sig = NULL,

    #' @field burn_in_iters number of iterations to perform with fixed densities.
    burn_in_iters = NULL,

    #' @description
    #' Initialize the strategy.
    #'
    #' @param init_gauss_std initial gwidth of the .
    #' @param blur_sig radius of gaussian smoothing kernel applied to the loss.
    #' @param burn_in_iters number of iterations to perform with fixed densities.
    initialize = function(init_gauss_std = 5, blur_sig = 5, burn_in_iters = 20) {
      self$init_gauss_std <- init_gauss_std
      self$blur_sig <- blur_sig
      self$burn_in_iters <- burn_in_iters
    },

    #' @description
    #' Calculate and return a new component (usually a voigt component).
    #'
    #' @param x a vector of x coordinates.
    #' @param p_y density being fit.
    #' @param components a list of AbstractComponents.
    #' @return new components list with the new peak added.
    place_new_peak = function(x, p_y, components) {
      # integrate(x, f) = dintegrate %*% f
      dintegrate <- sapply(
        seq_along(x),
        function(i) integrate(x, replace(numeric(length(x)), i, 1))
      )

      # calculate loss field
      d <- mixture_density(x, components)
      closs <- p_y * log(d)
      copt <- p_y * log(p_y)
      copt <- ifelse(is.na(copt), 0, copt)
      relative_entropy <- (copt - closs)
      relative_entropy_int <- relative_entropy * dintegrate

      # smooth loss field
      winsize <- base::ceiling(length(x) / 4)
      gaussian <- exp(-1 / 2 * (seq(-0, winsize) / self$blur_sig)^2)
      k <- kernel(gaussian / (2 * sum(gaussian) - gaussian[1]))
      smooth_lossfield_int_cor <- kernapply(relative_entropy_int, k, circular = TRUE)
      new_pos <- x[which.max(smooth_lossfield_int_cor)]

      # debugplot(x, p_y, d, smooth_lossfield_int_cor, lossfield_int_cor, new_pos)

      # new peak
      new_comp <- ComponentTrunctVoigt$new()$set_params(
        pos = new_pos,
        gwidth = self$init_gauss_std,
        lwidth = self$init_gauss_std / 10
      )$set_pi(1 / length(components))

      # optimize placed peak parameters
      fixed_comps <- lapply(components, function(comp) {
        comp$to_fixed_component(x)
      })
      fixed_comps <- c(fixed_comps, new_comp)
      for (i in 1:self$burn_in_iters) {
        resp <- calculate_responsibilities(x, fixed_comps)
        maximize_Qk(x, p_y, resp, fixed_comps)
        calc_and_assign_new_pi(x, resp * p_y, fixed_comps)
      }

      new_comps <- c(components, new_comp)
      mapply(function(fixed_comp, comp) {
        comp$set_pi(fixed_comp$pi)
      }, fixed_comps, new_comps)
      new_comps
    }
  )
)
