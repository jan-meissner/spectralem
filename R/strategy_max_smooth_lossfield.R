#' Max Smoothed Lossfield Placement Strategy
#'
#' @name StrategyMaxSmoothLossfield
#' @title Max Smoothed Lossfield Placement Strategy
#' @rdname StrategyMaxSmoothLossfield
#'
#' @description
#' Places a new Voigt peak where the lossfield is maximal.
#'
#' @keywords internal
StrategyMaxSmoothLossfield <- R6::R6Class(
  "StrategyMaxSmoothLossfield",
  inherit = AbstractStrategy,
  public = list(
    #' @field typical_width typical width of peaks, used as the width of the blur kernel and for the gwidth of the newly placed peak
    typical_width = NULL,

    #' @field min_width minimal width of a placed peak
    min_width = NULL,

    #' @field burn_in_iters number of iterations to perform with fixed densities
    burn_in_iters = NULL,

    #' @field possible_peak_positions range of possible peak positions
    possible_peak_positions = NULL,

    #' @description
    #' Initialize the strategy.
    #'
    #' @param typical_width typical width of peaks, used as the width of the blur kernel and for the gwidth of the newly placed peak
    #' @param min_width minimal gwidth and lwidth for the newly placed peak
    #' @param burn_in_iters number of iterations to perform with fixed densities.
    #' @param possible_peak_positions range of possible peak positions
    initialize = function(min_width, typical_width, burn_in_iters, possible_peak_positions) {
      self$typical_width <- typical_width
      self$burn_in_iters <- burn_in_iters
      self$min_width <- min_width
      self$possible_peak_positions <- possible_peak_positions
    },

    #' @description
    #' Calculate and return a new component (usually a voigt component).
    #'
    #' @param x a vector of x coordinates.
    #' @param p_y density being fit.
    #' @param components a list of AbstractComponents.
    #' @return new components list with the new peak added.
    place_new_peak = function(x, p_y, components) {
      # special case if there are no other components
      if (length(components) == 0) {
        new_comp <- list(private$get_voigt_at_pos(x[which.max(p_y)], 1))
        comps <- private$optimize_new_component_parameters(x, p_y, list(), new_comp)
        return(comps)
      }

      # compute smoothed lossfield
      smooth_lossfield_int_cor <- private$get_smoothed_lossfield(x, p_y, components)

      # new pos where lossfield is maximal
      new_pos <- x[which.max(smooth_lossfield_int_cor)]
      new_comp <- private$get_voigt_at_pos(new_pos, 1 / length(components))

      # optimize placed peak parameters, improves algorithm slightly, not strictly needed.
      private$optimize_new_component_parameters(x, p_y, components, new_comp)
    }
  ),
  private = list(
    get_smoothed_lossfield = function(x, p_y, components) {
      # integrate(x, f) = dintegrate %*% f
      dintegrate <- sapply(
        seq_along(x),
        function(i) integrate(x, replace(numeric(length(x)), i, 1))
      )

      # calculate integrand of relative entropy (kuberlack divergence)
      d <- mixture_density(x, components)
      closs <- p_y * log(d)
      copt <- p_y * log(p_y)
      copt <- ifelse(is.na(copt), 0, copt)
      relative_entropy <- (copt - closs)
      relative_entropy_int <- relative_entropy * dintegrate # correction to account for non-uniform spaced x

      # smooth integrand
      winsize <- base::ceiling(length(x) / 4)
      gaussian <- exp(-1 / 2 * (seq(-0, winsize) * mean(diff(x)) / self$typical_width)^2)
      k <- kernel(gaussian / (2 * sum(gaussian) - gaussian[1]))
      smooth_lossfield_int_cor <- kernapply(relative_entropy_int, k, circular = TRUE)

      # debugplot(x, p_y, d, smooth_lossfield_int_cor, relative_entropy_int, x[which.max(smooth_lossfield_int_cor)])

      smooth_lossfield_int_cor
    },

    get_voigt_at_pos = function(pos, pi) {
      ComponentTrunctVoigt$new(
        min_width = self$min_width,
        possible_peak_positions = self$possible_peak_positions
      )$set_params(
        pos = pos,
        gwidth = self$typical_width,
        lwidth = self$typical_width / 10 # min_width?
      )$set_pi(pi)
    },

    optimize_new_component_parameters = function(x, p_y, components, new_comp) {
      # fix all components, but the new one, and run the EM-algorthm for a few iterations in that setup

      fixed_comps <- lapply(components, function(comp) {
        comp$to_fixed_component(x)
      })

      # em algorithm
      fixed_comps <- c(fixed_comps, new_comp)
      for (i in 1:self$burn_in_iters) {
        resp <- calculate_responsibilities(x, fixed_comps)
        maximize_Qk(x, p_y, resp, fixed_comps)
        calc_and_assign_new_pi(x, resp * p_y, fixed_comps)
      }

      # set pi's from the ran em algorithm
      new_comps <- c(components, new_comp)
      mapply(function(fixed_comp, comp) {
        comp$set_pi(fixed_comp$pi)
      }, fixed_comps, new_comps)
      new_comps
    }
  )
)
