#' @name AbstractStrategy
#' @title Abstract Interface for a Peak Placement Strategy.
#' @rdname AbstractStrategy
#'
#' @description
#' AbstractStrategy is the abstract base class peak placement for strategies.
#' All subclasses must implement place_new_peak.
#'
#' @keywords internal
AbstractStrategy <- R6::R6Class(
  "AbstractStrategy",
  public = list(
    #' @description
    #' Abstract. Calculate and return a new component.
    #'
    #' @param x a vector of x coordinates.
    #' @param p_y density being fit.
    #' @param components list of AbstractComponent.
    #' @return new components list with the new peak added.
    place_new_peak = function(x, p_y, components) {
      stop("Not implemented!")
    }
  )
)
