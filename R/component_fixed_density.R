#' Component which models a fixed background function.
#'
#' @name ComponentFixedDensity
#' @title Component which models a fixed background function.
#' @rdname ComponentFixedDensity
#'
#' @description
#' Requires fixed background function.
#'
#' @export
ComponentFixedDensity <- R6::R6Class(
  "ComponentFixedDensity",
  inherit = AbstractComponent,
  public = list(
    #' @field f vector that defines density up to normalization.
    f = NULL,

    #' @description
    #' Initializes a component with constant density proportional to f. Assumes f positive.
    #'
    #' @param f a positive vector.
    initialize = function(f) {
      if (any(f < 0)) {
        stop("f must be positive")
      }
      self$f <- f
    },

    #' @description
    #' Calculate the density at the given x-coordinates.
    #'
    #' @param x a vector of points at which the density is evaluated.
    density = function(x) {
      private$norm <- integrate(x, self$f)
      return(self$f / private$norm)
    },

    #' @description
    #' Fits the component to the signal rf.
    #'
    #' @param x coordinates of points on the x-axis.
    #' @param rf coordinates of function values.
    fit = function(x, rf) {
      self
    },

    #' @description
    #' Gets amplitude of the voigt profile.
    #'
    get_amp = function() {
      self$pi / private$norm
    }
  ),
  list(
    norm = NULL
  )
)
