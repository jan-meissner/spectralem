#' @name AbstractComponent
#' @title Abstract Interface for a Component
#' @rdname AbstractComponent
#'
#' @description
#' AbstractComponent is the abstract base class for components.
#' All subclasses must implement fit and density.
#'
#' @export
AbstractComponent <- R6::R6Class(
  "AbstractComponent",
  public = list(
    #' @field pi the class probability of the component in the mixture model
    pi = 1,

    #' @description
    #' Calculate the density multiplied with pi at the given x-coordinates.
    #'
    #' @param x a vector of points at which the density is evaluated.
    pi_density = function(x) {
      self$pi * self$density(x)
    },

    #' @description
    #' Abstract. Calculate the density at the given x-coordinates.
    #'
    #' @param x a vector of points at which the density is evaluated.
    density = function(x) {
      stop("Not implemented!")
    },

    #' @description
    #' Abstract. Fits the component to the signal rf by maximizing Q.
    #' If performed approximately turns the algorithm into a ECM algorithm.
    #'
    #' @param x coordinates of points on the x-axis.
    #' @param rf coordinates of function values.
    fit = function(x, rf) {
      stop("Not implemented!")
    },

    #' @description
    #' Setter for pi.
    #
    #' @param pi the new value.
    set_pi = function(pi) {
      self$pi <- pi
      self
    },

    #' @description
    #' Calculates Q between rf and the density for the component as defined in !!!paper-doi!!!.
    #'
    #' @param x vector of x values.
    #' @param rf vector of the product between p_y and the responsiblity vector of the given component.
    Q = function(x, rf) {
      private$Q_calls <- private$Q_calls + 1
      q <- rf * log(self$density(x))
      if (!all(is.finite(q))) {
        stop("Integrand is not finite!")
      }
      integrate(x, q)
    },

    #' @description
    #' Returns a component with a fixed density identical to the density of the current component.
    #'
    #' @param x vector of x values.
    to_fixed_component = function(x) {
      ComponentFixedDensity$new(self$density(x))$set_pi(self$pi)
    }
  ),
  list(
    Q_calls = 0
  )
)
