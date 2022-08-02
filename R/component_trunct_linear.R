#' @name ComponentTrunctLinear
#' @title Component which models a linear background
#' @rdname ComponentTrunctLinear
#'
#' @description
#' Component which models a linear background a x + b.
#'
#' @keywords internal
ComponentTrunctLinear <- R6::R6Class(
  "ComponentTrunctLinear",
  inherit = AbstractComponent,
  public = list(
    #' @field a slope
    a = NULL,

    #' @field b offset
    b = NULL,

    #' @description
    #' Initializes the linear component.
    #'
    #' @param nu_max internal parameter, should be left as is.
    #' @param pi initial value for pi.
    initialize = function(pi = 1) {
      private$angle <- 0
      self$pi <- pi
    },

    #' @description
    #' Calculate the density at the given x-coordinates.
    #'
    #' @param x a vector of points at which the density is evaluated.
    density = function(x) {
      fx <- pmax(0, sin(private$angle) * x + cos(private$angle))
      ftrunct <- fx / integrate(x, fx)
      return(ftrunct)
    },

    #' @description
    #' Fits the component to the signal rf.
    #'
    #' @param x coordinates of points on the x-axis.
    #' @param rf coordinates of function values.
    fit = function(x, rf) {
      r <- range(x)
      # restrict range of angle such that sin(angle)*xrange+cos(angle) >= 0
      angle_range <- c( # kind of complicated
        if (max(r) >= 0) atan2(-1 / max(r), 1) else atan2(1 / max(r), -1),
        if (min(r) >= 0) atan2(1 / min(r), -1) else atan2(-1 / min(r), 1)
      )

      FUN <- function(angle) {
        private$angle <- angle
        res <- self$Q(x, rf)
        res
      }

      private$angle <- optimize(FUN, angle_range, maximum = TRUE, tol = private$tol)$maximum

      fx <- pmax(0, sin(private$angle) * x + cos(private$angle))
      norm <- integrate(x, fx)
      self$a <- sin(private$angle) / norm * self$pi
      self$b <- cos(private$angle) / norm * self$pi
      self
    },

    #' @description
    #' Setter for pi.
    #
    #' @param pi the new value.
    set_pi = function(pi) {
      private$rescale_a_b(pi)
      self$pi <- pi
    }
  ),
  private = list(
    tol = 1e-16, # tolerance of fit
    angle = NULL,
    rescale_a_b = function(new_pi) {
      self$a <- self$a * new_pi / self$pi
      self$b <- self$b * new_pi / self$pi
    }
  )
)
