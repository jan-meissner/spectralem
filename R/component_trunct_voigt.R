#' @name ComponentTrunctVoigt
#' @title Component which models a truncated voigt profile
#' @rdname ComponentTrunctVoigt
#'
#' @description
#' Component which models a truncated voigt profile. Models a voigt peak in the signal.
#' Truncation is needed to deal with peaks which have significant mass outside the x range.
#'
#' @keywords internal
#' @importFrom MASS ginv
#' @importFrom R6 is.R6
ComponentTrunctVoigt <- R6::R6Class(
  "ComponentTrunctVoigt",
  inherit = AbstractComponent,
  public = list(
    #' @field pos position of the voigt profile
    pos = NULL,

    #' @field gwidth gaussian width of the voigt profile
    gwidth = NULL,

    #' @field lwidth gaussian width of the voigt profile
    lwidth = NULL,

    #' @description
    #' Initialize the voigt component.
    #'
    #' @param pi initial value for pi.
    #' @param min_width minimal value for lwidth and gwidth
    #' @param possible_peak_positions range of possible peak positions
    initialize = function(min_width, possible_peak_positions, pi = 1) {
      private$min_width <- min_width
      private$possible_peak_positions <- possible_peak_positions
      self$pi <- pi
    },

    #' @description
    #' Calculate the density at the given x-coordinates.
    #'
    #' @param x a vector of points at which the density is evaluated.
    density = function(x) {
      private$cache_x <- x
      f <- RcppFaddeeva::Voigt(x, self$pos, self$gwidth, self$lwidth)
      ftrunct <- f / integrate(x, f)
      return(ftrunct)
    },

    #' @description
    #' Fits the component to the signal rf by maximizing Q.
    #' Ensures monotonicity by comparing with last fit.
    #'
    #' @param x coordinates of points on the x-axis.
    #' @param rf coordinates of function values.
    fit = function(x, rf) {
      minusQ <- function(cur) -private$Q_with_gradient(x, rf, cur[1], cur[2], cur[3])$Q
      minusQ_grad <- function(cur) -private$Q_with_gradient(x, rf, cur[1], cur[2], cur[3])$grad

      # plotdebug(x, rf)
      # guess initial values for optim with $guess blindly
      blind_guess <- private$guess(x, rf)

      # last fit as guess
      last_guess <- c(self$pos, self$gwidth, self$lwidth)

      # take better of the two guesses
      minusQ_blind_guess <- minusQ(blind_guess)
      minusQ_last_guess <- if (is.null(last_guess)) Inf else minusQ(last_guess)
      if (minusQ_blind_guess < minusQ_last_guess) {
        guess <- blind_guess
        guess_minusQ <- minusQ_blind_guess
      } else {
        guess <- last_guess
        guess_minusQ <- minusQ_last_guess
      }

      pos_upper <- max(private$possible_peak_positions)
      pos_lower <- min(private$possible_peak_positions)
      guess <- c(max(min(guess[1], pos_upper), pos_lower), guess[2], guess[3])


      # optim
      res <- optim(
        guess,
        lower = c(pos_lower, private$min_width, private$min_width),
        upper = c(pos_upper, Inf, Inf),
        minusQ,
        minusQ_grad,
        method = "L-BFGS-B"
      )

      # ensure that L-BFGS-B is at least as good as guess:
      est <- res$par
      if (guess_minusQ < res$value) {
        warning("L-BFGS-B failed, using guess instead")
        est <- guess
      }

      # cat('pos', est[[1]], 'widths', est[[2]], est[[3]], '\n')
      self$set_params(est[[1]], est[[2]], est[[3]])
      self
    },

    #' @description
    #' Setter for pos, gwidth and lwidth.
    #'
    #' @param pos the new pos
    #' @param gwidth the new gwidth
    #' @param lwidth the new lwidth
    set_params = function(pos = NULL, gwidth = NULL, lwidth = NULL) {
      if (!is.null(pos)) {
        self$pos <- pos
      }
      if (!is.null(gwidth)) {
        self$gwidth <- gwidth
      }
      if (!is.null(lwidth)) {
        self$lwidth <- lwidth
      }
      self
    },

    #' @description
    #' Gets amplitude of the voigt profile.
    #'
    get_amp = function() {
      x <- private$cache_x
      if (is.null(x)) {
        stop("Amp can only be calculated after self$density has been called.")
      }
      f <- RcppFaddeeva::Voigt(x, self$pos, self$gwidth, self$lwidth)
      self$pi / integrate(x, f)
    }
  ),
  private = list(
    min_width = NULL,
    possible_peak_positions = NULL,
    cache_x = NULL,
    cache_last_Q_gradient = NULL,
    cache_last_Q = NULL,
    cache_last_Q_signature = c(0, 0, 0),
    cache_dintegrate = numeric(0),
    Q_calls = 0, # diagnostic
    set_dintegrate = function(x) {
      if (length(x) != length(private$cache_dintegrate)) {
        private$cache_dintegrate <- sapply(
          seq_along(x),
          function(i) integrate(x, replace(numeric(length(x)), i, 1))
        )
      }
    },

    #' Calculate Q and it's gradient.
    #'
    Q_with_gradient = function(x, rf, pos, gwidth, lwidth) {
      if (all(private$cache_last_Q_signature == c(pos, gwidth, lwidth))) {
        return(list(grad = private$cache_last_Q_gradient, Q = private$cache_last_Q))
      }
      # forward
      n <- (gwidth * sqrt(2)) # 10 dn <-  sqrt(2)*gwidth %*% dgwidth
      zr <- (x - pos) / n # 9 dzr <- (-1/n) %*% dpos - zr/n %*% dn
      zi <- lwidth / n # 8 dzi <- 1/n %*% dlwidth - zi/n %*% dn
      w <- RcppFaddeeva::Faddeeva_w(zr + 1i * zi) # 7
      Rw <- Re(w) # 7 # dRw <- (2i/sqrt(pi) - 2 * z * w) * (dzr + 1i dzi)
      # dRw <- (- 2) * Re(z * w) %*% dzr + 2 * (Im(z*w) - 1/sqrt(pi)) %*% dzi
      f <- Rw / n / sqrt(pi) ## 5 # df <- diag(1/(n*sqrt(pi))) %*% Rw - diag(f/(n)) %*% dn
      o <- integrate(x, f) # 4 do = dintegrate %*% df
      b <- f / o # 3 # db = diag(1/o) %*% df - diag(b/o) %*% do
      q <- rf * log(b) # 2 dq = diag(rf/b) %*% db
      r <- integrate(x, q) # 1 dr = dintegrate %*% dq

      # helpers
      zw <- (zr + 1i * zi) * w

      # cache dintegrate
      private$set_dintegrate(x)

      # backwards
      drdq <- private$cache_dintegrate # 1
      drdb <- drdq * rf / b # 2
      drdo <- sum(drdb * (-b / o)) # 3
      drdf <- drdb / o + drdo * private$cache_dintegrate # 3 + 4
      drdRw <- drdf / (n * sqrt(pi)) # 5
      drdzr <- drdRw * (-2) * Re(zw) # 7
      drdzi <- drdRw * 2 * (Im(zw) - 1 / sqrt(pi)) # 7
      drdlwidth <- sum(drdzi / n) # 8
      drdpos <- sum(drdzr / (-n)) # 9
      drdn <- sum(drdf * (-f / (n)) +
        drdzi * (-zi / n) +
        drdzr * (-zr / n)) # 5 + 8 + 9
      drdgwidth <- drdn * sqrt(2) #

      gradient <- c(drdpos, drdgwidth, drdlwidth)

      # cache last call
      private$cache_last_Q_gradient <- gradient
      private$cache_last_Q <- r
      private$cache_last_Q_signature <- c(pos, gwidth, lwidth)
      private$Q_calls <- private$Q_calls + 1

      return(list(grad = gradient, Q = r))
    },

    #' Guess pos, gwidth, lwidth by fitting a pseudo voigt. To be treated as a blackbox.
    #'
    guess = function(x, y) {
      windowsize <- 5
      max_guess_width <- diff(range(x)) / 4
      #############################

      l2loss <- function(x, y, est, csigma, cmean) {
        z <- (x - cmean) / csigma
        fz <- exp(-log(2) * z^2) * sqrt(log(2) / base::pi)
        gz <- 1 / base::pi / (1 + z^2)
        f_est <- ((1 - est[2]) * fz + est[2] * gz) * est[1] / csigma
        return(sum((y - f_est)^2))
      }

      partial_fit <- function(x, y, csigma, cmean) {
        z <- (x - cmean) / csigma
        fz <- exp(-log(2) * z^2) * sqrt(log(2) / base::pi)
        gz <- 1 / base::pi / (1 + z^2)
        X <- matrix(0, length(x), 2)
        X[, ncol(X) - 1] <- fz
        X[, ncol(X)] <- gz - fz
        D <- t(X) %*% X
        d <- c(csigma * y %*% X)
        est <- MASS::ginv(D) %*% d
        if ((est[2] >= 0) && (est[1] >= est[2])) {
          est <- est
        } else {
          estA <- c(sum(d) / sum(D), sum(d) / sum(D))
          estB <- c(d[1] / D[1, 1], 0)
          if (t(estA) %*% D %*% estA - 2 * t(d) %*% estA >= t(estB) %*% D %*% estB - 2 * t(d) %*% estB) {
            est <- estB
          } else {
            est <- estA
          }
        }
        if (est[1] <= 0) {
          return(c(0, 0))
        } else {
          return(c(est[1], est[2] / est[1]))
        }
      }

      # fit quadratic model to the max and thus get a estimate for mean of pseudo voigt
      maskcenter <- min(max(which.max(y), windowsize + 1), length(y) - windowsize)
      mask <- seq(maskcenter - windowsize, maskcenter + windowsize)

      scale <- 1 / diff(range(x[mask]))
      offset <- min(x[mask])
      coef <- lm(y[mask] ~ poly((x[mask] - offset) * scale, 2, raw = TRUE))$coef
      mean <- min(max((-coef[[2]] / 2 / coef[[3]]) * scale + offset, range(x)[1]), range(x)[2])

      # incase quadratic fit fails fit simply take max
      r <- range(x[mask])
      if (!is.finite(mean) || !(mean >= r[1] && mean <= r[2])) {
        mean <- x[which.max(y)]
      }

      # fit pseudo voigt
      optim_func <- function(arg) {
        l2loss(x, y, partial_fit(x, y, arg, mean), arg, mean)
      }
      sigma <- optimize(optim_func, c(0, max_guess_width))$minimum
      res <- partial_fit(x, y, sigma, mean)
      eta <- res[2]

      # converting pseudo voigt to voigt
      ufl <- uniroot(function(arg) {
        1.36603 * arg -
          0.47719 * arg^2 +
          0.11116 * arg^3 - eta
      }, c(0, 1))$root
      ufg <- uniroot(function(ufg) {
        ufg^5 +
          2.69269 * ufg^4 * ufl +
          2.42843 * ufg^3 * ufl^2 +
          4.47163 * ufg^2 * ufl^3 +
          0.07842 * ufg * ufl^4 +
          ufl^5 - 1
      }, c(0, 1))$root
      c(
        pos = mean,
        gwidth = max(ufg * sigma, private$min_width),
        lwidth = max(ufl * sigma, private$min_width)
      )
    }
  )
)
