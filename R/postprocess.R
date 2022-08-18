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


