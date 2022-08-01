library(testthat)

check_if_fit_params_are_correct <- function(x, res = Null) {
  fit <- numeric(length(x))
  fp <- res$fit_params
  for (i in seq_along(fp$pos)) {
    fit <- fit + fp$amp[i] * RcppFaddeeva::Voigt(x, fp$pos[i], fp$gwidth[i], fp$lwidth[i])
  }

  fit <- fit + fp$a * x + fp$b

  j <- 1
  for (comp in res$components) {
    if ("ComponentFixedDensity" %in% class(comp)) {
      fit <- fit + fp$background_amps[j] * comp$f
      j <- j + 1
    }
  }
  # plot(x, fit)
  # lines(x, y)
  # lines(x, res$fit)
  expect_equal(sum((res$fit - fit)^2), 0, tol = 1e-10)
}

test_that("test em_algorithm basecase", {
  # library(spectralem)
  seed <- 213
  kp <- 4
  noise <- 0.001
  data <- data.synthetic.stormyclouds(seed, kp, noise)
  x <- data$x
  y <- data$y

  res <- spectralem(
    x,
    y,
    max_peaks = 4,
    max_iter = 15,
    add_component_every_iters = 2,
  )

  # plot_fit(x, y, res)
  # ok fit;
  expect_equal(sort(data$params$pos), sort(res$fit_params$pos), tol = 1e-3)
  expect_equal(sum((res$fit - y)^2), 0, tol = 1e-3)
  check_if_fit_params_are_correct(x, res)
})

test_that("test em_algorithm basecase running_start", {
  # library(spectralem)
  seed <- 213
  kp <- 4
  noise <- 0.001
  data <- data.synthetic.stormyclouds(seed, kp, noise)
  x <- data$x
  y <- data$y

  res <- spectralem(
    x,
    y,
    max_peaks = 4,
    max_iter = 6,
    add_component_every_iters = 2,
    start_peaks = list(
      pos = c(332, 577, 697, 813),
      gwidth = c(5, 5, 5, 5),
      lwidth = c(0.1, 0.1, 0.1, 0.1)
    )
  )

  # plot_fit(x, y, res)
  # ok fit;
  expect_equal(sort(data$params$pos), sort(res$fit_params$pos), tol = 1e-3)
  expect_equal(sum((res$fit - y)^2), 0, tol = 1e-3)
  check_if_fit_params_are_correct(x, res)
})


test_that("test em_algorithm trivial cover center x range", {
  # library(spectralem)
  x <- seq(0, 1, 0.001)-0.5
  y <- -3 * x + 5

  res <- expect_warning(spectralem(
    x,
    y,
    max_peaks = 0,
    max_iter = 2,
    add_component_every_iters = 1,
  ))
  # plot_fit(x, y, res)
  # ok fit;
  expect_equal(sum((res$fit - y)^2), 0, tol = 1e-3)
  check_if_fit_params_are_correct(x, res)
})


test_that("test em_algorithm trivial cover negative x range", {
  # library(spectralem)
  x <- seq(0, 1, 0.001)-2
  y <- -0.1 * x + 5

  res <- expect_warning(spectralem(
    x,
    y,
    max_peaks = 0,
    max_iter = 2,
    add_component_every_iters = 1,
  ))
  #plot_fit(x, y, res)
  # ok fit;
  expect_equal(sum((res$fit - y)^2), 0, tol = 1e-3)
  check_if_fit_params_are_correct(x, res)
})

test_that("test em_algorithm slope", {
  # library(spectralem)
  seed <- 4232
  kp <- 2
  noise <- 0.003
  data <- data.synthetic.stormyclouds(seed, kp, noise)
  x <- data$x
  y <- data$y + x * 0.001
  res <- expect_warning(
    spectralem(x, y, max_peaks = 2, max_iter = 6, add_component_every_iters = 2),
    "Strong background"
  )
  #plot_fit(x, y, res)
  # failure fit;
  check_if_fit_params_are_correct(x, res)
})

test_that("test em_algorithm cutoff peaks", {
  seed <- 12311258
  kp <- 9
  noise <- 0.00003
  data <- data.synthetic.stormyclouds(seed, kp, noise)
  x <- data$x[seq(100, 400)]
  y <- data$y[seq(100, 400)]

  res <- spectralem(x, y, max_peaks = 4, max_iter = 10, add_component_every_iters = 2)

  # plot_fit(x, y, res)
  expect_equal(sort(data$params$pos)[seq(1, 4)], sort(res$fit_params$pos), tol = 1e-2)
  check_if_fit_params_are_correct(x, res)
})
