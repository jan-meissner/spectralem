library(testthat)

check_if_fit_params_are_correct <- function(x, res = Null, background = list()) {
  fit <- voigt.model(x, res$fit_params)

  fp <- res$fit_params
  for (j in seq_along(background)) {
    fit <- fit + fp$background_amps[j] * background[[j]]
  }
  expect_equal(sum((res$fit - fit)^2), 0, tol = 1e-10)
}


test_that("test em_algorithm basecase", {
  data <- synthetic.signal(213, 4, 0.001)

  res <- spectralem(data$x, data$y, K = 4)

  # plot_fit(x, y, res)
  # ok fit;
  expect_equal(sort(data$params$pos), sort(res$fit_params$pos), tol = 1e-3)
  expect_equal(sum((res$fit - data$y)^2), 0, tol = 1e-3)
  check_if_fit_params_are_correct(data$x, res)
})


test_that("test em_algorithm basecase_small_x_scale_and_negative", {
  data <- synthetic.signal(213, 4, 0.001)
  scale <- 1e-8 / diff(range(data$x))
  x <- data$x * scale - 3
  y <- data$y

  res <- spectralem(x, y, K = 4)

  # plot_fit(x, y, res)
  # ok fit;
  expect_equal(sort(scale * data$params$pos - 3), sort(res$fit_params$pos), tol = 1e-3)
  expect_equal(sum((res$fit - y)^2), 0, tol = 1e-3)
  check_if_fit_params_are_correct(x, res)
})

test_that("test em_algorithm basecase_big_x_scale_and_center_x_negative_y", {
  data <- synthetic.signal(213, 4, 0.001)
  scale <- 2e10 / diff(range(data$x))
  offset <- 500
  x <- (data$x - offset) * scale
  y <- data$y - 1e5

  res <- spectralem(x, y, K = 4)

  # plot_fit(x, y, res)
  # ok fit;
  expect_equal(sort(scale * (data$params$pos - offset) - 3), sort(res$fit_params$pos), tol = 1e-3)
  expect_equal(sum((res$fit - y)^2), 0, tol = 1e-3)
  check_if_fit_params_are_correct(x, res)
})


test_that("test em_algorithm basecase_running_start", {
  data <- synthetic.signal(213, 4, 0.001)

  res <- spectralem(
    data$x, data$y,
    K = 4, max_iter = 10,
    start_peaks = list(
      pos = c(332, 577, 697, 813),
      gwidth = c(5, 5, 5, 5),
      lwidth = c(0.1, 0.1, 0.1, 0.1)
    )
  )

  # ok fit;
  expect_equal(sort(data$params$pos), sort(res$fit_params$pos), tol = 1e-3)
  expect_equal(sum((res$fit - data$y)^2), 0, tol = 1e-3)
  check_if_fit_params_are_correct(data$x, res)
})

test_that("test em_algorithm no_background", {
  data <- synthetic.signal(213, 4, 0.0000000001)
  y <- data$y + 0.00002 * data$x

  res <- spectralem(data$x, y, K = 4, background_model = list(linear = FALSE))

  # plot_fit(data$x, y, res)
  # ok fit;
  expect_equal(sum((res$fit - y)^2), 0, tol = 3e-1)
})

test_that("test em_algorithm trivial_cover_center_x_range", {
  x <- seq(0, 1, 0.001) - 0.5
  y <- -3 * x + 5

  res <- expect_warning(spectralem(x, y,
    K = 0,
    max_iter = 2,
    add_component_every_iters = 1,
  ))
  # ok fit;
  expect_equal(sum((res$fit - y)^2), 0, tol = 1e-3)
  check_if_fit_params_are_correct(x, res)
})


test_that("test em_algorithm trivial_cover_negative_x_range", {
  x <- seq(0, 1, 0.001) - 2
  y <- -0.1 * x + 5

  res <- expect_warning(spectralem(
    x,
    y,
    K = 0,
    max_iter = 2,
    add_component_every_iters = 1,
  ))
  # ok fit;
  expect_equal(sum((res$fit - y)^2), 0, tol = 1e-3)
  check_if_fit_params_are_correct(x, res)
})

test_that("test em_algorithm slope", {
  data <- synthetic.signal(4232, 2, 0.003)
  x <- data$x
  y <- data$y + x * 0.001
  res <- expect_warning(
    spectralem(x, y, K = 2, max_iter = 6, add_component_every_iters = 2),
    "Strong background"
  )
  # failure fit;
  check_if_fit_params_are_correct(x, res)
})

test_that("test em_algorithm cutoff_peaks", {
  data <- synthetic.signal(213, 4, 0.001)
  x <- data$x
  y <- data$y

  res <- spectralem(x, y, K = 4)

  # ok fit;
  expect_equal(sort(data$params$pos), sort(res$fit_params$pos), tol = 1e-3)
  expect_equal(sum((res$fit - y)^2), 0, tol = 1e-3)
  check_if_fit_params_are_correct(x, res)
})


test_that("test em_algorithm spline_background", {
  data <- synthetic.signal(213, 4, 0.001)
  x <- data$x
  y <- data$y

  xnorm <- x / max(x)
  f1 <- xnorm + 0.01 * xnorm^2 + 0.01 * xnorm^3 - 0.4 * xnorm^4 - 0.2 * xnorm^5
  f1 <- f1 - min(f1)

  f2 <- 0.01 * xnorm^2
  f2 <- f2 - min(f2)

  f3 <- 5 * xnorm + -8 * xnorm^2 + 0.5 * xnorm^3 + 2.9 * xnorm^4 - 0.2 * xnorm^5
  f3 <- f3 - min(f3)

  y <- y + 0.3 * (0.1 * f1 + 0.4 * f2 + 0.3 * f3)

  res <- spectralem(
    x,
    y,
    K = 4,
    background_model = list(
      linear = TRUE, f1, f2, f3
    )
  )

  # okay-ish fit; strong back ground is problematic for small peaks
  expect_equal(sort(data$params$pos), sort(res$fit_params$pos), tol = 1e-2)
  expect_equal(sum((res$fit - y)^2), 0, tol = 1e-2)
  check_if_fit_params_are_correct(x, res, list(f1, f2, f3))
})

test_that("test em_algorithm dense_x", {
  data <- synthetic.signal(555, 4, 0.01)

  x <- seq(min(data$x), max(data$x), 0.23) # super dense x
  y <- voigt.model(x, data$params)

  res <- spectralem(x, y, K = 4)

  expect_equal(sort(data$params$pos), sort(res$fit_params$pos), tol = 1e-3)
  expect_equal(sum((res$fit - y)^2), 0, tol = 1e-2)
  check_if_fit_params_are_correct(x, res)
})

test_that("test em_algorithm high_K", {
  data <- synthetic.signal(213, 30, 0.005)

  res <- spectralem(data$x, data$y, K = 30)

  expect_equal(sort(data$params$pos), sort(res$fit_params$pos), tol = 1e-1)
  expect_equal(sum((res$fit - data$y)^2), 0, tol = 1e-1)
  check_if_fit_params_are_correct(data$x, res)
})
