library(testthat)
test_that("test component_trunct_voigt1", {
  set.seed(1)
  pos <- 0.331
  gwidth <- 0.01
  lwidth <- 0.0001
  x <- seq(0, 1, 0.001)
  y <- RcppFaddeeva::Voigt(x, pos, gwidth, lwidth) + rnorm(length(x), sd = 1e-6)
  p_y <- y - min(y)
  p_y <- p_y / integrate(x, p_y)

  c <- ComponentTrunctVoigt$new()
  fit_p_y <- c$fit(x, p_y)$density(x)

  # ok accuracy
  expect_equal(c$pos / pos, 1, tolerance = 5e-2)
  expect_equal(c$gwidth / gwidth, 1, tolerance = 5e-2)
  expect_equal(c$lwidth / lwidth, 1, tolerance = 5e-2)
})

test_that("test component_trunct_voigt cutoff", {
  # peak pos is outside range(x)
  set.seed(1)
  pos <- 1.1
  gwidth <- 0.1
  lwidth <- 0.01
  x <- seq(0, 1, 0.001)
  y <- RcppFaddeeva::Voigt(x, pos, gwidth, lwidth) + rnorm(length(x), sd = 1e-6)
  p_y <- y - min(y)
  p_y <- p_y / integrate(x, p_y)

  c <- ComponentTrunctVoigt$new()
  fit_p_y <- c$fit(x, p_y)$density(x)

  # library("plotly")
  # fig <- plotly::plot_ly(x = x)
  # fig <- fig %>% plotly::add_lines(y = p_y, name = 'True p_y', type = 'scatter', mode = 'lines')
  # fig <- fig %>% plotly::add_lines(y = fit_p_y, name = 'fit_p_y', type = 'scatter', mode = 'lines')
  # fig <- fig %>% plotly::layout(xaxis = list(title = "", zeroline = FALSE),
  #                              yaxis = list(title = "", exponentformat = 'e', zeroline = FALSE))
  # fig
  # ok accuracy
  expect_equal(c$pos / pos, 1, tolerance = 30e-2)
  expect_equal(c$gwidth / gwidth, 1, tolerance = 30e-2)
  expect_equal(c$lwidth / lwidth, 1, tolerance = 1)
})
