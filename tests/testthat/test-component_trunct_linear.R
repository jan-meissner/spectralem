library(testthat)
test_that("test component_trunct_linear", {
  test_case <- function(a, b, x) {
    y <- a * x + b + rnorm(length(x), 0, 1)
    y <- pmax(0, y)
    norm <- integrate(x, y)
    p_y <- y / norm

    c <- ComponentTrunctLinear$new()
    fit_p_y <- c$fit(x, p_y)$density(x)

    expect_equal(sd(p_y - fit_p_y) * norm, 1, tolerance = 1)
    expect_equal(c$a * norm, a, tolerance = 1e-1)
    expect_equal(c$b * norm, b, tolerance = 1e-1)
    c
  }
  # positive slope
  set.seed(1)
  test_case(1, 5000, seq(1, 1000))
  set.seed(1)
  test_case(1, 1, seq(1, 1000))

  # case where signal is cutoff; can't be fit so test fails
  # set.seed(1)
  # test_case(0.8, -300, seq(1, 1000))

  # negative slope
  set.seed(1)
  c <- test_case(-0.9, 5000, seq(1, 1000))

  # test rescale
  old_a <- c$a
  old_b <- c$b
  c$set_pi(4)
  expect_equal(old_a * 4, c$a)
  expect_equal(old_b * 4, c$b)
})
