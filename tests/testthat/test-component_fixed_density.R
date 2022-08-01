library(testthat)
test_that("test component_fixed_density", {
  expect_error(ComponentFixedDensity$new(c(-1, 1, 1, 1)))
})
