library(testthat)
test_that("test component_abstract", {
  comp <- AbstractComponent$new()
  expect_error(comp$fit(1, 1))
  expect_error(comp$density(1))
})
