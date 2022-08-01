library(testthat)
library(spectralem)

if (Sys.getenv("NOT_CRAN") == "true") {
  test_check("spectralem")
}
