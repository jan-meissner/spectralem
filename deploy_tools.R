library(roxygen2)
library(devtools)

# document
devtools::document()
devtools::build_manual()

# stylr
library(styler)
styler::style_pkg(filetype = c(".R", ".Rmd", ".Rmarkdown", ".Rnw"))

# update R
install.packages("installr")
library(installr)
updateR()

#RcppFaddeeva
install.packages("RcppFaddeeva")
library(remotes)
remotes::install_version("RcppFaddeeva", "0.2.2")

#coverage
cov <- devtools::test_coverage()
covr::report(
  x = cov,
  browse = interactive(),
  file = "latest_coverage.html"
)

devtools::run_examples()

# building
devtools::build()
#devtools::release()



