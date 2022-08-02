# Based off https://github.com/mlverse/torch
# document
library(roxygen2)
library(devtools)

devtools::document()
devtools::build_manual()

# stylr
styler::style_pkg()

# update R
install.packages("installr")
library(installr)
updateR()

#RcppFaddeeva
install.packages("RcppFaddeeva")
library(remotes)
remotes::install_version("RcppFaddeeva", "0.2.2")

# building
devtools::build()
#https://builder.r-hub.io/

#coverage
cov <- devtools::test_coverage()
covr::report(
  x = cov,
  browse = interactive(),
  file = "latest_coverage.html"
)

devtools::run_examples()

use_github_action_check_standard(save_as = "R-CMD-check.yaml")


#devtools::release()



