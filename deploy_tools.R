library(devtools)
library(roxygen2)
main_dir <- 'C:/Users/Jan/Desktop/ISW_SpectraBayes/ISW_SpectraBayes/spectralem'
setwd(file.path(main_dir, "spectralem"))

#https://www.marinedatascience.co/blog/2020/01/09/checklist-for-r-package-re-submissions-on-cran/


# document
document()
devtools::build_manual()

# stylr
styler::style_pkg()

# website
install.packages("pkgdown")
usethis::use_pkgdown()
pkgdown::build_site()
#build_site_github_pages(
#  pkg = ".",
#  ...,
#  dest_dir = "docs",
#  clean = TRUE,
#  install = FALSE,
#  new_process = FALSE
#)

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
devtools::test_coverage()
devtools::run_examples()

# rhub
library(rhub)
#usethis::use_cran_comments()
rhub::check_on_linux()
#cran_prep <- rhub::check_for_cran()
rhub::cran_summary()


