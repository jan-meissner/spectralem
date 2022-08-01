

#####################
get_toluene_measured_data <- function(){
  path <- "C:\\Users\\Jan\\Desktop\\ISW_SpectraBayes\\ISW_SpectraBayes\\bayes_mcmc\\dev\\wip_gibbs\\golden_gibbs\\toluene_measured.csv"
  toluene_measured <- read.csv(path, header=FALSE, sep=";")
  x <- as.vector(toluene_measured[,1])
  ytrue <- as.vector(toluene_measured[,2])
  return(list(x=x,y=ytrue))
}

library(remotes)
remotes::install_version("RcppFaddeeva", "0.2.2")
library(RcppFaddeeva)
remotes::install_github("jan-meissner/spectralem")

library("spectralem")
data <- data.synthetic.stormyclouds(seed = 1, kp = 30, noise = 0.001)
x <- data$x
y <- data$y

res <- spectralem(x, y, max_peaks = 3, max_iter = 40)


install.packages('plotly')
library("plotly")
fig <- plotly::plot_ly(x = x, y = y, name = 'true', type = 'scatter', mode = 'lines')
fig <- fig %>% plotly::add_lines(y = res$fit, name = 'fit', type = 'scatter', mode = 'lines')
fig

#library(proftools)
#pd <- proftools::readProfileData("profile.out")
#head(funSummary(pd), 10)


library(babynames) # provide the dataset: a dataframe called babynames
library(dplyr)

# Keep only 3 names
don <- babynames %>%
filter(name %in% c("Ashley", "Patricia", "Helen")) %>%
filter(sex=="F")

library(ggplot2)
pd <- data.frame(y = c(res$fit, y), x=rep(x, 2), type = rep(c("fit","true"), each = length(x)))
pd %>% ggplot( aes(x=x, y=y, color=type)) + geom_line() + theme_bw()