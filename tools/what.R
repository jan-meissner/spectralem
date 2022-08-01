

#####################
get_toluene_measured_data <- function(){
  path <- "C:\\Users\\Jan\\Desktop\\ISW_SpectraBayes\\ISW_SpectraBayes\\bayes_mcmc\\dev\\wip_gibbs\\golden_gibbs\\toluene_measured.csv"
  toluene_measured <- read.csv(path, header=FALSE, sep=";")
  x <- as.vector(toluene_measured[,1])
  ytrue <- as.vector(toluene_measured[,2])
  return(list(x=x,y=ytrue))
}


library("plotly")
library("spectralem")
data <- get_toluene_measured_data()
x <- data$x
y <- data$y

res <- spectralem(
  x,
  y,
  max_peaks = 20,
  max_iter = 300,
  placement_strategy = StrategyMaxSmoothLossfield$new()
)


#fig <- plotly::plot_ly(x = x)
#fig <- fig %>% plotly::add_lines(y = y, name = 'true', type = 'scatter', mode = 'lines')
#fig <- fig %>% plotly::add_lines(y = res$fit, name = 'fit', type = 'scatter', mode = 'lines')
#fig <- fig %>% plotly::layout(xaxis = list(title = "", zeroline = FALSE),
#                              yaxis = list(title = "", exponentformat = 'e', zeroline = FALSE))
#fig %>% show()

#library(proftools)
#pd <- proftools::readProfileData("profile.out")
#head(funSummary(pd), 10)
