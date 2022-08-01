debugplot <- function(x, p_y, d, smooth_lossfield_int_cor, lossfield_int_cor, new_pos) {
  vline <- function(x = 0, color = "green") {
    list(
      type = "line",
      y0 = 0,
      y1 = 1,
      yref = "paper",
      x0 = x,
      x1 = x,
      line = list(color = color, dash = "dot")
    )
  }
  fig <- plotly::plot_ly(x = x)
  fig <- fig %>% plotly::add_lines(y = p_y, name = "y", type = "scatter", mode = "lines")
  # fig <- fig %>% plotly::add_lines(y = brrs, name = 'residual resp field', type = 'scatter', mode = 'lines', opacity = 0.5)
  fig <- fig %>% plotly::add_lines(y = d, name = "fit", type = "scatter", mode = "lines")
  fig <- fig %>% plotly::add_lines(y = smooth_lossfield_int_cor, name = "smoothloss", type = "smoothloss", mode = "lines", opacity = 0.5, line = list(color = "green"))
  fig <- fig %>% plotly::add_lines(y = lossfield_int_cor, name = "lossfield", type = "scatter", mode = "lines", opacity = 0.1, line = list(color = "green"))
  fig <- fig %>% plotly::layout(
    xaxis = list(title = "", zeroline = FALSE),
    yaxis = list(title = "", exponentformat = "e", zeroline = FALSE),
    shapes = list(vline(new_pos, color = "red"))
  )
  fig %>% show()
}

plotdebug <- function(x, y, res = NULL) {
  fig <- plotly::plot_ly(x = x)
  fig <- fig %>% plotly::add_lines(y = y, name = "y", type = "scatter", mode = "lines")
  fig <- fig %>% plotly::layout(
    xaxis = list(title = "", zeroline = FALSE),
    yaxis = list(title = "", exponentformat = "e", zeroline = FALSE)
  )
  fig %>% show()
}

plot_fit <- function(x, y, res = NULL) {
  fig <- plotly::plot_ly(x = x)
  fig <- fig %>% plotly::add_lines(y = y, name = "y", type = "scatter", mode = "lines")
  if (!is.null(res)) fig <- fig %>% plotly::add_lines(y = res$fit, name = "fit", type = "scatter", mode = "lines")
  fig <- fig %>% plotly::layout(
    xaxis = list(title = "", zeroline = FALSE),
    yaxis = list(title = "", exponentformat = "e", zeroline = FALSE)
  )
  fig %>% show()
}
