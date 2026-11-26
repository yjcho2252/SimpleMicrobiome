req_plot_legend_tab <- function(session, tab) {
  active_tab <- session$rootScope()$userData$active_plot_tab
  if (is.function(active_tab)) shiny::req(identical(active_tab(), tab), cancelOutput = TRUE)
}

plot_legend_output <- function(ns, id, result_ids, class = NULL) {
  shiny::uiOutput(
    ns(id),
    class = paste(c("app-plot-legend", class), collapse = " "),
    `data-legend-results` = paste(ns(result_ids), collapse = " ")
  )
}

plot_legend_box <- function(width, title, ...) {
  shiny::tags$div(
    class = "app-result-card",
    style = sprintf("width: %dpx;", as.integer(width)),
    shiny::tags$div(style = "font-weight: 600; margin-bottom: 4px;", title),
    shiny::tags$div(...)
  )
}
