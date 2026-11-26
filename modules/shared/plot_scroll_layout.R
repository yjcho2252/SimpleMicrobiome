plot_scroll_layout <- function(settings, results) {
  htmltools::tagAppendAttributes(
    shiny::sidebarLayout(
      htmltools::tagAppendAttributes(settings, tabindex = "0", role = "region", "aria-label" = "Plot settings"),
      htmltools::tagAppendAttributes(results, tabindex = "0", "aria-label" = "Analysis results")
    ),
    class = "plot-scroll-layout"
  )
}
