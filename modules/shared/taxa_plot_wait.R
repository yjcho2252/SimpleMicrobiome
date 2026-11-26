plot_data_wait_view <- function(ns, ...) {
  shiny::tagList(
    shiny::conditionalPanel(
      "output.plot_data_ready !== 'ready'", ns = ns,
      shiny::div(class = "taxa-plot-stage", shiny::uiOutput(ns("plot_data_wait")))
    ),
    shiny::conditionalPanel("output.plot_data_ready === 'ready'", ..., ns = ns)
  )
}

register_plot_data_wait <- function(output, ps_obj, plot_label) {
  output$plot_data_ready <- shiny::renderText({
    ps <- ps_obj()
    if (!is.null(ps) && phyloseq::nsamples(ps) >= 1L) "ready" else "waiting"
  })
  shiny::outputOptions(output, "plot_data_ready", suspendWhenHidden = FALSE)
  output$plot_data_wait <- shiny::renderUI({
    validate_taxa_plot_data(ps_obj(), plot_label)
    NULL
  })
}

taxa_plot_loading <- function(plot_label) {
  shiny::div(
    class = "taxa-plot-loading", role = "status",
    shiny::span(sprintf("Preparing the %s. Please wait...", plot_label))
  )
}

validate_taxa_plot_data <- function(ps, plot_label) {
  shiny::validate(
    shiny::need(!is.null(ps), sprintf("Load and preprocess data to view the %s.", plot_label)),
    shiny::need(is.null(ps) || phyloseq::nsamples(ps) >= 1L, "At least one sample is required."),
    errorClass = "taxa-plot-wait"
  )
}
