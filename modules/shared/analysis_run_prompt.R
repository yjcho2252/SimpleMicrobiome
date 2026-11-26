analysis_run_view <- function(ns, run_id, ...) {
  started <- sprintf("input.%s > 0", run_id)
  shiny::tagList(
    shiny::conditionalPanel(
      paste0("!(", started, ")"), ns = ns,
      shiny::tags$div(style = "max-width: 900px;", shiny::plotOutput(ns("run_prompt"), height = 700))
    ),
    shiny::conditionalPanel(started, ..., ns = ns)
  )
}

analysis_run_prompt <- function(button_label) {
  force(button_label)
  shiny::renderPlot({
    graphics::plot.new()
    graphics::text(0.5, 0.5, sprintf("Click '%s' to start analysis.", button_label), cex = 0.85)
  }, res = 110)
}

automatic_analysis_view <- function(ns, ..., prompt_overlay = NULL) {
  prompt <- shiny::plotOutput(ns("run_prompt"), height = 700)
  if (!is.null(prompt_overlay)) prompt <- prompt_overlay(prompt)
  shiny::tagList(
    shiny::conditionalPanel(
      "output.result_visibility !== 'ready'", ns = ns,
      shiny::tags$div(style = "max-width: 900px;", prompt)
    ),
    shiny::conditionalPanel("output.result_visibility === 'ready'", ..., ns = ns)
  )
}

req_active_analysis <- function(active) {
  # Successful computations do not depend on tab visibility.
  if (!isTRUE(shiny::isolate(active()))) shiny::req(isTRUE(active()))
}

automatic_analysis <- function(active, inputs, compute) {
  opened <- shiny::reactiveVal(FALSE)
  shiny::observeEvent(active(), {
    if (isTRUE(active())) opened(TRUE)
  })
  settled <- shiny::debounce(shiny::reactive({
    if (opened()) tryCatch(inputs(), shiny.silent.error = function(e) NULL)
  }), 350)
  last_request <- NULL
  last_result <- NULL
  has_result <- FALSE

  shiny::reactive({
    shiny::req(opened())
    requested <- inputs()
    if (has_result && identical(requested, last_request)) return(last_result)
    req_active_analysis(active)
    shiny::validate(shiny::need(identical(requested, settled()), "Updating analysis..."))
    result <- shiny::isolate(compute(requested))
    last_request <<- requested
    last_result <<- result
    has_result <<- TRUE
    result
  })
}

register_automatic_analysis_status <- function(output, active, result, initial_message,
                                               updating_message = "Updating analysis...") {
  state <- shiny::reactive({
    shiny::req(isTRUE(active()), cancelOutput = TRUE)
    tryCatch({
      result()
      list(ready = TRUE)
    }, shiny.silent.error = function(e) {
      detail <- conditionMessage(e)
      if (identical(detail, "Updating analysis...")) detail <- updating_message
      list(ready = FALSE, detail = if (nzchar(detail)) detail else initial_message)
    }, error = function(e) {
      message("Automatic analysis failed: ", conditionMessage(e))
      list(ready = FALSE, detail = "The analysis could not be generated. Check the analysis settings.")
    })
  })
  output$result_visibility <- shiny::renderText(if (state()$ready) "ready" else "waiting")
  shiny::outputOptions(output, "result_visibility", suspendWhenHidden = FALSE)
  output$run_prompt <- shiny::renderPlot({
    current <- state()
    shiny::req(!current$ready)
    graphics::plot.new()
    graphics::text(0.5, 0.5, paste(strwrap(current$detail, width = 75), collapse = "\n"), cex = 0.85)
  }, res = 110)
}
