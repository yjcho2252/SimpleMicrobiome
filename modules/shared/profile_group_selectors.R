profile_group_selectors <- function(ns, primary_id, primary_label, secondary_id, secondary_label) {
  shiny::tags$fieldset(
    id = ns("profile_group_selectors"), class = "profile-group-selectors", disabled = NA,
    htmltools::tagAppendAttributes(
      shiny::textOutput(ns("profile_group_selector_state")),
      class = "profile-group-selector-state", hidden = NA
    ),
    shiny::selectInput(ns(primary_id), primary_label, choices = NULL),
    shiny::selectInput(ns(secondary_id), secondary_label, choices = NULL)
  )
}

register_profile_group_selectors <- function(input, output, session, columns, primary_id, secondary_id,
                                             primary_extra = character(), secondary_none,
                                             exclude_sample_id = FALSE) {
  available_columns <- shiny::reactive({
    tryCatch(columns(), shiny.silent.error = function(e) character())
  })
  output$profile_group_selector_state <- shiny::renderText({
    if (length(available_columns())) "enabled" else "disabled"
  })
  shiny::outputOptions(output, "profile_group_selector_state", suspendWhenHidden = FALSE)

  initialized <- FALSE
  primary_selection <- function(available) {
    defaults <- setdiff(available, "SampleID")
    choices <- c(available, primary_extra)
    current <- input[[primary_id]]
    if (length(current) == 1L && current %in% choices &&
        (initialized || nzchar(current))) return(current)
    if (length(defaults)) return(defaults[1L])
    head(choices, 1L)
  }

  shiny::observeEvent(available_columns(), {
    available <- available_columns()
    choices <- if (length(available)) c(available, primary_extra) else character()
    selected <- if (length(available)) primary_selection(available) else character()
    shiny::updateSelectInput(session, primary_id, choices = choices, selected = selected)
    initialized <<- length(available) > 0L
  }, ignoreNULL = FALSE)

  shiny::observeEvent(list(available_columns(), input[[primary_id]]), {
    available <- available_columns()
    if (!length(available)) {
      shiny::updateSelectInput(session, secondary_id, choices = character(), selected = character())
      return()
    }
    primary <- primary_selection(available)
    excluded <- c(primary, if (exclude_sample_id) "SampleID")
    choices <- c(secondary_none, if (nzchar(primary)) setdiff(available, excluded))
    current <- input[[secondary_id]]
    selected <- if (length(current) == 1L && current %in% choices) current else unname(secondary_none)
    shiny::updateSelectInput(session, secondary_id, choices = choices, selected = selected)
  })
}
