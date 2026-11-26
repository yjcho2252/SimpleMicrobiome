parse_plot_order <- function(text, available) {
  entries <- trimws(strsplit(if (is.null(text)) "" else text, "\n", fixed = TRUE)[[1L]])
  entries <- entries[nzchar(entries)]
  unknown <- setdiff(entries, available)
  duplicates <- unique(entries[duplicated(entries)])
  problems <- c(
    if (length(unknown)) paste("Unknown:", paste(unknown, collapse = ", ")),
    if (length(duplicates)) paste("Repeated:", paste(duplicates, collapse = ", "))
  )
  list(values = entries, problems = problems)
}

resolve_sample_plot_order <- function(metadata, specification, groups) {
  ids <- rownames(metadata)
  if (identical(specification$mode, "default") || is.null(specification)) return(NULL)
  keys <- lapply(specification$values, function(variable) {
    values <- metadata[[variable]]
    if (length(groups[[variable]]) || !is.numeric(values)) {
      plot_order_factor(values, variable, list(groups = groups))
    } else {
      values
    }
  })
  ids[do.call(order, c(keys, list(ids, na.last = TRUE)))]
}

## UI
plot_order_editor <- function(input_id, values) {
  div(
    class = "plot-order-editor",
    div(hidden = NA, textAreaInput(input_id, NULL, value = paste(values, collapse = "\n"))),
    helpText("Drag group levels or use the Up and Down keys on a focused item."),
    tags$ol(
      class = "plot-order-items",
      "aria-label" = "Group levels",
      lapply(values, function(value) tags$li(draggable = "true", tabindex = "0", value))
    )
  )
}

mod_preprocessing_ui <- function(id) {
  ns <- NS(id)
  tagList(
    singleton(tags$head(tags$script(src = "plot-order.js", defer = NA))),
    tags$style(HTML("
      .well h4 { font-size: 16px; }
      .well h5 { font-size: 13px; }
      .well .control-label { font-size: 12px; }
      .well .checkbox label { font-size: 12px; }
      .well .form-control { font-size: 12px; }
      .well .btn { font-size: 11px; }
      .well .irs-grid-text { display: none !important; }
      .plot-order-settings { padding: 8px 12px; margin-bottom: 8px; }
      .plot-order-settings > summary { cursor: pointer; font-weight: 600; }
      .plot-order-settings[open] > summary { margin-bottom: 12px; }
      .plot-order-items { max-height: 12rem; overflow-y: auto; padding-left: 2rem; margin: 6px 0; }
      .plot-order-items li { border: 1px solid var(--bs-border-color, #ced4da); border-radius: 4px;
        margin-bottom: 3px; padding: 4px 8px; cursor: grab; overflow-wrap: anywhere; }
      .plot-order-items li[data-drop='before'] { box-shadow: inset 0 2px currentColor; }
      .plot-order-items li[data-drop='after'] { box-shadow: inset 0 -2px currentColor; }
      .plot-order-items:empty { display: none; }
    ")),
    fluidRow(
      column(
        width = 2,
        wellPanel(
          h4(icon("filter"), "Sample Selection"),
          hr(),
          h5(icon("filter"), "Read Counts Filtering"),
          uiOutput(ns("read_count_filter_ui")),
          div(style = "background: #f8f9fa; padding: 10px; border-radius: 5px;",
              verbatimTextOutput(ns("sample_count_info"))
          ),
          hr(),
          h5(icon("layer-group"), "Group-Based Toggle"),
          selectInput(ns("toggle_group_var"), "Group Variable", choices = NULL),
          selectizeInput(
            ns("toggle_group_levels"),
            "Group Level",
            choices = NULL,
            multiple = FALSE,
            options = list(
              placeholder = "Select one level",
              plugins = list("remove_button")
            )
          ),
          actionButton(
            ns("toggle_group_selection"),
            "Toggle Selected Group Rows",
            icon = icon("exchange-alt"),
            class = "btn-secondary btn-sm",
            style = "font-size: 12px;",
            width = "100%"
          ),
          hr(),
          actionButton(ns("reset_selection"), "Select All", 
                       icon = icon("sync"), class = "btn-secondary btn-sm", style = "font-size: 12px;", width = "100%"),
          div(style = "height: 4px;"),
          actionButton(
            ns("unselect_all"),
            "Unselect All",
            icon = icon("ban"),
            class = "btn-secondary btn-sm",
            style = "font-size: 12px;",
            width = "100%"
          )
        )
      ),
      column(
        width = 9,
        tags$details(
          class = "well plot-order-settings",
          tags$summary("Plot Order"),
          fluidRow(
            column(
              width = 6,
              selectInput(ns("sample_order_mode"), "Sample Order",
                          c("Module default" = "default", "By metadata" = "metadata")),
              conditionalPanel(
                condition = "input.sample_order_mode === 'metadata'", ns = ns,
                selectInput(ns("sample_order_first"), "First sort column", choices = NULL),
                selectInput(ns("sample_order_second"), "Second sort column", choices = NULL, selectize = FALSE)
              ),
              verbatimTextOutput(ns("sample_order_preview")),
              actionButton(ns("apply_sample_order"), "Apply Sample Order")
            ),
            column(
              width = 6,
              selectInput(ns("group_order_var"), "Group Order", choices = NULL),
              uiOutput(ns("group_order_editor")),
              actionButton(ns("apply_group_order"), "Apply Group Order"),
              actionButton(ns("reset_group_order"), "Reset Group Order")
            )
          )
        ),
        wellPanel(
          h4(icon("list-check"), "Individual Sample Selection (Click to Toggle)"),
          p("After finishing selection, click another analysis tab to continue.", style = "margin: 0 0 8px 0; color: #4b5563; font-size: 12px;"),
          div(
            style = "width: 100%; overflow-x: auto;",
            DTOutput(ns("sample_table"))
          )
        )
      )
    )
  )
}

## Server
mod_preprocessing_server <- function(id, ps_obj_initial, active_tab) {
  moduleServer(id, function(input, output, session) {
    resolve_meta_colname <- function(requested, available) {
      if (is.null(requested) || is.null(available) || length(available) == 0) {
        return(requested)
      }
      if (requested %in% available) {
        return(requested)
      }
      matched <- available[make.names(available) == make.names(requested)]
      if (length(matched) > 0) {
        return(matched[1])
      }
      requested
    }
    
    ps_filtered <- reactiveVal(NULL)
    ps_selected_raw <- reactiveVal(NULL)
    last_applied_selection <- reactiveVal(NULL)
    default_min_reads <- 2000

    order_metadata <- reactive({
      req(ps_obj_initial())
      as(phyloseq::sample_data(ps_obj_initial()), "data.frame")
    })
    sample_order <- reactiveVal(NULL)
    group_orders <- reactiveVal(list())

    observeEvent(ps_obj_initial(), {
      sample_order(NULL)
      group_orders(list())
      updateSelectInput(session, "sample_order_mode", selected = "default")
      if (is.null(ps_obj_initial())) return()
      columns <- setdiff(names(order_metadata()), "SampleID")
      updateSelectInput(session, "sample_order_first", choices = columns)
      updateSelectInput(session, "sample_order_second", choices = c("None" = "", columns), selected = "")
      updateSelectInput(session, "group_order_var", choices = columns)
    }, ignoreNULL = FALSE)

    group_order_levels <- reactive({
      metadata <- order_metadata()
      req(input$group_order_var %in% names(metadata))
      levels(factor(metadata[[input$group_order_var]]))
    })

    group_order_input_id <- reactive({
      req(group_order_levels())
      paste0("group_order_text_", match(input$group_order_var, names(order_metadata())))
    })

    output$group_order_editor <- renderUI({
      available <- group_order_levels()
      plot_order_editor(
        session$ns(group_order_input_id()),
        plot_order_levels(available, group_orders()[[input$group_order_var]])
      )
    })

    group_order_draft <- reactive({
      parse_plot_order(input[[group_order_input_id()]], group_order_levels())
    })

    sample_order_draft <- reactive({
      metadata <- order_metadata()
      mode <- input$sample_order_mode
      req(mode)
      draft <- list(mode = mode, values = NULL, problems = NULL)
      if (identical(mode, "metadata")) {
        variables <- unique(c(input$sample_order_first, input$sample_order_second))
        draft$values <- variables[nzchar(variables)]
        if (!length(draft$values) || any(!draft$values %in% setdiff(names(metadata), "SampleID"))) {
          draft$problems <- "Select a metadata column."
        }
      }
      draft
    })

    plot_order <- reactive({
      groups <- group_orders()
      list(samples = resolve_sample_plot_order(order_metadata(), sample_order(), groups), groups = groups)
    })

    format_order_preview <- function(values) {
      shown <- paste(utils::head(values, 12L), collapse = " -> ")
      if (length(values) > 12L) {
        shown <- paste0(shown, " ... (", length(values), " total)")
      }
      shown
    }

    output$sample_order_preview <- renderText({
      draft <- sample_order_draft()
      if (length(draft$problems)) return(paste(draft$problems, collapse = "\n"))
      ids <- resolve_sample_plot_order(order_metadata(), draft, group_orders())
      if (is.null(ids)) return("Each plot uses its default sample sorting.")
      selected <- meta_df_for_dt()$SampleID[selected_indices_for_apply()]
      format_order_preview(intersect(ids, selected))
    })

    observeEvent(input$apply_sample_order, {
      draft <- sample_order_draft()
      if (length(draft$problems)) return(showNotification(paste(draft$problems, collapse = " "), type = "error"))
      sample_order(draft)
      showNotification("Sample order applied.")
    })

    observeEvent(input$apply_group_order, {
      draft <- group_order_draft()
      if (length(draft$problems)) return(showNotification(paste(draft$problems, collapse = " "), type = "error"))
      groups <- group_orders()
      groups[[input$group_order_var]] <- draft$values
      group_orders(groups)
      showNotification("Group order applied.")
    })

    observeEvent(input$reset_group_order, {
      req(group_order_levels())
      groups <- group_orders()
      groups[[input$group_order_var]] <- NULL
      group_orders(groups)
      updateTextAreaInput(session, group_order_input_id(), value = paste(group_order_levels(), collapse = "\n"))
      showNotification("Group order reset.")
    })
    
    meta_df_for_dt <- reactive({
      req(ps_obj_initial())
      ps <- ps_obj_initial()
      
      meta_df <- data.frame(phyloseq::sample_data(ps))
      
      if ("SampleID" %in% colnames(meta_df)) {
        meta_df <- meta_df %>% dplyr::select(-SampleID)
      }
      
      read_counts <- phyloseq::sample_sums(ps)
      
      df_for_dt <- meta_df %>%
        tibble::rownames_to_column("SampleID") %>% 
        dplyr::mutate(`Read Count` = as.numeric(read_counts[SampleID])) %>%
        dplyr::rename_with(~gsub("_", " ", .x)) %>%
        dplyr::select(`SampleID`, `Read Count`, everything())
      
      return(df_for_dt)
    })
    
    output$sample_table <- renderDT({
      df <- meta_df_for_dt()
      req(df)
      min_reads <- input$min_read_count
      if (is.null(min_reads)) min_reads <- default_min_reads
      initial_selected_rows <- which(df$`Read Count` >= min_reads)
      
      datatable(
        df,
        rownames = FALSE, 
        selection = list(
          mode = 'multiple', 
          selected = initial_selected_rows,
          target = 'row'
        ), 
        options = list(
          dom = 'fti',
          paging = FALSE,
          columnDefs = list(list(className = 'dt-center', targets = "_all")),
          scrollX = TRUE,
          scrollY = "65vh"
        ),
        callback = JS(
          sprintf(
            "table.on('init.dt', function() {
               var container = $(table.table().container());
               var filterBox = container.find('div.dataTables_filter');
               if (filterBox.find('.toggle-selection-btn').length === 0) {
                 $('<button type=\"button\" class=\"btn btn-secondary btn-sm toggle-selection-btn\" style=\"margin-left:8px;\">Toggle Selection</button>')
                   .appendTo(filterBox)
                   .on('click', function() {
                     Shiny.setInputValue('%s', Date.now(), {priority: 'event'});
                   });
               }
             });",
            session$ns("toggle_selection_all")
          )
        )
      ) 
    }, server = TRUE)
    
    proxy <- DT::dataTableProxy('sample_table')
    
    output$read_count_filter_ui <- renderUI({
      req(meta_df_for_dt())
      df <- meta_df_for_dt()
      max_reads <- max(df$`Read Count`, na.rm = TRUE)
      if (!is.finite(max_reads)) max_reads <- default_min_reads
      sliderInput(
        session$ns("min_read_count"),
        "",
        min = 0,
        max = ceiling(max_reads),
        value = min(default_min_reads, ceiling(max_reads)),
        step = 100
      )
    })

    observeEvent(ps_obj_initial(), {
      last_applied_selection(NULL)
      ps_filtered(ps_obj_initial())
      ps_selected_raw(ps_obj_initial())
    }, ignoreInit = FALSE, ignoreNULL = FALSE)

    observeEvent(input$min_read_count, {
      req(meta_df_for_dt())
      df <- meta_df_for_dt()
      selected_rows <- which(df$`Read Count` >= input$min_read_count)
      DT::selectRows(proxy, selected_rows)
    }, ignoreInit = TRUE)
    
    observe({
      req(meta_df_for_dt())
      df <- meta_df_for_dt()
      group_vars <- setdiff(colnames(df), c("SampleID", "Read Count"))
      selected_var <- input$toggle_group_var
      if (is.null(selected_var) || !selected_var %in% group_vars) {
        selected_var <- if (length(group_vars) > 0) group_vars[1] else NULL
      }
      updateSelectInput(session, "toggle_group_var", choices = group_vars, selected = selected_var)
    })
    
    observeEvent(list(meta_df_for_dt(), input$toggle_group_var), {
      df <- meta_df_for_dt()
      toggle_group_var <- resolve_meta_colname(input$toggle_group_var, colnames(df))
      validate(need(toggle_group_var %in% colnames(df), "Selected group variable is not available."))
      
      level_choices <- sort(unique(as.character(df[[toggle_group_var]])))
      level_choices <- level_choices[!is.na(level_choices) & nzchar(level_choices)]
      selected_levels <- isolate(input$toggle_group_levels)
      if (is.null(selected_levels) || !nzchar(selected_levels)) selected_levels <- character(0)
      selected_level_safe <- NULL
      if (length(selected_levels) > 0 && selected_levels[1] %in% level_choices) {
        selected_level_safe <- selected_levels[1]
      }
      updateSelectizeInput(
        session,
        "toggle_group_levels",
        choices = level_choices,
        selected = selected_level_safe,
        server = TRUE
      )
    }, ignoreInit = FALSE)
    
    observeEvent(input$reset_selection, {
      req(meta_df_for_dt())
      DT::selectRows(proxy, seq_len(nrow(meta_df_for_dt())))
    })

    observeEvent(input$unselect_all, {
      req(meta_df_for_dt())
      DT::selectRows(proxy, integer(0))
    })
    
    observeEvent(input$toggle_selection_all, {
      req(meta_df_for_dt())
      target_rows <- input$sample_table_rows_all
      if (is.null(target_rows) || length(target_rows) == 0) {
        return()
      }
      current_selected <- input$sample_table_rows_selected
      if (is.null(current_selected)) {
        current_selected <- integer(0)
      }
      
      selected_in_target <- intersect(current_selected, target_rows)
      toggled_target <- setdiff(target_rows, selected_in_target)
      selected_outside_target <- setdiff(current_selected, target_rows)
      new_selected <- sort(unique(c(selected_outside_target, toggled_target)))
      
      DT::selectRows(proxy, new_selected)
    })
    
    observeEvent(input$toggle_group_selection, {
      req(meta_df_for_dt(), input$toggle_group_var)
      selected_level <- input$toggle_group_levels
      validate(need(!is.null(selected_level) && nzchar(selected_level),
                    "Select one group level to toggle."))
      
      df <- meta_df_for_dt()
      toggle_group_var <- resolve_meta_colname(input$toggle_group_var, colnames(df))
      validate(need(toggle_group_var %in% colnames(df), "Selected group variable is not available."))
      
      group_values <- as.character(df[[toggle_group_var]])
      target_rows <- which(group_values == selected_level)
      validate(need(length(target_rows) > 0, "No rows match the selected group level."))
      
      current_selected <- input$sample_table_rows_selected
      if (is.null(current_selected)) {
        current_selected <- integer(0)
      }
      
      selected_in_target <- intersect(current_selected, target_rows)
      toggled_target <- setdiff(target_rows, selected_in_target)
      selected_outside_target <- setdiff(current_selected, target_rows)
      new_selected <- sort(unique(c(selected_outside_target, toggled_target)))
      
      DT::selectRows(proxy, new_selected)
    })

    selected_indices_for_apply <- reactive({
      req(meta_df_for_dt())
      selected_indices <- input$sample_table_rows_selected
      if (is.null(selected_indices)) {
        min_reads <- input$min_read_count
        if (is.null(min_reads)) {
          min_reads <- default_min_reads
        }
        selected_indices <- which(meta_df_for_dt()$`Read Count` >= min_reads)
      }
      sort(unique(as.integer(selected_indices)))
    })
    
    apply_current_selection <- function(selected_indices) {
      req(ps_obj_initial(), meta_df_for_dt())
      if (is.null(selected_indices)) {
        selected_indices <- integer(0)
      } else {
        selected_indices <- sort(unique(as.integer(selected_indices)))
      }
      
      if (identical(selected_indices, last_applied_selection())) {
        return(invisible(NULL))
      }

      if (length(selected_indices) == 0) {
        ps_new <- NULL
      } else {
        selected_samples <- meta_df_for_dt()$SampleID[selected_indices]
        ps_new <- phyloseq::prune_samples(selected_samples, ps_obj_initial())
      }
      ps_selected_raw(ps_new)

      if (!is.null(ps_new) && phyloseq::nsamples(ps_new) > 0) {
        ps_new <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_new) > 0, ps_new)
      }
      
      ps_filtered(ps_new)
      last_applied_selection(selected_indices)
      invisible(NULL)
    }
    
    observeEvent(active_tab(), {
      current_tab <- active_tab()
      if (!is.null(current_tab) && !identical(current_tab, "Preprocessing")) {
        selected_indices <- selected_indices_for_apply()
        total_rows <- nrow(meta_df_for_dt())
        if (length(selected_indices) == 0 && total_rows > 0) {
          selected_indices <- seq_len(total_rows)
          DT::selectRows(proxy, selected_indices)
        }
        apply_current_selection(selected_indices)
      }
    }, ignoreInit = TRUE)
    
    output$sample_count_info <- renderText({
      req(ps_obj_initial(), meta_df_for_dt())
      total_samples <- nrow(meta_df_for_dt())
      selected_count <- length(selected_indices_for_apply())
      paste0(
        sprintf("Total: %d samples\n", total_samples),
        sprintf("Selected: %d samples", selected_count)
      )
    })

    return(list(
      ps_filtered_raw = ps_filtered,
      ps_selected_raw = ps_selected_raw,
      plot_order = plot_order
    ))
  })
}
