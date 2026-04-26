library(shiny)
library(phyloseq)
library(dplyr)
library(DT)

## UI
mod_preprocessing_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML("
      .well h4 { font-size: 16px; }
      .well h5 { font-size: 13px; }
      .well .control-label { font-size: 12px; }
      .well .checkbox label { font-size: 12px; }
      .well .form-control { font-size: 12px; }
      .well .btn { font-size: 11px; }
      .well .irs-grid-text { display: none !important; }
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
          h5(icon("layer-group"), "Group-based Toggle"),
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
          br(),
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
    last_applied_indices <- reactiveVal(NULL)
    default_min_reads <- 2000
    
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
      last_applied_indices(NULL)
      ps_filtered(ps_obj_initial())
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
      
      last_indices <- last_applied_indices()
      if (is.null(last_indices)) {
        last_indices <- integer(0)
      } else {
        last_indices <- sort(unique(as.integer(last_indices)))
      }
      if (identical(selected_indices, last_indices)) {
        return(invisible(NULL))
      }
      
      if (length(selected_indices) == 0) {
        ps_new <- NULL
      } else {
        selected_samples <- meta_df_for_dt()$SampleID[selected_indices]
        ps_new <- phyloseq::prune_samples(selected_samples, ps_obj_initial())
      }
      
      if (!is.null(ps_new) && phyloseq::nsamples(ps_new) > 0) {
        ps_new <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_new) > 0, ps_new)
      }
      
      ps_filtered(ps_new)
      last_applied_indices(selected_indices)
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
        sprintf("Total: %d\n", total_samples),
        sprintf("Selected: %d", selected_count)
      )
    })
    
    return(list(
      ps_filtered_raw = ps_filtered
    ))
  })
}
