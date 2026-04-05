library(shiny)
library(phyloseq)
library(dplyr)
library(DT)

## UI
mod_preprocessing_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        width = 3,
        wellPanel(
          h4(icon("filter"), "Sample Selection Overview"),
          p("All samples over 2000 reads are selected by default."),
          p("Unselect rows in the table to EXCLUDE samples.", style = "margin-top: -8px;"),
          hr(),
          div(style = "background: #f8f9fa; padding: 10px; border-radius: 5px;",
              verbatimTextOutput(ns("sample_count_info"))
          ),
          hr(),
          h5(icon("layer-group"), "Group-based Toggle"),
          selectInput(ns("toggle_group_var"), "Group variable", choices = NULL),
          selectizeInput(
            ns("toggle_group_levels"),
            "Group level",
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
          actionButton(ns("reset_selection"), "Reset: Select All Samples", 
                       icon = icon("sync"), class = "btn-secondary btn-sm", style = "font-size: 12px;", width = "100%")
        )
      ),
      column(
        width = 9,
        wellPanel(
          h4(icon("list-check"), "Individual Sample Selection (Click to Toggle)"),
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
    
    ps_filtered <- reactiveVal(NULL) 
    ps_normalized <- reactiveVal(NULL) 
    
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
      initial_selected_rows <- which(df$`Read Count` >= 2000)
      
      datatable(
        df,
        rownames = FALSE, 
        selection = list(
          mode = 'multiple', 
          selected = initial_selected_rows,
          target = 'row'
        ), 
        options = list(
          dom = 'ftip', 
          pageLength = 20,
          columnDefs = list(list(className = 'dt-center', targets = "_all")),
          scrollX = TRUE
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
      validate(need(input$toggle_group_var %in% colnames(df), "Selected group variable is not available."))
      
      level_choices <- sort(unique(as.character(df[[input$toggle_group_var]])))
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
      validate(need(input$toggle_group_var %in% colnames(df), "Selected group variable is not available."))
      
      group_values <- as.character(df[[input$toggle_group_var]])
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
    
    observe({      
      req(ps_obj_initial(), meta_df_for_dt())
      
      selected_indices <- input$sample_table_rows_selected
      
      if (is.null(selected_indices)) {
        ps_new <- ps_obj_initial()
      } else if (length(selected_indices) == 0) {
        ps_new <- phyloseq::prune_samples(character(0), ps_obj_initial())
      } else {
        selected_samples <- meta_df_for_dt()$SampleID[selected_indices]
        ps_new <- phyloseq::prune_samples(selected_samples, ps_obj_initial())
      }
      
      if (phyloseq::nsamples(ps_new) > 0) {
        ps_new <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_new) > 0, ps_new)
      }
      
      ps_filtered(ps_new)
    })
    
    observe({
      ps_obj <- ps_filtered()
      req(ps_obj)
      
      if (phyloseq::nsamples(ps_obj) == 0) {
        ps_normalized(ps_obj)
        return()
      }
      if (phyloseq::ntaxa(ps_obj) == 0) {
        ps_normalized(ps_obj)
        return()
      }
      
      ps_res <- tryCatch({
        phyloseq::transform_sample_counts(ps_obj, function(x) (x / sum(x)) * 100)
      }, error = function(e) {
        showNotification(
          paste0(
            "Transformation Error: ", e$message,
            " [method=tss",
            ", nsamples=", phyloseq::nsamples(ps_obj),
            ", ntaxa=", phyloseq::ntaxa(ps_obj), "]"
          ),
          type = "error"
        )
        ps_obj
      })
      
      ps_normalized(ps_res)
    })

    output$sample_count_info <- renderText({
      req(ps_obj_initial())
      ps_now <- ps_normalized()
      
      res <- sprintf("Initial Total: %d samples\n", phyloseq::nsamples(ps_obj_initial()))
      if (!is.null(ps_now)) {
        res <- paste0(res, sprintf("Currently Selected: %d samples\n", phyloseq::nsamples(ps_now)))
        res <- paste0(res, sprintf("Active Taxa (ASVs): %d", phyloseq::ntaxa(ps_now)))
      } else {
        res <- paste0(res, "Initializing selection...")
      }
      res
    })
    
    return(list(
      ps_filtered_raw = ps_filtered,
      ps_normalized = ps_normalized
    ))
  })
}
