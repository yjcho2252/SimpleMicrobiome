## UI
mod_spieceasi_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML("
      .tabbable > .nav-tabs {
        flex-wrap: nowrap;
        overflow-x: auto;
        overflow-y: hidden;
      }
      .tabbable > .nav-tabs .nav-link {
        font-size: 12px;
        padding: 5px 8px;
        white-space: nowrap;
      }
      .spieceasi-dt-panel {
        margin-top: 10px;
        margin-bottom: 36px;
        overflow: visible;
      }
      .spieceasi-dt-panel::after {
        content: '';
        display: block;
        clear: both;
      }
      .spieceasi-dt-panel .dataTables_wrapper {
        padding-bottom: 20px;
      }
      .well h4 { font-size: 16px; }
      .well h5 { font-size: 13px; }
      .well .control-label { font-size: 12px; }
      .well .checkbox label { font-size: 12px; }
      .well .form-control { font-size: 12px; }
      .well .btn { font-size: 11px; }
      .spieceasi-run-btn .spieceasi-run-icon {
        display: none;
        margin-right: 6px;
      }
      .spieceasi-run-btn .spieceasi-run-icon i {
        display: inline-block;
      }
      .spieceasi-run-btn.is-running .spieceasi-run-icon {
        display: inline-block;
        animation: spieceasi-run-spin 0.9s linear infinite;
      }
      .spieceasi-run-btn.is-running .spieceasi-run-label {
        opacity: 0.85;
      }
      @keyframes spieceasi-run-spin {
        from { transform: rotate(0deg); }
        to { transform: rotate(360deg); }
      }
    ")),
    plot_scroll_layout(
      sidebarPanel(
        width = 2,
        h4(icon("diagram-project"), "SPIEC-EASI"),
        hr(),
        selectInput(ns("group_var"), "1. Primary Grouping Variable", choices = NULL),
        checkboxInput(ns("use_subgroup"), "Select Subgroup", value = FALSE),
        conditionalPanel(
          condition = "input.use_subgroup",
          ns = ns,
          selectInput(ns("primary_level"), "Primary Level to Include", choices = NULL),
          selectInput(ns("secondary_var"), "Secondary Grouping Variable", choices = c("None"), selected = "None")
        ),
        selectizeInput(
          ns("group_levels"),
          "2. Group Levels",
          choices = NULL,
          selected = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Select one or more levels",
            plugins = list("remove_button")
          )
        ),
        selectInput(
          ns("analysis_mode"),
          "3. Analysis Mode",
          choices = c("Single network", "Compare two groups"),
          selected = "Single network"
        ),
        selectInput(ns("tax_level"), "4. Taxonomic Level", choices = c("ASV", "Genus", "Species", "Strain"), selected = "Genus"),
        numericInput(ns("prevalence_filter_pct"), "5. Prevalence Filter (%)", value = 10, min = 0, max = 100, step = 1),
        selectInput(ns("node_size_by"), "6. Node Size", choices = c("Connectivity", "Abundance"), selected = "Connectivity"),
        selectInput(ns("node_color_by"), "7. Node Color", choices = c("None"), selected = "None"),
        tags$details(
          style = "margin-bottom: 8px;",
          tags$summary("Advanced Options"),
          numericInput(ns("min_edge_weight"), "8. Minimum Absolute Edge Weight", value = 0.1, min = 0, max = 1, step = 0.01),
          numericInput(ns("max_taxa"), "9. Max Taxa for Network", value = 200, min = 20, max = 2000, step = 10),
          numericInput(ns("seed"), "10. Seed", value = 1001, min = 1, step = 1)
        ),
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Display Options"),
        numericInput(ns("plot_width"), "Plot Width (px)", value = 800, min = 400, max = 2400, step = 50),
        numericInput(ns("plot_height"), "Plot Height (px)", value = 550, min = 300, max = 2400, step = 50),
        tags$div(
          style = "display: flex; align-items: center; gap: 8px; flex-wrap: wrap;",
          tags$button(
            id = ns("run_network_btn"),
            type = "button",
            class = "action-button btn btn-danger spieceasi-run-btn",
            `data-val` = 0,
            style = "font-size: 12px;",
            tags$span(class = "spieceasi-run-icon", icon("spinner", class = "fa-spin")),
            tags$span(class = "spieceasi-run-label", "Run SPIEC-EASI")
          ),
          tags$span("May take a long time.", style = "font-size: 11px; color: #b94a48;")
        ),
        tags$script(HTML(
          "Shiny.addCustomMessageHandler('toggle-spieceasi-run-btn', function(msg) {
             var btn = document.getElementById(msg.id);
             if (!btn) return;
             btn.disabled = !!msg.disabled;
             btn.classList.toggle('is-running', !!msg.running);
             var label = btn.querySelector('.spieceasi-run-label');
             if (label && msg.label) label.textContent = msg.label;
           });
           Shiny.addCustomMessageHandler('set-tab-container-width', function(msg) {
             var el = document.getElementById(msg.id);
             if (!el) return;
             el.style.width = msg.width;
             el.style.maxWidth = '100%';
           });"
        ))
      ),
      mainPanel(
        h4("SPIEC-EASI Network"),
        analysis_run_view(ns, "run_network_btn",
        tags$div(
          id = ns("spieceasi_tab_container"),
          style = "max-width: 100%;",
          tabsetPanel(
            id = ns("spieceasi_active_tab"),
          tabPanel(
            "All",
            downloadButton(ns("download_network_plot_all"), "Download Plot (PNG)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
            tags$div(
              style = "margin-top: 8px;",
              plotOutput(ns("network_plot_all"), height = "auto")
            )
          ),
          tabPanel(
            "Connected",
            downloadButton(ns("download_network_plot_connected"), "Download Plot (PNG)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
            tags$div(
              style = "margin-top: 8px;",
              plotOutput(ns("network_plot_connected"), height = "auto")
            )
          ),
          tabPanel(
            "Table",
            downloadButton(ns("download_edge_table_all"), "Download Table (TSV)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
            tags$div(
              class = "spieceasi-dt-panel",
              DTOutput(ns("edge_table_all"))
            )
          ),
          tabPanel("Summary", verbatimTextOutput(ns("network_summary"))),
          tabPanel(
            "Comparison Network",
            downloadButton(ns("download_comparison_network_plot"), "Download Plot (PNG)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
            tags$div(
              style = "margin-top: 8px;",
              plotOutput(ns("comparison_network_plot"), height = "auto")
            )
          ),
          tabPanel(
            "Differential Edges",
            downloadButton(ns("download_comparison_edge_table"), "Download Table (TSV)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
            tags$div(
              class = "spieceasi-dt-panel",
              DTOutput(ns("comparison_edge_table"))
            )
          ),
          tabPanel("Comparison Summary", verbatimTextOutput(ns("comparison_summary"))),
          tabPanel(
            "Hub Table",
            downloadButton(ns("download_hub_table"), "Download Table (TSV)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
            tags$div(
              class = "spieceasi-dt-panel",
              DTOutput(ns("hub_table"))
            )
          )
          )
        ),
        plot_legend_output(ns, "spieceasi_legend_box", c("network_plot_all", "network_plot_connected", "comparison_network_plot", "edge_table_all", "network_summary", "comparison_edge_table", "comparison_summary", "hub_table")),
        uiOutput(ns("spieceasi_status_separator")),
        h5(icon("circle-info"), "SPIEC-EASI Status"),
        uiOutput(ns("spieceasi_status_box"))
        )
      )
    )
  )
}

## Server
mod_spieceasi_server <- function(id, ps_obj) {
  moduleServer(id, function(input, output, session) {
    output$run_prompt <- analysis_run_prompt("Run SPIEC-EASI")
    supported_taxonomic_ranks <- c("Genus", "Species", "Strain")

    observeEvent(ps_obj(), {
      ps <- ps_obj()
      available_ranks <- tryCatch(phyloseq::rank_names(ps), error = function(e) character(0))
      rank_choices <- c("ASV", supported_taxonomic_ranks[supported_taxonomic_ranks %in% available_ranks])
      current_rank <- isolate(input$tax_level)
      selected_rank <- if (!is.null(current_rank) && current_rank %in% rank_choices) {
        current_rank
      } else if ("Genus" %in% rank_choices) {
        "Genus"
      } else {
        "ASV"
      }
      updateSelectInput(session, "tax_level", choices = rank_choices, selected = selected_rank)
    }, ignoreInit = FALSE)

    draw_wait_message <- function(message_text) {
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      graphics::par(mar = c(0, 0, 0, 0))
      graphics::plot.new()
      graphics::text(0.5, 0.5, message_text, cex = 0.85)
    }

    observe({
      width_px <- suppressWarnings(as.integer(input$plot_width))
      if (!is.finite(width_px) || is.na(width_px) || width_px <= 0) width_px <- 800L
      session$sendCustomMessage(
        "set-tab-container-width",
        list(id = session$ns("spieceasi_tab_container"), width = paste0(width_px, "px"))
      )
    })

    apply_disambiguated_taxrank <- function(ps, tax_level) {
      if (is.null(ps) || identical(tax_level, "ASV")) {
        return(ps)
      }
      tt_obj <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
      if (is.null(tt_obj)) {
        return(ps)
      }
      tt <- as.data.frame(tt_obj, stringsAsFactors = FALSE)
      tax_cols <- colnames(tt)
      if (!tax_level %in% tax_cols) {
        return(ps)
      }
      target_raw <- as.character(tt[[tax_level]])
      target_norm <- tolower(trimws(target_raw))
      target_norm <- gsub("^[a-z]__", "", target_norm)
      idx_placeholder <- !is.na(target_norm) & grepl("(uncultured|unassigned)", target_norm)
      if (!any(idx_placeholder)) {
        return(ps)
      }
      tax_ranks <- phyloseq::rank_names(ps)
      rank_pos <- match(tax_level, tax_ranks)
      parent_rank <- NULL
      if (!is.na(rank_pos) && rank_pos > 1) {
        parent_candidates <- rev(tax_ranks[seq_len(rank_pos - 1)])
        parent_candidates <- parent_candidates[parent_candidates %in% tax_cols]
        if (length(parent_candidates) > 0) parent_rank <- parent_candidates[1]
      }
      if (is.null(parent_rank)) {
        parent_val <- rep("UnclassifiedParent", nrow(tt))
      } else {
        parent_val <- as.character(tt[[parent_rank]])
        parent_val[is.na(parent_val) | !nzchar(parent_val)] <- "UnclassifiedParent"
      }
      tt[[tax_level]][idx_placeholder] <- paste0(parent_val[idx_placeholder], "|", target_raw[idx_placeholder])
      phyloseq::tax_table(ps) <- phyloseq::tax_table(as.matrix(tt))
      ps
    }

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

    group_var_resolved <- reactive({
      req(ps_obj(), input$group_var)
      meta_df <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
      resolve_meta_colname(input$group_var, colnames(meta_df))
    })

    secondary_var_selected <- reactive({
      req(ps_obj())
      if (!isTRUE(input$use_subgroup)) {
        return(NULL)
      }
      req(input$secondary_var)
      if (identical(input$secondary_var, "None")) {
        return(NULL)
      }
      meta_df <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
      resolve_meta_colname(input$secondary_var, colnames(meta_df))
    })

    observeEvent(input$use_subgroup, {
      if (!isTRUE(input$use_subgroup)) {
        if (!is.null(input$primary_level) && !identical(as.character(input$primary_level), "All")) {
          updateSelectInput(session, "primary_level", selected = "All")
        }
        if (!is.null(input$secondary_var) && !identical(input$secondary_var, "None")) {
          updateSelectInput(session, "secondary_var", selected = "None")
        }
      }
    }, ignoreInit = TRUE)

    observeEvent(ps_obj(), {
      req(ps_obj())
      meta_df <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
      meta_cols <- colnames(meta_df)
      candidate_groups <- setdiff(meta_cols, "SampleID")
      group_choices <- candidate_groups[vapply(candidate_groups, function(col) {
        level_values <- as.character(meta_df[[col]])
        level_values <- level_values[!is.na(level_values) & nzchar(level_values)]
        length(unique(level_values)) <= 10
      }, logical(1))]
      selected_group <- if (length(group_choices) > 0) group_choices[1] else NULL
      updateSelectInput(session, "group_var", choices = group_choices, selected = selected_group)
    }, ignoreNULL = FALSE)

    observeEvent(list(ps_obj(), input$group_var), {
      req(ps_obj(), input$group_var)
      meta_df <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
      primary_var <- group_var_resolved()
      validate(
        need(primary_var %in% colnames(meta_df), paste0("ERROR: Metadata variable '", input$group_var, "' not found."))
      )

      level_choices <- sort(unique(as.character(meta_df[[primary_var]])))
      level_choices <- level_choices[!is.na(level_choices) & nzchar(level_choices)]
      level_choices <- c("All", level_choices)
      selected_primary_level <- isolate(input$primary_level)
      if (is.null(selected_primary_level) || !selected_primary_level %in% level_choices) {
        selected_primary_level <- "All"
      }
      updateSelectInput(session, "primary_level", choices = level_choices, selected = selected_primary_level)

      secondary_choices <- setdiff(colnames(meta_df), c("SampleID", primary_var))
      selected_secondary <- isolate(input$secondary_var)
      if (is.null(selected_secondary) || (!selected_secondary %in% c("None", secondary_choices))) {
        selected_secondary <- "None"
      }
      updateSelectInput(
        session,
        "secondary_var",
        choices = c("None" = "None", secondary_choices),
        selected = selected_secondary
      )
    }, ignoreNULL = FALSE)

    observeEvent(list(ps_obj(), input$group_var, input$primary_level, input$secondary_var), {
      req(ps_obj(), input$group_var, input$primary_level, input$secondary_var)
      meta_df <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
      primary_var <- group_var_resolved()
      secondary_var <- secondary_var_selected()
      validate(
        need(primary_var %in% colnames(meta_df), paste0("ERROR: Metadata variable '", input$group_var, "' not found."))
      )
      primary_values <- as.character(meta_df[[primary_var]])
      if (is.null(secondary_var)) {
        level_choices <- sort(unique(primary_values))
        level_choices <- level_choices[!is.na(level_choices) & nzchar(level_choices)]
      } else {
        validate(
          need(secondary_var %in% colnames(meta_df), paste0("ERROR: Metadata variable '", input$secondary_var, "' not found."))
        )
        secondary_values <- as.character(meta_df[[secondary_var]])
        if (!is.null(input$primary_level) && identical(as.character(input$primary_level), "All")) {
          subset_idx <- !is.na(primary_values)
        } else {
          subset_idx <- !is.na(primary_values) & (primary_values == input$primary_level)
        }
        level_choices <- sort(unique(secondary_values[subset_idx]))
        level_choices <- level_choices[!is.na(level_choices) & nzchar(level_choices)]
      }

      selected_levels <- isolate(input$group_levels)
      if (is.null(selected_levels)) selected_levels <- character(0)
      selected_levels <- intersect(selected_levels, level_choices)
      if (length(selected_levels) == 0) selected_levels <- head(level_choices, 2)
      updateSelectizeInput(session, "group_levels", choices = level_choices, selected = selected_levels, server = TRUE)
    }, ignoreNULL = FALSE)

    observe({
      req(ps_obj())
      ps <- ps_obj()
      color_choices <- "None"
      tt <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
      if (!is.null(tt)) {
        tax_ranks <- setdiff(colnames(tt), "Kingdom")
        if (!is.null(tax_ranks) && length(tax_ranks) > 0) {
          color_choices <- c("None", tax_ranks)
        }
      }
      selected_color <- input$node_color_by
      if (is.null(selected_color) || !selected_color %in% color_choices) {
        selected_color <- "None"
      }
      updateSelectInput(session, "node_color_by", choices = color_choices, selected = selected_color)
    })

    build_network_inputs <- reactive({
      req(ps_obj())
      ps <- ps_obj()
      validate(
        need(phyloseq::nsamples(ps) > 2, "At least 3 samples are required for network analysis."),
        need(phyloseq::ntaxa(ps) > 2, "At least 3 taxa are required for network analysis.")
      )

      ps_use <- ps
      if (input$tax_level != "ASV") {
        validate(
          need(!is.null(phyloseq::tax_table(ps_use, errorIfNULL = FALSE)), "Taxonomy table is required to aggregate by taxonomic level.")
        )
        tax_cols <- colnames(phyloseq::tax_table(ps_use))
        validate(
          need(input$tax_level %in% tax_cols, paste("Taxonomic rank", input$tax_level, "not found in taxonomy table."))
        )
        ps_use <- apply_disambiguated_taxrank(ps_use, input$tax_level)
        ps_use <- phyloseq::tax_glom(ps_use, taxrank = input$tax_level, NArm = TRUE)
        tax_df <- as.data.frame(phyloseq::tax_table(ps_use), stringsAsFactors = FALSE)
        rank_values <- as.character(tax_df[[input$tax_level]])
        rank_values[is.na(rank_values) | rank_values == ""] <- phyloseq::taxa_names(ps_use)[is.na(rank_values) | rank_values == ""]
        phyloseq::taxa_names(ps_use) <- make.unique(rank_values)
      }

      otu_mat <- as.matrix(phyloseq::otu_table(ps_use))
      if (!phyloseq::taxa_are_rows(phyloseq::otu_table(ps_use))) otu_mat <- t(otu_mat)

      prevalence_cutoff <- as.numeric(input$prevalence_filter_pct)
      if (is.na(prevalence_cutoff)) prevalence_cutoff <- 10
      prevalence_cutoff <- max(0, min(100, prevalence_cutoff))
      keep_prev <- rowMeans(otu_mat > 0) * 100 >= prevalence_cutoff
      otu_mat <- otu_mat[keep_prev, , drop = FALSE]
      validate(need(nrow(otu_mat) >= 3, "Too few taxa remain after prevalence filtering. Lower the cutoff."))

      max_taxa <- as.integer(input$max_taxa)
      if (is.na(max_taxa) || max_taxa < 3) max_taxa <- 200
      if (nrow(otu_mat) > max_taxa) {
        otu_mat <- otu_mat[order(rowSums(otu_mat), decreasing = TRUE)[seq_len(max_taxa)], , drop = FALSE]
      }

      sample_ids <- colnames(otu_mat)
      sd <- phyloseq::sample_data(ps_use, errorIfNULL = FALSE)
      if (is.null(sd)) {
        sample_df <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
        rownames(sample_df) <- sample_ids
      } else {
        sample_df <- as.data.frame(sd, stringsAsFactors = FALSE)
        sample_df$sample_id <- rownames(sample_df)
        sample_df <- sample_df[sample_ids, , drop = FALSE]
      }

      taxa_annotation <- data.frame(node_name = rownames(otu_mat), stringsAsFactors = FALSE)
      tt_use <- phyloseq::tax_table(ps_use, errorIfNULL = FALSE)
      if (!is.null(tt_use)) {
        tax_df <- as.data.frame(tt_use, stringsAsFactors = FALSE)
        tax_df$node_name <- rownames(tax_df)
        taxa_annotation <- dplyr::left_join(taxa_annotation, tax_df, by = "node_name")
      }

      list(otu_mat = otu_mat, sample_df = sample_df, taxa_annotation = taxa_annotation)
    })

    network_running <- reactiveVal(FALSE)

    network_run_nonce <- reactiveVal(0L)
    network_next_allowed_at <- reactiveVal(as.POSIXct(NA))
    last_run_settings <- reactiveVal(NULL)

    normalize_setting <- function(x) {
      if (is.null(x)) return(character(0))
      sort(unique(as.character(x)))
    }

    current_analysis_settings <- reactive({
      ps <- ps_obj()
      req(ps)
      otu_signature_matrix <- as.matrix(phyloseq::otu_table(ps))
      otu_signature_values <- as.numeric(otu_signature_matrix)
      otu_signature_values[!is.finite(otu_signature_values)] <- 0
      taxonomy_obj <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
      taxonomy_signature <- if (is.null(taxonomy_obj)) {
        character(0)
      } else {
        taxonomy_matrix <- as.matrix(taxonomy_obj)
        c(colnames(taxonomy_matrix), as.character(taxonomy_matrix))
      }
      metadata_values <- list()
      sample_data_obj <- phyloseq::sample_data(ps, errorIfNULL = FALSE)
      if (!is.null(sample_data_obj)) {
        metadata_df <- as.data.frame(sample_data_obj, stringsAsFactors = FALSE)
        selected_metadata <- setdiff(
          unique(c(as.character(input$group_var), as.character(input$secondary_var))),
          c(NA_character_, "", "None", "All")
        )
        for (column_name in intersect(selected_metadata, colnames(metadata_df))) {
          metadata_values[[column_name]] <- as.character(metadata_df[[column_name]])
        }
      }
      list(
        dataset = list(
          sample_names = sort(phyloseq::sample_names(ps)),
          taxa_names = sort(phyloseq::taxa_names(ps)),
          rank_names = tryCatch(phyloseq::rank_names(ps), error = function(e) character(0)),
          taxonomy_signature = taxonomy_signature,
          otu_signature = c(
            dim(otu_signature_matrix),
            sum(otu_signature_values),
            sum(otu_signature_values^2),
            sum(otu_signature_values * ((seq_along(otu_signature_values) %% 1009L) + 1L))
          ),
          metadata_values = metadata_values
        ),
        tax_level = normalize_setting(input$tax_level),
        group_var = normalize_setting(input$group_var),
        primary_level = normalize_setting(input$primary_level),
        secondary_var = normalize_setting(input$secondary_var),
        group_levels = normalize_setting(input$group_levels),
        analysis_mode = normalize_setting(input$analysis_mode),
        prevalence_filter_pct = suppressWarnings(as.numeric(input$prevalence_filter_pct)),
        max_taxa = suppressWarnings(as.integer(input$max_taxa)),
        seed = suppressWarnings(as.integer(input$seed))
      )
    })

    results_stale <- reactive({
      previous <- last_run_settings()
      if (is.null(previous)) return(FALSE)
      !identical(previous, current_analysis_settings())
    })

    observeEvent(input$run_network_btn, {
      now_ts <- Sys.time()
      next_allowed <- network_next_allowed_at()
      if (isTRUE(network_running())) {
        showNotification("SPIEC-EASI is already running. Please wait for completion.", type = "message", duration = 3)
        return(invisible(NULL))
      }
      if (!is.na(next_allowed) && now_ts < next_allowed) {
        wait_sec <- ceiling(as.numeric(difftime(next_allowed, now_ts, units = "secs")))
        showNotification(paste0("Please wait ", wait_sec, " second(s) before running again."), type = "message", duration = 2)
        return(invisible(NULL))
      }
      network_next_allowed_at(now_ts + 10)
      network_run_nonce(network_run_nonce() + 1L)
    }, ignoreInit = TRUE)

    network_result <- eventReactive(network_run_nonce(), {
      req(!is.null(input$run_network_btn) && input$run_network_btn > 0)
      req(build_network_inputs())
      validate(
        need(requireNamespace("NetCoMi", quietly = TRUE), "NetCoMi package is not installed. Please install 'NetCoMi'."),
        need(requireNamespace("igraph", quietly = TRUE), "igraph package is not installed. Please install 'igraph'.")
      )
      network_running(TRUE)
      run_settings <- isolate(current_analysis_settings())

      session$sendCustomMessage("toggle-spieceasi-run-btn", list(id = session$ns("run_network_btn"), disabled = TRUE, running = TRUE, label = "Running..."))
      on.exit({
        network_running(FALSE)
        session$sendCustomMessage("toggle-spieceasi-run-btn", list(id = session$ns("run_network_btn"), disabled = FALSE, running = FALSE, label = "Run SPIEC-EASI"))
      }, add = TRUE)

      set.seed(as.integer(input$seed))
      built <- build_network_inputs()
      otu_mat <- built$otu_mat
      sample_df <- built$sample_df
      primary_var_resolved <- group_var_resolved()
      secondary_var_resolved <- secondary_var_selected()
      effective_group_var <- if (is.null(secondary_var_resolved)) primary_var_resolved else secondary_var_resolved

      result <- tryCatch({
        withProgress(message = "Running NetCoMi SPIEC-EASI...", value = 0, {
          sample_ids <- colnames(otu_mat)
          if (!is.null(primary_var_resolved) && primary_var_resolved %in% colnames(sample_df)) {
            primary_values <- as.character(sample_df[[primary_var_resolved]])
            names(primary_values) <- rownames(sample_df)
            primary_values <- primary_values[sample_ids]
          } else {
            primary_values <- rep("All", length(sample_ids))
            names(primary_values) <- sample_ids
          }
          if (!is.null(secondary_var_resolved) && secondary_var_resolved %in% colnames(sample_df)) {
            selected_primary_level <- input$primary_level
            if (is.null(selected_primary_level) || identical(as.character(selected_primary_level), "All")) {
              sample_ids <- sample_ids[!is.na(primary_values)]
            } else {
              sample_ids <- sample_ids[!is.na(primary_values) & primary_values == selected_primary_level]
            }
          }

          if (!is.null(effective_group_var) && effective_group_var != "All" && effective_group_var %in% colnames(sample_df)) {
            group_values <- as.character(sample_df[[effective_group_var]])
            names(group_values) <- rownames(sample_df)
            group_values <- group_values[sample_ids]
            group_values <- trimws(group_values)
            group_values[is.na(group_values) | group_values == ""] <- "(Missing)"
            group_samples <- split(sample_ids, group_values)
          } else {
            group_samples <- list(All = sample_ids)
          }

          selected_levels <- input$group_levels
          analysis_mode <- if (is.null(input$analysis_mode)) "Single network" else input$analysis_mode
          if (!is.null(effective_group_var) && effective_group_var != "All") {
            if (is.null(selected_levels)) selected_levels <- character(0)
            selected_levels <- trimws(as.character(selected_levels))
            selected_levels <- selected_levels[nzchar(selected_levels)]
            validate(need(length(selected_levels) > 0, "Select one or more group levels."))
            if (identical(analysis_mode, "Compare two groups")) {
              validate(need(length(selected_levels) == 2, "Compare mode requires selecting exactly 2 group levels."))
            }
            available_levels <- names(group_samples)
            available_levels_key <- trimws(as.character(available_levels))
            selected_levels_key <- trimws(as.character(selected_levels))
            matched_levels <- available_levels[available_levels_key %in% selected_levels_key]
            debug_msg <- paste0(
              "No selected group levels are available for network estimation.\n",
              "Selected (raw): ", paste(selected_levels, collapse = ", "), "\n",
              "Available (raw): ", paste(available_levels, collapse = ", "), "\n",
              "Selected (normalized): ", paste(selected_levels_key, collapse = ", "), "\n",
              "Available (normalized): ", paste(available_levels_key, collapse = ", "), "\n",
              "Group variable (effective): ", as.character(effective_group_var)
            )
            validate(need(length(matched_levels) > 0, debug_msg))
            group_samples <- group_samples[matched_levels]
          }
          validate(need(length(group_samples) > 0, "No group level is available for network estimation."))
          if (identical(analysis_mode, "Compare two groups")) {
            validate(need(length(group_samples) == 2, "Compare mode requires exactly 2 eligible groups."))
          }

          get_association_matrix <- function(x) {
            if (is.null(x)) return(NULL)
            if (is.list(x) && !is.data.frame(x)) {
              if (!is.null(x$refit)) return(as.matrix(x$refit))
              if (length(x) >= 1) return(as.matrix(x[[1]]))
            }
            as.matrix(x)
          }

          build_group_network <- function(group_name, ids) {
            sub_otu <- otu_mat[, ids, drop = FALSE]
            run_spieceasi <- function(ncores) {
              NetCoMi::netConstruct(
                data = t(sub_otu),
                dataType = "counts",
                measure = "spieceasi",
                measurePar = list(
                  pulsar.params = list(ncores = ncores)
                ),
                verbose = 0
              )
            }

            default_ncores <- if (.Platform$OS.type == "windows") 1L else 4L
            net_obj <- tryCatch(
              run_spieceasi(default_ncores),
              error = function(e) {
                msg <- conditionMessage(e)
                is_parallel_worker_error <- grepl("jobs failed", msg, ignore.case = TRUE) && grepl("sparseiCov", msg, ignore.case = TRUE)
                if (!is_parallel_worker_error) stop(e)
                showNotification("SPIEC-EASI parallel run failed due to worker export issue. Retrying with 1 core.", type = "warning", duration = 8)
                run_spieceasi(1)
              }
            )

            adja <- net_obj$adjaMat1
            if (is.null(adja)) adja <- net_obj$adjaMat
            validate(need(!is.null(adja), paste("NetCoMi returned no adjacency matrix for SPIEC-EASI group:", group_name)))
            adja_mat <- as.matrix(adja)

            assoc <- get_association_matrix(net_obj$assoMat1)
            if (is.null(assoc)) assoc <- get_association_matrix(net_obj$assoMat)
            if (is.null(assoc)) assoc <- get_association_matrix(net_obj$assoEst1)
            if (is.null(assoc)) assoc <- get_association_matrix(net_obj$assoEst)

            weight_mat <- adja_mat
            if (!is.null(assoc) &&
                is.matrix(assoc) &&
                nrow(assoc) == nrow(adja_mat) &&
                ncol(assoc) == ncol(adja_mat)) {
              if (!is.null(rownames(adja_mat)) &&
                  !is.null(colnames(adja_mat)) &&
                  !is.null(rownames(assoc)) &&
                  !is.null(colnames(assoc)) &&
                  all(rownames(adja_mat) %in% rownames(assoc)) &&
                  all(colnames(adja_mat) %in% colnames(assoc))) {
                assoc <- assoc[rownames(adja_mat), colnames(adja_mat), drop = FALSE]
              }
              weight_mat[adja_mat != 0] <- assoc[adja_mat != 0]
            }

            graph <- igraph::graph_from_adjacency_matrix(weight_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
            edge_df <- igraph::as_data_frame(graph, what = "edges")
            if (nrow(edge_df) > 0 && "weight" %in% colnames(edge_df)) {
              edge_df$weight <- as.numeric(edge_df$weight)
              edge_df$sign <- ifelse(edge_df$weight > 0, "Positive", ifelse(edge_df$weight < 0, "Negative", "Zero"))
            } else {
              edge_df$weight <- numeric(0)
              edge_df$sign <- character(0)
            }
            edge_df$group <- rep(group_name, nrow(edge_df))
            sample_totals <- colSums(sub_otu, na.rm = TRUE)
            sample_totals[sample_totals <= 0 | is.na(sample_totals)] <- 1
            sub_otu_rel <- sweep(sub_otu, 2, sample_totals, "/")
            abundance <- rowMeans(sub_otu_rel, na.rm = TRUE)
            list(graph = graph, edge_table = edge_df, abundance = abundance, net_obj = net_obj)
          }

          networks <- list()
          group_status <- list()
          for (gname in names(group_samples)) {
            n_samp <- length(group_samples[[gname]])
            if (n_samp < 3) {
              group_status[[gname]] <- data.frame(
                group = gname,
                n_samples = n_samp,
                status = "Skipped (<3 samples)",
                nodes = NA_integer_,
                edges = NA_integer_,
                message = "At least 3 samples are required.",
                stringsAsFactors = FALSE
              )
              next
            }

            g_res <- tryCatch(
              build_group_network(gname, group_samples[[gname]]),
              error = function(e) e
            )
            if (inherits(g_res, "error")) {
              group_status[[gname]] <- data.frame(
                group = gname,
                n_samples = n_samp,
                status = "Failed",
                nodes = NA_integer_,
                edges = NA_integer_,
                message = conditionMessage(g_res),
                stringsAsFactors = FALSE
              )
            } else {
              networks[[gname]] <- g_res
              group_status[[gname]] <- data.frame(
                group = gname,
                n_samples = n_samp,
                status = "Success",
                nodes = igraph::vcount(g_res$graph),
                edges = igraph::ecount(g_res$graph),
                message = "",
                stringsAsFactors = FALSE
              )
            }
          }

          validate(need(length(networks) > 0, "No group network was successfully estimated."))
          edge_table <- do.call(rbind, lapply(networks, `[[`, "edge_table"))
          if (is.null(edge_table)) edge_table <- data.frame()
          status_df <- do.call(rbind, group_status)
          if (is.null(status_df)) status_df <- data.frame()

          list(
            networks = networks,
            edge_table = edge_table,
            n_taxa = nrow(otu_mat),
            n_samples = ncol(otu_mat),
            group_var = if (is.null(effective_group_var) || effective_group_var == "All") "All" else effective_group_var,
            group_status = status_df,
            taxa_annotation = built$taxa_annotation
          )
        })
      }, error = function(e) {
        showNotification(paste("SPIEC-EASI execution failed:", e$message), type = "error", duration = NULL)
        NULL
      })

      validate(need(!is.null(result), "Network analysis failed. Check the error message above."))
      result$run_settings <- run_settings
      last_run_settings(run_settings)
      result
    })

    display_network_result <- reactive({
      res <- network_result()
      validate(
        need(
          !results_stale(),
          "Analysis settings changed after the last run. Click 'Run SPIEC-EASI' to update all result tabs."
        )
      )
      res
    })

    get_connected_graph <- function(g) {
      connected_nodes <- which(igraph::degree(g) > 0)
      igraph::induced_subgraph(g, vids = connected_nodes)
    }

    filter_graph_by_weight <- function(g, cutoff) {
      if (igraph::ecount(g) == 0) return(g)
      edge_weight <- igraph::E(g)$weight
      edge_weight <- as.numeric(edge_weight)
      edge_weight[!is.finite(edge_weight)] <- 0
      keep_edges <- abs(edge_weight) >= cutoff
      if (!any(keep_edges)) {
        return(igraph::delete_edges(g, igraph::E(g)))
      }
      igraph::subgraph.edges(g, eids = igraph::E(g)[keep_edges], delete.vertices = FALSE)
    }

    current_min_edge_weight <- function() {
      min_edge_weight <- suppressWarnings(as.numeric(input$min_edge_weight))
      if (is.na(min_edge_weight) || min_edge_weight < 0) min_edge_weight <- 0.1
      min_edge_weight
    }

    get_current_eligible_groups <- function(built, group_var_value, group_levels_value) {
      if (is.null(built)) return(character(0))
      if (is.null(group_var_value) || group_var_value == "All" || !group_var_value %in% colnames(built$sample_df)) {
        return("All")
      }
      gvals <- as.character(built$sample_df[[group_var_value]])
      gvals <- trimws(gvals)
      gvals[is.na(gvals) | gvals == ""] <- "(Missing)"
      cnt <- table(gvals)
      groups <- names(cnt[cnt >= 3])
      if (!is.null(group_levels_value) && length(group_levels_value) > 0) {
        groups <- intersect(groups, trimws(as.character(group_levels_value)))
      }
      groups
    }

    output$comparison_status <- renderText({
      if (isTRUE(results_stale())) {
        return("Status: Analysis settings changed after the last run. Click 'Run SPIEC-EASI' to update all result tabs.")
      }
      selected_group_var <- if (!is.null(input$secondary_var) && !identical(input$secondary_var, "None")) {
        input$secondary_var
      } else if (is.null(input$group_var)) {
        "All"
      } else {
        input$group_var
      }
      analysis_mode <- if (is.null(input$analysis_mode)) "Single network" else input$analysis_mode
      built <- tryCatch(build_network_inputs(), error = function(e) NULL)
      sample_count_line <- "Sample counts by level: (run data not ready)"
      if (!is.null(built) && selected_group_var != "All" && selected_group_var %in% colnames(built$sample_df)) {
        gvals <- as.character(built$sample_df[[selected_group_var]])
        gvals[is.na(gvals) | gvals == ""] <- "(Missing)"
        cnt <- sort(table(gvals), decreasing = TRUE)
        selected_levels <- input$group_levels
        if (!is.null(selected_levels) && length(selected_levels) > 0) {
          cnt <- cnt[names(cnt) %in% selected_levels]
        }
        if (length(cnt) == 0) {
          sample_count_line <- "Sample counts by level: None matched selection"
        } else {
          count_labels <- paste0(names(cnt), "=", as.integer(cnt), ifelse(as.integer(cnt) >= 3, " (ok)", " (<3)"))
          sample_count_line <- paste("Sample counts by level:", paste(count_labels, collapse = ", "))
        }
      } else if (!is.null(built) && selected_group_var == "All") {
        sample_count_line <- paste0("Sample counts by level: All=", ncol(built$otu_mat), " (ok)")
      }

      res <- tryCatch(network_result(), error = function(e) NULL)

      if (is.null(res)) {
        return(paste(
          paste0("Mode: ", analysis_mode),
          paste0("Group variable: ", selected_group_var),
          sample_count_line,
          "Status: Run network to compute available groups.",
          sep = "\n"
        ))
      }

      groups <- get_current_eligible_groups(built, selected_group_var, input$group_levels)
      run_groups <- names(res$networks)
      lines <- c(
        paste0("Mode: ", analysis_mode),
        paste0("Group variable: ", selected_group_var),
        paste0("Available groups (n=", length(groups), "): ", if (length(groups) == 0) "None" else paste(groups, collapse = ", ")),
        paste0("Selected levels: ", if (is.null(input$group_levels) || length(input$group_levels) == 0) "All" else paste(input$group_levels, collapse = ", ")),
        sample_count_line
      )
      if (identical(analysis_mode, "Compare two groups")) {
        if (length(run_groups) == 2) {
          lines <- c(lines, paste0("Compare target: ", paste(run_groups, collapse = " vs ")))
        } else {
          lines <- c(lines, paste0("Compare target: Need exactly 2 successful groups (current=", length(run_groups), ")"))
        }
      }

      if (!is.null(res$group_status) && nrow(res$group_status) > 0) {
        lines <- c(lines, "Run status by group:")
        for (i in seq_len(nrow(res$group_status))) {
          st <- res$group_status[i, , drop = FALSE]
          lines <- c(
            lines,
            paste0(
              "- ", st$group,
              " | samples=", st$n_samples,
              " | status=", st$status,
              ifelse(is.na(st$nodes), "", paste0(" | nodes=", st$nodes)),
              ifelse(is.na(st$edges), "", paste0(" | edges=", st$edges)),
              ifelse(is.na(st$message) || st$message == "", "", paste0(" | msg=", st$message))
            )
          )
        }
      }
      lines <- c(lines, "Status: Facet by available groups.")

      paste(lines, collapse = "\n")
    })

    summarize_network_structure <- function(net_res) {
      top_n <- 5
      g <- net_res$graph
      if (is.null(g) || igraph::vcount(g) == 0) {
        return(list(
          hub_taxa = character(0),
          hub_scores = numeric(0),
          modularity_netanalyze = NA_real_,
          modularity_source = "Unavailable",
          netanalyze_error = NULL
        ))
      }

      extract_modularity <- function(obj) {
        candidates <- numeric(0)
        walk_obj <- function(x, path = "") {
          if (is.null(x)) return(invisible(NULL))
          if (is.numeric(x) && length(x) == 1 && is.finite(x) && grepl("mod", tolower(path))) {
            candidates <<- c(candidates, as.numeric(x))
            return(invisible(NULL))
          }
          if (!is.list(x)) return(invisible(NULL))
          nms <- names(x)
          if (is.null(nms)) nms <- rep("", length(x))
          for (i in seq_along(x)) {
            child_path <- paste(c(path, nms[[i]]), collapse = "/")
            walk_obj(x[[i]], child_path)
          }
        }
        walk_obj(obj)
        if (length(candidates) == 0) NA_real_ else candidates[[1]]
      }

      mod_netanalyze <- NA_real_
      modularity_source <- "Unavailable"
      netanalyze_error <- NULL
      hub_taxa <- character(0)
      hub_scores <- numeric(0)
      if (!is.null(net_res$net_obj)) {
        ana <- tryCatch(
          suppressWarnings(NetCoMi::netAnalyze(net_res$net_obj, verbose = 0)),
          error = function(e) {
            netanalyze_error <<- conditionMessage(e)
            NULL
          }
        )
        mod_netanalyze <- extract_modularity(ana)
        if (is.finite(mod_netanalyze)) modularity_source <- "NetCoMi::netAnalyze"
      }

      if (igraph::ecount(g) > 0) {
        vertex_names <- igraph::V(g)$name
        edge_strength <- abs(as.numeric(igraph::E(g)$weight))
        edge_strength[!is.finite(edge_strength)] <- 0
        eigen_weights <- if (length(edge_strength) > 0 && any(edge_strength > 0)) edge_strength else NULL

        degree_values <- as.numeric(igraph::degree(g))
        betweenness_values <- as.numeric(igraph::betweenness(g, directed = FALSE, weights = NA, normalized = TRUE))
        eigen_values <- tryCatch(
          as.numeric(igraph::eigen_centrality(g, directed = FALSE, weights = eigen_weights, scale = TRUE)$vector),
          error = function(e) rep(0, length(vertex_names))
        )
        safe_z <- function(x) {
          sx <- stats::sd(x, na.rm = TRUE)
          if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
          as.numeric((x - mean(x, na.rm = TRUE)) / sx)
        }
        score <- 0.4 * safe_z(degree_values) + 0.3 * safe_z(betweenness_values) + 0.3 * safe_z(eigen_values)
        names(score) <- vertex_names
        hub_scores <- sort(score, decreasing = TRUE)
        hub_taxa <- head(names(hub_scores), top_n)
        hub_scores <- hub_scores[hub_taxa]

        if (!is.finite(mod_netanalyze)) {
          mod_netanalyze <- tryCatch({
            community <- igraph::cluster_louvain(g, weights = eigen_weights)
            as.numeric(igraph::modularity(community, weights = eigen_weights))
          }, error = function(e) NA_real_)
          if (is.finite(mod_netanalyze)) modularity_source <- "igraph fallback"
        }
      }

      list(
        hub_taxa = hub_taxa,
        hub_scores = hub_scores,
        modularity_netanalyze = mod_netanalyze,
        modularity_source = modularity_source,
        netanalyze_error = netanalyze_error
      )
    }

    network_topology <- reactive({
      res <- display_network_result()
      lapply(res$networks, summarize_network_structure)
    })

    comparison_metrics <- reactive({
      res <- display_network_result()
      if (!identical(res$run_settings$analysis_mode, "Compare two groups")) return(NULL)
      if (length(res$networks) != 2) return(NULL)
      groups <- names(res$networks)
      min_edge_weight <- current_min_edge_weight()
      g1 <- filter_graph_by_weight(res$networks[[groups[1]]]$graph, min_edge_weight)
      g2 <- filter_graph_by_weight(res$networks[[groups[2]]]$graph, min_edge_weight)
      edge1 <- igraph::as_data_frame(g1, what = "edges")
      edge2 <- igraph::as_data_frame(g2, what = "edges")
      key1 <- if (nrow(edge1) > 0) paste(pmin(edge1$from, edge1$to), pmax(edge1$from, edge1$to), sep = "||") else character(0)
      key2 <- if (nrow(edge2) > 0) paste(pmin(edge2$from, edge2$to), pmax(edge2$from, edge2$to), sep = "||") else character(0)
      common <- intersect(key1, key2)
      only1 <- setdiff(key1, key2)
      only2 <- setdiff(key2, key1)
      data.frame(
        metric = c("Nodes", "Edges", "Connected components", "Shared edges", "Unique edges"),
        group_1 = c(igraph::vcount(g1), igraph::ecount(g1), igraph::components(g1)$no, length(common), length(only1)),
        group_2 = c(igraph::vcount(g2), igraph::ecount(g2), igraph::components(g2)$no, length(common), length(only2)),
        delta = c(
          igraph::vcount(g1) - igraph::vcount(g2),
          igraph::ecount(g1) - igraph::ecount(g2),
          igraph::components(g1)$no - igraph::components(g2)$no,
          0,
          length(only1) - length(only2)
        ),
        stringsAsFactors = FALSE
      )
    })

    comparison_edge_table <- reactive({
      res <- display_network_result()
      if (!identical(res$run_settings$analysis_mode, "Compare two groups")) return(data.frame())
      if (length(res$networks) != 2) return(data.frame())
      groups <- names(res$networks)
      g1_name <- groups[1]
      g2_name <- groups[2]
      min_edge_weight <- current_min_edge_weight()
      g1 <- filter_graph_by_weight(res$networks[[g1_name]]$graph, min_edge_weight)
      g2 <- filter_graph_by_weight(res$networks[[g2_name]]$graph, min_edge_weight)
      e1 <- igraph::as_data_frame(g1, what = "edges")
      e2 <- igraph::as_data_frame(g2, what = "edges")
      if (nrow(e1) > 0) {
        e1$key <- paste(pmin(e1$from, e1$to), pmax(e1$from, e1$to), sep = "||")
        e1 <- e1[, c("key", "from", "to", "weight"), drop = FALSE]
        names(e1)[names(e1) == "weight"] <- "weight_group_1"
      } else {
        e1 <- data.frame(key = character(0), from = character(0), to = character(0), weight_group_1 = numeric(0), stringsAsFactors = FALSE)
      }
      if (nrow(e2) > 0) {
        e2$key <- paste(pmin(e2$from, e2$to), pmax(e2$from, e2$to), sep = "||")
        e2 <- e2[, c("key", "from", "to", "weight"), drop = FALSE]
        names(e2)[names(e2) == "weight"] <- "weight_group_2"
      } else {
        e2 <- data.frame(key = character(0), from = character(0), to = character(0), weight_group_2 = numeric(0), stringsAsFactors = FALSE)
      }
      out <- merge(e1, e2, by = "key", all = TRUE, suffixes = c("_1", "_2"), sort = FALSE)
      out$from <- ifelse(is.na(out$from_1), out$from_2, out$from_1)
      out$to <- ifelse(is.na(out$to_1), out$to_2, out$to_1)
      out$from_1 <- NULL
      out$to_1 <- NULL
      out$from_2 <- NULL
      out$to_2 <- NULL
      out$status <- ifelse(
        !is.na(out$weight_group_1) & !is.na(out$weight_group_2),
        "Shared",
        ifelse(!is.na(out$weight_group_1), paste0("Only ", g1_name), paste0("Only ", g2_name))
      )
      out$weight_delta <- as.numeric(out$weight_group_1) - as.numeric(out$weight_group_2)
      out <- out[, c("from", "to", "status", "weight_group_1", "weight_group_2", "weight_delta"), drop = FALSE]
      names(out)[names(out) == "weight_group_1"] <- paste0("weight_", g1_name)
      names(out)[names(out) == "weight_group_2"] <- paste0("weight_", g2_name)
      out
    })

    output$comparison_summary <- renderText({
      res <- display_network_result()
      if (!identical(res$run_settings$analysis_mode, "Compare two groups")) {
        return("Comparison mode is disabled. Select 'Compare two groups' to view results.")
      }
      if (length(res$networks) != 2) {
        return("Comparison requires exactly 2 successfully estimated groups.")
      }
      groups <- names(res$networks)
      metrics <- comparison_metrics()
      topology <- network_topology()
      s1 <- topology[[groups[1]]]
      s2 <- topology[[groups[2]]]
      delta_modularity_netanalyze <- s1$modularity_netanalyze - s2$modularity_netanalyze
      score_text <- function(s) {
        if (length(s$hub_scores) == 0) return("NA")
        vals <- s$hub_scores
        if (length(vals) > 5) vals <- vals[seq_len(5)]
        paste(paste0(names(vals), "=", sprintf("%.2f", as.numeric(vals))), collapse = ", ")
      }
      lines <- c(
        "Group Comparison Summary",
        paste0("Group 1: ", groups[1]),
        paste0("Group 2: ", groups[2]),
        "Note: This summary is based on per-group estimated networks."
      )
      if (!is.null(metrics) && nrow(metrics) > 0) {
        lines <- c(lines, "", "Metrics:")
        for (i in seq_len(nrow(metrics))) {
          lines <- c(
            lines,
            paste0(
              "- ", metrics$metric[i], ": ",
              groups[1], "=", metrics$group_1[i], ", ",
              groups[2], "=", metrics$group_2[i], ", delta=", metrics$delta[i]
            )
          )
        }
      }
      lines <- c(
        lines,
        "",
        "Network Topology (Comparison):",
        "- Hub score method: 0.4*z(degree) + 0.3*z(betweenness) + 0.3*z(eigenvector)",
        paste0("- Hub taxa (", groups[1], "): ", if (length(s1$hub_taxa) > 0) paste(s1$hub_taxa, collapse = ", ") else "NA"),
        paste0("- Hub taxa (", groups[2], "): ", if (length(s2$hub_taxa) > 0) paste(s2$hub_taxa, collapse = ", ") else "NA"),
        paste0("- Hub score (", groups[1], "): ", score_text(s1)),
        paste0("- Hub score (", groups[2], "): ", score_text(s2)),
        paste0("- Modularity (", s1$modularity_source, ", ", groups[1], "): ", if (is.finite(s1$modularity_netanalyze)) sprintf("%.3f", s1$modularity_netanalyze) else "NA"),
        paste0("- Modularity (", s2$modularity_source, ", ", groups[2], "): ", if (is.finite(s2$modularity_netanalyze)) sprintf("%.3f", s2$modularity_netanalyze) else "NA"),
        paste0("- Delta modularity (", groups[1], " - ", groups[2], "): ", if (is.finite(delta_modularity_netanalyze)) sprintf("%.3f", delta_modularity_netanalyze) else "NA")
      )
      paste(lines, collapse = "\n")
    })

    hub_table <- reactive({
      res <- display_network_result()
      topology <- network_topology()
      rows <- list()
      for (gname in names(res$networks)) {
        s <- topology[[gname]]
        taxa <- s$hub_taxa
        if (length(taxa) == 0) next
        scores <- s$hub_scores
        score_vals <- rep(NA_real_, length(taxa))
        if (length(scores) > 0) {
          score_vals <- as.numeric(scores[taxa])
        }
        rows[[length(rows) + 1]] <- data.frame(
          group = gname,
          rank = seq_along(taxa),
          taxa = taxa,
          hub_score = round(score_vals, 2),
          stringsAsFactors = FALSE
        )
      }
      if (length(rows) == 0) {
        return(data.frame(group = character(0), rank = integer(0), taxa = character(0), hub_score = numeric(0)))
      }
      do.call(rbind, rows)
    })

    output$hub_table <- renderDT({
      table <- hub_table()
      DT::datatable(
        table,
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE)
      ) %>%
        DT::formatRound(names(table)[vapply(table, is.numeric, logical(1))], digits = 2)
    })

    draw_comparison_network_plot <- function(edge_tbl, edge_palette, run_settings) {
      validate(
        need(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 package is required for network plotting."),
        need(requireNamespace("igraph", quietly = TRUE), "igraph package is required for network plotting."),
        need(nrow(edge_tbl) > 0, "No differential edge is available for comparison plotting.")
      )
      g <- igraph::graph_from_data_frame(edge_tbl[, c("from", "to"), drop = FALSE], directed = FALSE)
      if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
        validate(need(FALSE, "No edge is available to build the comparison network plot."))
      }
      lay <- igraph::layout_with_fr(g)
      node_df <- data.frame(
        name = igraph::V(g)$name,
        x = lay[, 1],
        y = lay[, 2],
        degree = igraph::degree(g),
        stringsAsFactors = FALSE
      )

      edge_df <- igraph::as_data_frame(g, what = "edges")
      edge_df <- dplyr::left_join(edge_df, edge_tbl, by = c("from", "to"))
      edge_df <- dplyr::left_join(edge_df, node_df[, c("name", "x", "y")], by = c("from" = "name"))
      edge_df <- dplyr::left_join(edge_df, node_df[, c("name", "x", "y")], by = c("to" = "name"), suffix = c("", "_end"))
      names(edge_df)[names(edge_df) == "x"] <- "x"
      names(edge_df)[names(edge_df) == "y"] <- "y"
      names(edge_df)[names(edge_df) == "x_end"] <- "xend"
      names(edge_df)[names(edge_df) == "y_end"] <- "yend"
      edge_df$weight_delta <- as.numeric(edge_df$weight_delta)
      edge_df$weight_delta[!is.finite(edge_df$weight_delta)] <- 0
      edge_df$edge_width <- pmax(abs(edge_df$weight_delta), 0.05)

      ggplot2::ggplot() +
        ggplot2::geom_segment(
          data = edge_df,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = status, linewidth = edge_width),
          alpha = 0.85
        ) +
        ggplot2::geom_point(
          data = node_df,
          ggplot2::aes(x = x, y = y, size = degree),
          shape = 21,
          fill = "#F7F7F7",
          color = "#222222",
          stroke = 0.25
        ) +
        ggrepel::geom_text_repel(
          data = node_df,
          ggplot2::aes(x = x, y = y, label = name),
          size = 2.8,
          max.overlaps = Inf,
          seed = as.integer(run_settings$seed)
        ) +
        ggplot2::scale_color_manual(
          values = edge_palette,
          drop = FALSE,
          name = "Edge class",
          na.translate = FALSE
        ) +
        ggplot2::scale_linewidth(name = "|Weight delta|", range = c(0.3, 2.1)) +
        ggplot2::scale_size(name = "Node degree", range = c(2.5, 7)) +
        ggplot2::theme_void() +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(color = "grey70", fill = NA, linewidth = 0.4),
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          legend.title = ggplot2::element_text(size = 9),
          legend.text = ggplot2::element_text(size = 8)
        ) +
        ggplot2::ggtitle(paste0(run_settings$tax_level, "-Level SPIEC-EASI Differential Network (Two-Group Comparison)"))
    }

    output$comparison_network_plot <- renderPlot({
      if (is.null(input$run_network_btn) || input$run_network_btn < 1) {
        draw_wait_message("Click 'Run SPIEC-EASI' to start analysis.")
        return(invisible(NULL))
      }
      if (isTRUE(network_running())) {
        draw_wait_message("SPIEC-EASI analysis is running. Please wait...")
        return(invisible(NULL))
      }
      res <- display_network_result()
      validate(need(identical(res$run_settings$analysis_mode, "Compare two groups"), "Enable 'Compare two groups' to view differential network."))
      edge_tbl <- comparison_edge_table()
      validate(need(nrow(edge_tbl) > 0, "No differential edge is available for comparison plotting."))
      status_vals <- unique(edge_tbl$status)
      only_status <- status_vals[grepl("^Only ", status_vals)]
      palette <- c("Shared" = "#7F8C8D")
      if (length(only_status) >= 1) palette[only_status[1]] <- "#1F78B4"
      if (length(only_status) >= 2) palette[only_status[2]] <- "#E31A1C"
      p <- draw_comparison_network_plot(edge_tbl, palette, res$run_settings)
      print(p)
    },
    height = function() {
      if (is.null(input$run_network_btn) || input$run_network_btn < 1) {
        200
      } else {
        h <- suppressWarnings(as.numeric(input$plot_height))
        if (!is.finite(h) || is.na(h) || h <= 0) 550 else h
      }
    },
    width = function() {
      w <- suppressWarnings(as.numeric(input$plot_width))
      if (!is.finite(w) || is.na(w) || w <= 0) 800 else w
    })

    output$download_comparison_network_plot <- downloadHandler(
      filename = function() paste0("spieceasi_differential_network_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) {
        res <- display_network_result()
        validate(need(identical(res$run_settings$analysis_mode, "Compare two groups"), "Enable 'Compare two groups' to download differential network."))
        edge_tbl <- comparison_edge_table()
        validate(need(nrow(edge_tbl) > 0, "No differential edge is available for download."))
        status_vals <- unique(edge_tbl$status)
        only_status <- status_vals[grepl("^Only ", status_vals)]
        palette <- c("Shared" = "#7F8C8D")
        if (length(only_status) >= 1) palette[only_status[1]] <- "#1F78B4"
        if (length(only_status) >= 2) palette[only_status[2]] <- "#E31A1C"

        plot_width <- as.integer(input$plot_width)
        plot_height <- as.integer(input$plot_height)
        if (is.na(plot_width) || plot_width <= 0) plot_width <- 800L
        if (is.na(plot_height) || plot_height <= 0) plot_height <- 550L

        dpi_val <- 300
        grDevices::png(
          file,
          width = plot_width / 72,
          height = plot_height / 72,
          units = "in",
          res = dpi_val
        )
        on.exit(grDevices::dev.off(), add = TRUE)
        p <- draw_comparison_network_plot(edge_tbl, palette, res$run_settings)
        print(p)
      }
    )

    output$network_summary <- renderText({
      res <- display_network_result()
      lines <- c(
        "SPIEC-EASI Network Summary",
        paste0("Taxonomic level: ", res$run_settings$tax_level),
        paste0("Analysis mode: ", res$run_settings$analysis_mode),
        paste0("Prevalence cutoff: ", res$run_settings$prevalence_filter_pct, "%"),
        paste0("Maximum taxa: ", res$run_settings$max_taxa),
        paste0("Seed: ", res$run_settings$seed),
        paste0("Samples used: ", res$n_samples),
        paste0("Taxa used: ", res$n_taxa),
        paste0("Group variable: ", res$group_var),
        ""
      )

      for (gname in names(res$networks)) {
        g_all <- res$networks[[gname]]$graph
        g_connected <- get_connected_graph(g_all)
        lines <- c(
          lines,
          paste0("[", gname, "] All"),
          paste0("Nodes: ", igraph::vcount(g_all)),
          paste0("Edges: ", igraph::ecount(g_all)),
          paste0("Connected components: ", if (igraph::vcount(g_all) == 0) 0 else igraph::components(g_all)$no),
          paste0("[", gname, "] Connected"),
          paste0("Nodes: ", igraph::vcount(g_connected)),
          paste0("Edges: ", igraph::ecount(g_connected)),
          paste0("Connected components: ", if (igraph::vcount(g_connected) == 0) 0 else igraph::components(g_connected)$no),
          ""
        )
      }
      paste(lines, collapse = "\n")
    })

    draw_network_plot <- function(res, connected_only = FALSE) {
      validate(
        need(requireNamespace("ggplot2", quietly = TRUE), "ggplot2 package is required for network plotting."),
        need(requireNamespace("ggrepel", quietly = TRUE), "ggrepel package is required for network label repulsion."),
        need(length(res$networks) > 0, "No group network is available for plotting.")
      )

      size_var <- input$node_size_by
      color_var <- input$node_color_by
      if (is.null(color_var)) color_var <- "None"
      min_edge_weight <- suppressWarnings(as.numeric(input$min_edge_weight))
      if (is.na(min_edge_weight) || min_edge_weight < 0) min_edge_weight <- 0.1
      node_rows <- list()
      edge_rows <- list()
      idx <- 1L

      for (gname in names(res$networks)) {
        g <- res$networks[[gname]]$graph
        g <- filter_graph_by_weight(g, min_edge_weight)
        if (connected_only) g <- get_connected_graph(g)
        if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) next

        edge_weight <- igraph::E(g)$weight
        edge_abs <- abs(edge_weight)
        edge_abs[!is.finite(edge_abs)] <- 0

        layout_mat <- igraph::layout_with_fr(g, weights = pmax(edge_abs, .Machine$double.eps)) * 0.65
        node_df <- data.frame(
          name = igraph::V(g)$name,
          x = layout_mat[, 1],
          y = layout_mat[, 2],
          group = gname,
          degree = igraph::degree(g),
          abundance = as.numeric(res$networks[[gname]]$abundance[igraph::V(g)$name]),
          stringsAsFactors = FALSE
        )
        node_df$abundance[is.na(node_df$abundance)] <- 0
        node_df$node_size <- if (identical(size_var, "Abundance")) node_df$abundance else as.numeric(node_df$degree)

        edge_df <- igraph::as_data_frame(g, what = "edges")
        edge_df$edge_sign <- ifelse(edge_df$weight > 0, "Positive", ifelse(edge_df$weight < 0, "Negative", "Zero"))
        edge_df$edge_strength <- abs(edge_df$weight)
        edge_df$group <- gname

        node_rows[[idx]] <- node_df
        edge_rows[[idx]] <- edge_df
        idx <- idx + 1L
      }

      validate(need(length(node_rows) > 0, "No network panels available after filtering."))
      node_plot_df <- do.call(rbind, node_rows)
      edge_plot_df <- do.call(rbind, edge_rows)
      node_plot_df$node_color_group <- "All Nodes"
      if (!identical(color_var, "None") && !is.null(res$taxa_annotation) && color_var %in% colnames(res$taxa_annotation)) {
        color_map <- res$taxa_annotation[, c("node_name", color_var), drop = FALSE]
        color_map[[color_var]] <- as.character(color_map[[color_var]])
        color_map[[color_var]][is.na(color_map[[color_var]]) | !nzchar(color_map[[color_var]])] <- paste0(color_var, "__Unassigned")
        color_map <- color_map[!duplicated(color_map$node_name), , drop = FALSE]
        node_plot_df <- dplyr::left_join(node_plot_df, color_map, by = c("name" = "node_name"))
        node_plot_df$node_color_group <- as.character(node_plot_df[[color_var]])
        node_plot_df$node_color_group[is.na(node_plot_df$node_color_group) | !nzchar(node_plot_df$node_color_group)] <- paste0(color_var, "__Unassigned")
      }
      node_plot_df$group <- factor(node_plot_df$group, levels = names(res$networks))
      edge_plot_df$group <- factor(edge_plot_df$group, levels = names(res$networks))

      edge_plot_df <- merge(edge_plot_df, node_plot_df[, c("name", "x", "y", "group")], by.x = c("from", "group"), by.y = c("name", "group"), all.x = TRUE, sort = FALSE)
      edge_plot_df <- merge(edge_plot_df, node_plot_df[, c("name", "x", "y", "group")], by.x = c("to", "group"), by.y = c("name", "group"), all.x = TRUE, sort = FALSE, suffixes = c("", "_end"))
      names(edge_plot_df)[names(edge_plot_df) == "x"] <- "x"
      names(edge_plot_df)[names(edge_plot_df) == "y"] <- "y"
      names(edge_plot_df)[names(edge_plot_df) == "x_end"] <- "xend"
      names(edge_plot_df)[names(edge_plot_df) == "y_end"] <- "yend"

      p <- ggplot2::ggplot() +
        ggplot2::geom_segment(
          data = edge_plot_df,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = edge_sign, linewidth = edge_strength),
          alpha = 0.7
        ) +
        ggplot2::geom_point(
          data = node_plot_df,
          ggplot2::aes(x = x, y = y, size = node_size, fill = node_color_group),
          shape = 21,
          color = "#222222",
          stroke = 0.2
        ) +
        ggrepel::geom_text_repel(
          data = node_plot_df,
          ggplot2::aes(x = x, y = y, label = name),
          size = 2.8,
          max.overlaps = Inf,
          seed = as.integer(res$run_settings$seed)
        ) +
        ggplot2::facet_wrap(~group, scales = "free") +
        ggplot2::scale_color_manual(values = c("Positive" = "#4575B4", "Negative" = "#D73027", "Zero" = "#999999"), drop = FALSE, name = "Edge sign") +
        ggplot2::scale_linewidth(name = "Edge strength (|w|)", range = c(0.2, 1.6)) +
        ggplot2::scale_size(name = if (identical(size_var, "Abundance")) "Node size: Relative Abundance" else "Node size: Connectivity", range = c(2.4, 7.2)) +
        ggplot2::theme_void() +
        ggplot2::theme(
          strip.text = ggplot2::element_text(face = "bold", size = 10, margin = ggplot2::margin(t = 4, b = 4)),
          strip.background = ggplot2::element_rect(fill = "grey95", color = "grey65", linewidth = 0.4),
          panel.spacing = grid::unit(6, "mm"),
          panel.border = ggplot2::element_rect(color = "grey70", fill = NA, linewidth = 0.4),
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12, margin = ggplot2::margin(b = 10)),
          legend.title = ggplot2::element_text(size = 9),
          legend.text = ggplot2::element_text(size = 8)
        ) +
        ggplot2::ggtitle(
          if (connected_only) {
            paste0("NetCoMi SPIEC-EASI ", res$run_settings$tax_level, "-Level Network (Connected)")
          } else {
            paste0("NetCoMi SPIEC-EASI ", res$run_settings$tax_level, "-Level Network (All)")
          }
        )

      if (identical(color_var, "None")) {
        p <- p + ggplot2::scale_fill_manual(values = c("All Nodes" = "#F1C40F"), guide = "none")
      } else {
        color_levels <- sort(unique(node_plot_df$node_color_group))
        palette <- stats::setNames(grDevices::hcl.colors(length(color_levels), "Dark 3"), color_levels)
        p <- p + ggplot2::scale_fill_manual(
          name = paste0("Node color: ", color_var),
          values = palette,
          drop = FALSE
        )
      }

      print(p)
    }

    edge_table_all <- reactive({
      edge_tbl <- display_network_result()$edge_table
      min_edge_weight <- suppressWarnings(as.numeric(input$min_edge_weight))
      if (is.na(min_edge_weight) || min_edge_weight < 0) min_edge_weight <- 0.1
      if (!"weight" %in% colnames(edge_tbl)) return(edge_tbl)
      edge_tbl <- edge_tbl[is.finite(edge_tbl$weight) & abs(edge_tbl$weight) >= min_edge_weight, , drop = FALSE]
      edge_tbl
    })

    output$network_plot_all <- renderPlot({
      if (is.null(input$run_network_btn) || input$run_network_btn < 1) {
        draw_wait_message("Click 'Run SPIEC-EASI' to start analysis.")
        return(invisible(NULL))
      }
      if (isTRUE(network_running())) {
        draw_wait_message("SPIEC-EASI analysis is running. Please wait...")
        return(invisible(NULL))
      }
      res <- display_network_result()
      draw_network_plot(res, connected_only = FALSE)
    },
    height = function() {
      if (is.null(input$run_network_btn) || input$run_network_btn < 1) {
        200
      } else {
        h <- suppressWarnings(as.numeric(input$plot_height))
        if (!is.finite(h) || is.na(h) || h <= 0) 550 else h
      }
    },
    width = function() {
      w <- suppressWarnings(as.numeric(input$plot_width))
      if (!is.finite(w) || is.na(w) || w <= 0) 800 else w
    })

    output$network_plot_connected <- renderPlot({
      if (is.null(input$run_network_btn) || input$run_network_btn < 1) {
        draw_wait_message("Click 'Run SPIEC-EASI' to start analysis.")
        return(invisible(NULL))
      }
      if (isTRUE(network_running())) {
        draw_wait_message("SPIEC-EASI analysis is running. Please wait...")
        return(invisible(NULL))
      }
      res <- display_network_result()
      draw_network_plot(res, connected_only = TRUE)
    },
    height = function() {
      if (is.null(input$run_network_btn) || input$run_network_btn < 1) {
        200
      } else {
        h <- suppressWarnings(as.numeric(input$plot_height))
        if (!is.finite(h) || is.na(h) || h <= 0) 550 else h
      }
    },
    width = function() {
      w <- suppressWarnings(as.numeric(input$plot_width))
      if (!is.finite(w) || is.na(w) || w <= 0) 800 else w
    })
    output$edge_table_all <- renderDT({
      table <- edge_table_all()
      datatable(table, options = list(scrollX = TRUE)) %>%
        DT::formatRound(names(table)[vapply(table, is.numeric, logical(1))], digits = 2)
    })
    output$comparison_edge_table <- renderDT({
      table <- comparison_edge_table()
      datatable(table, options = list(scrollX = TRUE, pageLength = 10)) %>%
        DT::formatRound(names(table)[vapply(table, is.numeric, logical(1))], digits = 2)
    })

    output$download_network_plot_all <- downloadHandler(
      filename = function() paste0("spieceasi_network_all_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) {
        res <- display_network_result()
        plot_width <- as.integer(input$plot_width)
        plot_height <- as.integer(input$plot_height)
        if (is.na(plot_width) || plot_width <= 0) plot_width <- 800L
        if (is.na(plot_height) || plot_height <= 0) plot_height <- 550L

        dpi_val <- 300
        grDevices::png(
          file,
          width = plot_width / 72,
          height = plot_height / 72,
          units = "in",
          res = dpi_val
        )
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_network_plot(res, connected_only = FALSE)
      }
    )

    output$download_network_plot_connected <- downloadHandler(
      filename = function() paste0("spieceasi_network_connected_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) {
        res <- display_network_result()
        plot_width <- as.integer(input$plot_width)
        plot_height <- as.integer(input$plot_height)
        if (is.na(plot_width) || plot_width <= 0) plot_width <- 800L
        if (is.na(plot_height) || plot_height <= 0) plot_height <- 550L

        dpi_val <- 300
        grDevices::png(
          file,
          width = plot_width / 72,
          height = plot_height / 72,
          units = "in",
          res = dpi_val
        )
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_network_plot(res, connected_only = TRUE)
      }
    )

    output$download_edge_table_all <- downloadHandler(
      filename = function() paste0("spieceasi_edges_all_", Sys.Date(), ".tsv"),
      content = function(file) readr::write_tsv(edge_table_all(), file)
    )

    output$download_comparison_edge_table <- downloadHandler(
      filename = function() paste0("spieceasi_differential_edges_", Sys.Date(), ".tsv"),
      content = function(file) readr::write_tsv(comparison_edge_table(), file)
    )

    output$download_hub_table <- downloadHandler(
      filename = function() paste0("spieceasi_hub_table_", Sys.Date(), ".tsv"),
      content = function(file) readr::write_tsv(hub_table(), file)
    )

    output$spieceasi_legend_box <- renderUI({
      req_plot_legend_tab(session, "SPIEC-EASI")
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 800 else input$plot_width
      tags$div(
        style = sprintf("width: %spx;", box_width),
        tags$div(
          style = "color: #1f2937; font-size: 12.5px; line-height: 1.55;",
          spieceasi_figure_legend()
        )
      )
    })
    outputOptions(output, "spieceasi_legend_box", suspendWhenHidden = FALSE, priority = -10)

    output$spieceasi_status_box <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 800 else input$plot_width
      tags$div(
        class = "app-result-card app-status-card",
        style = paste0("width: ", box_width, "px; max-width: 100%;"),
        verbatimTextOutput(session$ns("comparison_status"))
      )
    })

    output$spieceasi_status_separator <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 800 else input$plot_width
      tags$hr(
        style = paste(
          sprintf("width: %spx;", box_width),
          "max-width: 100%;",
          "margin: 14px 0 12px 0;",
          "border-top: 1px solid #d1d5db;"
        )
      )
    })

    spieceasi_figure_legend <- function() {
      res <- display_network_result()
      run_settings <- res$run_settings
      active_tab <- if (is.null(input$spieceasi_active_tab) || !nzchar(input$spieceasi_active_tab)) "Network Plot (All)" else input$spieceasi_active_tab
      title_text <- paste0("SPIEC-EASI ", active_tab)
      body_text <- if (identical(active_tab, "Network Plot (All)")) {
        "This panel shows all nodes and edges passing the minimum absolute edge-weight threshold. Edge color indicates sign (positive/negative), edge width indicates absolute weight, and node size reflects the selected node-size metric."
      } else if (identical(active_tab, "Network Plot (Connected)")) {
        "Only taxa with at least one retained edge are shown. Edge color indicates the sign of the association, and thicker edges indicate larger absolute weights."
      } else if (identical(active_tab, "Comparison Network")) {
        "This panel visualizes differential edges between two selected groups. It is available only when analysis mode is set to Compare two groups."
      } else if (identical(active_tab, "Edge Table")) {
        "This table lists filtered edges and associated statistics for all estimated group networks."
      } else if (identical(active_tab, "Differential Edges")) {
        "This table reports edge-level differences between two groups in comparison mode."
      } else if (identical(active_tab, "Hub Table")) {
        "This table summarizes hub candidates based on node centrality metrics."
      } else {
        "This panel summarizes network-level topology and run metadata."
      }
      tags$div(
        tags$div(
          style = "font-weight: 600; margin-bottom: 4px;",
          title_text
        ),
        tags$div(
          paste0(
            body_text,
            " Settings: taxonomic level = ",
            run_settings$tax_level,
            ", analysis mode = ",
            run_settings$analysis_mode,
            ", prevalence cutoff = ",
            run_settings$prevalence_filter_pct,
            "%, minimum |edge weight| = ",
            input$min_edge_weight,
            "."
          )
        )
      )
    }

  })
}
