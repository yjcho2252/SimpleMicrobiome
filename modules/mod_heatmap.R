## UI
mod_heatmap_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML("
      .simple-result-card {
        width: 600px;
        max-width: 100%;
        height: 130px;
        overflow-y: auto;
        font-size: 12px;
        line-height: 1.45;
        background: #f8fafc;
        border: 1px solid #d9e2ec;
        border-radius: 10px;
        padding: 10px 12px;
        box-shadow: 0 1px 2px rgba(16, 24, 40, 0.04);
        margin-bottom: 10px;
      }
      .simple-result-card pre {
        margin: 0;
        padding: 0;
        border: 0;
        background: transparent;
        font-size: 12px;
        line-height: 1.45;
        white-space: pre-wrap;
      }
      .well h4 { font-size: 16px; }
      .well h5 { font-size: 13px; }
      .well .control-label { font-size: 12px; }
      .well .checkbox label { font-size: 12px; }
      .well .form-control { font-size: 12px; }
      .well .btn { font-size: 11px; }
      .heatmap-run-btn .heatmap-run-icon {
        display: none;
        margin-right: 6px;
      }
      .heatmap-run-btn .heatmap-run-icon i {
        display: inline-block;
      }
      .heatmap-run-btn.is-running .heatmap-run-icon {
        display: inline-block;
        animation: heatmap-run-spin 0.9s linear infinite;
      }
      .heatmap-run-btn.is-running .heatmap-run-label {
        opacity: 0.85;
      }
      @keyframes heatmap-run-spin {
        from { transform: rotate(0deg); }
        to { transform: rotate(360deg); }
      }
    ")),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        h4(icon("table-cells"), "Association Heatmap"),
        hr(),
        selectInput(
          ns("group_var"),
          "1. Primary grouping variable",
          choices = character(0),
          selected = NULL
        ),
        checkboxInput(ns("use_subgroup"), "Select subgroup", value = FALSE),
        conditionalPanel(
          condition = "input.use_subgroup",
          ns = ns,
          selectInput(
            ns("primary_level"),
            "Primary level to include",
            choices = character(0),
            selected = NULL
          ),
          selectInput(
            ns("secondary_var"),
            "Secondary grouping variable",
            choices = c("None"),
            selected = "None"
          )
        ),
        selectInput(
          ns("group_type"),
          "2. Group type",
          choices = c("Auto" = "auto", "Continuous" = "continuous", "Categorical" = "categorical"),
          selected = "auto"
        ),
        selectInput(
          ns("tax_level"),
          "3. Taxonomic level",
          choices = c("ASV", "Genus", "Species", "Strain"),
          selected = "Genus"
        ),
        selectInput(
          ns("association_mode"),
          "4. Association mode",
          choices = c("Abundance-based" = "abundance", "Prevalence-based (binary)" = "prevalence"),
          selected = "abundance"
        ),
        selectInput(
          ns("value_scale"),
          "5. Heatmap value scale",
          choices = c("Raw" = "raw", "Z-score (by taxa)" = "zscore"),
          selected = "raw"
        ),
        tags$details(
          tags$summary("Advanced options"),
          selectInput(
            ns("transform_method"),
            "Data transform",
            choices = c("TSS" = "tss", "CLR" = "clr"),
            selected = "clr"
          ),
          selectInput(
            ns("corr_method"),
            "Correlation method",
            choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
            selected = "pearson"
          ),
          numericInput(
            ns("prevalence_filter_pct"),
            "Prevalence filter (%)",
            value = 10,
            min = 0,
            max = 100,
            step = 1
          ),
          numericInput(
            ns("max_taxa"),
            "Max taxa",
            value = 30,
            min = 10,
            max = 300,
            step = 5
          ),
          selectizeInput(
            ns("selected_taxa"),
            "Selected taxa (optional)",
            choices = NULL,
            selected = NULL,
            multiple = TRUE,
            options = list(placeholder = "If selected, only these taxa are shown")
          ),
          checkboxInput(ns("show_row_dendrogram"), "Show row dendrogram (cluster row)", value = TRUE),
          checkboxInput(ns("show_col_dendrogram"), "Show column dendrogram (cluster column)", value = TRUE),
          selectInput(
            ns("color_palette"),
            "Color palette",
            choices = c(
              "Blue - White - Red" = "blue_red",
              "Blue - White - Orange (colorblind friendly)" = "blue_orange",
              "Purple - White - Green" = "purple_green",
              "Teal - White - Brown" = "teal_brown"
            ),
            selected = "blue_red"
          ),
          numericInput(ns("abs_limit"), "Color scale abs limit (0 = auto)", value = 0, min = 0, max = 100, step = 0.1),
          checkboxInput(ns("apply_sig_filter"), "Mask non-significant cells (FDR)", value = FALSE),
          numericInput(ns("q_cutoff"), "FDR cutoff (q)", value = 0.05, min = 0.0001, max = 1, step = 0.01)
        ),
        hr(),
        h4(icon("sliders"), "Display Options"),
        selectInput(ns("plot_size_mode"), "Plot Size", choices = c("Auto" = "auto", "Manual" = "manual")),
        conditionalPanel(
          "input.plot_size_mode === 'manual'", ns = ns,
          numericInput(ns("plot_width"), "Plot Width (px)", value = 900, min = 400, max = 2400, step = 50),
          numericInput(ns("plot_height"), "Plot Height (px)", value = 550, min = 300, max = 2400, step = 50)
        ),
        numericInput(ns("base_size"), "Base Font Size:", value = 11, min = 6, max = 30, step = 1),
        tags$button(
          id = ns("run_heatmap"),
          type = "button",
          class = "action-button btn btn-danger heatmap-run-btn",
          `data-val` = 0,
          style = "font-size: 12px;",
          tags$span(class = "heatmap-run-icon", icon("spinner", class = "fa-spin")),
          tags$span(class = "heatmap-run-label", "Run Heatmap")
        ),
        hr(),
        h5(icon("download"), "Download"),
        tags$div(
          style = "display: flex; gap: 4px; align-items: center; flex-wrap: nowrap;",
          downloadButton(
            ns("download_heatmap_png"),
            "Download Plot (PNG)",
            style = "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
          ),
          downloadButton(
            ns("download_corr_tsv"),
            "Download Matrix (TSV)",
            style = "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
          )
        ),
        tags$script(HTML(
          "Shiny.addCustomMessageHandler('toggle-heatmap-run-btn', function(msg) {
             var btn = document.getElementById(msg.id);
             if (!btn) return;
             btn.disabled = !!msg.disabled;
             btn.classList.toggle('is-running', !!msg.running);
             var label = btn.querySelector('.heatmap-run-label');
             if (label && msg.label) label.textContent = msg.label;
            });
           Shiny.addCustomMessageHandler('toggle-input-disabled', function(msg) {
             var el = document.getElementById(msg.id);
             if (!el) return;
             el.disabled = !!msg.disabled;
             if (el.selectize) {
               if (msg.disabled) {
                 el.selectize.disable();
               } else {
                 el.selectize.enable();
               }
             }
           });"
        ))
      ),
      mainPanel(
        h4("Association Heatmap"),
        plotOutput(ns("corr_heatmap_plot"), height = "auto"),
        uiOutput(ns("heatmap_legend_box")),
        uiOutput(ns("heatmap_status_separator")),
        h5(icon("circle-info"), "Heatmap Status"),
        uiOutput(ns("heatmap_status_box"))
      )
    )
  )
}

## Server
mod_heatmap_server <- function(id, ps_obj, meta_vars = NULL, plot_order = function() NULL) {
  moduleServer(id, function(input, output, session) {
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

    analysis_target <- reactive("taxa_group")
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

    observe({
      disable_prevalence_irrelevant <- identical(input$association_mode, "prevalence")
      session$sendCustomMessage(
        "toggle-input-disabled",
        list(id = session$ns("transform_method"), disabled = disable_prevalence_irrelevant)
      )
      session$sendCustomMessage(
        "toggle-input-disabled",
        list(id = session$ns("corr_method"), disabled = disable_prevalence_irrelevant)
      )
      session$sendCustomMessage(
        "toggle-input-disabled",
        list(id = session$ns("value_scale"), disabled = disable_prevalence_irrelevant)
      )
    })

    observeEvent(input$association_mode, {
      if (identical(input$association_mode, "prevalence")) {
        updateSelectInput(session, "value_scale", selected = "raw")
      }
    })

    sanitize_taxa_names <- function(ps_data, rank_name) {
      tt <- phyloseq::tax_table(ps_data, errorIfNULL = FALSE)
      if (is.null(tt) || !rank_name %in% colnames(tt)) {
        return(ps_data)
      }
      tax_df <- as.data.frame(tt, stringsAsFactors = FALSE)
      rank_values <- as.character(tax_df[[rank_name]])
      rank_norm <- tolower(trimws(rank_values))
      rank_norm <- gsub("^[a-z]__", "", rank_norm)
      idx_placeholder <- !is.na(rank_norm) & grepl("(uncultured|unassigned)", rank_norm)

      if (any(idx_placeholder)) {
        tax_ranks <- phyloseq::rank_names(ps_data)
        rank_pos <- match(rank_name, tax_ranks)
        parent_candidates <- character(0)
        if (!is.na(rank_pos) && rank_pos > 1) {
          parent_candidates <- rev(tax_ranks[seq_len(rank_pos - 1)])
          parent_candidates <- parent_candidates[parent_candidates %in% colnames(tax_df)]
        }
        parent_val <- rep("UnclassifiedParent", nrow(tax_df))
        if (length(parent_candidates) > 0) {
          for (parent_rank in parent_candidates) {
            candidate_val <- as.character(tax_df[[parent_rank]])
            candidate_norm <- tolower(trimws(candidate_val))
            candidate_norm <- gsub("^[a-z]__", "", candidate_norm)
            candidate_is_placeholder <- is.na(candidate_norm) | !nzchar(candidate_norm) | grepl("(uncultured|unassigned)", candidate_norm)
            use_idx <- idx_placeholder & parent_val == "UnclassifiedParent" & !candidate_is_placeholder
            parent_val[use_idx] <- candidate_val[use_idx]
          }
        }
        rank_values[idx_placeholder] <- parent_val[idx_placeholder]
      }

      missing_idx <- is.na(rank_values) | !nzchar(rank_values)
      rank_values[missing_idx] <- phyloseq::taxa_names(ps_data)[missing_idx]
      phyloseq::taxa_names(ps_data) <- make.unique(rank_values)
      ps_data
    }

    metadata_df <- reactive({
      req(ps_obj())
      sd <- phyloseq::sample_data(ps_obj(), errorIfNULL = FALSE)
      if (is.null(sd)) {
        return(data.frame())
      }
      as.data.frame(sd, stringsAsFactors = FALSE)
    })

    secondary_var_resolved <- reactive({
      meta_df <- metadata_df()
      if (!isTRUE(input$use_subgroup)) {
        return(NULL)
      }
      if (is.null(input$secondary_var) || identical(input$secondary_var, "None")) {
        return(NULL)
      }
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

    observe({
      meta_df <- metadata_df()
      candidate_cols <- setdiff(colnames(meta_df), "SampleID")
      current_selected <- isolate(input$group_var)
      selected_val <- if (!is.null(current_selected) && current_selected %in% candidate_cols) {
        current_selected
      } else if (length(candidate_cols) > 0) {
        candidate_cols[1]
      } else {
        character(0)
      }
      updateSelectInput(session, "group_var", choices = candidate_cols, selected = selected_val)

      level_choices <- character(0)
      if (!is.null(selected_val) && nzchar(selected_val) && selected_val %in% colnames(meta_df)) {
        level_vec <- as.character(meta_df[[selected_val]])
        level_vec <- level_vec[!is.na(level_vec) & nzchar(level_vec)]
        if (length(level_vec) > 0) {
          level_counts <- table(level_vec)
          if (length(level_counts) < 10) {
            level_choices <- sort(names(level_counts))
          }
        }
      }
      level_choices <- c("All", level_choices)
      selected_primary_level <- isolate(input$primary_level)
      if (is.null(selected_primary_level) || !selected_primary_level %in% level_choices) {
        selected_primary_level <- "All"
      }
      updateSelectInput(session, "primary_level", choices = level_choices, selected = selected_primary_level)

      secondary_choices <- if (!is.null(selected_val) && nzchar(selected_val)) setdiff(candidate_cols, selected_val) else candidate_cols
      selected_secondary <- isolate(input$secondary_var)
      if (!isTRUE(input$use_subgroup)) {
        selected_secondary <- "None"
      } else if (is.null(selected_secondary) || !selected_secondary %in% c("None", secondary_choices)) {
        selected_secondary <- "None"
      }
      updateSelectInput(session, "secondary_var", choices = c("None" = "None", secondary_choices), selected = selected_secondary)
    })

    observeEvent(list(metadata_df(), input$group_var), {
      meta_df <- metadata_df()
      req(nrow(meta_df) >= 0)
      group_col <- input$group_var
      if (is.null(group_col) || !nzchar(group_col)) {
        updateSelectInput(session, "primary_level", choices = character(0), selected = character(0))
        return()
      }
      resolved_group_col <- resolve_meta_colname(group_col, colnames(meta_df))
      if (is.null(resolved_group_col) || !resolved_group_col %in% colnames(meta_df)) {
        updateSelectInput(session, "primary_level", choices = character(0), selected = character(0))
        return()
      }
      level_vec <- as.character(meta_df[[resolved_group_col]])
      level_vec <- level_vec[!is.na(level_vec) & nzchar(level_vec)]
      if (length(level_vec) > 0) {
        level_counts <- table(level_vec)
        if (length(level_counts) < 10) {
          level_choices <- sort(names(level_counts))
        } else {
          level_choices <- character(0)
        }
      } else {
        level_choices <- character(0)
      }
      level_choices <- c("All", level_choices)
      selected_primary_level <- isolate(input$primary_level)
      if (is.null(selected_primary_level) || !selected_primary_level %in% level_choices) {
        selected_primary_level <- "All"
      }
      updateSelectInput(session, "primary_level", choices = level_choices, selected = selected_primary_level)
    }, ignoreInit = FALSE)

    heatmap_rank_matrix <- reactive({
      req(ps_obj(), input$tax_level)
      ps_data <- ps_obj()
      validate(
        need(phyloseq::ntaxa(ps_data) >= 1, "No taxa available.")
      )

      if (!identical(input$tax_level, "ASV")) {
        tt <- phyloseq::tax_table(ps_data, errorIfNULL = FALSE)
        validate(
          need(!is.null(tt), "Taxonomy table is required for taxonomic aggregation.")
        )
        tax_cols <- colnames(tt)
        validate(
          need(input$tax_level %in% tax_cols, paste0("Taxonomic rank '", input$tax_level, "' is not available."))
        )
        ps_data <- phyloseq::tax_glom(ps_data, taxrank = input$tax_level, NArm = FALSE)
        ps_data <- sanitize_taxa_names(ps_data, input$tax_level)
      }

      otu_mat <- as.matrix(phyloseq::otu_table(ps_data))
      if (!phyloseq::taxa_are_rows(phyloseq::otu_table(ps_data))) {
        otu_mat <- t(otu_mat)
      }
      otu_mat
    })

    taxa_choices_preview <- reactive({
      otu_mat <- heatmap_rank_matrix()
      prevalence_cutoff <- as.numeric(input$prevalence_filter_pct)
      if (is.na(prevalence_cutoff)) prevalence_cutoff <- 10
      prevalence_cutoff <- max(0, min(100, prevalence_cutoff))
      keep_prev <- rowMeans(otu_mat > 0, na.rm = TRUE) * 100 >= prevalence_cutoff
      otu_mat <- otu_mat[keep_prev, , drop = FALSE]

      max_taxa <- as.integer(input$max_taxa)
      if (is.na(max_taxa) || max_taxa < 3) max_taxa <- 30
      if (nrow(otu_mat) > max_taxa) {
        otu_mat <- otu_mat[order(rowSums(otu_mat, na.rm = TRUE), decreasing = TRUE)[seq_len(max_taxa)], , drop = FALSE]
      }

      rownames(otu_mat)
    })

    observe({
      taxa_choices <- taxa_choices_preview()
      current_selected <- isolate(input$selected_taxa)
      selected_keep <- if (length(current_selected) > 0) intersect(current_selected, taxa_choices) else character(0)
      updateSelectizeInput(session, "selected_taxa", choices = taxa_choices, selected = selected_keep, server = TRUE)
    })

    heatmap_running <- reactiveVal(FALSE)
    heatmap_run_nonce <- reactiveVal(NULL)
    heatmap_next_allowed_at <- reactiveVal(as.POSIXct(NA))

    phyloseq_signature <- function(ps) {
      otu <- as(phyloseq::otu_table(ps), "matrix")
      if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
      values <- as.numeric(otu)
      metadata <- as.data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE)
      metadata[] <- lapply(metadata, as.character)
      tax_obj <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
      taxonomy <- if (is.null(tax_obj)) NULL else as.data.frame(tax_obj, stringsAsFactors = FALSE)
      if (!is.null(taxonomy)) taxonomy[] <- lapply(taxonomy, as.character)
      list(
        samples = sort(phyloseq::sample_names(ps)),
        taxa = sort(phyloseq::taxa_names(ps)),
        ranks = if (is.null(tax_obj)) character(0) else colnames(tax_obj),
        abundance = c(
          total = sum(values, na.rm = TRUE),
          squares = sum(values^2, na.rm = TRUE),
          weighted = sum(values * seq_along(values), na.rm = TRUE)
        ),
        metadata = metadata,
        taxonomy = taxonomy
      )
    }

    current_analysis_settings <- reactive({
      req(ps_obj())
      list(
        data = phyloseq_signature(ps_obj()),
        group_var = input$group_var,
        use_subgroup = isTRUE(input$use_subgroup),
        primary_level = input$primary_level,
        secondary_var = input$secondary_var,
        group_type = input$group_type,
        tax_level = input$tax_level,
        transform_method = input$transform_method,
        association_mode = input$association_mode,
        corr_method = input$corr_method,
        value_scale = input$value_scale,
        prevalence_filter_pct = as.numeric(input$prevalence_filter_pct),
        max_taxa = as.integer(input$max_taxa),
        apply_sig_filter = isTRUE(input$apply_sig_filter),
        q_cutoff = as.numeric(input$q_cutoff),
        show_row_dendrogram = isTRUE(input$show_row_dendrogram),
        show_col_dendrogram = isTRUE(input$show_col_dendrogram)
      )
    })

    observeEvent(input$run_heatmap, {
      ps_data <- tryCatch(ps_obj(), shiny.silent.error = function(e) NULL)
      if (is.null(ps_data)) {
        showNotification("Load data and check the preprocessing settings before running the heatmap.", type = "warning")
        return(invisible(NULL))
      }
      if (identical(analysis_target(), "taxa_group") && !isTruthy(input$group_var)) {
        showNotification("Select a primary group variable before running the heatmap.", type = "warning")
        return(invisible(NULL))
      }
      now_ts <- Sys.time()
      next_allowed <- heatmap_next_allowed_at()
      if (isTRUE(heatmap_running())) {
        showNotification("Heatmap is already running. Please wait for completion.", type = "message", duration = 3)
        return(invisible(NULL))
      }
      if (!is.na(next_allowed) && now_ts < next_allowed) {
        wait_sec <- ceiling(as.numeric(difftime(next_allowed, now_ts, units = "secs")))
        showNotification(paste0("Please wait ", wait_sec, " second(s) before running again."), type = "message", duration = 2)
        return(invisible(NULL))
      }
      heatmap_next_allowed_at(now_ts + 10)
      heatmap_run_nonce(input$run_heatmap)
    }, ignoreInit = TRUE)

    corr_payload <- eventReactive(heatmap_run_nonce(), {
      heatmap_running(TRUE)
      session$sendCustomMessage(
        "toggle-heatmap-run-btn",
        list(id = session$ns("run_heatmap"), disabled = TRUE, running = TRUE, label = "Running...")
      )
      on.exit(
        {
          heatmap_running(FALSE)
          session$sendCustomMessage(
            "toggle-heatmap-run-btn",
            list(id = session$ns("run_heatmap"), disabled = FALSE, running = FALSE, label = "Run Heatmap")
          )
        },
        add = TRUE
      )
      validate(need(
        all(vapply(list(input$tax_level, input$transform_method, input$association_mode,
                        input$corr_method, input$value_scale), isTruthy, logical(1))),
        "Heatmap settings are incomplete. Check the analysis options and run again."
      ))
      ps_data <- ps_obj()
      validate(
        need(phyloseq::nsamples(ps_data) >= 3, "At least 3 samples are required."),
        need(phyloseq::ntaxa(ps_data) >= 3, "At least 3 taxa are required.")
      )

      otu_mat <- heatmap_rank_matrix()

      prevalence_cutoff <- as.numeric(input$prevalence_filter_pct)
      if (is.na(prevalence_cutoff)) prevalence_cutoff <- 10
      prevalence_cutoff <- max(0, min(100, prevalence_cutoff))
      keep_prev <- rowMeans(otu_mat > 0, na.rm = TRUE) * 100 >= prevalence_cutoff
      otu_mat <- otu_mat[keep_prev, , drop = FALSE]
      validate(
        need(nrow(otu_mat) >= 3, "Too few taxa remain after prevalence filtering.")
      )

      max_taxa <- as.integer(input$max_taxa)
      if (is.na(max_taxa) || max_taxa < 3) max_taxa <- 30
      if (nrow(otu_mat) > max_taxa) {
        otu_mat <- otu_mat[order(rowSums(otu_mat, na.rm = TRUE), decreasing = TRUE)[seq_len(max_taxa)], , drop = FALSE]
      }

      transformed <- if (identical(input$transform_method, "clr")) {
        log_mat <- log(otu_mat + 1)
        out <- apply(log_mat, 2, function(col_x) col_x - mean(col_x, na.rm = TRUE))
        if (!is.matrix(out)) out <- matrix(out, nrow = nrow(log_mat))
        rownames(out) <- rownames(log_mat)
        colnames(out) <- colnames(log_mat)
        out
      } else {
        col_totals <- colSums(otu_mat, na.rm = TRUE)
        col_totals[col_totals <= 0 | is.na(col_totals)] <- 1
        sweep(otu_mat, 2, col_totals, "/")
      }

      analysis_mat <- if (identical(input$association_mode, "prevalence")) {
        (otu_mat > 0) * 1
      } else {
        transformed
      }

      taxa_var <- apply(analysis_mat, 1, stats::var, na.rm = TRUE)
      keep_var <- is.finite(taxa_var) & taxa_var > 0
      analysis_mat <- analysis_mat[keep_var, , drop = FALSE]
      validate(
        need(nrow(analysis_mat) >= 3, "Too few taxa remain after variance filtering.")
      )

      meta_df <- metadata_df()
      sample_ids <- colnames(analysis_mat)
      if (nrow(meta_df) > 0) {
        meta_df <- meta_df[sample_ids, , drop = FALSE]
      }

      p_mat <- NULL
      q_mat <- NULL
      sig_supported <- FALSE
      sig_note <- "Not computed."
      column_group <- NULL
      column_panels <- NULL
      row_hclust <- NULL
      col_hclust <- NULL

      if (identical(analysis_target(), "taxa_taxa")) {
        corr_method_eff <- if (identical(input$association_mode, "prevalence")) "pearson" else input$corr_method
        assoc_mat <- stats::cor(t(analysis_mat), method = corr_method_eff, use = "pairwise.complete.obs")
        assoc_type <- if (identical(input$association_mode, "prevalence")) {
          "Taxa vs Taxa prevalence association (phi)"
        } else {
          "Taxa vs Taxa correlation"
        }
        sig_supported <- TRUE
        sig_note <- if (identical(input$association_mode, "prevalence")) {
          "P-values come from pairwise phi (binary Pearson) tests."
        } else {
          "P-values come from pairwise correlation tests."
        }

        n_tax <- nrow(analysis_mat)
        p_mat <- matrix(NA_real_, nrow = n_tax, ncol = n_tax, dimnames = list(rownames(analysis_mat), rownames(analysis_mat)))
        diag(p_mat) <- 0
        if (n_tax >= 2) {
          for (i in seq_len(n_tax - 1)) {
            xi <- analysis_mat[i, ]
            for (j in (i + 1):n_tax) {
              xj <- analysis_mat[j, ]
              p_val <- tryCatch(
                suppressWarnings(stats::cor.test(xi, xj, method = corr_method_eff, exact = FALSE)$p.value),
                error = function(e) NA_real_
              )
              p_mat[i, j] <- p_val
              p_mat[j, i] <- p_val
            }
          }
          upper_idx <- upper.tri(p_mat)
          p_vals <- p_mat[upper_idx]
          q_vals <- rep(NA_real_, length(p_vals))
          ok <- is.finite(p_vals)
          if (any(ok)) {
            q_vals[ok] <- stats::p.adjust(p_vals[ok], method = "BH")
          }
          q_mat <- matrix(NA_real_, nrow = n_tax, ncol = n_tax, dimnames = dimnames(p_mat))
          q_mat[upper_idx] <- q_vals
          q_mat[lower.tri(q_mat)] <- t(q_mat)[lower.tri(q_mat)]
          diag(q_mat) <- 0
        }
      } else {
        req(input$group_var)
        validate(
          need(input$group_var %in% colnames(meta_df), "Selected group variable is not available in metadata.")
        )
        current_secondary_var <- secondary_var_resolved()
        effective_group_var <- if (!is.null(current_secondary_var)) current_secondary_var else input$group_var
        validate(
          need(effective_group_var %in% colnames(meta_df), "Selected secondary/primary grouping variable is not available in metadata.")
        )
        group_vec <- meta_df[[effective_group_var]]
        names(group_vec) <- rownames(meta_df)
        group_vec <- group_vec[sample_ids]
        primary_vec <- NULL
        secondary_levels <- NULL
        if (!is.null(current_secondary_var)) {
          primary_vec <- meta_df[[input$group_var]]
          names(primary_vec) <- rownames(meta_df)
          primary_vec <- primary_vec[sample_ids]
          secondary_levels <- unique(as.character(group_vec[!is.na(primary_vec)]))
          secondary_levels[is.na(secondary_levels) | !nzchar(secondary_levels)] <- "(Missing)"
          secondary_levels <- unique(secondary_levels)
          if (!isTruthy(input$primary_level) || identical(as.character(input$primary_level), "All")) {
            keep_idx <- !is.na(primary_vec)
          } else {
            keep_idx <- !is.na(primary_vec) & (as.character(primary_vec) == as.character(input$primary_level))
          }
          group_vec <- group_vec[keep_idx]
          analysis_mat <- analysis_mat[, keep_idx, drop = FALSE]
          sample_ids <- sample_ids[keep_idx]
          primary_vec <- as.character(primary_vec[keep_idx])
        }

        group_type_selected <- input$group_type
        if (is.null(group_type_selected) || !group_type_selected %in% c("auto", "continuous", "categorical")) {
          group_type_selected <- "auto"
        }
        group_num_auto <- suppressWarnings(as.numeric(group_vec))
        auto_numeric_like <- sum(is.finite(group_num_auto)) >= 3
        n_unique_numeric <- length(unique(group_num_auto[is.finite(group_num_auto)]))
        force_categorical_numeric <- n_unique_numeric > 0 && n_unique_numeric <= 5
        is_numeric_group <- if (identical(group_type_selected, "continuous")) {
          TRUE
        } else if (identical(group_type_selected, "categorical")) {
          FALSE
        } else {
          (is.numeric(group_vec) || is.integer(group_vec) || auto_numeric_like)
        }
        if (isTRUE(force_categorical_numeric)) {
          is_numeric_group <- FALSE
        }

        if (is_numeric_group) {
          group_num <- suppressWarnings(as.numeric(group_vec))
          validate(
            need(sum(is.finite(group_num)) >= 3, "Continuous group requires at least 3 numeric values.")
          )
          corr_method_eff <- if (identical(input$association_mode, "prevalence")) "pearson" else input$corr_method
          assoc_vals <- apply(analysis_mat, 1, function(taxa_x) {
            suppressWarnings(stats::cor(taxa_x, group_num, method = corr_method_eff, use = "pairwise.complete.obs"))
          })
          assoc_mat <- matrix(as.numeric(assoc_vals), ncol = 1)
          rownames(assoc_mat) <- names(assoc_vals)
          colnames(assoc_mat) <- effective_group_var
          assoc_type <- if (identical(input$association_mode, "prevalence")) {
            "Taxa prevalence vs continuous group association (point-biserial)"
          } else {
            "Taxa vs continuous group correlation"
          }
          sig_supported <- TRUE
          sig_note <- if (identical(input$association_mode, "prevalence")) {
            "P-values come from point-biserial correlation tests of taxon presence against the continuous variable."
          } else {
            "P-values come from correlation tests of taxon abundance against the continuous variable."
          }

          p_vals <- apply(analysis_mat, 1, function(taxa_x) {
            tryCatch(
              suppressWarnings(stats::cor.test(taxa_x, group_num, method = corr_method_eff, exact = FALSE)$p.value),
              error = function(e) NA_real_
            )
          })
          q_vals <- rep(NA_real_, length(p_vals))
          ok <- is.finite(p_vals)
          if (any(ok)) {
            q_vals[ok] <- stats::p.adjust(p_vals[ok], method = "BH")
          }
          p_mat <- matrix(as.numeric(p_vals), ncol = 1, dimnames = list(names(p_vals), effective_group_var))
          q_mat <- matrix(as.numeric(q_vals), ncol = 1, dimnames = list(names(q_vals), effective_group_var))
        } else {
          group_chr <- as.character(group_vec)
          group_chr[is.na(group_chr) | !nzchar(group_chr)] <- "(Missing)"
          levels_keep <- if (is.null(primary_vec)) unique(group_chr) else secondary_levels
          column_group <- effective_group_var
          validate(
            need(length(levels_keep) >= 2, "Categorical group needs at least 2 levels.")
          )
          if (is.null(primary_vec)) {
            column_levels <- levels_keep
            column_ids <- column_levels
          } else {
            primary_levels <- unique(primary_vec)
            validate(need(length(primary_levels) > 0L, "No samples are available for the selected primary level."))
            column_panels <- expand.grid(
              secondary = levels_keep, primary = primary_levels,
              stringsAsFactors = FALSE
            )
            column_levels <- column_panels$secondary
            column_ids <- make.unique(paste(column_panels$primary, column_levels, sep = " / "))
            rownames(column_panels) <- column_ids
          }
          assoc_mat <- matrix(
            NA_real_, nrow = nrow(analysis_mat), ncol = length(column_ids),
            dimnames = list(rownames(analysis_mat), column_ids)
          )
          assoc_type <- if (identical(input$association_mode, "prevalence")) {
            "Taxa prevalence difference by categorical group level (in-group minus out-group)"
          } else {
            "Taxa vs categorical group level mean"
          }
          sig_supported <- TRUE
          sig_note <- if (identical(input$association_mode, "prevalence")) {
            "For each taxon, Fisher's exact test compares presence in the indicated group with presence in all other groups combined."
          } else {
            "For each taxon, a Wilcoxon rank-sum test compares the indicated group with all other groups combined."
          }

          if (!is.null(column_panels)) {
            sig_note <- paste(sig_note, "Comparisons stay within each primary group, and BH adjustment covers all valid tests across the displayed panels.")
          }
          p_mat <- matrix(
            NA_real_,
            nrow = nrow(assoc_mat),
            ncol = ncol(assoc_mat),
            dimnames = list(rownames(assoc_mat), colnames(assoc_mat))
          )
          for (j in seq_along(column_levels)) {
            in_panel <- if (is.null(column_panels)) rep(TRUE, length(group_chr)) else
              primary_vec == column_panels$primary[j]
            in_group <- which(in_panel & group_chr == column_levels[j])
            out_group <- which(in_panel & group_chr != column_levels[j])
            if (length(in_group)) {
              in_mean <- rowMeans(analysis_mat[, in_group, drop = FALSE], na.rm = TRUE)
              if (!identical(input$association_mode, "prevalence")) {
                assoc_mat[, j] <- in_mean
              } else if (length(out_group)) {
                assoc_mat[, j] <- in_mean - rowMeans(analysis_mat[, out_group, drop = FALSE], na.rm = TRUE)
              }
            }
            if (length(in_group) < 2 || length(out_group) < 2) {
              next
            }
            p_vec <- apply(analysis_mat, 1, function(taxa_x) {
              if (identical(input$association_mode, "prevalence")) {
                in_present <- sum(taxa_x[in_group] > 0, na.rm = TRUE)
                in_absent <- sum(taxa_x[in_group] <= 0, na.rm = TRUE)
                out_present <- sum(taxa_x[out_group] > 0, na.rm = TRUE)
                out_absent <- sum(taxa_x[out_group] <= 0, na.rm = TRUE)
                contingency <- matrix(c(in_present, in_absent, out_present, out_absent), nrow = 2, byrow = TRUE)
                if (any(contingency < 0) || any(rowSums(contingency) == 0) || any(colSums(contingency) == 0)) {
                  return(NA_real_)
                }
                tryCatch(
                  stats::fisher.test(contingency)$p.value,
                  error = function(e) NA_real_
                )
              } else {
                tryCatch(
                  stats::wilcox.test(
                    x = as.numeric(taxa_x[in_group]),
                    y = as.numeric(taxa_x[out_group]),
                    exact = FALSE
                  )$p.value,
                  error = function(e) NA_real_
                )
              }
            })
            p_mat[, j] <- as.numeric(p_vec)
          }

          q_mat <- matrix(
            NA_real_,
            nrow = nrow(p_mat),
            ncol = ncol(p_mat),
            dimnames = dimnames(p_mat)
          )
          p_vals_all <- as.numeric(p_mat)
          ok <- is.finite(p_vals_all)
          if (any(ok)) {
            q_vals_all <- rep(NA_real_, length(p_vals_all))
            q_vals_all[ok] <- stats::p.adjust(p_vals_all[ok], method = "BH")
            q_mat[] <- matrix(q_vals_all, nrow = nrow(p_mat), ncol = ncol(p_mat))
          }
        }
      }

      assoc_mat[!is.finite(assoc_mat)] <- if (is.null(column_panels)) 0 else NA_real_
      raw_assoc_mat <- assoc_mat

      if (identical(analysis_target(), "taxa_taxa")) {
        assoc_mat <- pmax(-1, pmin(1, assoc_mat))
        if (is.null(dim(assoc_mat))) {
          assoc_mat <- matrix(
            assoc_mat,
            nrow = nrow(raw_assoc_mat),
            ncol = ncol(raw_assoc_mat),
            dimnames = dimnames(raw_assoc_mat)
          )
        }
        diag(assoc_mat) <- 1
      }

      if (identical(input$value_scale, "zscore") && !identical(input$association_mode, "prevalence")) {
        z_by_row <- t(apply(assoc_mat, 1, function(row_x) {
          row_sd <- stats::sd(row_x, na.rm = TRUE)
          if (!is.finite(row_sd) || row_sd <= 0) {
            return(ifelse(is.finite(row_x), 0, NA_real_))
          }
          (row_x - mean(row_x, na.rm = TRUE)) / row_sd
        }))
        if (is.matrix(z_by_row)) {
          rownames(z_by_row) <- rownames(assoc_mat)
          colnames(z_by_row) <- colnames(assoc_mat)
          assoc_mat <- z_by_row
        }
      }

      if (isTRUE(input$show_row_dendrogram) && nrow(assoc_mat) >= 3) {
        row_dist <- stats::dist(assoc_mat, method = "euclidean")
        row_hclust <- tryCatch(
          stats::hclust(row_dist, method = "average"),
          error = function(e) NULL
        )
        row_ord <- if (!is.null(row_hclust)) row_hclust$order else seq_len(nrow(assoc_mat))
        assoc_mat <- assoc_mat[row_ord, , drop = FALSE]
      }

      if (is.null(column_panels) && isTRUE(input$show_col_dendrogram) && ncol(assoc_mat) >= 3) {
        col_dist <- stats::dist(t(assoc_mat), method = "euclidean")
        col_hclust <- tryCatch(
          stats::hclust(col_dist, method = "average"),
          error = function(e) NULL
        )
        col_ord <- if (!is.null(col_hclust)) col_hclust$order else seq_len(ncol(assoc_mat))
        assoc_mat <- assoc_mat[, col_ord, drop = FALSE]
      }

      if (!is.null(p_mat)) {
        p_mat <- p_mat[rownames(assoc_mat), colnames(assoc_mat), drop = FALSE]
      }
      if (!is.null(q_mat)) {
        q_mat <- q_mat[rownames(assoc_mat), colnames(assoc_mat), drop = FALSE]
      }
      raw_assoc_mat <- raw_assoc_mat[rownames(assoc_mat), colnames(assoc_mat), drop = FALSE]

      heat_df <- as.data.frame(as.table(assoc_mat), stringsAsFactors = FALSE)
      colnames(heat_df) <- c("AxisX", "AxisY", "Association")
      raw_df <- as.data.frame(as.table(raw_assoc_mat), stringsAsFactors = FALSE)
      colnames(raw_df) <- c("AxisX", "AxisY", "AssociationRaw")
      heat_df <- dplyr::left_join(heat_df, raw_df, by = c("AxisX", "AxisY"))
      if (!is.null(p_mat)) {
        p_df <- as.data.frame(as.table(p_mat), stringsAsFactors = FALSE)
        colnames(p_df) <- c("AxisX", "AxisY", "PValue")
        heat_df <- dplyr::left_join(heat_df, p_df, by = c("AxisX", "AxisY"))
      } else {
        heat_df$PValue <- NA_real_
      }
      if (!is.null(q_mat)) {
        q_df <- as.data.frame(as.table(q_mat), stringsAsFactors = FALSE)
        colnames(q_df) <- c("AxisX", "AxisY", "QValue")
        heat_df <- dplyr::left_join(heat_df, q_df, by = c("AxisX", "AxisY"))
      } else {
        heat_df$QValue <- NA_real_
      }
      q_cutoff <- as.numeric(input$q_cutoff)
      if (is.na(q_cutoff)) q_cutoff <- 0.05
      q_cutoff <- max(0, min(1, q_cutoff))

      if (isTRUE(input$apply_sig_filter) && sig_supported) {
        sig_keep <- is.finite(heat_df$QValue) & heat_df$QValue <= q_cutoff
        heat_df$Association[!sig_keep] <- NA_real_
      }
      heat_df$SigMark <- ""
      if (sig_supported) {
        sig_cells <- is.finite(heat_df$QValue) & heat_df$QValue <= q_cutoff
        heat_df$SigMark[sig_cells] <- "*"
      }

      list(
        assoc_mat = assoc_mat,
        column_group = column_group,
        column_panels = column_panels,
        heat_df = heat_df,
        n_samples = ncol(analysis_mat),
        n_taxa = nrow(analysis_mat),
        assoc_type = assoc_type,
        sig_supported = sig_supported,
        sig_note = sig_note,
        q_cutoff = q_cutoff,
        row_hclust = row_hclust,
        col_hclust = col_hclust,
        n_tests = sum(is.finite(heat_df$QValue)),
        n_sig = sum(is.finite(heat_df$QValue) & heat_df$QValue <= q_cutoff, na.rm = TRUE),
        run_settings = isolate(current_analysis_settings())
      )
    }, ignoreNULL = TRUE)

    display_corr_payload <- reactive({
      req(heatmap_run_nonce())
      payload <- corr_payload()
      validate(
        need(
          identical(payload$run_settings, current_analysis_settings()),
          "Analysis settings changed after the last run. Click 'Run Heatmap' to update all result tabs."
        )
      )
      payload
    })

    observeEvent(corr_payload(), {
      taxa_choices <- rownames(corr_payload()$assoc_mat)
      current_selected <- isolate(input$selected_taxa)
      selected_keep <- if (length(current_selected) > 0) intersect(current_selected, taxa_choices) else character(0)
      updateSelectizeInput(session, "selected_taxa", choices = taxa_choices, selected = selected_keep, server = TRUE)
    }, ignoreInit = FALSE)

    cell_dims_px <- reactive({
      list(cell_width = 30, cell_height = 10)
    })

    heatmap_palette <- reactive({
      switch(
        input$color_palette,
        blue_orange = list(colors = c("#3B4CC0", "#F7F7F7", "#E66101"), label = "Blue - White - Orange"),
        purple_green = list(colors = c("#762A83", "#F7F7F7", "#1B7837"), label = "Purple - White - Green"),
        teal_brown = list(colors = c("#018571", "#F5F5F5", "#A6611A"), label = "Teal - White - Brown"),
        list(colors = c("#2166AC", "#F7F7F7", "#B2182B"), label = "Blue - White - Red")
      )
    })

    heatmap_plot_obj <- reactive({
      payload <- display_corr_payload()
      heat_df <- payload$heat_df
      assoc_mat <- payload$assoc_mat
      panels <- payload$column_panels
      if (!is.null(panels)) {
        primary_order <- plot_order_factor(panels$primary, input$group_var, plot_order(), unique(panels$primary))
        secondary_order <- plot_order_factor(panels$secondary, payload$column_group, plot_order(), unique(panels$secondary))
        panel_order <- order(primary_order, secondary_order)
        panels <- panels[panel_order, , drop = FALSE]
        assoc_mat <- assoc_mat[, rownames(panels), drop = FALSE]
        panel_split <- droplevels(primary_order[panel_order])
      } else if (!is.null(payload$column_group) && !isTRUE(input$show_col_dendrogram)) {
        columns <- colnames(assoc_mat)
        ordered <- levels(plot_order_factor(columns, payload$column_group, plot_order(), columns))
        assoc_mat <- assoc_mat[, ordered, drop = FALSE]
      }
      if (length(input$selected_taxa) > 0) {
        taxa_keep <- intersect(input$selected_taxa, rownames(assoc_mat))
        validate(need(length(taxa_keep) > 0, "None of the selected taxa are available in the current matrix."))
        assoc_mat <- assoc_mat[taxa_keep, , drop = FALSE]
        heat_df <- heat_df[as.character(heat_df$AxisX) %in% taxa_keep, , drop = FALSE]
      }
      validate(
        need(nrow(heat_df) > 0, "No cells to plot with current options.")
      )
      validate(
        need(requireNamespace("ComplexHeatmap", quietly = TRUE), "Package 'ComplexHeatmap' is required."),
        need(requireNamespace("circlize", quietly = TRUE), "Package 'circlize' is required.")
      )

      finite_associations <- abs(heat_df$Association[is.finite(heat_df$Association)])
      fill_limit <- if (length(finite_associations) > 0) {
        as.numeric(stats::quantile(finite_associations, probs = 0.98, names = FALSE, type = 8))
      } else {
        1
      }
      if (!is.finite(fill_limit) || fill_limit <= 0) {
        fill_limit <- max(finite_associations, na.rm = TRUE)
      }
      if (!is.finite(fill_limit) || fill_limit <= 0) fill_limit <- 1
      abs_limit_input <- as.numeric(input$abs_limit)
      if (is.finite(abs_limit_input) && abs_limit_input > 0) {
        fill_limit <- abs_limit_input
      }

      base_size <- as.numeric(input$base_size)
      if (!is.finite(base_size)) base_size <- 11
      base_size <- max(6, min(30, base_size))
      sig_fontsize_pt <- max(9, base_size * 1.15)

      # Let the manual layout allocate the remaining canvas space to the body.
      body_width <- body_height <- NULL
      if (!identical(input$plot_size_mode, "manual")) {
        dims_cell <- cell_dims_px()
        body_width <- grid::unit(ncol(assoc_mat) * dims_cell$cell_width * 0.75, "pt")
        body_height <- grid::unit(nrow(assoc_mat) * dims_cell$cell_height * 0.75, "pt")
      }

      sig_mat <- matrix("", nrow = nrow(assoc_mat), ncol = ncol(assoc_mat), dimnames = dimnames(assoc_mat))
      if ("SigMark" %in% colnames(heat_df) && nrow(heat_df) > 0) {
        for (i in seq_len(nrow(heat_df))) {
          rx <- as.character(heat_df$AxisX[i])
          cy <- as.character(heat_df$AxisY[i])
          if (!is.na(rx) && !is.na(cy) && rx %in% rownames(sig_mat) && cy %in% colnames(sig_mat)) {
            sig_mat[rx, cy] <- heat_df$SigMark[i]
          }
        }
      }

      col_fun <- circlize::colorRamp2(
        c(-fill_limit, 0, fill_limit),
        heatmap_palette()$colors
      )

      legend_title <- if (identical(input$association_mode, "prevalence")) {
        if (identical(analysis_target(), "taxa_group")) {
          "Delta prevalence\n(in vs out)"
        } else {
          "Phi coefficient\n(binary association)"
        }
      } else if (identical(input$value_scale, "zscore")) {
        "Z-score\n(by taxa)"
      } else {
        if (identical(input$transform_method, "clr")) {
          "CLR\nvalue"
        } else {
          "TSS\nvalue"
        }
      }

      show_row_dend_opt <- isTRUE(input$show_row_dendrogram)
      show_col_dend_opt <- is.null(panels) && isTRUE(input$show_col_dendrogram)

      # Recompute clustering on the currently displayed matrix (after taxa filtering)
      # so dendrogram size always matches matrix dimensions.
      row_hc_current <- NULL
      if (show_row_dend_opt && nrow(assoc_mat) >= 3) {
        row_hc_current <- tryCatch(
          stats::hclust(stats::dist(assoc_mat, method = "euclidean"), method = "average"),
          error = function(e) NULL
        )
      }
      col_hc_current <- NULL
      if (show_col_dend_opt && ncol(assoc_mat) >= 3) {
        col_hc_current <- tryCatch(
          stats::hclust(stats::dist(t(assoc_mat), method = "euclidean"), method = "average"),
          error = function(e) NULL
        )
      }

      cluster_rows_opt <- if (show_row_dend_opt && !is.null(row_hc_current)) row_hc_current else FALSE
      cluster_cols_opt <- if (show_col_dend_opt && !is.null(col_hc_current)) col_hc_current else FALSE

      if (isTRUE(input$apply_sig_filter) && isTRUE(payload$sig_supported)) {
        assoc_mat[sig_mat != "*"] <- NA_real_
      }

      ComplexHeatmap::Heatmap(
        assoc_mat,
        name = legend_title,
        col = col_fun,
        na_col = "#D9D9D9",
        rect_gp = grid::gpar(col = "#FFFFFF", lwd = 0.4),
        border = TRUE,
        border_gp = grid::gpar(col = "#595959", lwd = 1.2),
        cluster_rows = cluster_rows_opt,
        cluster_columns = cluster_cols_opt,
        column_split = if (is.null(panels)) NULL else panel_split,
        cluster_column_slices = FALSE,
        column_labels = if (is.null(panels)) colnames(assoc_mat) else panels$secondary,
        show_row_dend = show_row_dend_opt,
        show_column_dend = show_col_dend_opt,
        row_dend_side = "left",
        row_names_side = if (show_row_dend_opt) "right" else "left",
        column_names_side = "bottom",
        column_names_rot = 90,
        row_names_gp = grid::gpar(fontsize = max(6, base_size - 2)),
        column_names_gp = grid::gpar(fontsize = max(6, base_size - 2)),
        heatmap_legend_param = list(
          title = legend_title,
          title_position = "topcenter",
          title_gp = grid::gpar(fontsize = max(6, base_size - 2)),
          labels_gp = grid::gpar(fontsize = max(5, base_size - 3))
        ),
        width = body_width,
        height = body_height,
        column_title = if (is.null(panels)) NULL else "%s",
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (identical(sig_mat[i, j], "*")) {
            grid::grid.text(
              "*", x = x, y = y - grid::unit(2.5, "pt"),
              gp = grid::gpar(fontsize = sig_fontsize_pt, fontface = "bold")
            )
          }
        }
      )
    })

    plot_dims_px <- reactive({
      if (identical(input$plot_size_mode, "manual")) {
        clamp_dimension <- function(value, fallback, minimum) {
          if (is.null(value) || !is.finite(value)) value <- fallback
          as.integer(max(minimum, min(2400, value)))
        }
        return(list(
          width_px = clamp_dimension(input$plot_width, 900L, 400L),
          height_px = clamp_dimension(input$plot_height, 550L, 300L)
        ))
      }
      payload <- tryCatch(display_corr_payload(), error = function(e) NULL)
      if (is.null(payload)) {
        return(list(width_px = 900L, height_px = 500L))
      }
      assoc_mat <- payload$assoc_mat
      if (length(input$selected_taxa) > 0) {
        taxa_keep <- intersect(input$selected_taxa, rownames(assoc_mat))
        if (length(taxa_keep) > 0) {
          assoc_mat <- assoc_mat[taxa_keep, , drop = FALSE]
        }
      }
      n_cols <- ncol(assoc_mat)
      n_rows <- nrow(assoc_mat)
      dims_cell <- cell_dims_px()
      cell_width <- dims_cell$cell_width
      cell_height <- dims_cell$cell_height

      # Reserve margins for axis labels and legend, then scale panel by cell size.
      left_extra_px <- if (isTRUE(input$show_row_dendrogram)) 120 else 0
      panel_gap_px <- if (is.null(payload$column_panels)) 0 else
        8 * (length(unique(payload$column_panels$primary)) - 1L)
      width_px <- as.integer((n_cols * cell_width) + 340 + left_extra_px + panel_gap_px)
      height_px <- as.integer((n_rows * cell_height) + 230)

      list(width_px = width_px, height_px = height_px)
    })

    output$corr_heatmap_plot <- renderPlot(
      {
        if (is.null(heatmap_run_nonce())) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Click 'Run Heatmap' to start analysis.", cex = 0.85)
          return(invisible(NULL))
        }
        ht <- heatmap_plot_obj()
        grid::grid.newpage()
        left_pad_mm <- if (isTRUE(input$show_row_dendrogram)) 16 else 10
        top_pad_mm <- 10
        ComplexHeatmap::draw(
          ht,
          heatmap_legend_side = "right",
          padding = grid::unit(c(top_pad_mm, 2, 6, left_pad_mm), "mm")
        )
      },
      width = function() plot_dims_px()$width_px,
      height = function() plot_dims_px()$height_px,
      res = 110
    )

    output$heatmap_status <- renderText({
      req(display_corr_payload())
      payload <- display_corr_payload()
      primary_group_label <- if (is.null(input$group_var) || !nzchar(input$group_var)) "(Not selected)" else input$group_var
      primary_level_label <- if (is.null(input$primary_level) || !nzchar(input$primary_level)) "(All / Not selected)" else as.character(input$primary_level)
      current_secondary_var <- secondary_var_resolved()
      secondary_group_label <- if (is.null(current_secondary_var)) "(None)" else current_secondary_var
      analysis_scope_label <- if (is.null(current_secondary_var)) {
        "Primary group is used directly (no primary-level prefilter)."
      } else if (!is.null(payload$column_panels)) {
        "Secondary groups are analyzed separately within each displayed primary group."
      } else if (!isTruthy(input$primary_level) || identical(input$primary_level, "All")) {
        "All primary levels are pooled for the continuous secondary variable."
      } else {
        paste0(
          "Filtered to Primary level '", primary_level_label,
          "' and compared across Secondary group levels."
        )
      }
      analyzed_group_levels <- if (ncol(payload$assoc_mat) > 1) {
        paste(colnames(payload$assoc_mat), collapse = ", ")
      } else {
        colnames(payload$assoc_mat)[1]
      }
      test_method_label <- if (identical(input$group_type, "continuous")) {
        if (identical(input$association_mode, "prevalence")) {
          "Point-biserial correlation test (binary taxa vs continuous group) + BH correction"
        } else {
          paste0("Correlation test (", tools::toTitleCase(input$corr_method), ") + BH correction")
        }
      } else if (identical(input$group_type, "categorical")) {
        if (identical(input$association_mode, "prevalence")) {
          "One-vs-rest Fisher exact test (taxa presence x group) + BH correction"
        } else {
          "One-vs-rest Wilcoxon rank-sum test (taxa x group) + BH correction"
        }
      } else {
        if (identical(input$association_mode, "prevalence")) {
          "Auto: continuous uses point-biserial test, categorical uses one-vs-rest Fisher exact + BH correction"
        } else {
          "Auto: continuous uses correlation test, categorical uses one-vs-rest Wilcoxon + BH correction"
        }
      }
      n_taxa_display <- if (length(input$selected_taxa) > 0) {
        length(intersect(input$selected_taxa, rownames(payload$assoc_mat)))
      } else {
        payload$n_taxa
      }
      paste(
        c(
          paste0("Analysis target: ", if (identical(analysis_target(), "taxa_taxa")) "Taxa vs Taxa" else "Taxa vs Group"),
          paste0("Association type: ", payload$assoc_type),
          paste0("1) Primary group variable: ", primary_group_label),
          paste0("2) Primary level filter: ", primary_level_label),
          "Primary-level dropdown is shown only when the number of distinct levels is fewer than 10.",
          paste0("3) Secondary group variable: ", secondary_group_label),
          paste0("Current analysis scope: ", analysis_scope_label),
          paste0("Analyzed group levels/axis: ", analyzed_group_levels),
          paste0("Samples used: ", payload$n_samples),
          paste0("Taxa shown: ", n_taxa_display),
          paste0("Association mode: ", if (identical(input$association_mode, "prevalence")) "Prevalence-based (binary)" else "Abundance-based"),
          paste0("Transform: ", if (identical(input$association_mode, "prevalence")) "Not applied (binary prevalence)" else if (identical(input$transform_method, "clr")) "CLR" else "TSS"),
          paste0("Method: ", if (identical(input$association_mode, "prevalence")) "Pearson/phi for correlation context" else tools::toTitleCase(input$corr_method)),
          paste0(
            "Heatmap scale: ",
            if (identical(input$association_mode, "prevalence")) {
              "Raw (delta prevalence for categorical groups)"
            } else if (identical(input$value_scale, "zscore")) {
              "Z-score (by taxa)"
            } else {
              "Raw"
            }
          ),
          paste0("Color palette: ", heatmap_palette()$label),
          paste0("Color scale abs limit: ", if (is.finite(as.numeric(input$abs_limit)) && as.numeric(input$abs_limit) > 0) format(as.numeric(input$abs_limit), digits = 3) else "Auto (98th percentile)"),
          paste0("Plot size mode: ", if (identical(input$plot_size_mode, "manual")) "Manual" else "Auto"),
          if (identical(input$plot_size_mode, "manual")) paste0("Plot width (px): ", plot_dims_px()$width_px) else NULL,
          if (identical(input$plot_size_mode, "manual")) paste0("Plot height (px): ", plot_dims_px()$height_px) else NULL,
          paste0("Significance mark size: linked to base font size (", format(as.numeric(input$base_size), digits = 3), ")"),
          paste0("Significance filter: ", if (isTRUE(input$apply_sig_filter)) "On" else "Off"),
          paste0("Significance test method: ", test_method_label),
          paste0("FDR cutoff (q): ", format(as.numeric(input$q_cutoff), digits = 3)),
          paste0("Significance note: ", payload$sig_note),
          if (isTRUE(payload$sig_supported)) paste0("Significance tests: ", payload$n_tests) else NULL,
          if (isTRUE(payload$sig_supported)) paste0("Significant cells (*): ", payload$n_sig) else NULL,
          paste0("Taxonomic level: ", input$tax_level),
          if (identical(analysis_target(), "taxa_group")) paste0("Group type: ", tools::toTitleCase(input$group_type)) else NULL,
          if (identical(analysis_target(), "taxa_group")) {
            paste0(
              "Group variable: ",
              if (!is.null(current_secondary_var)) current_secondary_var else input$group_var
            )
          } else NULL
        ),
        collapse = "\n"
      )
    })

    output$heatmap_status_box <- renderUI({
      req(plot_dims_px())
      box_width <- plot_dims_px()$width_px
      if (is.null(box_width) || !is.finite(box_width)) {
        box_width <- 600
      }
      tags$div(
        class = "simple-result-card",
        style = paste0("width: ", box_width, "px; max-width: 100%;"),
        verbatimTextOutput(session$ns("heatmap_status"))
      )
    })

    output$heatmap_legend_box <- renderUI({
      req(plot_dims_px())
      box_width <- plot_dims_px()$width_px
      if (is.null(box_width) || !is.finite(box_width)) {
        box_width <- 600
      }
      tags$div(
        style = paste(
          "margin-top: 12px;",
          sprintf("width: %spx;", box_width),
          "max-width: 100%;",
          "padding: 12px 14px;",
          "border: 1px solid #e5e7eb;",
          "border-left: 4px solid #6b7280;",
          "border-radius: 8px;",
          "background: linear-gradient(180deg, #fcfcfd 0%, #f7f8fa 100%);",
          "box-shadow: 0 1px 2px rgba(0,0,0,0.04);",
          "box-sizing: border-box;"
        ),
        tags$div(
          style = "color: #1f2937; font-size: 12.5px; line-height: 1.55;",
          uiOutput(session$ns("heatmap_figure_legend"))
        )
      )
    })

    output$heatmap_status_separator <- renderUI({
      req(plot_dims_px())
      box_width <- plot_dims_px()$width_px
      if (is.null(box_width) || !is.finite(box_width)) {
        box_width <- 600
      }
      tags$hr(
        style = paste(
          sprintf("width: %spx;", box_width),
          "max-width: 100%;",
          "margin: 14px 0 12px 0;",
          "border-top: 1px solid #d1d5db;"
        )
      )
    })

    output$heatmap_figure_legend <- renderUI({
      req(display_corr_payload())
      req(input$tax_level, input$transform_method, input$association_mode, input$corr_method, input$value_scale)
      payload <- tryCatch(display_corr_payload(), error = function(e) NULL)
      assoc_type_label <- if (!is.null(payload) && "assoc_type" %in% names(payload)) payload$assoc_type else "N/A"
      significance_text <- "No significance tests were performed."
      if (isTRUE(payload$sig_supported)) {
        significance_text <- paste0(
          "An asterisk (*) marks a cell with a Benjamini–Hochberg-adjusted p-value (q) ≤ ",
          format(payload$q_cutoff), ". ", payload$sig_note
        )
      }
      if (isTRUE(payload$sig_supported) && isTRUE(input$apply_sig_filter)) {
        significance_text <- paste0(
          significance_text,
          " Cells above this cutoff and cells without a valid q-value are masked."
        )
      }
      current_secondary_var <- secondary_var_resolved()
      group_label <- if (!is.null(current_secondary_var)) {
        current_secondary_var
      } else if (is.null(input$group_var) || !nzchar(input$group_var)) {
        "selected group variable"
      } else {
        input$group_var
      }
      title_text <- "Association heatmap"
      body_text <- if (!is.null(payload$column_panels)) {
        cell_summary <- if (identical(input$association_mode, "prevalence")) {
          "Each cell shows the difference in prevalence between the indicated secondary group and the other secondary groups in that panel."
        } else {
          paste0(
            "Each cell summarizes the mean ", toupper(input$transform_method),
            "-transformed abundance for that secondary group.",
            if (identical(input$value_scale, "zscore")) " Values are standardized by taxon across the displayed panels." else ""
          )
        }
        paste0(
          "Panels represent ", input$group_var, ", and columns within each panel represent ",
          group_label, ". ", cell_summary,
          " Panels share the same color scale. Gray marks unavailable or masked values."
        )
      } else if (identical(input$association_mode, "prevalence")) {
        paste0(
          "Cells represent presence/absence-based association values at ",
          tolower(input$tax_level),
          " level. Taxa abundance is binarized to prevalence (0/1). For categorical groups, cell values are delta prevalence (in-group minus out-group). Color indicates direction and magnitude, with value scale fixed to raw values",
          ". Group variable is ",
          group_label,
          ". Association type: ",
          assoc_type_label,
          "."
        )
      } else {
        transform_desc <- if (identical(input$value_scale, "zscore")) {
          if (identical(input$transform_method, "clr")) {
            "CLR-transformed abundance with row-wise z-score scaling"
          } else {
            "TSS-transformed abundance with row-wise z-score scaling"
          }
        } else {
          paste0(toupper(input$transform_method), " transformed abundance")
        }
        paste0(
          "Cells represent taxa-to-group association values at ",
          tolower(input$tax_level),
          " level. The analysis uses ",
          transform_desc,
          " and ",
          tools::toTitleCase(input$corr_method),
          " statistics. Color indicates association direction and magnitude, with value scale set to ",
          if (identical(input$value_scale, "zscore")) "row-wise z-score" else "raw association value",
          ". Group variable is ",
          group_label,
          ". Association type: ",
          assoc_type_label,
          "."
        )
      }
      tags$div(
        tags$div(style = "font-weight: 600; margin-bottom: 4px;", title_text),
        tags$div(
          style = "white-space: normal; overflow-wrap: anywhere;",
          tags$div(body_text),
          tags$div(style = "margin-top: 4px;", significance_text)
        )
      )
    })

    output$download_heatmap_png <- downloadHandler(
      filename = function() {
        paste0(
          "association_heatmap_",
          analysis_target(),
          "_",
          if (identical(input$transform_method, "clr")) "clr" else "tss",
          "_",
          input$corr_method,
          "_",
          Sys.Date(),
          ".png"
        )
      },
      content = function(file) {
        req(display_corr_payload())
        dims <- plot_dims_px()
        dpi_val <- 300
        grDevices::png(
          filename = file,
          width = dims$width_px / 110,
          height = dims$height_px / 110,
          units = "in",
          res = dpi_val,
          bg = "white"
        )
        on.exit(grDevices::dev.off(), add = TRUE)
        ht <- heatmap_plot_obj()
        grid::grid.newpage()
        left_pad_mm <- if (isTRUE(input$show_row_dendrogram)) 16 else 10
        top_pad_mm <- 10
        ComplexHeatmap::draw(
          ht,
          heatmap_legend_side = "right",
          padding = grid::unit(c(top_pad_mm, 2, 6, left_pad_mm), "mm")
        )
      }
    )

    output$download_corr_tsv <- downloadHandler(
      filename = function() paste0("association_matrix_", analysis_target(), "_", Sys.Date(), ".tsv"),
      content = function(file) {
        req(display_corr_payload())
        mat_df <- as.data.frame(display_corr_payload()$assoc_mat, stringsAsFactors = FALSE)
        mat_df <- cbind(Taxon = rownames(mat_df), mat_df, stringsAsFactors = FALSE)
        utils::write.table(mat_df, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    )

  })
}
