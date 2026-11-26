## UI
mod_biplot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML("
      .well h4 { font-size: 16px; }
      .well h5 { font-size: 13px; }
      .well .control-label { font-size: 12px; }
      .well .checkbox label { font-size: 12px; }
      .well .form-control { font-size: 12px; }
      .well .btn { font-size: 11px; }
      .biplot-run-btn .biplot-run-icon {
        display: none;
        margin-right: 6px;
      }
      .biplot-run-btn .biplot-run-icon i {
        display: inline-block;
      }
      .biplot-run-btn.is-running .biplot-run-icon {
        display: inline-block;
        animation: biplot-run-spin 0.9s linear infinite;
      }
      .biplot-run-btn.is-running .biplot-run-label {
        opacity: 0.85;
      }
      @keyframes biplot-run-spin {
        from { transform: rotate(0deg); }
        to { transform: rotate(360deg); }
      }
    ")),
    plot_scroll_layout(
      sidebarPanel(
        width = 2,
        h4(icon("arrows-up-down-left-right"), "Constrained Ordination"),
        hr(),
        selectInput(ns("group_var"), "Primary Explanatory Variable", choices = NULL),
        hr(),
        h4(icon("sliders"), "Analysis Settings"),
        selectInput(
          ns("analysis_method"),
          "Analysis Method",
          choices = c("dbRDA (capscale)" = "dbrda", "CCA" = "cca"),
          selected = "dbrda"
        ),
        selectInput(
          ns("distance_metric"),
          "Distance Metric",
          choices = c(
            "Bray-Curtis" = "bray",
            "Jaccard" = "jaccard",
            "Jensen-Shannon" = "jsd",
            "Aitchison" = "aitchison"
          ),
          selected = "bray"
        ),
        selectInput(
          ns("tax_level"),
          "Taxonomic Level",
          choices = c("ASV", "Genus", "Species", "Strain"),
          selected = "Genus"
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'jaccard'", ns("distance_metric")),
          textOutput(ns("jaccard_rarefaction_depth"))
        ),
        tags$details(style = "margin-top: 8px; margin-bottom: 8px;",
          tags$summary(strong("Advanced Model Settings")),
          br(),
          selectizeInput(
            ns("additional_vars"),
            "Additional Explanatory Variables",
            choices = NULL,
            multiple = TRUE,
            options = list(plugins = list("remove_button"))
          ),
          selectizeInput(
            ns("conditioning_vars"),
            "Conditioning Variables",
            choices = NULL,
            multiple = TRUE,
            options = list(plugins = list("remove_button"))
          ),
          checkboxInput(ns("standardize_continuous"), "Standardize Continuous Variables", value = TRUE),
          uiOutput(ns("model_formula_preview"))
        ),
        tags$details(style = "margin-top: 8px; margin-bottom: 8px;",
          tags$summary(strong("Taxa Settings")),
          br(),
          checkboxInput(ns("show_taxa_vectors"), "Show Taxa Scores", value = FALSE),
          numericInput(
            ns("prevalence_filter_pct"),
            "Prevalence Filter (%)",
            value = 10,
            min = 0,
            max = 100,
            step = 1
          ),
          numericInput(
            ns("max_taxa"),
            "Max Taxa",
            value = 100,
            min = 10,
            max = 500,
            step = 10
          ),
          numericInput(
            ns("top_taxa_vectors"),
            "Top Taxa Scores",
            value = 15,
            min = 0,
            max = 100,
            step = 1
          ),
          selectizeInput(
            ns("selected_taxa"),
            "Selected Taxa (Optional)",
            choices = NULL,
            selected = NULL,
            multiple = TRUE,
            options = list(
              placeholder = "Selected taxa are drawn instead of top taxa scores",
              plugins = list("remove_button")
            )
          )
        ),
        tags$details(style = "margin-top: 8px; margin-bottom: 8px;",
          tags$summary(strong("Explanatory Vectors")),
          br(),
          checkboxInput(ns("show_environment_vectors"), "Show Explanatory Vectors", value = TRUE),
          selectizeInput(
            ns("environment_vector_vars"),
            "Explanatory Vectors to Display",
            choices = NULL,
            multiple = TRUE,
            options = list(plugins = list("remove_button"))
          ),
          checkboxInput(ns("environment_only_significant"), "Show Only Significant Variables", value = FALSE),
          numericInput(ns("environment_p_cutoff"), "p-value Cutoff", value = 0.05, min = 0.001, max = 1, step = 0.01),
          checkboxInput(ns("show_environment_pvalues"), "Show p-values in Labels", value = FALSE),
          numericInput(ns("max_environment_vectors"), "Maximum Vectors", value = 10, min = 1, max = 50, step = 1)
        ),
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Display Options"),
        numericInput(ns("plot_width"), "Plot Width (px)", value = 700, min = 400, max = 2400, step = 50),
        numericInput(ns("plot_height"), "Plot Height (px)", value = 500, min = 300, max = 2400, step = 50),
        numericInput(ns("base_size"), "Base Font Size", value = 14, min = 6, max = 30, step = 1),
        numericInput(ns("dot_size"), "Dot Size (Point Size)", value = 5, min = 0.5, max = 10, step = 0.5),
        checkboxInput(ns("show_dot_outline"), "Show Dot Outline", value = TRUE),
        checkboxInput(ns("show_ellipses"), "Show Group Ellipses", value = FALSE),
        checkboxInput(ns("show_group_hulls"), "Show Group Hulls", value = FALSE),
        checkboxInput(ns("show_centroid_lines"), "Show Sample-to-Centroid Lines", value = FALSE),
        checkboxInput(ns("connect_subject_samples"), "Connect Samples within Subjects", value = FALSE),
        uiOutput(ns("subject_id_selector")),
        checkboxInput(ns("show_group_centroid"), "Show Primary-Group Centroids", value = TRUE),
        checkboxInput(ns("show_additional_centroids"), "Show Centroids for Additional Categorical Variables", value = FALSE),
        checkboxInput(ns("show_sample_names"), "Show Sample Names", value = FALSE),
        tags$button(
          id = ns("run_biplot"),
          type = "button",
          class = "action-button btn btn-danger biplot-run-btn",
          `data-val` = 0,
          style = "font-size: 12px;",
          tags$span(class = "biplot-run-icon", icon("spinner", class = "fa-spin")),
          tags$span(class = "biplot-run-label", "Run Ordination")
        ),
        hr(),
        h5(icon("download"), "Download"),
        tags$div(
          style = "display: flex; gap: 4px; align-items: center; flex-wrap: nowrap;",
          downloadButton(
            ns("download_biplot_png"),
            "Download Plot (PNG)",
            style = "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
          ),
          downloadButton(
            ns("download_biplot_scores"),
            "Download Scores (TSV)",
            style = "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
          )
        ),
        tags$script(HTML(
          "Shiny.addCustomMessageHandler('toggle-biplot-run-btn', function(msg) {
             var btn = document.getElementById(msg.id);
             if (!btn) return;
             btn.disabled = !!msg.disabled;
             btn.classList.toggle('is-running', !!msg.running);
             var label = btn.querySelector('.biplot-run-label');
             if (label && msg.label) label.textContent = msg.label;
            });
           $(document).on('change', 'select[id$=\"analysis_method\"]', function() {
             var method = $(this).val();
             var panel = $(this).closest('.well, .sidebar, .form-group, body');
             var distanceSelect = panel.find('select[id$=\"distance_metric\"]');
             var disable = (method === 'cca');
             distanceSelect.prop('disabled', disable);
             distanceSelect.closest('.form-group').find('.control-label').toggleClass('text-muted', disable);
           });
           $(document).on('shiny:connected', function() {
             $('select[id$=\"analysis_method\"]').trigger('change');
           });"
        ))
      ),
      mainPanel(
        h4("Constrained Ordination"),
        analysis_run_view(ns, "run_biplot",
        plotOutput(ns("biplot_plot"), height = "auto"),
        plot_legend_output(ns, "biplot_legend_box", c("biplot_plot")),
        uiOutput(ns("biplot_status_separator")),
        tabsetPanel(
          id = ns("biplot_results_tab"),
          tabPanel("Model Summary", uiOutput(ns("biplot_status_box"))),
          tabPanel("Term Tests", tableOutput(ns("biplot_term_tests"))),
          tabPanel("Collinearity", tableOutput(ns("biplot_vif_table"))),
          tabPanel("Sample Exclusions", tableOutput(ns("biplot_exclusions")))
        )
        )
      )
    )
  )
}

## Server
mod_biplot_server <- function(id, ps_obj, meta_vars = NULL) {
  moduleServer(id, function(input, output, session) {
    output$run_prompt <- analysis_run_prompt("Run Ordination")
    output$subject_id_selector <- renderUI({
      req(isTRUE(input$connect_subject_samples), ps_obj())
      metadata <- phyloseq::sample_data(ps_obj(), errorIfNULL = FALSE)
      req(!is.null(metadata))
      choices <- setdiff(colnames(metadata), input$group_var)
      choices <- choices[gsub("[^a-z0-9]", "", tolower(choices)) != "sampleid"]
      if (!length(choices)) return(helpText("No metadata columns are available for Subject ID."))
      selectInput(session$ns("subject_id_var"), "Subject ID Variable", choices = choices)
    })
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

    jaccard_rarefaction <- reactive({
      req(identical(input$distance_metric, "jaccard"), ps_obj())
      reads <- phyloseq::sample_sums(ps_obj())
      validate(
        need(
          length(reads) > 0L && all(is.finite(reads) & reads >= 1),
          "Jaccard rarefaction requires selected samples with positive read counts."
        )
      )
      list(depth = min(reads), subsample = any(reads > min(reads)))
    })

    output$jaccard_rarefaction_depth <- renderText({
      info <- jaccard_rarefaction()
      paste0("Jaccard rarefaction depth: ", scales::comma(info$depth), " reads per sample.")
    })

    observeEvent(ps_obj(), {
      req(ps_obj())
      sd <- phyloseq::sample_data(ps_obj(), errorIfNULL = FALSE)
      if (is.null(sd)) {
        updateSelectInput(session, "group_var", choices = character(0), selected = NULL)
        return()
      }
      meta_df <- as.data.frame(sd, stringsAsFactors = FALSE)
      group_choices <- setdiff(colnames(meta_df), "SampleID")
      sel <- if (!is.null(input$group_var) && input$group_var %in% group_choices) input$group_var else if (length(group_choices) > 0) group_choices[1] else NULL
      updateSelectInput(session, "group_var", choices = group_choices, selected = sel)
      additional_selected <- intersect(isolate(input$additional_vars), setdiff(group_choices, sel))
      conditioning_selected <- intersect(isolate(input$conditioning_vars), setdiff(group_choices, c(sel, additional_selected)))
      updateSelectizeInput(
        session, "additional_vars",
        choices = setdiff(group_choices, sel),
        selected = additional_selected,
        server = TRUE
      )
      updateSelectizeInput(
        session, "conditioning_vars",
        choices = setdiff(group_choices, sel),
        selected = conditioning_selected,
        server = TRUE
      )
    }, ignoreInit = FALSE)

    observeEvent(input$group_var, {
      req(ps_obj())
      sd <- phyloseq::sample_data(ps_obj(), errorIfNULL = FALSE)
      if (is.null(sd)) return()
      meta_df <- as.data.frame(sd, stringsAsFactors = FALSE)
      all_choices <- setdiff(colnames(meta_df), "SampleID")
      primary <- input$group_var
      additional <- setdiff(intersect(input$additional_vars, all_choices), primary)
      conditioning <- setdiff(intersect(input$conditioning_vars, all_choices), c(primary, additional))
      updateSelectizeInput(
        session, "additional_vars",
        choices = setdiff(all_choices, primary),
        selected = additional,
        server = TRUE
      )
      updateSelectizeInput(
        session, "conditioning_vars",
        choices = setdiff(all_choices, primary),
        selected = conditioning,
        server = TRUE
      )
    }, ignoreInit = TRUE)

    observeEvent(list(input$additional_vars, input$conditioning_vars), {
      additional <- unique(as.character(input$additional_vars))
      conditioning <- unique(as.character(input$conditioning_vars))
      conditioning_clean <- setdiff(conditioning, additional)
      if (!identical(conditioning, conditioning_clean)) {
        updateSelectizeInput(session, "conditioning_vars", selected = conditioning_clean, server = TRUE)
      }
    }, ignoreInit = TRUE)

    observeEvent(list(input$group_var, input$additional_vars), {
      req(ps_obj())
      sd <- phyloseq::sample_data(ps_obj(), errorIfNULL = FALSE)
      if (is.null(sd)) return()
      meta_df <- as.data.frame(sd, stringsAsFactors = FALSE)
      all_choices <- setdiff(colnames(meta_df), "SampleID")
      primary <- input$group_var
      additional <- setdiff(intersect(input$additional_vars, all_choices), primary)
      explanatory <- unique(c(primary, additional))
      continuous <- explanatory[vapply(meta_df[explanatory], function(x) is.numeric(x) || is.integer(x), logical(1))]
      selected_vectors <- intersect(isolate(input$environment_vector_vars), continuous)
      if (length(selected_vectors) == 0L) selected_vectors <- continuous
      updateSelectizeInput(
        session, "environment_vector_vars",
        choices = continuous,
        selected = selected_vectors,
        server = TRUE
      )
    }, ignoreInit = FALSE)

    selected_model_variables <- reactive({
      primary <- input$group_var
      additional <- setdiff(unique(input$additional_vars), primary)
      conditioning <- setdiff(unique(input$conditioning_vars), c(primary, additional))
      list(primary = primary, explanatory = unique(c(primary, additional)), conditioning = conditioning)
    })

    output$model_formula_preview <- renderUI({
      vars <- selected_model_variables()
      req(length(vars$primary) == 1L)
      rhs <- paste(vars$explanatory, collapse = " + ")
      if (length(vars$conditioning) > 0L) {
        rhs <- paste(c(rhs, paste0("Condition(", vars$conditioning, ")")), collapse = " + ")
      }
      sample_line <- NULL
      if (!is.null(ps_obj())) {
        meta <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
        used <- unique(c(vars$explanatory, vars$conditioning))
        if (length(used) > 0L && all(used %in% colnames(meta))) {
          complete <- stats::complete.cases(meta[, used, drop = FALSE])
          sample_line <- paste0("Complete samples: ", sum(complete), " / ", nrow(meta))
        }
      }
      tags$div(
        style = "font-size: 11px; padding: 8px; background: #f8fafc; border: 1px solid #d9e2ec; border-radius: 6px; margin-bottom: 8px;",
        tags$strong("Model Formula"),
        tags$div(paste0("Community ~ ", rhs), style = "overflow-wrap: anywhere;"),
        if (!is.null(sample_line)) tags$div(sample_line, style = "margin-top: 4px; color: #4b5563;")
      )
    })

    biplot_running <- reactiveVal(FALSE)
    biplot_run_nonce <- reactiveVal(0L)
    biplot_next_allowed_at <- reactiveVal(as.POSIXct(NA))

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
        ranks = phyloseq::rank_names(ps),
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
        additional_vars = sort(as.character(input$additional_vars)),
        conditioning_vars = sort(as.character(input$conditioning_vars)),
        standardize_continuous = isTRUE(input$standardize_continuous),
        tax_level = input$tax_level,
        analysis_method = input$analysis_method,
        distance_metric = input$distance_metric,
        prevalence_filter_pct = as.numeric(input$prevalence_filter_pct),
        max_taxa = as.integer(input$max_taxa),
        top_taxa_vectors = as.integer(input$top_taxa_vectors)
      )
    })

    observeEvent(input$run_biplot, {
      now_ts <- Sys.time()
      next_allowed <- biplot_next_allowed_at()
      if (isTRUE(biplot_running())) {
        showNotification("Constrained ordination is already running. Please wait for completion.", type = "message", duration = 3)
        return(invisible(NULL))
      }
      if (!is.na(next_allowed) && now_ts < next_allowed) {
        wait_sec <- ceiling(as.numeric(difftime(next_allowed, now_ts, units = "secs")))
        showNotification(paste0("Please wait ", wait_sec, " second(s) before running again."), type = "message", duration = 2)
        return(invisible(NULL))
      }
      biplot_next_allowed_at(now_ts + 10)
      biplot_run_nonce(biplot_run_nonce() + 1L)
    }, ignoreInit = TRUE)

    biplot_payload <- eventReactive(biplot_run_nonce(), {
      biplot_running(TRUE)
      session$sendCustomMessage(
        "toggle-biplot-run-btn",
        list(id = session$ns("run_biplot"), disabled = TRUE, running = TRUE, label = "Running...")
      )
      on.exit(
        {
          biplot_running(FALSE)
          session$sendCustomMessage(
            "toggle-biplot-run-btn",
            list(id = session$ns("run_biplot"), disabled = FALSE, running = FALSE, label = "Run Ordination")
          )
        },
        add = TRUE
      )
      req(ps_obj(), input$group_var, input$tax_level, input$analysis_method, input$distance_metric)
      validate(
        need(requireNamespace("vegan", quietly = TRUE), "Package 'vegan' is required for ordination.")
      )
      ps_data <- ps_obj()
      validate(
        need(phyloseq::nsamples(ps_data) >= 4, "At least 4 samples are required."),
        need(phyloseq::ntaxa(ps_data) >= 3, "At least 3 taxa are required.")
      )

      jaccard_info <- NULL
      if (identical(input$analysis_method, "dbrda") && identical(input$distance_metric, "jaccard")) {
        jaccard_info <- jaccard_rarefaction()
        if (isTRUE(jaccard_info$subsample)) {
          ps_data <- withr::with_preserve_seed(
            phyloseq::rarefy_even_depth(
              ps_data,
              sample.size = jaccard_info$depth,
              rngseed = 42,
              replace = FALSE,
              trimOTUs = TRUE,
              verbose = FALSE
            )
          )
        }
      }

      if (!identical(input$tax_level, "ASV")) {
        validate(
          need(!is.null(phyloseq::tax_table(ps_data, errorIfNULL = FALSE)), "Taxonomy table is required for aggregation.")
        )
        tax_cols <- colnames(phyloseq::tax_table(ps_data))
        validate(
          need(input$tax_level %in% tax_cols, paste0("Taxonomic rank '", input$tax_level, "' is not available."))
        )
        ps_data <- phyloseq::tax_glom(ps_data, taxrank = input$tax_level, NArm = FALSE)
        ps_data <- sanitize_taxa_names(ps_data, input$tax_level)
      }

      otu <- as.matrix(phyloseq::otu_table(ps_data))
      if (!phyloseq::taxa_are_rows(phyloseq::otu_table(ps_data))) {
        otu <- t(otu)
      }
      prev <- as.numeric(input$prevalence_filter_pct)
      if (is.na(prev)) prev <- 10
      prev <- max(0, min(100, prev))
      keep_prev <- rowMeans(otu > 0, na.rm = TRUE) * 100 >= prev
      otu <- otu[keep_prev, , drop = FALSE]
      validate(need(nrow(otu) >= 3, "Too few taxa remain after prevalence filtering."))

      max_taxa <- as.integer(input$max_taxa)
      if (is.na(max_taxa) || max_taxa < 3) max_taxa <- 100
      if (nrow(otu) > max_taxa) {
        otu <- otu[order(rowSums(otu, na.rm = TRUE), decreasing = TRUE)[seq_len(max_taxa)], , drop = FALSE]
      }

      if (identical(input$analysis_method, "dbrda") && identical(input$distance_metric, "aitchison")) {
        log_mat <- log(otu + 1)
        otu_model <- apply(log_mat, 2, function(x) x - mean(x, na.rm = TRUE))
        if (!is.matrix(otu_model)) otu_model <- matrix(otu_model, nrow = nrow(log_mat))
        rownames(otu_model) <- rownames(log_mat)
        colnames(otu_model) <- colnames(log_mat)
      } else if (identical(input$analysis_method, "dbrda") && identical(input$distance_metric, "jaccard")) {
        otu_model <- ifelse(otu > 0, 1, 0)
        dimnames(otu_model) <- dimnames(otu)
      } else {
        totals <- colSums(otu, na.rm = TRUE)
        totals[totals <= 0 | is.na(totals)] <- 1
        otu_model <- sweep(otu, 2, totals, "/")
      }

      taxa_var <- apply(otu_model, 1, stats::var, na.rm = TRUE)
      otu_model <- otu_model[is.finite(taxa_var) & taxa_var > 0, , drop = FALSE]
      validate(need(nrow(otu_model) >= 3, "Too few taxa remain after variance filtering."))

      sd <- phyloseq::sample_data(ps_data, errorIfNULL = FALSE)
      validate(need(!is.null(sd), "Metadata is required for constrained ordination."))
      meta_df <- as.data.frame(sd, stringsAsFactors = FALSE)
      selected_vars <- selected_model_variables()
      explanatory_vars <- selected_vars$explanatory
      conditioning_vars <- selected_vars$conditioning
      all_model_vars <- unique(c(explanatory_vars, conditioning_vars))
      validate(
        need(length(explanatory_vars) > 0L, "Select at least one explanatory variable."),
        need(all(all_model_vars %in% colnames(meta_df)), "One or more selected model variables are unavailable in metadata.")
      )

      sample_ids <- colnames(otu_model)
      meta_df <- meta_df[sample_ids, , drop = FALSE]
      model_source <- meta_df[, all_model_vars, drop = FALSE]
      nonempty <- vapply(model_source, function(x) {
        if (is.character(x) || is.factor(x)) !is.na(x) & nzchar(trimws(as.character(x))) else !is.na(x)
      }, logical(nrow(model_source)))
      if (is.null(dim(nonempty))) {
        nonempty <- matrix(nonempty, ncol = 1L, dimnames = list(rownames(model_source), all_model_vars))
      }
      keep_samples <- if (is.null(dim(nonempty))) nonempty else apply(nonempty, 1, all)
      excluded_samples <- data.frame(SampleID = character(0), MissingVariables = character(0))
      if (any(!keep_samples)) {
        excluded_samples <- data.frame(
          SampleID = rownames(model_source)[!keep_samples],
          MissingVariables = apply(nonempty[!keep_samples, , drop = FALSE], 1, function(x) {
            paste(all_model_vars[!x], collapse = ", ")
          }),
          stringsAsFactors = FALSE
        )
      }
      validate(need(sum(keep_samples) >= 4, "At least 4 complete samples are required for the selected model variables."))
      otu_model <- otu_model[, keep_samples, drop = FALSE]
      meta_df <- meta_df[keep_samples, , drop = FALSE]

      internal_names <- c(
        paste0("term_", seq_along(explanatory_vars)),
        paste0("condition_", seq_along(conditioning_vars))
      )
      names(internal_names) <- all_model_vars
      model_df <- data.frame(row.names = rownames(meta_df))
      variable_types <- character(length(all_model_vars))
      names(variable_types) <- all_model_vars
      for (variable in all_model_vars) {
        value <- meta_df[[variable]]
        internal <- internal_names[[variable]]
        if (is.numeric(value) || is.integer(value)) {
          value <- as.numeric(value)
          validate(need(is.finite(stats::sd(value)) && stats::sd(value) > 0, paste0("Variable '", variable, "' has no usable variation.")))
          if (isTRUE(input$standardize_continuous)) value <- as.numeric(scale(value))
          model_df[[internal]] <- value
          variable_types[[variable]] <- "Continuous"
        } else {
          value <- droplevels(as.factor(as.character(value)))
          validate(need(nlevels(value) >= 2L, paste0("Variable '", variable, "' must contain at least two levels.")))
          model_df[[internal]] <- value
          variable_types[[variable]] <- "Categorical"
        }
      }
      comm_df <- as.data.frame(t(otu_model), stringsAsFactors = FALSE)
      validate(
        need(nrow(model_df) == nrow(comm_df), "Model data and community matrix row counts do not match."),
        need(all(rownames(model_df) == rownames(comm_df)), "Sample IDs between model data and community matrix do not match.")
      )

      model_type <- input$analysis_method
      metric_label <- switch(
        input$distance_metric,
        aitchison = "Aitchison",
        jaccard = "Jaccard",
        jsd = "Jensen-Shannon",
        "Bray-Curtis"
      )
      transformation_label <- if (identical(model_type, "dbrda")) {
        switch(
          input$distance_metric,
          aitchison = "CLR (log(count + 1), centered within each sample)",
          jaccard = paste0(
            "Rarefaction to ", scales::comma(jaccard_info$depth),
            " reads per sample, followed by presence/absence"
          ),
          "TSS (relative abundance)"
        )
      } else {
        "TSS (relative abundance)"
      }
      explanatory_terms <- unname(internal_names[explanatory_vars])
      condition_terms <- unname(internal_names[conditioning_vars])
      rhs_terms <- explanatory_terms
      if (length(condition_terms) > 0L) rhs_terms <- c(rhs_terms, paste0("Condition(", condition_terms, ")"))
      response_name <- if (identical(model_type, "dbrda")) "dist_obj" else "comm_df"
      model_formula <- stats::as.formula(paste(response_name, "~", paste(rhs_terms, collapse = " + ")))
      set.seed(42L)
      if (identical(model_type, "dbrda")) {
        dist_obj <- tryCatch(
          {
            if (input$distance_metric %in% c("jaccard", "jsd")) {
              distance_ps <- phyloseq::phyloseq(
                phyloseq::otu_table(as.matrix(t(comm_df)), taxa_are_rows = TRUE)
              )
              phyloseq::distance(distance_ps, method = input$distance_metric)
            } else {
              dist_method <- switch(
                input$distance_metric,
                aitchison = "euclidean",
                "bray"
              )
              vegan::vegdist(comm_df, method = dist_method)
            }
          },
          error = function(e) {
            validate(need(FALSE, paste0("Distance calculation failed: ", conditionMessage(e))))
          }
        )
        model <- tryCatch(
          vegan::capscale(model_formula, data = model_df, comm = comm_df),
          error = function(e) {
            validate(need(FALSE, paste0("dbRDA model fitting failed: ", conditionMessage(e))))
          }
        )
      } else {
        model <- tryCatch(
          vegan::cca(model_formula, data = model_df),
          error = function(e) {
            validate(need(FALSE, paste0("CCA model fitting failed: ", conditionMessage(e))))
          }
        )
        metric_label <- "Not used (CCA)"
      }

      anova_all <- tryCatch(vegan::anova.cca(model, permutations = 999), error = function(e) NULL)
      p_all <- if (!is.null(anova_all) && nrow(as.data.frame(anova_all)) >= 1) as.data.frame(anova_all)[1, "Pr(>F)"] else NA_real_
      term_tests <- tryCatch(as.data.frame(vegan::anova.cca(model, by = "margin", permutations = 999)), error = function(e) data.frame())
      if (nrow(term_tests) > 0L) {
        term_tests$InternalTerm <- rownames(term_tests)
        reverse_names <- stats::setNames(names(internal_names), internal_names)
        term_tests$Variable <- unname(reverse_names[term_tests$InternalTerm])
        term_tests$Variable[is.na(term_tests$Variable)] <- term_tests$InternalTerm[is.na(term_tests$Variable)]
        term_tests$Type <- unname(variable_types[term_tests$Variable])
        rownames(term_tests) <- NULL
      }
      vif_values <- tryCatch(vegan::vif.cca(model), error = function(e) numeric(0))
      if (length(vif_values) > 0L) {
        vif_labels <- names(vif_values)
        for (variable_index in seq_along(internal_names)) {
          variable <- names(internal_names)[variable_index]
          internal <- unname(internal_names[variable_index])
          if (is.na(variable) || !nzchar(variable) || is.na(internal) || !nzchar(internal)) next
          matched <- startsWith(vif_labels, internal)
          vif_labels[matched] <- paste0(variable, substring(vif_labels[matched], nchar(internal) + 1L))
        }
        names(vif_values) <- vif_labels
      }

      site_scores <- vegan::scores(model, display = "sites", choices = 1:2, scaling = 2)
      species_scores <- tryCatch(vegan::scores(model, display = "species", choices = 1:2, scaling = 2), error = function(e) NULL)
      biplot_scores <- tryCatch(vegan::scores(model, display = "bp", choices = 1:2, scaling = 2), error = function(e) NULL)

      site_df <- as.data.frame(site_scores, stringsAsFactors = FALSE)
      validate(
        need(ncol(site_df) >= 2, "Ordination did not return 2D site scores.")
      )
      colnames(site_df)[1:2] <- c("Axis1", "Axis2")
      site_df$SampleID <- rownames(site_df)
      primary_internal <- internal_names[[input$group_var]]
      site_df$Group <- model_df[[primary_internal]]
      additional_centroid_df <- data.frame()
      categorical_additional <- setdiff(explanatory_vars, input$group_var)
      categorical_additional <- categorical_additional[variable_types[categorical_additional] == "Categorical"]
      if (length(categorical_additional) > 0L) {
        centroid_parts <- lapply(categorical_additional, function(variable) {
          level_values <- model_df[[internal_names[[variable]]]]
          data.frame(Axis1 = site_df$Axis1, Axis2 = site_df$Axis2, Level = level_values) %>%
            dplyr::group_by(Level) %>%
            dplyr::summarise(Axis1 = mean(Axis1), Axis2 = mean(Axis2), .groups = "drop") %>%
            dplyr::mutate(Variable = variable, Label = paste0(variable, ": ", Level))
        })
        additional_centroid_df <- dplyr::bind_rows(centroid_parts)
      }

      species_df <- data.frame()
      species_all_df <- data.frame()
      if (!is.null(species_scores)) {
        species_all_df <- as.data.frame(species_scores, stringsAsFactors = FALSE)
        if (ncol(species_all_df) >= 2) {
          colnames(species_all_df)[1:2] <- c("Axis1", "Axis2")
          species_all_df$Taxon <- rownames(species_all_df)
          species_all_df$Length <- sqrt(species_all_df$Axis1^2 + species_all_df$Axis2^2)
          species_all_df <- species_all_df[order(species_all_df$Length, decreasing = TRUE), , drop = FALSE]
          species_df <- species_all_df
          n_top <- as.integer(input$top_taxa_vectors)
          if (is.na(n_top) || n_top < 0) n_top <- 15
          if (n_top == 0) {
            species_df <- species_df[0, , drop = FALSE]
          } else if (nrow(species_df) > n_top) {
            species_df <- species_df[seq_len(n_top), , drop = FALSE]
          }
        } else {
          species_all_df <- data.frame()
          species_df <- data.frame()
        }
      }

      bp_df <- data.frame()
      if (!is.null(biplot_scores)) {
        bp_df <- as.data.frame(biplot_scores, stringsAsFactors = FALSE)
        if (ncol(bp_df) >= 2) {
          colnames(bp_df)[1:2] <- c("Axis1", "Axis2")
          bp_df$InternalTerm <- rownames(bp_df)
          reverse_names <- stats::setNames(names(internal_names), internal_names)
          bp_df$Variable <- unname(reverse_names[bp_df$InternalTerm])
          bp_df$Variable[is.na(bp_df$Variable)] <- bp_df$InternalTerm[is.na(bp_df$Variable)]
        } else {
          bp_df <- data.frame()
        }
      }

      environment_df <- data.frame()
      continuous_explanatory <- explanatory_vars[variable_types[explanatory_vars] == "Continuous"]
      if (nrow(bp_df) > 0L && length(continuous_explanatory) > 0L) {
        environment_df <- bp_df[bp_df$Variable %in% continuous_explanatory, , drop = FALSE]
        p_column <- intersect(c("Pr(>F)", "Pr(>Chi)"), colnames(term_tests))
        environment_df$p_value <- NA_real_
        if (length(p_column) > 0L && nrow(term_tests) > 0L) {
          p_lookup <- stats::setNames(term_tests[[p_column[1]]], term_tests$Variable)
          environment_df$p_value <- unname(p_lookup[environment_df$Variable])
        }
      }

      eig <- vegan::eigenvals(model)
      axis1_var <- if (length(eig) >= 1) round(100 * eig[1] / sum(eig), 1) else NA_real_
      axis2_var <- if (length(eig) >= 2) round(100 * eig[2] / sum(eig), 1) else NA_real_

      list(
        site_df = site_df,
        species_df = species_df,
        species_all_df = species_all_df,
        bp_df = bp_df,
        environment_df = environment_df,
        additional_centroid_df = additional_centroid_df,
        primary_is_categorical = identical(variable_types[[input$group_var]], "Categorical"),
        explanatory_vars = explanatory_vars,
        conditioning_vars = conditioning_vars,
        variable_types = variable_types,
        term_tests = term_tests,
        vif = vif_values,
        model_formula_label = paste(
          "Community ~",
          paste(c(explanatory_vars, paste0("Condition(", conditioning_vars, ")")), collapse = " + ")
        ),
        model_label = if (identical(model_type, "cca")) "CCA" else "dbRDA",
        model_type = model_type,
        metric_label = metric_label,
        transformation_label = transformation_label,
        permutation_count = 999L,
        permutation_seed = 42L,
        axis1_var = axis1_var,
        axis2_var = axis2_var,
        model_p = as.numeric(p_all),
        n_samples = nrow(site_df),
        n_samples_total = nrow(model_source),
        n_samples_excluded = sum(!keep_samples),
        excluded_samples = excluded_samples,
        n_taxa = nrow(otu_model),
        run_settings = isolate(current_analysis_settings())
      )
    }, ignoreNULL = FALSE)

    display_biplot_payload <- reactive({
      req(!is.null(input$run_biplot) && input$run_biplot > 0)
      payload <- biplot_payload()
      validate(
        need(
          identical(payload$run_settings, current_analysis_settings()),
          "Analysis settings changed after the last run. Click 'Run Ordination' to update all result tabs."
        )
      )
      payload
    })

    observeEvent(biplot_payload(), {
      payload <- biplot_payload()
      all_taxa <- if (nrow(payload$species_all_df) > 0) payload$species_all_df$Taxon else character(0)
      current_selected <- isolate(input$selected_taxa)
      selected_keep <- if (length(current_selected) > 0) intersect(current_selected, all_taxa) else character(0)
      updateSelectizeInput(session, "selected_taxa", choices = all_taxa, selected = selected_keep, server = TRUE)
    }, ignoreInit = FALSE)

    biplot_plot <- reactive({
      req(display_biplot_payload())
      payload <- display_biplot_payload()
      base_size <- input$base_size
      if (is.null(base_size) || !is.finite(base_size)) {
        base_size <- 14
      }
      dot_size <- input$dot_size
      if (is.null(dot_size) || !is.finite(dot_size)) {
        dot_size <- 5
      }
      text_size_taxa <- max(2.2, base_size / 3.6)
      text_size_group_vector <- max(2.3, base_size / 3.4)
      text_size_sample <- max(2.1, base_size / 4.0)
      site_df <- payload$site_df
      taxa_df <- payload$species_df
      if (length(input$selected_taxa) > 0 && nrow(payload$species_all_df) > 0) {
        taxa_df <- payload$species_all_df[payload$species_all_df$Taxon %in% input$selected_taxa, , drop = FALSE]
      }
      bp_df <- payload$bp_df
      environment_df <- payload$environment_df

      p <- ggplot2::ggplot(site_df, ggplot2::aes(x = Axis1, y = Axis2, color = Group)) +
        ggplot2::theme_bw(base_size = base_size) +
          plot_axis_label_theme(base_size, 14) +
        ggplot2::labs(
          title = if (identical(payload$model_type, "cca")) {
            paste0(input$tax_level, "-Level Constrained Ordination (CCA)")
          } else {
            paste0(input$tax_level, "-Level Constrained Ordination (dbRDA, ", payload$metric_label, ")")
          },
          x = paste0(if (identical(payload$model_type, "cca")) "CCA1" else "CAP1", if (is.finite(payload$axis1_var)) paste0(" (", payload$axis1_var, "%)") else ""),
          y = paste0(if (identical(payload$model_type, "cca")) "CCA2" else "CAP2", if (is.finite(payload$axis2_var)) paste0(" (", payload$axis2_var, "%)") else ""),
          color = input$group_var
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = base_size + 3, face = "bold"),
          axis.title.x = ggplot2::element_text(face = "bold", size = base_size + 1),
          axis.title.y = ggplot2::element_text(face = "bold", size = base_size + 1)
        )

      if (isTRUE(input$connect_subject_samples) && isTruthy(input$subject_id_var)) {
        metadata <- as.data.frame(phyloseq::sample_data(ps_obj()))
        p <- add_subject_pair_lines(p, site_df, metadata, input$subject_id_var, "Axis1", "Axis2")
      }

      centroid_df <- data.frame()
      if ((isTRUE(input$show_centroid_lines) || isTRUE(input$show_group_centroid)) && isTRUE(payload$primary_is_categorical)) {
        centroid_df <- site_df %>%
          dplyr::group_by(Group) %>%
          dplyr::summarise(
            CentroidAxis1 = mean(Axis1, na.rm = TRUE),
            CentroidAxis2 = mean(Axis2, na.rm = TRUE),
            .groups = "drop"
          )
      }
      if (isTRUE(input$show_group_hulls) && isTRUE(payload$primary_is_categorical)) {
        hull_parts <- lapply(split(site_df, as.character(site_df$Group)), function(group_df) {
          group_df <- group_df[
            is.finite(group_df$Axis1) & is.finite(group_df$Axis2),
            ,
            drop = FALSE
          ]
          if (nrow(group_df) < 3L || nrow(unique(group_df[, c("Axis1", "Axis2"), drop = FALSE])) < 3L) {
            return(NULL)
          }
          group_df[grDevices::chull(group_df$Axis1, group_df$Axis2), , drop = FALSE]
        })
        hull_parts <- hull_parts[!vapply(hull_parts, is.null, logical(1))]
        if (length(hull_parts) > 0L) {
          hull_df <- dplyr::bind_rows(hull_parts)
          p <- p + ggplot2::geom_polygon(
            data = hull_df,
            ggplot2::aes(x = Axis1, y = Axis2, group = Group, fill = Group, colour = Group),
            inherit.aes = FALSE,
            alpha = 0.14,
            linewidth = 0.5,
            show.legend = FALSE
          )
        }
      }
      if (isTRUE(input$show_ellipses) && isTRUE(payload$primary_is_categorical)) {
        p <- p + ggplot2::stat_ellipse(ggplot2::aes(group = Group), level = 0.80)
      }
      if (isTRUE(input$show_centroid_lines) && nrow(centroid_df) > 0L) {
        centroid_line_df <- dplyr::left_join(site_df, centroid_df, by = "Group")
        p <- p + ggplot2::geom_segment(
          data = centroid_line_df,
          ggplot2::aes(
            x = Axis1,
            y = Axis2,
            xend = CentroidAxis1,
            yend = CentroidAxis2,
            colour = Group
          ),
          inherit.aes = FALSE,
          linewidth = 0.3,
          alpha = 0.45,
          show.legend = FALSE
        )
      }
      if (isTRUE(input$show_dot_outline)) {
        p <- p + ggplot2::geom_point(
          ggplot2::aes(
            colour = ggplot2::stage(start = Group, after_scale = "#262626"),
            fill = ggplot2::after_scale(colour)
          ),
          shape = 21, size = dot_size, stroke = 0.7
        )
      } else {
        p <- p + ggplot2::geom_point(size = dot_size, alpha = 0.9)
      }

      if (isTRUE(input$show_group_centroid) && isTRUE(payload$primary_is_categorical)) {
        if (nrow(centroid_df) > 0) {
          p <- p + ggplot2::geom_label(
            data = centroid_df,
            mapping = ggplot2::aes(x = CentroidAxis1, y = CentroidAxis2, label = Group, colour = Group),
            fill = "white",
            text.colour = "#262626",
            fontface = "bold",
            size = text_size_sample,
            linewidth = 0.45,
            label.padding = grid::unit(0.18, "lines"),
            show.legend = FALSE,
            inherit.aes = FALSE
          )
        }
      }
      if (isTRUE(input$show_additional_centroids) && nrow(payload$additional_centroid_df) > 0L) {
        p <- p + ggplot2::geom_label(
          data = payload$additional_centroid_df,
          mapping = ggplot2::aes(x = Axis1, y = Axis2, label = Label),
          fill = "white",
          colour = "#0F766E",
          text.colour = "#262626",
          fontface = "bold",
          size = text_size_sample,
          linewidth = 0.45,
          label.padding = grid::unit(0.18, "lines"),
          show.legend = FALSE,
          inherit.aes = FALSE
        )
      }

      if (isTRUE(input$show_taxa_vectors) && nrow(taxa_df) > 0) {
        scale_mult <- 0.8
        taxa_repel_df <- data.frame(
          Axis1 = taxa_df$Axis1 * scale_mult,
          Axis2 = taxa_df$Axis2 * scale_mult,
          Taxon = taxa_df$Taxon,
          repel_point_size = 1,
          stringsAsFactors = FALSE
        )
        centroid_obstacles <- list()
        if (isTRUE(input$show_group_centroid) && nrow(centroid_df) > 0L) {
          centroid_obstacles[[length(centroid_obstacles) + 1L]] <- data.frame(
            Axis1 = centroid_df$CentroidAxis1,
            Axis2 = centroid_df$CentroidAxis2,
            Taxon = "",
            repel_point_size = 12,
            stringsAsFactors = FALSE
          )
        }
        if (isTRUE(input$show_additional_centroids) && nrow(payload$additional_centroid_df) > 0L) {
          centroid_obstacles[[length(centroid_obstacles) + 1L]] <- data.frame(
            Axis1 = payload$additional_centroid_df$Axis1,
            Axis2 = payload$additional_centroid_df$Axis2,
            Taxon = "",
            repel_point_size = 12,
            stringsAsFactors = FALSE
          )
        }
        if (length(centroid_obstacles) > 0L) {
          taxa_repel_df <- dplyr::bind_rows(taxa_repel_df, dplyr::bind_rows(centroid_obstacles))
        }
        p <- p +
          ggplot2::geom_segment(
            data = taxa_df,
            ggplot2::aes(x = 0, y = 0, xend = Axis1 * scale_mult, yend = Axis2 * scale_mult),
            inherit.aes = FALSE,
            color = "#7B1FA2",
            linewidth = 0.9,
            arrow = grid::arrow(length = grid::unit(0.22, "cm"))
          )
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_label_repel(
            data = taxa_repel_df,
            ggplot2::aes(
              x = Axis1,
              y = Axis2,
              label = Taxon,
              point.size = repel_point_size
            ),
            inherit.aes = FALSE,
            show.legend = FALSE,
            color = "#4A148C",
            fill = "white",
            size = text_size_taxa,
            fontface = "bold",
            label.size = 0.25,
            min.segment.length = 0,
            box.padding = 0.3,
            point.padding = 0.25,
            max.overlaps = Inf,
            seed = 1002
          )
        } else {
          p <- p + ggplot2::geom_label(
            data = taxa_df,
            ggplot2::aes(x = Axis1 * scale_mult, y = Axis2 * scale_mult, label = Taxon),
            inherit.aes = FALSE,
            color = "#4A148C",
            fill = "white",
            size = text_size_taxa,
            fontface = "bold",
            label.size = 0.25
          )
        }
      }

      if (isTRUE(input$show_environment_vectors) && nrow(environment_df) > 0L) {
        selected_environment <- input$environment_vector_vars
        if (length(selected_environment) > 0L) {
          environment_df <- environment_df[environment_df$Variable %in% selected_environment, , drop = FALSE]
        }
        if (isTRUE(input$environment_only_significant)) {
          cutoff <- suppressWarnings(as.numeric(input$environment_p_cutoff))
          if (!is.finite(cutoff)) cutoff <- 0.05
          environment_df <- environment_df[is.finite(environment_df$p_value) & environment_df$p_value <= cutoff, , drop = FALSE]
        }
        max_vectors <- suppressWarnings(as.integer(input$max_environment_vectors))
        if (!is.finite(max_vectors) || max_vectors < 1L) max_vectors <- 10L
        if (nrow(environment_df) > max_vectors) {
          environment_df <- environment_df[order(environment_df$p_value, na.last = TRUE), , drop = FALSE]
          environment_df <- environment_df[seq_len(max_vectors), , drop = FALSE]
        }
        if (nrow(environment_df) > 0L) {
          environment_df$Label <- environment_df$Variable
          if (isTRUE(input$show_environment_pvalues)) {
            p_text <- ifelse(
              is.finite(environment_df$p_value),
              paste0("p=", format(environment_df$p_value, digits = 3)),
              "p=NA"
            )
            environment_df$Label <- paste0(environment_df$Variable, " (", p_text, ")")
          }
          p <- p +
            ggplot2::geom_segment(
              data = environment_df,
              ggplot2::aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
              inherit.aes = FALSE,
              color = "#1B5E20",
              linewidth = 0.7,
              linetype = "dashed",
              arrow = grid::arrow(length = grid::unit(0.17, "cm"))
            ) +
            ggplot2::geom_text(
              data = environment_df,
              ggplot2::aes(x = Axis1, y = Axis2, label = Label),
              inherit.aes = FALSE,
              color = "#1B5E20",
              size = text_size_group_vector,
              vjust = -0.4
            )
        }
      }

      if (isTRUE(input$show_sample_names)) {
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_text_repel(
            data = site_df,
            ggplot2::aes(x = Axis1, y = Axis2, label = SampleID, color = Group),
            size = text_size_sample,
            show.legend = FALSE,
            box.padding = 0.2,
            point.padding = 0.15
          )
        } else {
          p <- p + ggplot2::geom_text(
            data = site_df,
            ggplot2::aes(x = Axis1, y = Axis2, label = SampleID, color = Group),
            size = text_size_sample,
            hjust = -0.1,
            show.legend = FALSE
          )
        }
      }

      p
    })

    output$biplot_plot <- renderPlot(
      {
        if (is.null(input$run_biplot) || input$run_biplot < 1) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Click 'Run Ordination' to start analysis.", cex = 0.85)
          return(invisible(NULL))
        }
        if (isTRUE(biplot_running())) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Constrained ordination is running. Please wait...", cex = 0.85)
          return(invisible(NULL))
        }
        biplot_plot()
      },
      width = function() input$plot_width,
      height = function() input$plot_height
    )

    output$biplot_status <- renderText({
      req(display_biplot_payload())
      payload <- display_biplot_payload()
      pval_txt <- if (is.finite(payload$model_p)) format(payload$model_p, digits = 3) else "NA"
      term_lines <- character(0)
      if (nrow(payload$term_tests) > 0L) {
        p_column <- intersect(c("Pr(>F)", "Pr(>Chi)"), colnames(payload$term_tests))
        if (length(p_column) > 0L) {
          keep <- payload$term_tests$Variable %in% payload$explanatory_vars
          term_rows <- payload$term_tests[keep, , drop = FALSE]
          term_lines <- vapply(seq_len(nrow(term_rows)), function(i) {
            paste0(
              "  ", term_rows$Variable[i], " [", term_rows$Type[i], "]: p=",
              format(term_rows[[p_column[1]]][i], digits = 3)
            )
          }, character(1))
        }
      }
      vif_lines <- character(0)
      if (length(payload$vif) > 0L) {
        vif_lines <- paste0("  ", names(payload$vif), ": ", format(payload$vif, digits = 3))
      }
      taxa_vectors_shown <- if (length(input$selected_taxa) > 0 && nrow(payload$species_all_df) > 0) {
        sum(payload$species_all_df$Taxon %in% input$selected_taxa)
      } else {
        nrow(payload$species_df)
      }
      continuous_model_vars <- names(payload$variable_types)[payload$variable_types == "Continuous"]
      standardization_line <- if (
        isTRUE(input$standardize_continuous) && length(continuous_model_vars) > 0L
      ) {
        paste0(
          "Continuous-variable standardization: mean 0 and standard deviation 1 (",
          paste(continuous_model_vars, collapse = ", "),
          ")"
        )
      } else {
        NULL
      }
      paste(
        c(
          paste0("Model: ", payload$model_label),
          paste0("Distance metric: ", payload$metric_label),
          paste0("Transformation: ", payload$transformation_label),
          standardization_line,
          paste0("Formula: ", payload$model_formula_label),
          paste0("Primary variable: ", input$group_var),
          paste0("Samples used: ", payload$n_samples, " / ", payload$n_samples_total),
          paste0("Samples excluded for missing values: ", payload$n_samples_excluded),
          paste0("Taxa used: ", payload$n_taxa),
          paste0("Taxa scores shown: ", taxa_vectors_shown),
          paste0(
            "Permutation setup: ",
            payload$permutation_count,
            " permutations (initial seed ",
            payload$permutation_seed,
            ")"
          ),
          paste0(payload$model_label, " model p-value (perm): ", pval_txt),
          if (length(term_lines) > 0L) c("Marginal term tests:", term_lines) else NULL,
          if (length(vif_lines) > 0L) c("Collinearity (VIF):", vif_lines) else NULL
        ),
        collapse = "\n"
      )
    })

    output$biplot_status_box <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 700 else input$plot_width
      tags$div(
        class = "app-result-card app-status-card",
        style = paste0("width: ", box_width, "px; max-width: 100%; height: 280px;"),
        verbatimTextOutput(session$ns("biplot_status"))
      )
    })

    output$biplot_term_tests <- renderTable({
      payload <- display_biplot_payload()
      tests <- payload$term_tests
      validate(need(nrow(tests) > 0L, "Marginal term tests are unavailable for this model."))
      tests <- tests[tests$Variable %in% payload$explanatory_vars, , drop = FALSE]
      statistic_columns <- intersect(c("Df", "Variance", "SumOfSqs", "ChiSquare", "F", "Pr(>F)", "Pr(>Chi)"), colnames(tests))
      tests[, c("Variable", "Type", statistic_columns), drop = FALSE]
    }, striped = TRUE, bordered = TRUE, spacing = "s", digits = 4)

    output$biplot_vif_table <- renderTable({
      payload <- display_biplot_payload()
      validate(need(length(payload$vif) > 0L, "VIF values are unavailable for this model."))
      data.frame(
        ModelTerm = names(payload$vif),
        VIF = as.numeric(payload$vif),
        Assessment = ifelse(payload$vif >= 10, "High", ifelse(payload$vif >= 5, "Elevated", "Acceptable")),
        check.names = FALSE
      )
    }, striped = TRUE, bordered = TRUE, spacing = "s", digits = 3)

    output$biplot_exclusions <- renderTable({
      payload <- display_biplot_payload()
      validate(need(nrow(payload$excluded_samples) > 0L, "No samples were excluded for missing model variables."))
      payload$excluded_samples
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$biplot_legend_box <- renderUI({
      req_plot_legend_tab(session, "Constrained Ordination")
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 700 else input$plot_width
      tags$div(
        style = sprintf("width: %spx;", box_width),
        tags$div(
          style = "color: #1f2937; font-size: 12.5px; line-height: 1.55;",
          biplot_figure_legend()
        )
      )
    })
    outputOptions(output, "biplot_legend_box", suspendWhenHidden = FALSE, priority = -10)

    output$biplot_status_separator <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 700 else input$plot_width
      tags$hr(
        style = paste(
          sprintf("width: %spx;", box_width),
          "max-width: 100%;",
          "margin: 14px 0 12px 0;",
          "border-top: 1px solid #d1d5db;"
        )
      )
    })

    biplot_figure_legend <- function() {
      payload <- display_biplot_payload()
      req(input$group_var, input$tax_level, input$analysis_method, input$distance_metric)
      distance_label <- switch(
        input$distance_metric,
        aitchison = "Aitchison",
        jaccard = "Jaccard",
        jsd = "Jensen-Shannon",
        "Bray-Curtis"
      )
      taxa_scores_shown <- isTRUE(input$show_taxa_vectors) && if (
        length(input$selected_taxa) > 0L && nrow(payload$species_all_df) > 0L
      ) {
        any(payload$species_all_df$Taxon %in% input$selected_taxa)
      } else {
        nrow(payload$species_df) > 0L
      }
      group_hulls_shown <- isTRUE(input$show_group_hulls) &&
        isTRUE(payload$primary_is_categorical) &&
        any(table(payload$site_df$Group) >= 3L)
      continuous_model_vars <- names(payload$variable_types)[payload$variable_types == "Continuous"]
      standardization_sentence <- if (
        isTRUE(input$standardize_continuous) && length(continuous_model_vars) > 0L
      ) {
        paste0(
          " Continuous model variables (",
          paste(continuous_model_vars, collapse = ", "),
          ") were standardized to mean 0 and standard deviation 1 before model fitting; ",
          "categorical variables were not standardized."
        )
      } else {
        ""
      }
      body_text <- if (identical(input$analysis_method, "cca")) {
        paste0(
          "CCA constrained ordination projects sample relationships on CCA1 and CCA2 at ",
          tolower(input$tax_level),
          " level. Distance metric is not used in CCA. Input abundances are transformed using ",
          payload$transformation_label,
          ". Points are samples colored by ",
          input$group_var,
          ". ",
          if (taxa_scores_shown) {
            paste0(
              "Purple arrows represent taxa scores from the fitted constrained ordination; ",
              "their directions and relative lengths describe associations with the displayed axes, ",
              "and individual taxa p-values are not calculated. "
            )
          } else "",
          if (isTRUE(input$show_environment_vectors)) "Green dashed arrows show fitted continuous environmental variables. " else "",
          if (isTRUE(input$show_ellipses) && isTRUE(payload$primary_is_categorical)) "Ellipses show the fitted 80% contours for each group. " else "",
          if (group_hulls_shown) "Shaded polygons show convex hulls for primary-group samples. " else "",
          if (isTRUE(input$show_centroid_lines) && isTRUE(payload$primary_is_categorical)) "Lines connect samples to their primary-group centroids. " else "",
          if (isTRUE(input$connect_subject_samples) && isTruthy(input$subject_id_var)) "Thin gray lines connect pairs of samples from the same subject. " else "",
          if (isTRUE(input$show_group_centroid) && isTRUE(payload$primary_is_categorical)) "Boxed labels indicate primary-group centroids. " else "",
          if (isTRUE(input$show_sample_names)) "Sample labels are displayed near points." else ""
        )
      } else {
        paste0(
          "dbRDA constrained ordination projects sample relationships on CAP1 and CAP2 using ",
          distance_label,
          " distance at ",
          tolower(input$tax_level),
          " level. Input abundances are transformed using ",
          payload$transformation_label,
          ". Points are samples colored by ",
          input$group_var,
          ". ",
          if (taxa_scores_shown) {
            paste0(
              "Purple arrows represent taxa scores from the fitted constrained ordination; ",
              "their directions and relative lengths describe associations with the displayed axes, ",
              "and individual taxa p-values are not calculated. "
            )
          } else "",
          if (isTRUE(input$show_environment_vectors)) "Green dashed arrows show fitted continuous environmental variables. " else "",
          if (isTRUE(input$show_ellipses) && isTRUE(payload$primary_is_categorical)) "Ellipses show the fitted 80% contours for each group. " else "",
          if (group_hulls_shown) "Shaded polygons show convex hulls for primary-group samples. " else "",
          if (isTRUE(input$show_centroid_lines) && isTRUE(payload$primary_is_categorical)) "Lines connect samples to their primary-group centroids. " else "",
          if (isTRUE(input$connect_subject_samples) && isTruthy(input$subject_id_var)) "Thin gray lines connect pairs of samples from the same subject. " else "",
          if (isTRUE(input$show_group_centroid) && isTRUE(payload$primary_is_categorical)) "Boxed labels indicate primary-group centroids. " else "",
          if (isTRUE(input$show_sample_names)) "Sample labels are displayed near points." else ""
        )
      }
      body_text <- paste0(
        body_text,
        standardization_sentence,
        " Permutation-based tests use ",
        payload$permutation_count,
        " permutations with initial seed ",
        payload$permutation_seed,
        "."
      )
      tags$div(
        tags$div(style = "font-weight: 600; margin-bottom: 4px;", "Constrained ordination"),
        tags$div(body_text)
      )
    }

    output$download_biplot_png <- downloadHandler(
      filename = function() paste0("constrained_ordination_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      contentType = "image/png",
      content = function(file) {
        req(!is.null(input$run_biplot) && input$run_biplot > 0)
        req(display_biplot_payload())
        dpi_val <- 300
        plot_width <- suppressWarnings(as.numeric(input$plot_width))
        plot_height <- suppressWarnings(as.numeric(input$plot_height))
        if (!is.finite(plot_width) || is.na(plot_width) || plot_width <= 0) plot_width <- 700
        if (!is.finite(plot_height) || is.na(plot_height) || plot_height <= 0) plot_height <- 500
        width_in <- plot_width / 72
        height_in <- plot_height / 72
        p <- biplot_plot()
        ggplot2::ggsave(file, plot = p, device = "png", width = width_in, height = height_in, units = "in", dpi = dpi_val)
      }
    )

    output$download_biplot_scores <- downloadHandler(
      filename = function() paste0("constrained_ordination_scores_", Sys.Date(), ".tsv"),
      content = function(file) {
        req(display_biplot_payload())
        payload <- display_biplot_payload()
        site_out <- payload$site_df %>%
          dplyr::mutate(Type = "Sample", Label = SampleID, Group = as.character(Group)) %>%
          dplyr::select(Type, Label, Axis1, Axis2, Group)
        taxa_out <- if (nrow(payload$species_df) > 0) {
          payload$species_df %>% dplyr::mutate(Type = "Taxa", Label = Taxon, Group = NA_character_) %>% dplyr::select(Type, Label, Axis1, Axis2, Group)
        } else {
          data.frame(Type = character(0), Label = character(0), Axis1 = numeric(0), Axis2 = numeric(0), Group = character(0))
        }
        environment_out <- if (nrow(payload$environment_df) > 0) {
          payload$environment_df %>%
            dplyr::mutate(Type = "Environment", Label = Variable, Group = NA_character_) %>%
            dplyr::select(Type, Label, Axis1, Axis2, Group)
        } else {
          data.frame(Type = character(0), Label = character(0), Axis1 = numeric(0), Axis2 = numeric(0), Group = character(0))
        }
        out_df <- dplyr::bind_rows(site_out, taxa_out, environment_out)
        utils::write.table(out_df, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    )

  })
}
