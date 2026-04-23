## UI
mod_sparcc_ui <- function(id) {
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
      .well h4 { font-size: 16px; }
      .well h5 { font-size: 13px; }
      .well .control-label { font-size: 12px; }
      .well .checkbox label { font-size: 12px; }
      .well .form-control { font-size: 12px; }
      .well .btn { font-size: 11px; }
    ")),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        h4(icon("diagram-project"), "SparCC"),
        hr(),
        selectInput(ns("group_var"), "1. Group panel", choices = c("All"), selected = "All"),
        selectizeInput(
          ns("group_levels"),
          "2. Group levels",
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
          "3. Analysis mode",
          choices = c("Single network", "Compare two groups"),
          selected = "Single network"
        ),
        selectInput(ns("tax_level"), "4. Taxonomic level", choices = c("ASV", "Genus", "Species"), selected = "Genus"),
        numericInput(ns("prevalence_filter_pct"), "5. Prevalence filter (%)", value = 10, min = 0, max = 100, step = 1),
        selectInput(ns("node_size_by"), "6. Node size", choices = c("Connectivity", "Abundance"), selected = "Connectivity"),
        selectInput(ns("node_color_by"), "7. Node color", choices = c("None"), selected = "None"),
        tags$details(
          style = "margin-bottom: 8px;",
          tags$summary("Advanced Options"),
          numericInput(ns("min_edge_weight"), "8. Minimum absolute edge weight", value = 0.1, min = 0, max = 1, step = 0.01),
          numericInput(ns("max_taxa"), "9. Max taxa for network", value = 200, min = 20, max = 2000, step = 10),
          numericInput(ns("seed"), "10. Seed", value = 1001, min = 1, step = 1)
        ),
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Plot Dimensions"),
        numericInput(ns("plot_width"), "Plot width (px)", value = 800, min = 400, max = 2400, step = 50),
        numericInput(ns("plot_height"), "Plot height (px)", value = 550, min = 300, max = 2400, step = 50),
        tags$div(
          style = "display: flex; align-items: center; gap: 8px; flex-wrap: wrap;",
          actionButton(ns("run_network_btn"), "Run SparCC", class = "btn-danger", style = "font-size: 12px;"),
          tags$span("May take a long time.", style = "font-size: 11px; color: #b94a48;")
        ),
        tags$script(HTML(
          "Shiny.addCustomMessageHandler('toggle-sparcc-run-btn', function(msg) {
             var btn = document.getElementById(msg.id);
             if (!btn) return;
             btn.disabled = !!msg.disabled;
             if (msg.label) btn.textContent = msg.label;
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
        h4("SparCC Network"),
        tags$div(
          id = ns("sparcc_tab_container"),
          style = "max-width: 100%;",
          tabsetPanel(
            id = ns("sparcc_active_tab"),
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
            DTOutput(ns("edge_table_all"))
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
          tabPanel("Differential Edges", DTOutput(ns("comparison_edge_table"))),
          tabPanel("Comparison Summary", verbatimTextOutput(ns("comparison_summary"))),
          tabPanel("Hub Table", DTOutput(ns("hub_table")))
          )
        ),
        uiOutput(ns("sparcc_legend_box")),
        uiOutput(ns("sparcc_status_separator")),
        h5(icon("circle-info"), "SparCC Status"),
        uiOutput(ns("sparcc_status_box"))
      )
    )
  )
}

## Server
mod_sparcc_server <- function(id, ps_obj) {
  moduleServer(id, function(input, output, session) {
    draw_wait_message <- function(message_text) {
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      graphics::par(mar = c(0, 0, 0, 0))
      graphics::plot.new()
      graphics::text(0.5, 0.5, message_text)
    }

    observe({
      width_px <- suppressWarnings(as.integer(input$plot_width))
      if (!is.finite(width_px) || is.na(width_px) || width_px <= 0) width_px <- 800L
      session$sendCustomMessage(
        "set-tab-container-width",
        list(id = session$ns("sparcc_tab_container"), width = paste0(width_px, "px"))
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
      parent_candidates <- character(0)
      if (!is.na(rank_pos) && rank_pos > 1) {
        parent_candidates <- rev(tax_ranks[seq_len(rank_pos - 1)])
        parent_candidates <- parent_candidates[parent_candidates %in% tax_cols]
      }
      parent_val <- rep("UnclassifiedParent", nrow(tt))
      if (length(parent_candidates) > 0) {
        for (parent_rank in parent_candidates) {
          candidate_val <- as.character(tt[[parent_rank]])
          candidate_norm <- tolower(trimws(candidate_val))
          candidate_norm <- gsub("^[a-z]__", "", candidate_norm)
          candidate_is_placeholder <- is.na(candidate_norm) | !nzchar(candidate_norm) | grepl("(uncultured|unassigned)", candidate_norm)
          use_idx <- idx_placeholder & parent_val == "UnclassifiedParent" & !candidate_is_placeholder
          parent_val[use_idx] <- candidate_val[use_idx]
        }
      }
      tt[[tax_level]][idx_placeholder] <- parent_val[idx_placeholder]
      phyloseq::tax_table(ps) <- phyloseq::tax_table(as.matrix(tt))
      ps
    }

    group_levels_choices_cache <- reactiveVal(character(0))

    observeEvent(list(ps_obj(), input$group_var), {
      req(ps_obj())
      ps <- ps_obj()

      group_choices <- "All"
      if (!is.null(phyloseq::sample_data(ps, errorIfNULL = FALSE))) {
        group_choices <- c("All", colnames(as.data.frame(phyloseq::sample_data(ps))))
      }
      selected_group <- if (!is.null(input$group_var) && input$group_var %in% group_choices) input$group_var else "All"
      updateSelectInput(session, "group_var", choices = group_choices, selected = selected_group)

      level_choices <- character(0)
      if (selected_group != "All" && !is.null(phyloseq::sample_data(ps, errorIfNULL = FALSE))) {
        sd_df <- as.data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE)
        if (selected_group %in% colnames(sd_df)) {
          level_choices <- unique(as.character(sd_df[[selected_group]]))
          level_choices <- sort(level_choices[!is.na(level_choices) & level_choices != ""])
        }
      }
      selected_levels <- input$group_levels
      if (is.null(selected_levels)) selected_levels <- character(0)
      selected_levels <- intersect(selected_levels, level_choices)

      old_choices <- group_levels_choices_cache()
      old_selected <- input$group_levels
      if (is.null(old_selected)) old_selected <- character(0)
      choices_changed <- !identical(old_choices, level_choices)
      selection_changed <- !setequal(selected_levels, old_selected)
      if (choices_changed || selection_changed) {
        freezeReactiveValue(input, "group_levels")
        updateSelectizeInput(session, "group_levels", choices = level_choices, selected = selected_levels, server = TRUE)
        group_levels_choices_cache(level_choices)
      }
    }, ignoreInit = FALSE)

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

    network_result <- eventReactive(input$run_network_btn, {
      req(build_network_inputs())
      validate(
        need(requireNamespace("NetCoMi", quietly = TRUE), "NetCoMi package is not installed. Please install 'NetCoMi'."),
        need(requireNamespace("igraph", quietly = TRUE), "igraph package is not installed. Please install 'igraph'.")
      )
      network_running(TRUE)

      session$sendCustomMessage("toggle-sparcc-run-btn", list(id = session$ns("run_network_btn"), disabled = TRUE, label = "Running..."))
      on.exit({
        network_running(FALSE)
        session$sendCustomMessage("toggle-sparcc-run-btn", list(id = session$ns("run_network_btn"), disabled = FALSE, label = "Run SparCC"))
      }, add = TRUE)

      set.seed(as.integer(input$seed))
      built <- build_network_inputs()
      otu_mat <- built$otu_mat
      sample_df <- built$sample_df
      group_var <- input$group_var
      group_var_resolved <- group_var
      if (!is.null(group_var) && group_var != "All" && !group_var %in% colnames(sample_df)) {
        matched_col <- colnames(sample_df)[make.names(colnames(sample_df)) == make.names(group_var)]
        if (length(matched_col) >= 1) {
          group_var_resolved <- matched_col[1]
        } else {
          group_var_resolved <- "All"
        }
      }

      result <- tryCatch({
        withProgress(message = "Running NetCoMi SparCC...", value = 0, {
          sample_ids <- colnames(otu_mat)
          if (!is.null(group_var_resolved) && group_var_resolved != "All" && group_var_resolved %in% colnames(sample_df)) {
            group_values <- as.character(sample_df[[group_var_resolved]])
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
          if (!is.null(group_var_resolved) && group_var_resolved != "All") {
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
              "Group variable (input/resolved): ", as.character(group_var), " / ", as.character(group_var_resolved)
            )
            validate(need(length(matched_levels) > 0, debug_msg))
            group_samples <- group_samples[matched_levels]
          }
          validate(need(length(group_samples) > 0, "No group level is available for network estimation."))
          if (identical(analysis_mode, "Compare two groups")) {
            validate(need(length(group_samples) == 2, "Compare mode requires exactly 2 eligible groups."))
          }

          build_group_network <- function(group_name, ids) {
            sub_otu <- otu_mat[, ids, drop = FALSE]
            net_obj <- NetCoMi::netConstruct(
              data = t(sub_otu),
              dataType = "counts",
              measure = "sparcc",
              measurePar = list(ncpus = 4),
              verbose = 0
            )
            adja <- net_obj$adjaMat1
            if (is.null(adja)) adja <- net_obj$adjaMat
            validate(need(!is.null(adja), paste("NetCoMi returned no adjacency matrix for SparCC group:", group_name)))

            adja_mat <- as.matrix(adja)
            assoc <- net_obj$assoMat1
            if (is.null(assoc)) assoc <- net_obj$assoMat
            if (!is.null(assoc)) assoc <- as.matrix(assoc)

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
            edge_df$group <- group_name

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
            group_var = if (is.null(group_var) || group_var == "All") "All" else group_var,
            group_status = status_df,
            taxa_annotation = built$taxa_annotation
          )
        })
      }, error = function(e) {
        showNotification(paste("SparCC execution failed:", e$message), type = "error", duration = NULL)
        NULL
      })

      validate(need(!is.null(result), "Network analysis failed. Check the error message above."))
      result
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
      selected_group_var <- if (is.null(input$group_var)) "All" else input$group_var
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
        return(list(hub_taxa = character(0), hub_scores = numeric(0), modularity_netanalyze = NA_real_))
      }

      extract_hubs <- function(obj) {
        if (is.null(obj)) return(character(0))
        if (is.character(obj)) return(obj[nzchar(obj)])
        if (is.factor(obj)) return(as.character(obj))
        if (is.list(obj)) {
          nms <- names(obj)
          if (!is.null(nms)) {
            idx <- which(grepl("hub", tolower(nms)))
            for (i in idx) {
              v <- extract_hubs(obj[[i]])
              if (length(v) > 0) return(v)
            }
          }
          for (el in obj) {
            v <- extract_hubs(el)
            if (length(v) > 0) return(v)
          }
        }
        character(0)
      }

      extract_modularity <- function(obj) {
        if (is.null(obj)) return(NA_real_)
        if (is.numeric(obj) && length(obj) == 1 && is.finite(obj)) return(as.numeric(obj))
        if (is.list(obj)) {
          nms <- names(obj)
          if (!is.null(nms)) {
            idx <- which(grepl("mod", tolower(nms)))
            for (i in idx) {
              v <- extract_modularity(obj[[i]])
              if (is.finite(v)) return(v)
            }
          }
          for (el in obj) {
            v <- extract_modularity(el)
            if (is.finite(v)) return(v)
          }
        }
        NA_real_
      }

      extract_named_numeric <- function(obj, key_pattern = NULL) {
        out <- list()
        walk_obj <- function(x, nm = "") {
          if (is.null(x)) return()
          if (is.numeric(x) && !is.null(names(x))) {
            if (is.null(key_pattern) || grepl(key_pattern, tolower(nm))) {
              out[[length(out) + 1]] <<- x
            }
            return()
          }
          if (is.list(x)) {
            nms <- names(x)
            if (is.null(nms)) nms <- rep("", length(x))
            for (i in seq_along(x)) walk_obj(x[[i]], nms[[i]])
          }
        }
        walk_obj(obj)
        out
      }

      mod_netanalyze <- NA_real_
      hub_taxa <- character(0)
      hub_scores <- numeric(0)
      if (!is.null(net_res$net_obj)) {
        ana <- tryCatch(
          suppressWarnings(NetCoMi::netAnalyze(net_res$net_obj, verbose = 0)),
          error = function(e) NULL
        )
        hub_taxa <- extract_hubs(ana)
        mod_netanalyze <- extract_modularity(ana)

        if (!is.null(ana)) {
          vertex_names <- igraph::V(g)$name
          deg_cands <- extract_named_numeric(ana, "degree")
          btw_cands <- extract_named_numeric(ana, "betwe")
          eig_cands <- extract_named_numeric(ana, "eigen")
          pick_vec <- function(cands) {
            if (length(cands) == 0) return(NULL)
            for (v in cands) {
              nm <- names(v)
              if (!is.null(nm) && all(vertex_names %in% nm)) return(as.numeric(v[vertex_names]))
            }
            NULL
          }
          deg_v <- pick_vec(deg_cands)
          btw_v <- pick_vec(btw_cands)
          eig_v <- pick_vec(eig_cands)
          if (!is.null(deg_v) && !is.null(btw_v) && !is.null(eig_v)) {
            safe_z <- function(x) {
              sx <- stats::sd(x, na.rm = TRUE)
              if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
              as.numeric((x - mean(x, na.rm = TRUE)) / sx)
            }
            score <- 0.4 * safe_z(deg_v) + 0.3 * safe_z(btw_v) + 0.3 * safe_z(eig_v)
            names(score) <- vertex_names
            hub_scores <- sort(score, decreasing = TRUE)
          }
        }
      }

      if (length(hub_taxa) == 0) {
        deg <- igraph::degree(g)
        hub_taxa <- names(sort(as.numeric(deg), decreasing = TRUE))
      }
      hub_taxa <- unique(as.character(hub_taxa))
      hub_taxa <- hub_taxa[nzchar(hub_taxa)]
      hub_taxa <- hub_taxa[hub_taxa %in% igraph::V(g)$name]
      if (length(hub_scores) > 0) {
        hub_taxa <- names(hub_scores)
      }
      if (length(hub_taxa) > top_n) hub_taxa <- hub_taxa[seq_len(top_n)]
      if (length(hub_scores) > 0) {
        hub_scores <- hub_scores[hub_taxa]
      }

      list(hub_taxa = hub_taxa, hub_scores = hub_scores, modularity_netanalyze = mod_netanalyze)
    }

    comparison_metrics <- reactive({
      req(network_result())
      res <- network_result()
      if (!identical(input$analysis_mode, "Compare two groups")) return(NULL)
      if (length(res$networks) != 2) return(NULL)
      groups <- names(res$networks)
      g1 <- res$networks[[groups[1]]]$graph
      g2 <- res$networks[[groups[2]]]$graph
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
      req(network_result())
      res <- network_result()
      if (!identical(input$analysis_mode, "Compare two groups")) return(data.frame())
      if (length(res$networks) != 2) return(data.frame())
      groups <- names(res$networks)
      g1_name <- groups[1]
      g2_name <- groups[2]
      e1 <- igraph::as_data_frame(res$networks[[g1_name]]$graph, what = "edges")
      e2 <- igraph::as_data_frame(res$networks[[g2_name]]$graph, what = "edges")
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
      req(network_result())
      res <- network_result()
      if (!identical(input$analysis_mode, "Compare two groups")) {
        return("Comparison mode is disabled. Select 'Compare two groups' to view results.")
      }
      if (length(res$networks) != 2) {
        return("Comparison requires exactly 2 successfully estimated groups.")
      }
      groups <- names(res$networks)
      metrics <- comparison_metrics()
      s1 <- summarize_network_structure(res$networks[[groups[1]]])
      s2 <- summarize_network_structure(res$networks[[groups[2]]])
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
        paste0("- Modularity (netAnalyze, ", groups[1], "): ", if (is.finite(s1$modularity_netanalyze)) sprintf("%.3f", s1$modularity_netanalyze) else "NA"),
        paste0("- Modularity (netAnalyze, ", groups[2], "): ", if (is.finite(s2$modularity_netanalyze)) sprintf("%.3f", s2$modularity_netanalyze) else "NA"),
        paste0("- Delta modularity (netAnalyze, ", groups[1], " - ", groups[2], "): ", if (is.finite(delta_modularity_netanalyze)) sprintf("%.3f", delta_modularity_netanalyze) else "NA")
      )
      paste(lines, collapse = "\n")
    })

    hub_table <- reactive({
      req(network_result())
      res <- network_result()
      rows <- list()
      for (gname in names(res$networks)) {
        s <- summarize_network_structure(res$networks[[gname]])
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
      DT::datatable(
        hub_table(),
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE)
      )
    })

    draw_comparison_network_plot <- function(edge_tbl, edge_palette) {
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
          seed = as.integer(input$seed)
        ) +
        ggplot2::scale_color_manual(
          values = edge_palette,
          breaks = names(edge_palette),
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
        ggplot2::ggtitle(paste0(input$tax_level, "-Level SparCC Differential Network (Two-Group Comparison)"))
    }

    output$comparison_network_plot <- renderPlot({
      if (is.null(input$run_network_btn) || input$run_network_btn < 1) {
        draw_wait_message("Click 'Run SparCC' to start analysis.")
        return(invisible(NULL))
      }
      if (isTRUE(network_running())) {
        draw_wait_message("SparCC analysis is running. Please wait...")
        return(invisible(NULL))
      }
      req(network_result())
      validate(need(identical(input$analysis_mode, "Compare two groups"), "Enable 'Compare two groups' to view differential network."))
      edge_tbl <- comparison_edge_table()
      validate(need(nrow(edge_tbl) > 0, "No differential edge is available for comparison plotting."))
      status_vals <- unique(edge_tbl$status)
      only_status <- status_vals[grepl("^Only ", status_vals)]
      palette <- c("Shared" = "#7F8C8D")
      if (length(only_status) >= 1) palette[only_status[1]] <- "#1F78B4"
      if (length(only_status) >= 2) palette[only_status[2]] <- "#E31A1C"
      draw_comparison_network_plot(edge_tbl, palette)
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
      filename = function() paste0("sparcc_differential_network_", Sys.Date(), ".png"),
      content = function(file) {
        req(network_result())
        validate(need(identical(input$analysis_mode, "Compare two groups"), "Enable 'Compare two groups' to download differential network."))
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

        grDevices::png(file, width = plot_width * 2L, height = plot_height * 2L, res = 120)
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_comparison_network_plot(edge_tbl, palette)
      }
    )
    outputOptions(output, "comparison_network_plot", suspendWhenHidden = FALSE)
    outputOptions(output, "download_comparison_network_plot", suspendWhenHidden = FALSE)


    output$network_summary <- renderText({
      req(network_result())
      res <- network_result()
      lines <- c(
        "SparCC Network Summary",
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
        rank_prefix_map <- c(
          Kingdom = "k__",
          Phylum = "p__",
          Class = "c__",
          Order = "o__",
          Family = "f__",
          Genus = "g__",
          Species = "s__"
        )
        get_rank_unassigned_label <- function(rank_name) {
          prefix <- rank_prefix_map[[rank_name]]
          if (is.null(prefix) || !nzchar(prefix)) {
            return(paste0(rank_name, "__Unassigned"))
          }
          paste0(prefix, "Unassigned")
        }
        color_map <- res$taxa_annotation[, c("node_name", color_var), drop = FALSE]
        color_map[[color_var]] <- as.character(color_map[[color_var]])
        color_map[[color_var]][is.na(color_map[[color_var]]) | !nzchar(color_map[[color_var]])] <- get_rank_unassigned_label(color_var)
        color_map <- color_map[!duplicated(color_map$node_name), , drop = FALSE]
        node_plot_df <- dplyr::left_join(node_plot_df, color_map, by = c("name" = "node_name"))
        node_plot_df$node_color_group <- as.character(node_plot_df[[color_var]])
        node_plot_df$node_color_group[is.na(node_plot_df$node_color_group) | !nzchar(node_plot_df$node_color_group)] <- get_rank_unassigned_label(color_var)
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
          seed = as.integer(input$seed)
        ) +
        ggplot2::facet_wrap(~group, scales = "free") +
        ggplot2::scale_color_manual(
          values = c("Positive" = "#4575B4", "Negative" = "#D73027", "Zero" = "#999999"),
          breaks = c("Positive", "Negative", "Zero"),
          drop = FALSE,
          name = "Edge sign"
        ) +
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
            paste0("NetCoMi SparCC ", input$tax_level, "-Level Network (Connected)")
          } else {
            paste0("NetCoMi SparCC ", input$tax_level, "-Level Network (All)")
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
          drop = FALSE,
          guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 6))
        )
      }

      print(p)
    }

    edge_table_all <- reactive({
      req(network_result())
      edge_tbl <- network_result()$edge_table
      min_edge_weight <- suppressWarnings(as.numeric(input$min_edge_weight))
      if (is.na(min_edge_weight) || min_edge_weight < 0) min_edge_weight <- 0.1
      if (!"weight" %in% colnames(edge_tbl)) return(edge_tbl)
      edge_tbl <- edge_tbl[is.finite(edge_tbl$weight) & abs(edge_tbl$weight) >= min_edge_weight, , drop = FALSE]
      edge_tbl
    })

    output$network_plot_all <- renderPlot({
      if (is.null(input$run_network_btn) || input$run_network_btn < 1) {
        draw_wait_message("Click 'Run SparCC' to start analysis.")
        return(invisible(NULL))
      }
      if (isTRUE(network_running())) {
        draw_wait_message("SparCC analysis is running. Please wait...")
        return(invisible(NULL))
      }
      req(network_result())
      draw_network_plot(network_result(), connected_only = FALSE)
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
        draw_wait_message("Click 'Run SparCC' to start analysis.")
        return(invisible(NULL))
      }
      if (isTRUE(network_running())) {
        draw_wait_message("SparCC analysis is running. Please wait...")
        return(invisible(NULL))
      }
      req(network_result())
      draw_network_plot(network_result(), connected_only = TRUE)
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
    output$edge_table_all <- renderDT({ datatable(edge_table_all(), options = list(scrollX = TRUE)) })
    output$comparison_edge_table <- renderDT({
      datatable(comparison_edge_table(), options = list(scrollX = TRUE, pageLength = 15))
    })

    output$download_network_plot_all <- downloadHandler(
      filename = function() paste0("sparcc_network_all_", Sys.Date(), ".png"),
      content = function(file) {
        req(network_result())
        plot_width <- as.integer(input$plot_width)
        plot_height <- as.integer(input$plot_height)
        if (is.na(plot_width) || plot_width <= 0) plot_width <- 800L
        if (is.na(plot_height) || plot_height <= 0) plot_height <- 550L

        grDevices::png(file, width = plot_width * 2L, height = plot_height * 2L, res = 120)
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_network_plot(network_result(), connected_only = FALSE)
      }
    )

    output$download_network_plot_connected <- downloadHandler(
      filename = function() paste0("sparcc_network_connected_", Sys.Date(), ".png"),
      content = function(file) {
        req(network_result())
        plot_width <- as.integer(input$plot_width)
        plot_height <- as.integer(input$plot_height)
        if (is.na(plot_width) || plot_width <= 0) plot_width <- 800L
        if (is.na(plot_height) || plot_height <= 0) plot_height <- 550L

        grDevices::png(file, width = plot_width * 2L, height = plot_height * 2L, res = 120)
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_network_plot(network_result(), connected_only = TRUE)
      }
    )

    output$download_edge_table_all <- downloadHandler(
      filename = function() paste0("sparcc_edges_all_", Sys.Date(), ".tsv"),
      content = function(file) readr::write_tsv(edge_table_all(), file)
    )

    output$sparcc_legend_box <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 800 else input$plot_width
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
          uiOutput(session$ns("sparcc_figure_legend"))
        )
      )
    })

    output$sparcc_status_box <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 800 else input$plot_width
      tags$div(
        class = "simple-result-card",
        style = paste0("width: ", box_width, "px; max-width: 100%;"),
        verbatimTextOutput(session$ns("comparison_status"))
      )
    })

    output$sparcc_status_separator <- renderUI({
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

    output$sparcc_figure_legend <- renderUI({
      req(input$tax_level, input$analysis_mode)
      active_tab <- if (is.null(input$sparcc_active_tab) || !nzchar(input$sparcc_active_tab)) "Network Plot (All)" else input$sparcc_active_tab
      title_text <- paste0("SparCC ", active_tab)
      body_text <- if (identical(active_tab, "Network Plot (All)")) {
        "This panel shows all nodes and edges passing the minimum absolute edge-weight threshold. Edge color indicates sign (positive/negative), edge width indicates absolute weight, and node size reflects the selected node-size metric."
      } else if (identical(active_tab, "Network Plot (Connected)")) {
        "This panel shows the connected component view after edge filtering, excluding isolated nodes to emphasize connected structure."
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
            input$tax_level,
            ", analysis mode = ",
            input$analysis_mode,
            ", prevalence cutoff = ",
            input$prevalence_filter_pct,
            "%, minimum |edge weight| = ",
            input$min_edge_weight,
            "."
          )
        )
      )
    })
  })
}

