mod_biplot_ui <- function(id) {
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
    ")),
    tags$style(HTML(paste0(
      "#", ns("selected_taxa"), "-selectized::placeholder { font-size: 10px; }",
      "#", ns("selected_taxa"), "-selectized::-webkit-input-placeholder { font-size: 10px; }",
      "#", ns("selected_taxa"), "-selectized::-moz-placeholder { font-size: 10px; }",
      "#", ns("selected_taxa"), "-selectized:-ms-input-placeholder { font-size: 10px; }"
    ))),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        h4(icon("arrows-up-down-left-right"), "Association Biplot"),
        hr(),
        selectInput(ns("group_var"), "1. Group variable", choices = NULL),
        selectInput(
          ns("tax_level"),
          "2. Taxonomic level",
          choices = c("ASV", "Genus", "Species"),
          selected = "Genus"
        ),
        selectInput(
          ns("analysis_method"),
          "3. Analysis method",
          choices = c("dbRDA (capscale)" = "dbrda", "CCA" = "cca"),
          selected = "dbrda"
        ),
        selectInput(
          ns("distance_metric"),
          "4. Distance metric",
          choices = c("Bray-Curtis" = "bray", "Aitchison" = "aitchison"),
          selected = "bray"
        ),
        hr(),
        h4(icon("dna"), "Taxa Settings"),
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
          value = 100,
          min = 10,
          max = 500,
          step = 10
        ),
        numericInput(
          ns("top_taxa_vectors"),
          "Top taxa vectors",
          value = 15,
          min = 0,
          max = 100,
          step = 1
        ),
        selectizeInput(
          ns("selected_taxa"),
          "Selected taxa (optional)",
          choices = NULL,
          selected = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Selected taxa are drawn instead of top taxa vectors",
            plugins = list("remove_button")
          )
        ),
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Plot Settings"),
        numericInput(ns("plot_width"), "Plot width (px)", value = 700, min = 400, max = 2400, step = 50),
        numericInput(ns("plot_height"), "Plot height (px)", value = 500, min = 300, max = 2400, step = 50),
        numericInput(ns("base_size"), "Base Font Size:", value = 11, min = 6, max = 30, step = 1),
        numericInput(ns("dot_size"), "Dot Size (point size):", value = 3, min = 0.5, max = 10, step = 0.5),
        checkboxInput(ns("show_dot_outline"), "Show Dot Outline", value = TRUE),
        checkboxInput(ns("show_taxa_vectors"), "Show taxa vectors", value = TRUE),
        checkboxInput(ns("show_group_vectors"), "Show group vectors", value = TRUE),
        checkboxInput(ns("show_group_centroid"), "Show group centroids", value = TRUE),
        checkboxInput(ns("show_sample_names"), "Show sample names", value = FALSE),
        actionButton(ns("run_biplot"), "Run Biplot", class = "btn-danger", style = "font-size: 12px;"),
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
             if (msg.label) btn.textContent = msg.label;
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
        h4("Association Biplot"),
        plotOutput(ns("biplot_plot"), height = "auto"),
        uiOutput(ns("biplot_legend_box")),
        uiOutput(ns("biplot_status_separator")),
        h5(icon("circle-info"), "Biplot Status"),
        uiOutput(ns("biplot_status_box"))
      )
    )
  )
}

mod_biplot_server <- function(id, ps_obj, meta_vars = NULL) {
  moduleServer(id, function(input, output, session) {
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
    }, ignoreInit = FALSE)

    biplot_running <- reactiveVal(FALSE)

    biplot_payload <- eventReactive(input$run_biplot, {
      biplot_running(TRUE)
      session$sendCustomMessage(
        "toggle-biplot-run-btn",
        list(id = session$ns("run_biplot"), disabled = TRUE, label = "Running...")
      )
      on.exit(
        {
          biplot_running(FALSE)
          session$sendCustomMessage(
            "toggle-biplot-run-btn",
            list(id = session$ns("run_biplot"), disabled = FALSE, label = "Run Biplot")
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
      } else {
        totals <- colSums(otu, na.rm = TRUE)
        totals[totals <= 0 | is.na(totals)] <- 1
        otu_model <- sweep(otu, 2, totals, "/")
      }

      taxa_var <- apply(otu_model, 1, stats::var, na.rm = TRUE)
      otu_model <- otu_model[is.finite(taxa_var) & taxa_var > 0, , drop = FALSE]
      validate(need(nrow(otu_model) >= 3, "Too few taxa remain after variance filtering."))

      sd <- phyloseq::sample_data(ps_data, errorIfNULL = FALSE)
      validate(need(!is.null(sd), "Metadata is required for biplot."))
      meta_df <- as.data.frame(sd, stringsAsFactors = FALSE)
      validate(need(input$group_var %in% colnames(meta_df), "Selected group variable is not available in metadata."))

      sample_ids <- colnames(otu_model)
      meta_df <- meta_df[sample_ids, , drop = FALSE]
      group_raw <- meta_df[[input$group_var]]
      group_chr <- as.character(group_raw)
      keep_samples <- !is.na(group_raw) & nzchar(group_chr)
      validate(need(sum(keep_samples) >= 4, "At least 4 samples with valid group values are required."))
      otu_model <- otu_model[, keep_samples, drop = FALSE]
      meta_df <- meta_df[keep_samples, , drop = FALSE]
      group_raw <- meta_df[[input$group_var]]
      group_model <- if (is.numeric(group_raw) || is.integer(group_raw)) {
        as.numeric(group_raw)
      } else {
        group_chr <- as.character(group_raw)
        group_chr[is.na(group_chr) | !nzchar(group_chr)] <- "(Missing)"
        droplevels(as.factor(group_chr))
      }
      model_df <- data.frame(group_model = group_model, stringsAsFactors = FALSE)
      rownames(model_df) <- rownames(meta_df)
      comm_df <- as.data.frame(t(otu_model), stringsAsFactors = FALSE)
      validate(
        need(nrow(model_df) == nrow(comm_df), "Model data and community matrix row counts do not match."),
        need(all(rownames(model_df) == rownames(comm_df)), "Sample IDs between model data and community matrix do not match.")
      )

      model_type <- input$analysis_method
      permanova_p <- NA_real_
      permanova_r2 <- NA_real_
      metric_label <- if (identical(input$distance_metric, "aitchison")) "Aitchison" else "Bray-Curtis"
      if (identical(model_type, "dbrda")) {
        dist_method <- if (identical(input$distance_metric, "aitchison")) "euclidean" else "bray"
        dist_obj <- tryCatch(
          vegan::vegdist(comm_df, method = dist_method),
          error = function(e) {
            validate(need(FALSE, paste0("Distance calculation failed: ", conditionMessage(e))))
          }
        )
        model <- tryCatch(
          vegan::capscale(stats::as.formula("dist_obj ~ group_model"), data = model_df, comm = comm_df),
          error = function(e) {
            validate(need(FALSE, paste0("dbRDA model fitting failed: ", conditionMessage(e))))
          }
        )

        adonis_tbl <- tryCatch(
          vegan::adonis2(dist_obj ~ group_model, data = model_df, permutations = 999),
          error = function(e) NULL
        )
        if (!is.null(adonis_tbl)) {
          adf <- as.data.frame(adonis_tbl)
          if ("group_model" %in% rownames(adf)) {
            if ("Pr(>F)" %in% colnames(adf)) permanova_p <- as.numeric(adf["group_model", "Pr(>F)"])
            if ("R2" %in% colnames(adf)) permanova_r2 <- as.numeric(adf["group_model", "R2"])
          } else if (nrow(adf) >= 1) {
            if ("Pr(>F)" %in% colnames(adf)) permanova_p <- as.numeric(adf[1, "Pr(>F)"])
            if ("R2" %in% colnames(adf)) permanova_r2 <- as.numeric(adf[1, "R2"])
          }
        }
      } else {
        model <- tryCatch(
          vegan::cca(stats::as.formula("comm_df ~ group_model"), data = model_df),
          error = function(e) {
            validate(need(FALSE, paste0("CCA model fitting failed: ", conditionMessage(e))))
          }
        )
        metric_label <- "Not used (CCA)"
      }

      anova_all <- tryCatch(vegan::anova.cca(model, permutations = 999), error = function(e) NULL)
      p_all <- if (!is.null(anova_all) && nrow(as.data.frame(anova_all)) >= 1) as.data.frame(anova_all)[1, "Pr(>F)"] else NA_real_

      site_scores <- vegan::scores(model, display = "sites", choices = 1:2, scaling = 2)
      species_scores <- tryCatch(vegan::scores(model, display = "species", choices = 1:2, scaling = 2), error = function(e) NULL)
      biplot_scores <- tryCatch(vegan::scores(model, display = "bp", choices = 1:2, scaling = 2), error = function(e) NULL)

      site_df <- as.data.frame(site_scores, stringsAsFactors = FALSE)
      validate(
        need(ncol(site_df) >= 2, "Ordination did not return 2D site scores.")
      )
      colnames(site_df)[1:2] <- c("Axis1", "Axis2")
      site_df$SampleID <- rownames(site_df)
      site_df$Group <- as.character(model_df$group_model)

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
          bp_df$Variable <- rownames(bp_df)
        } else {
          bp_df <- data.frame()
        }
      }

      # For categorical groups, bp scores from model contrasts may omit the reference level.
      # Build explicit level-wise arrows from sample centroids so every group level is shown.
      group_arrow_df <- data.frame()
      if (is.factor(model_df$group_model)) {
        group_arrow_df <- site_df %>%
          dplyr::group_by(Group) %>%
          dplyr::summarise(
            Axis1 = mean(Axis1, na.rm = TRUE),
            Axis2 = mean(Axis2, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          dplyr::mutate(Variable = paste0("", Group))
      }

      eig <- vegan::eigenvals(model)
      axis1_var <- if (length(eig) >= 1) round(100 * eig[1] / sum(eig), 1) else NA_real_
      axis2_var <- if (length(eig) >= 2) round(100 * eig[2] / sum(eig), 1) else NA_real_

      list(
        site_df = site_df,
        species_df = species_df,
        species_all_df = species_all_df,
        bp_df = bp_df,
        group_arrow_df = group_arrow_df,
        model_label = if (identical(model_type, "cca")) "CCA" else "dbRDA",
        model_type = model_type,
        metric_label = metric_label,
        axis1_var = axis1_var,
        axis2_var = axis2_var,
        model_p = as.numeric(p_all),
        permanova_p = permanova_p,
        permanova_r2 = permanova_r2,
        n_samples = nrow(site_df),
        n_taxa = nrow(otu_model)
      )
    }, ignoreNULL = FALSE)

    observeEvent(biplot_payload(), {
      payload <- biplot_payload()
      all_taxa <- if (nrow(payload$species_all_df) > 0) payload$species_all_df$Taxon else character(0)
      current_selected <- isolate(input$selected_taxa)
      selected_keep <- if (length(current_selected) > 0) intersect(current_selected, all_taxa) else character(0)
      updateSelectizeInput(session, "selected_taxa", choices = all_taxa, selected = selected_keep, server = TRUE)
    }, ignoreInit = FALSE)

    biplot_plot <- reactive({
      req(biplot_payload())
      payload <- biplot_payload()
      base_size <- input$base_size
      if (is.null(base_size) || !is.finite(base_size)) {
        base_size <- 11
      }
      dot_size <- input$dot_size
      if (is.null(dot_size) || !is.finite(dot_size)) {
        dot_size <- 3
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
      group_arrow_df <- payload$group_arrow_df

      p <- ggplot2::ggplot(site_df, ggplot2::aes(x = Axis1, y = Axis2, color = Group)) +
        ggplot2::theme_bw(base_size = base_size) +
        ggplot2::labs(
          title = if (identical(payload$model_type, "cca")) {
            paste0(input$tax_level, "-Level Association Biplot (CCA)")
          } else {
            paste0(input$tax_level, "-Level Association Biplot (dbRDA, ", payload$metric_label, ")")
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
      if (isTRUE(input$show_dot_outline)) {
        outline_size <- dot_size + 0.4
        p <- p + ggplot2::geom_point(
          color = "black",
          size = outline_size,
          alpha = 0.9,
          show.legend = FALSE
        )
      }
      p <- p + ggplot2::geom_point(size = dot_size, alpha = 0.9)

      if (isTRUE(input$show_group_centroid)) {
        centroid_df <- site_df %>%
          dplyr::group_by(Group) %>%
          dplyr::summarise(Axis1 = mean(Axis1, na.rm = TRUE), Axis2 = mean(Axis2, na.rm = TRUE), .groups = "drop")
        if (nrow(centroid_df) > 0) {
          p <- p + ggplot2::geom_point(
            data = centroid_df,
            mapping = ggplot2::aes(x = Axis1, y = Axis2, color = Group),
            shape = 4,
            size = 4,
            stroke = 1.2,
            show.legend = FALSE
          )
        }
      }

      if (isTRUE(input$show_taxa_vectors) && nrow(taxa_df) > 0) {
        scale_mult <- 0.8
        p <- p +
          ggplot2::geom_segment(
            data = taxa_df,
            ggplot2::aes(x = 0, y = 0, xend = Axis1 * scale_mult, yend = Axis2 * scale_mult),
            inherit.aes = FALSE,
            color = "#8C2D04",
            linewidth = 0.5,
            arrow = grid::arrow(length = grid::unit(0.15, "cm"))
          )
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_text_repel(
            data = taxa_df,
            ggplot2::aes(x = Axis1 * scale_mult, y = Axis2 * scale_mult, label = Taxon),
            inherit.aes = FALSE,
            color = "#6A1B9A",
            size = text_size_taxa,
            min.segment.length = 0,
            box.padding = 0.2,
            point.padding = 0.1,
            max.overlaps = Inf
          )
        } else {
          p <- p + ggplot2::geom_text(
            data = taxa_df,
            ggplot2::aes(x = Axis1 * scale_mult, y = Axis2 * scale_mult, label = Taxon),
            inherit.aes = FALSE,
            color = "#6A1B9A",
            size = text_size_taxa,
            check_overlap = TRUE
          )
        }
      }

      if (isTRUE(input$show_group_vectors)) {
        if (nrow(group_arrow_df) > 0) {
          p <- p +
            ggplot2::geom_segment(
              data = group_arrow_df,
              ggplot2::aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
              inherit.aes = FALSE,
              color = "#1B5E20",
              linewidth = 0.7,
              linetype = "dashed",
              arrow = grid::arrow(length = grid::unit(0.17, "cm"))
            ) +
            ggplot2::geom_text(
              data = group_arrow_df,
              ggplot2::aes(x = Axis1, y = Axis2, label = Variable),
              inherit.aes = FALSE,
              color = "#1B5E20",
              size = text_size_group_vector,
              vjust = -0.4
            )
        } else if (nrow(bp_df) > 0) {
          p <- p +
            ggplot2::geom_segment(
              data = bp_df,
              ggplot2::aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
              inherit.aes = FALSE,
              color = "#1B5E20",
              linewidth = 0.7,
              linetype = "dashed",
              arrow = grid::arrow(length = grid::unit(0.17, "cm"))
            ) +
            ggplot2::geom_text(
              data = bp_df,
              ggplot2::aes(x = Axis1, y = Axis2, label = Variable),
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
          graphics::text(0.5, 0.5, "Click 'Run Biplot' to start analysis.", cex = 0.85)
          return(invisible(NULL))
        }
        if (isTRUE(biplot_running())) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Biplot analysis is running. Please wait...", cex = 0.85)
          return(invisible(NULL))
        }
        biplot_plot()
      },
      width = function() input$plot_width,
      height = function() input$plot_height
    )

    output$biplot_status <- renderText({
      req(biplot_payload())
      payload <- biplot_payload()
      pval_txt <- if (is.finite(payload$model_p)) format(payload$model_p, digits = 3) else "NA"
      permanova_p_txt <- if (is.finite(payload$permanova_p)) format(payload$permanova_p, digits = 3) else "NA"
      permanova_r2_txt <- if (is.finite(payload$permanova_r2)) format(payload$permanova_r2, digits = 3) else "NA"
      taxa_vectors_shown <- if (length(input$selected_taxa) > 0 && nrow(payload$species_all_df) > 0) {
        sum(payload$species_all_df$Taxon %in% input$selected_taxa)
      } else {
        nrow(payload$species_df)
      }
      paste(
        c(
          paste0("Model: ", payload$model_label),
          paste0("Distance metric: ", payload$metric_label),
          paste0("Group variable: ", input$group_var),
          paste0("Samples used: ", payload$n_samples),
          paste0("Taxa used: ", payload$n_taxa),
          paste0("Taxa vectors shown: ", taxa_vectors_shown),
          paste0(payload$model_label, " model p-value (perm): ", pval_txt),
          if (identical(payload$model_type, "dbrda")) paste0("PERMANOVA p-value (adonis2): ", permanova_p_txt) else NULL,
          if (identical(payload$model_type, "dbrda")) paste0("PERMANOVA R2 (adonis2): ", permanova_r2_txt) else NULL
        ),
        collapse = "\n"
      )
    })

    output$biplot_status_box <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 700 else input$plot_width
      tags$div(
        class = "simple-result-card",
        style = paste0("width: ", box_width, "px; max-width: 100%;"),
        verbatimTextOutput(session$ns("biplot_status"))
      )
    })

    output$biplot_legend_box <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 700 else input$plot_width
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
          uiOutput(session$ns("biplot_figure_legend"))
        )
      )
    })

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

    output$biplot_figure_legend <- renderUI({
      req(input$group_var, input$tax_level, input$analysis_method, input$distance_metric)
      distance_label <- if (identical(input$distance_metric, "aitchison")) "Aitchison" else "Bray-Curtis"
      body_text <- if (identical(input$analysis_method, "cca")) {
        paste0(
          "The CCA association biplot projects sample relationships on CCA1 and CCA2 at ",
          tolower(input$tax_level),
          " level. Distance metric is not used in CCA. Points are samples colored by ",
          input$group_var,
          ". ",
          if (isTRUE(input$show_taxa_vectors)) "Taxa arrows indicate taxa loadings. " else "",
          if (isTRUE(input$show_group_vectors)) "Group vectors show fitted group-direction effects. " else "",
          if (isTRUE(input$show_group_centroid)) "Cross markers indicate group centroids. " else "",
          if (isTRUE(input$show_sample_names)) "Sample labels are displayed near points." else ""
        )
      } else {
        paste0(
          "The dbRDA association biplot projects sample relationships on CAP1 and CAP2 using ",
          distance_label,
          " distance at ",
          tolower(input$tax_level),
          " level. Points are samples colored by ",
          input$group_var,
          ". ",
          if (isTRUE(input$show_taxa_vectors)) "Taxa arrows indicate taxa loadings. " else "",
          if (isTRUE(input$show_group_vectors)) "Group vectors show fitted group-direction effects. " else "",
          if (isTRUE(input$show_group_centroid)) "Cross markers indicate group centroids. " else "",
          if (isTRUE(input$show_sample_names)) "Sample labels are displayed near points." else ""
        )
      }
      tags$div(
        tags$div(style = "font-weight: 600; margin-bottom: 4px;", "Association biplot"),
        tags$div(body_text)
      )
    })

    output$download_biplot_png <- downloadHandler(
      filename = function() paste0("association_biplot_", Sys.Date(), ".png"),
      content = function(file) {
        dpi_val <- 300
        width_in <- input$plot_width / dpi_val
        height_in <- input$plot_height / dpi_val
        ggplot2::ggsave(file, plot = biplot_plot(), device = "png", width = width_in, height = height_in, units = "in", dpi = dpi_val)
      }
    )

    output$download_biplot_scores <- downloadHandler(
      filename = function() paste0("association_biplot_scores_", Sys.Date(), ".tsv"),
      content = function(file) {
        req(biplot_payload())
        payload <- biplot_payload()
        site_out <- payload$site_df %>% dplyr::mutate(Type = "Sample", Label = SampleID) %>% dplyr::select(Type, Label, Axis1, Axis2, Group)
        taxa_out <- if (nrow(payload$species_df) > 0) {
          payload$species_df %>% dplyr::mutate(Type = "Taxa", Label = Taxon, Group = NA_character_) %>% dplyr::select(Type, Label, Axis1, Axis2, Group)
        } else {
          data.frame(Type = character(0), Label = character(0), Axis1 = numeric(0), Axis2 = numeric(0), Group = character(0))
        }
        out_df <- dplyr::bind_rows(site_out, taxa_out)
        utils::write.table(out_df, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    )
  })
}
