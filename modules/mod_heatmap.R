heatmap_adjusted_pvalues <- function(plan) {
  values <- plan$values
  adjusted <- function(p) {
    valid <- is.finite(p)
    p[!valid] <- NA_real_
    p[valid] <- stats::p.adjust(p[valid], method = "BH")
    p
  }
  correlation_p <- function(x, y) {
    tryCatch(
      suppressWarnings(stats::cor.test(x, y, method = plan$method, exact = FALSE)$p.value),
      error = function(e) NA_real_
    )
  }
  if (identical(plan$type, "pairwise")) {
    q <- matrix(NA_real_, nrow(values), nrow(values), dimnames = rep(list(rownames(values)), 2L))
    pairs <- which(upper.tri(q), arr.ind = TRUE)
    p <- vapply(seq_len(nrow(pairs)), function(k) {
      correlation_p(values[pairs[k, 1L], ], values[pairs[k, 2L], ])
    }, numeric(1))
    q[pairs] <- adjusted(p)
    q[lower.tri(q)] <- t(q)[lower.tri(q)]
    return(q)
  }
  if (identical(plan$type, "continuous")) {
    p <- apply(values, 1, correlation_p, y = plan$group)
    return(matrix(adjusted(p), ncol = 1L, dimnames = list(rownames(values), plan$columns)))
  }
  p <- matrix(NA_real_, nrow(values), length(plan$columns),
              dimnames = list(rownames(values), plan$columns))
  for (j in seq_along(plan$groups)) {
    inside <- plan$groups[[j]]$inside
    outside <- plan$groups[[j]]$outside
    if (length(inside) < 2L || length(outside) < 2L) next
    p[, j] <- apply(values, 1, function(taxon) {
      tryCatch({
        if (plan$prevalence) {
          detected <- c(sum(taxon[inside] > 0, na.rm = TRUE), sum(taxon[outside] > 0, na.rm = TRUE))
          measured <- c(sum(is.finite(taxon[inside])), sum(is.finite(taxon[outside])))
          contingency <- cbind(detected, measured - detected)
          if (any(rowSums(contingency) == 0) || any(colSums(contingency) == 0)) return(NA_real_)
          stats::fisher.test(contingency)$p.value
        } else {
          stats::wilcox.test(taxon[inside], taxon[outside], exact = FALSE)$p.value
        }
      }, error = function(e) NA_real_)
    })
  }
  adjusted(p)
}

## UI
mod_heatmap_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$script(src = "heatmap-loading.js", defer = NA),
    tags$style(HTML("
      .well h4 { font-size: 16px; }
      .well h5 { font-size: 13px; }
      .well .control-label { font-size: 12px; }
      .well .checkbox label { font-size: 12px; }
      .well .form-control { font-size: 12px; }
      .well .btn { font-size: 11px; }
      .heatmap-plot-stage { position: relative; max-width: 100%; min-width: 0; overflow-x: auto; }
      .heatmap-plot-stage[data-fit-width=true] > .shiny-plot-output img {
        max-width: 100%;
        height: auto;
      }
      .heatmap-plot-stage[data-heatmap-loading=true] > .shiny-plot-output:empty {
        min-height: 700px;
      }
      .heatmap-wait-overlay {
        position: absolute;
        inset: 0;
        display: flex;
        align-items: center;
        justify-content: center;
        pointer-events: none;
        z-index: 1;
      }
      .heatmap-wait-overlay span {
        padding: 8px 12px;
        border-radius: 6px;
        background: rgba(255, 255, 255, 0.94);
        color: #334155;
      }
      .heatmap-plot-stage[data-heatmap-loading=false] .heatmap-drawing-overlay {
        display: none;
      }
    ")),
    plot_scroll_layout(
      sidebarPanel(
        width = 2,
        h4(icon("table-cells"), "Association Heatmap"),
        hr(),
        selectInput(
          ns("group_var"),
          "Primary Grouping Variable",
          choices = character(0),
          selected = NULL
        ),
        selectInput(
          ns("secondary_var"),
          "Secondary Grouping Variable",
          choices = "None",
          selected = "None"
        ),
        selectInput(ns("column_view"), "Column View", choices = c(
          "Group Summary" = "group", "Individual Samples" = "samples"
        )),
        conditionalPanel(
          "input.column_view !== 'samples'", ns = ns,
          selectInput(
            ns("group_type"),
            "Group Type",
            choices = c("Auto" = "auto", "Continuous" = "continuous", "Categorical" = "categorical"),
            selected = "auto"
          )
        ),
        hr(),
        h4(icon("sliders"), "Analysis & Plot Settings"),
        selectInput(
          ns("tax_level"),
          "Taxonomic Level",
          choices = c("ASV", "Genus", "Species", "Strain"),
          selected = "Genus"
        ),
        selectInput(
          ns("association_mode"),
          "Association Mode",
          choices = c("Abundance-based" = "abundance", "Prevalence-based (binary)" = "prevalence"),
          selected = "abundance"
        ),
        selectInput(
          ns("transform_method"),
          "Data Transform",
          choices = c("TSS" = "tss", "CLR" = "clr"),
          selected = "clr"
        ),
        selectInput(
          ns("value_scale"),
          "Heatmap Value Scale",
          choices = c("Raw" = "raw", "Z-score (by taxa)" = "zscore"),
          selected = "raw"
        ),
        tags$details(
          tags$summary("Advanced Options"),
          conditionalPanel(
            "input.column_view !== 'samples'", ns = ns,
            selectInput(
              ns("corr_method"),
              "Correlation Method",
              choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
              selected = "pearson"
            )
          ),
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
            value = 30,
            min = 10,
            max = 300,
            step = 5
          ),
          selectizeInput(
            ns("selected_taxa"),
            "Selected Taxa (Optional)",
            choices = NULL,
            selected = NULL,
            multiple = TRUE,
            options = list(placeholder = "If selected, only these taxa are shown")
          ),
          checkboxInput(ns("show_row_dendrogram"), "Show Row Dendrogram (Cluster Row)", value = TRUE),
          conditionalPanel(
            "input.column_view === 'samples'", ns = ns,
            checkboxInput(ns("cluster_sample_columns"), "Cluster Samples Within Groups", value = FALSE)
          ),
          conditionalPanel(
            "input.column_view !== 'samples'", ns = ns,
            checkboxInput(ns("show_col_dendrogram"), "Show Column Dendrogram (Cluster Column)", value = TRUE)
          ),
          numericInput(ns("abs_limit"), "Signed Color Scale Abs Limit (0 = Auto)", value = 0, min = 0, max = 100, step = 0.1),
          conditionalPanel(
            "input.column_view !== 'samples'", ns = ns,
            checkboxInput(ns("show_significance_marks"), "Show Significance Marks (*)", value = FALSE),
            checkboxInput(ns("apply_sig_filter"), "Mask Non-Significant Cells (FDR)", value = FALSE),
            numericInput(ns("q_cutoff"), "FDR Cutoff (q)", value = 0.05, min = 0.0001, max = 1, step = 0.01)
          )
        ),
        hr(),
        h4(icon("sliders"), "Display Options"),
        selectInput(ns("plot_size_mode"), "Plot Size", choices = c("Auto" = "auto", "Manual" = "manual")),
        conditionalPanel(
          "input.plot_size_mode === 'manual'", ns = ns,
          numericInput(ns("plot_width"), "Plot Width (px)", value = 900, min = 400, max = 2400, step = 50),
          numericInput(ns("plot_height"), "Plot Height (px)", value = 700, min = 300, max = 2400, step = 50)
        ),
        numericInput(ns("base_size"), "Base Font Size", value = 10, min = 6, max = 30, step = 1),
        selectInput(ns("x_label_angle"), "X-Axis Label Angle", choices = c(
          "0\u00b0 Horizontal" = "0", "45\u00b0 Diagonal" = "45", "90\u00b0 Vertical" = "90"
        ), selected = "90"),
        conditionalPanel(
          "output.palette_mode === 'diverging'", ns = ns,
          selectInput(ns("color_palette"), "Color Palette", diverging_plot_palettes, selected = "blue_red")
        ),
        conditionalPanel(
          "output.palette_mode === 'sequential'", ns = ns,
          selectInput(ns("sequential_color_palette"), "Color Palette", sequential_plot_palettes, selected = "plasma")
        ),
        checkboxInput(ns("show_cell_values"), "Show Cell Values", value = FALSE),
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
          "Shiny.addCustomMessageHandler('toggle-input-disabled', function(msg) {
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
        automatic_analysis_view(ns,
          div(
            class = "heatmap-plot-stage",
            `data-fit-width` = "false",
            `data-heatmap-loading` = "true",
            tags$script(HTML(sprintf(
              "$(document).on('shiny:inputchanged', function(event) {
                if (event.name !== '%s') return;
                const plot = document.getElementById('%s');
                if (!plot) return;
                const stage = plot.parentElement;
                stage.dataset.fitWidth = event.value === 'samples';
              });",
              ns("column_view"), ns("corr_heatmap_plot")
            ))),
            plotOutput(ns("corr_heatmap_plot"), height = "auto"),
            conditionalPanel(
              "input.column_view !== 'samples' && (input.show_significance_marks || input.apply_sig_filter) && (output.heatmap_test_state === 'off' || output.heatmap_test_state === 'waiting')",
              ns = ns,
              div(class = "heatmap-wait-overlay", role = "status", span("Calculating significance tests..."))
            ),
            conditionalPanel(
              "input.column_view === 'samples' || !(input.show_significance_marks || input.apply_sig_filter) || (output.heatmap_test_state !== 'off' && output.heatmap_test_state !== 'waiting')",
              ns = ns,
              div(class = "heatmap-wait-overlay heatmap-drawing-overlay", role = "status",
                  span("Preparing the heatmap. Please wait..."))
            )
          ),
          plot_legend_output(ns, "heatmap_legend_box", c("corr_heatmap_plot")),
          uiOutput(ns("heatmap_status_separator")),
          h5(icon("circle-info"), "Heatmap Status"),
          uiOutput(ns("heatmap_status_box")),
          prompt_overlay = function(prompt) {
            div(
              class = "heatmap-plot-stage",
              `data-heatmap-loading` = "true",
              prompt,
              div(class = "heatmap-wait-overlay heatmap-drawing-overlay", role = "status",
                  span("Preparing the heatmap. Please wait..."))
            )
          }
        )
      )
    )
  )
}

## Server
mod_heatmap_server <- function(id, ps_obj, meta_vars = NULL, plot_order = function() NULL, active = function() TRUE) {
  moduleServer(id, function(input, output, session) {
    supported_taxonomic_ranks <- c("Genus", "Species", "Strain")
    tss_value_scales <- c("Square Root (Hellinger)" = "hellinger", "Log TSS (Log10(% + 1))" = "log_tss")

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
      if (is.null(input$secondary_var) || identical(input$secondary_var, "None")) {
        return(NULL)
      }
      resolve_meta_colname(input$secondary_var, colnames(meta_df))
    })

    group_columns <- reactive(setdiff(colnames(metadata_df()), "SampleID"))

    observeEvent(group_columns(), {
      columns <- group_columns()
      current <- isolate(input$group_var)
      selected <- if (isTruthy(current) && current %in% columns) current else head(columns, 1L)
      updateSelectInput(session, "group_var", choices = columns, selected = selected)
    })

    observeEvent(list(group_columns(), input$group_var), {
      columns <- group_columns()
      primary <- resolve_meta_colname(input$group_var, columns)
      valid_primary <- isTruthy(primary) && primary %in% columns

      secondary_choices <- if (valid_primary) setdiff(columns, primary) else columns
      current_secondary <- isolate(input$secondary_var)
      selected_secondary <- if (isTruthy(current_secondary) &&
                                current_secondary %in% secondary_choices) current_secondary else "None"
      updateSelectInput(session, "secondary_var",
                        choices = c("None" = "None", secondary_choices), selected = selected_secondary)
    })

    heatmap_opened <- reactiveVal(FALSE)
    observeEvent(active(), {
      if (isTRUE(active())) heatmap_opened(TRUE)
    })

    heatmap_rank_matrix <- reactive({
      req(heatmap_opened())
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

      otu_mat <- as(phyloseq::otu_table(ps_data), "matrix")
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


    samples_view <- reactive(identical(input$column_view, "samples"))

    current_analysis_settings <- reactive({
      list(
        samples = samples_view(),
        group_var = input$group_var,
        secondary_var = input$secondary_var,
        group_type = if (!samples_view()) input$group_type else NULL,
        tax_level = input$tax_level,
        transform_method = input$transform_method,
        association_mode = input$association_mode,
        corr_method = if (!samples_view()) input$corr_method else NULL,
        value_scale = if (input$value_scale %in% tss_value_scales) "raw" else input$value_scale,
        prevalence_filter_pct = as.numeric(input$prevalence_filter_pct),
        max_taxa = as.integer(input$max_taxa)
      )
    })

    analysis_inputs <- reactive({
      ps <- ps_obj()
      validate(need(!is.null(ps), "Load and preprocess data to view the heatmap."))
      list(ps = ps, config = current_analysis_settings())
    })

    corr_payload <- automatic_analysis(active, analysis_inputs, function(request) {
      by_sample <- request$config$samples
      validate(need(
        all(vapply(list(input$tax_level, input$transform_method, input$association_mode,
                        input$value_scale), isTruthy, logical(1))) &&
          (by_sample || isTruthy(input$corr_method)),
        "Choose the required heatmap settings."
      ))
      ps_data <- request$ps
      minimum <- if (by_sample) 1L else 3L
      validate(
        need(phyloseq::nsamples(ps_data) >= minimum, if (by_sample) "At least one sample is required." else "At least 3 samples are required."),
        need(phyloseq::ntaxa(ps_data) >= minimum, if (by_sample) "At least one taxon is required." else "At least 3 taxa are required.")
      )

      otu_mat <- heatmap_rank_matrix()
      column_panels <- NULL
      column_group <- NULL
      sample_totals <- NULL
      if (by_sample) {
        meta <- metadata_df()
        primary <- resolve_meta_colname(input$group_var, names(meta))
        validate(need(isTruthy(primary) && primary %in% names(meta), "Select a primary grouping variable."))
        ids <- colnames(otu_mat)
        metadata_rows <- match(ids, rownames(meta))
        group_values <- function(column) {
          values <- as.character(meta[[column]][metadata_rows])
          values[is.na(values) | !nzchar(values)] <- "(Missing)"
          values
        }
        primary_values <- group_values(primary)
        column_group <- secondary_var_resolved()
        secondary <- if (isTruthy(column_group) && column_group %in% names(meta) && column_group != primary) {
          group_values(column_group)
        } else {
          column_group <- NULL
          rep("", length(ids))
        }
        sample_totals <- colSums(otu_mat, na.rm = TRUE)
        column_panels <- data.frame(
          primary = primary_values, secondary = secondary, sample = ids,
          row.names = ids, stringsAsFactors = FALSE
        )
      }

      prevalence_cutoff <- as.numeric(input$prevalence_filter_pct)
      if (is.na(prevalence_cutoff)) prevalence_cutoff <- 10
      prevalence_cutoff <- max(0, min(100, prevalence_cutoff))
      keep_prev <- rowMeans(otu_mat > 0, na.rm = TRUE) * 100 >= prevalence_cutoff
      otu_mat <- otu_mat[keep_prev, , drop = FALSE]
      validate(
        need(nrow(otu_mat) >= minimum, "Too few taxa remain after prevalence filtering.")
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
        col_totals <- if (by_sample) sample_totals else colSums(otu_mat, na.rm = TRUE)
        col_totals[col_totals <= 0 | is.na(col_totals)] <- 1
        sweep(otu_mat, 2, col_totals, "/")
      }

      analysis_mat <- if (identical(input$association_mode, "prevalence")) {
        (otu_mat > 0) * 1
      } else {
        transformed
      }

      if (!by_sample) {
        taxa_var <- apply(analysis_mat, 1, stats::var, na.rm = TRUE)
        keep_var <- is.finite(taxa_var) & taxa_var > 0
        analysis_mat <- analysis_mat[keep_var, , drop = FALSE]
        validate(need(nrow(analysis_mat) >= 3, "Too few taxa remain after variance filtering."))
      }

      meta_df <- metadata_df()
      sample_ids <- colnames(analysis_mat)
      if (nrow(meta_df) > 0) {
        meta_df <- meta_df[sample_ids, , drop = FALSE]
      }

      sig_note <- "Not computed."
      test_plan <- NULL

      if (by_sample) {
        assoc_mat <- analysis_mat
        assoc_type <- if (identical(input$association_mode, "prevalence")) "Sample-level detection" else
          "Sample-level abundance"
      } else if (identical(analysis_target(), "taxa_taxa")) {
        corr_method_eff <- if (identical(input$association_mode, "prevalence")) "pearson" else input$corr_method
        assoc_mat <- stats::cor(t(analysis_mat), method = corr_method_eff, use = "pairwise.complete.obs")
        assoc_type <- if (identical(input$association_mode, "prevalence")) {
          "Taxa vs Taxa prevalence association (phi)"
        } else {
          "Taxa vs Taxa correlation"
        }
        sig_note <- if (identical(input$association_mode, "prevalence")) {
          "P-values come from pairwise phi (binary Pearson) tests."
        } else {
          "P-values come from pairwise correlation tests."
        }

        test_plan <- list(type = "pairwise", method = corr_method_eff)
      } else {
        validate(need(isTruthy(input$group_var), "Select a primary grouping variable."))
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
          keep_idx <- !is.na(primary_vec)
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
          sig_note <- if (identical(input$association_mode, "prevalence")) {
            "P-values come from point-biserial correlation tests of taxon presence against the continuous variable."
          } else {
            "P-values come from correlation tests of taxon abundance against the continuous variable."
          }

          test_plan <- list(type = "continuous", method = corr_method_eff,
                            group = group_num, columns = effective_group_var)
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
            validate(need(length(primary_levels) > 0L, "No samples have a valid primary level."))
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
          sig_note <- if (identical(input$association_mode, "prevalence")) {
            "For each taxon, Fisher's exact test compares presence in the indicated group with presence in all other groups combined."
          } else {
            "For each taxon, a Wilcoxon rank-sum test compares the indicated group with all other groups combined."
          }

          if (!is.null(column_panels)) {
            sig_note <- paste(sig_note, "Comparisons stay within each primary group, and BH adjustment covers all valid tests across the displayed panels.")
          }
          group_indices <- vector("list", length(column_levels))
          for (j in seq_along(column_levels)) {
            in_panel <- if (is.null(column_panels)) rep(TRUE, length(group_chr)) else
              primary_vec == column_panels$primary[j]
            in_group <- which(in_panel & group_chr == column_levels[j])
            out_group <- which(in_panel & group_chr != column_levels[j])
            group_indices[[j]] <- list(inside = in_group, outside = out_group)
            if (length(in_group)) {
              in_mean <- rowMeans(analysis_mat[, in_group, drop = FALSE], na.rm = TRUE)
              if (!identical(input$association_mode, "prevalence")) {
                assoc_mat[, j] <- in_mean
              } else if (length(out_group)) {
                assoc_mat[, j] <- in_mean - rowMeans(analysis_mat[, out_group, drop = FALSE], na.rm = TRUE)
              }
            }
          }
          test_plan <- list(type = "categorical", groups = group_indices, columns = column_ids,
                            prevalence = identical(input$association_mode, "prevalence"))
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
        z_by_row <- vapply(seq_len(nrow(assoc_mat)), function(i) {
          row_x <- assoc_mat[i, ]
          row_sd <- stats::sd(row_x, na.rm = TRUE)
          if (!is.finite(row_sd) || row_sd <= 0) {
            return(ifelse(is.finite(row_x), 0, NA_real_))
          }
          (row_x - mean(row_x, na.rm = TRUE)) / row_sd
        }, numeric(ncol(assoc_mat)))
        assoc_mat[] <- t(matrix(z_by_row, nrow = ncol(assoc_mat)))
      }


      raw_assoc_mat <- raw_assoc_mat[rownames(assoc_mat), colnames(assoc_mat), drop = FALSE]

      heat_df <- as.data.frame(as.table(assoc_mat), stringsAsFactors = FALSE)
      colnames(heat_df) <- c("AxisX", "AxisY", "Association")
      raw_df <- as.data.frame(as.table(raw_assoc_mat), stringsAsFactors = FALSE)
      colnames(raw_df) <- c("AxisX", "AxisY", "AssociationRaw")
      heat_df <- dplyr::left_join(heat_df, raw_df, by = c("AxisX", "AxisY"))
      if (!is.null(test_plan)) test_plan$values <- analysis_mat

      list(
        assoc_mat = assoc_mat,
        column_group = column_group,
        column_panels = column_panels,
        heat_df = heat_df,
        n_samples = ncol(analysis_mat),
        n_taxa = nrow(analysis_mat),
        assoc_type = assoc_type,
        sig_note = sig_note,
        test_plan = test_plan,
        run_settings = request$config
      )
    })

    is_tss_abundance <- function(payload) {
      settings <- payload$run_settings
      (isTRUE(settings$samples) || !is.null(payload$column_group)) &&
        identical(settings$association_mode, "abundance") &&
        identical(settings$transform_method, "tss")
    }

    observe({
      req(heatmap_opened())
      choices <- c("Raw" = "raw", "Z-score (by taxa)" = "zscore")
      if (identical(input$association_mode, "abundance") && identical(input$transform_method, "tss") &&
          (samples_view() || is_tss_abundance(corr_payload()))) {
        choices <- c(choices, tss_value_scales)
      }
      selected <- isolate(input$value_scale)
      if (identical(input$association_mode, "prevalence") || !isTruthy(selected) || !selected %in% choices) selected <- "raw"
      updateSelectInput(session, "value_scale", choices = choices, selected = selected)
    })

    significance_requested <- reactive(
      !samples_view() && (isTRUE(input$show_significance_marks) || isTRUE(input$apply_sig_filter))
    )
    significance_enabled <- debounce(significance_requested, 150)
    significance_q <- reactive(heatmap_adjusted_pvalues(corr_payload()$test_plan))

    output$heatmap_test_state <- renderText({
      req(heatmap_opened())
      if (!significance_requested()) return("off")
      if (!isTRUE(significance_enabled())) return("waiting")
      tryCatch({
        significance_q()
        "ready"
      }, error = function(e) "unavailable")
    })
    outputOptions(output, "heatmap_test_state", suspendWhenHidden = FALSE)

    display_corr_payload <- reactive({
      payload <- corr_payload()
      scale <- input$value_scale
      if (is_tss_abundance(payload) && scale %in% tss_value_scales) {
        # Transform sample values or completed group means without changing the test inputs.
        transform <- if (identical(scale, "hellinger")) sqrt else function(x) log10(100 * x + 1)
        payload$assoc_mat <- transform(payload$assoc_mat)
        payload$heat_df$Association <- transform(payload$heat_df$Association)
        payload$run_settings$value_scale <- scale
      }
      payload$sig_tested <- significance_requested()
      payload$heat_df$QValue <- NA_real_
      if (payload$sig_tested) {
        req(isTRUE(significance_enabled()), cancelOutput = TRUE)
        q <- significance_q()
        payload$heat_df$QValue <- q[cbind(
          match(as.character(payload$heat_df$AxisX), rownames(q)),
          match(as.character(payload$heat_df$AxisY), colnames(q))
        )]
      }
      cutoff <- as.numeric(input$q_cutoff)
      if (length(cutoff) != 1L || !is.finite(cutoff)) cutoff <- 0.05
      payload$q_cutoff <- max(0, min(1, cutoff))
      significant <- is.finite(payload$heat_df$QValue) & payload$heat_df$QValue <= payload$q_cutoff
      payload$heat_df$SigMark <- ifelse(significant, "*", "")
      payload$n_tests <- sum(is.finite(payload$heat_df$QValue))
      payload$n_sig <- sum(significant)
      if (isTRUE(input$apply_sig_filter) && payload$sig_tested) {
        payload$heat_df$Association[!significant] <- NA_real_
      }
      payload
    })
    register_automatic_analysis_status(output, active, corr_payload,
                                       "Choose data and a primary group to view the heatmap.",
                                       updating_message = "Preparing the heatmap. Please wait...")

    observeEvent(corr_payload(), {
      taxa_choices <- rownames(corr_payload()$assoc_mat)
      current_selected <- isolate(input$selected_taxa)
      selected_keep <- if (length(current_selected) > 0) intersect(current_selected, taxa_choices) else character(0)
      updateSelectizeInput(session, "selected_taxa", choices = taxa_choices, selected = selected_keep, server = TRUE)
    }, ignoreInit = FALSE)

    heatmap_base_size <- reactive({
      size <- as.numeric(input$base_size)
      if (!isTruthy(size) || !is.finite(size)) size <- 10
      max(6, min(30, size))
    })

    heatmap_x_label_angle <- reactive({
      if (isTruthy(input$x_label_angle) && input$x_label_angle %in% c("0", "45", "90")) {
        as.numeric(input$x_label_angle)
      } else {
        90
      }
    })

    auto_cell_height <- function(n_rows) {
      max(14, min(26, floor((700 - 230) / max(1, n_rows))))
    }

    palette_mode <- reactive({
      payload <- display_corr_payload()
      req(payload)
      if (isTRUE(payload$run_settings$samples) && identical(input$association_mode, "prevalence")) {
        return("detection")
      }
      if (cell_values_are_percent(payload) || length(tss_value_scale_label(payload)) > 0L) {
        return("sequential")
      }
      "diverging"
    })

    output$palette_mode <- renderText(palette_mode())
    outputOptions(output, "palette_mode", suspendWhenHidden = FALSE)

    heatmap_palette <- reactive({
      mode <- palette_mode()
      if (mode == "detection") {
        return(list(colors = c("#F7F7F7", "#B2182B"), label = "Fixed (Absent / Present)"))
      }
      if (mode == "sequential") {
        choices <- sequential_plot_palettes
        palette <- input$sequential_color_palette
        default <- "plasma"
        color_function <- sequential_plot_colors
      } else {
        choices <- diverging_plot_palettes
        palette <- input$color_palette
        default <- "blue_red"
        color_function <- diverging_plot_colors
      }
      if (length(palette) != 1L || !palette %in% choices) palette <- default
      list(colors = color_function(palette), label = names(choices)[match(palette, choices)])
    })

    cell_values_are_percent <- function(payload) {
      is_tss_abundance(payload) && identical(payload$run_settings$value_scale, "raw")
    }

    tss_value_scale_label <- function(payload) {
      names(tss_value_scales)[match(payload$run_settings$value_scale, tss_value_scales, nomatch = 0L)]
    }

    heatmap_legend_title <- function(payload) {
      settings <- payload$run_settings
      if (cell_values_are_percent(payload)) {
        if (isTRUE(settings$samples)) "Relative\nabundance (%)" else "Mean relative\nabundance (%)"
      } else if (identical(settings$value_scale, "hellinger")) {
        if (isTRUE(settings$samples)) "Square root of\nrelative abundance" else "Square root of mean\nrelative abundance"
      } else if (identical(settings$value_scale, "log_tss")) {
        if (isTRUE(settings$samples)) "Log10(relative\nabundance % + 1)" else "Log10(mean relative\nabundance % + 1)"
      } else if (identical(settings$association_mode, "prevalence")) {
        if (isTRUE(settings$samples)) "Detection" else if (identical(analysis_target(), "taxa_group")) "Delta prevalence\n(in vs out)" else
          "Phi coefficient\n(binary association)"
      } else if (identical(settings$value_scale, "zscore")) {
        "Z-score\n(by taxa)"
      } else if (identical(settings$transform_method, "clr")) {
        "CLR\nvalue"
      } else {
        "TSS\nvalue"
      }
    }

    sample_tree_height_pt <- function(base_size) max(36, base_size * 4)

    heatmap_plot_spec <- function(canvas_width_mm) {
      payload <- display_corr_payload()
      heat_df <- payload$heat_df
      assoc_mat <- payload$assoc_mat
      panels <- payload$column_panels
      by_sample <- isTRUE(payload$run_settings$samples)
      show_subfacets <- by_sample && !is.null(panels) && !is.null(payload$column_group)
      if (!is.null(panels)) {
        order_config <- plot_order()
        primary_order <- plot_order_factor(panels$primary, input$group_var, order_config, unique(panels$primary))
        secondary_order <- plot_order_factor(panels$secondary, payload$column_group, order_config, unique(panels$secondary))
        panel_order <- if (by_sample) {
          sample_order <- plot_order_factor(panels$sample, "SampleID", order_config, panels$sample)
          order(primary_order, secondary_order, sample_order)
        } else {
          order(primary_order, secondary_order)
        }
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

      column_trees <- list()
      if (by_sample && !is.null(panels) && isTRUE(input$cluster_sample_columns)) {
        previous <- seq_len(nrow(panels) - 1L)
        starts <- c(TRUE, panels$primary[previous] != panels$primary[previous + 1L] |
          panels$secondary[previous] != panels$secondary[previous + 1L])
        groups <- split(seq_len(ncol(assoc_mat)), cumsum(starts))
        column_order <- seq_len(ncol(assoc_mat))
        for (columns in groups) {
          if (length(columns) < 2L) next
          tree <- tryCatch(
            stats::hclust(stats::dist(t(assoc_mat[, columns, drop = FALSE])), method = "average"),
            error = function(e) NULL
          )
          if (is.null(tree)) next
          column_order[columns] <- columns[tree$order]
          column_trees[[length(column_trees) + 1L]] <- list(columns = columns, tree = tree)
        }
        assoc_mat <- assoc_mat[, column_order, drop = FALSE]
        panels <- panels[column_order, , drop = FALSE]
      }

      raw_tss_mean <- cell_values_are_percent(payload)
      sequential_tss <- raw_tss_mean || length(tss_value_scale_label(payload)) > 0L
      sample_detection <- by_sample && identical(input$association_mode, "prevalence")
      finite_associations <- heat_df$Association[is.finite(heat_df$Association)]
      if (sample_detection) {
        fill_limit <- 1
      } else if (sequential_tss) {
        fill_limit <- max(c(0, finite_associations))
        if (fill_limit <= 0) fill_limit <- max(c(0, assoc_mat[is.finite(assoc_mat)]))
        if (fill_limit <= 0) fill_limit <- 0.01
      } else {
        absolute_values <- abs(finite_associations)
        fill_limit <- if (length(absolute_values)) {
          as.numeric(stats::quantile(absolute_values, probs = 0.98, names = FALSE, type = 8))
        } else {
          1
        }
        if (!is.finite(fill_limit) || fill_limit <= 0) fill_limit <- max(c(0, absolute_values))
        if (fill_limit <= 0) fill_limit <- 1
        abs_limit_input <- as.numeric(input$abs_limit)
        if (isTruthy(abs_limit_input) && is.finite(abs_limit_input) && abs_limit_input > 0) {
          fill_limit <- abs_limit_input
        }
      }

      base_size <- heatmap_base_size()
      font_family <- "sans"
      sig_fontsize_pt <- max(9, base_size * 1.15)

      # Let ComplexHeatmap allocate horizontal space after laying out labels and the legend.
      body_height <- NULL
      if (!identical(input$plot_size_mode, "manual")) {
        body_height <- grid::unit(nrow(assoc_mat) * auto_cell_height(nrow(assoc_mat)) / 110, "inches")
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

      colors <- heatmap_palette()$colors
      legend_ticks <- list()
      if (sample_detection) {
        col_fun <- circlize::colorRamp2(c(0, 1), colors)
        legend_ticks <- list(at = c(0, 1), labels = c("Absent", "Present"))
      } else if (sequential_tss) {
        col_fun <- circlize::colorRamp2(seq(0, fill_limit, length.out = length(colors)), colors)
        breaks <- pretty(c(0, fill_limit), n = 4)
        breaks <- c(0, breaks[breaks > 0 & breaks < fill_limit * 0.95], fill_limit)
        legend_ticks <- list(at = breaks, labels = scales::label_number(scale = if (raw_tss_mean) 100 else 1)(breaks))
      } else {
        col_fun <- circlize::colorRamp2(c(-fill_limit, 0, fill_limit), colors)
      }

      legend_title <- heatmap_legend_title(payload)

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

      if (isTRUE(input$apply_sig_filter) && payload$sig_tested) {
        assoc_mat[sig_mat != "*"] <- NA_real_
      }

      show_values <- isTRUE(input$show_cell_values)
      show_marks <- !by_sample && isTRUE(input$show_significance_marks)
      labels <- matrix("", nrow = nrow(assoc_mat), ncol = ncol(assoc_mat))
      if (show_values) {
        finite <- is.finite(assoc_mat)
        values <- assoc_mat[finite]
        labels[finite] <- if (raw_tss_mean) {
          sprintf("%.1f", 100 * values)
        } else {
          formatted <- sprintf("%.1f", values)
          replace(formatted, formatted == "-0.0", "0.0")
        }
      }
      if (show_marks) {
        marked <- sig_mat == "*" & is.finite(assoc_mat)
        labels[marked] <- paste0(labels[marked], "*")
      }

      legend_fonts <- heatmap_legend_font_sizes()
      legend_params <- c(list(
        title = legend_title,
        title_position = "topcenter",
        title_gp = grid::gpar(fontfamily = font_family, fontsize = legend_fonts[["title"]]),
        labels_gp = grid::gpar(fontfamily = font_family, fontsize = legend_fonts[["labels"]])
      ), legend_ticks)
      panel_titles <- NULL
      secondary_breaks <- integer()
      if (!is.null(panels)) {
        # Match the relative text size and vertical margins of the Taxa Summary facet strips.
        facet_labels <- levels(panel_split)
        facet_gp <- grid::gpar(fontfamily = font_family, fontsize = base_size * 0.8, fontface = "bold")
        label_height <- max(grid::grobHeight(grid::textGrob(facet_labels, gp = facet_gp)))
        primary_strip <- ComplexHeatmap::anno_block(
            labels = facet_labels,
            gp = grid::gpar(fill = "#F1F3F5", col = "#595959", lwd = 0.8),
            labels_gp = facet_gp,
            labels_offset = grid::unit(0.5, "npc"), labels_just = "center",
            height = label_height + grid::unit(base_size * 0.8, "pt")
        )
        annotations <- list()
        if (length(column_trees)) {
          annotations$trees <- ComplexHeatmap::AnnotationFunction(
            which = "column",
            n = ncol(assoc_mat),
            height = grid::unit(sample_tree_height_pt(base_size), "pt"),
            fun = function(index) {
              for (entry in column_trees) {
                positions <- match(entry$columns, index)
                if (anyNA(positions)) next
                tree <- entry$tree
                leaf_x <- (positions[order(tree$order)] - 0.5) / length(index)
                node_x <- numeric(nrow(tree$merge))
                max_height <- max(tree$height)
                node_y <- if (max_height > 0) 0.92 * tree$height / max_height else rep(0.08, length(node_x))
                left_x <- right_x <- left_y <- right_y <- node_x
                for (k in seq_along(node_x)) {
                  left <- tree$merge[k, 1L]
                  right <- tree$merge[k, 2L]
                  left_x[k] <- if (left < 0L) leaf_x[-left] else node_x[left]
                  right_x[k] <- if (right < 0L) leaf_x[-right] else node_x[right]
                  left_y[k] <- if (left < 0L) 0 else node_y[left]
                  right_y[k] <- if (right < 0L) 0 else node_y[right]
                  node_x[k] <- (left_x[k] + right_x[k]) / 2
                }
                grid::grid.segments(
                  x0 = c(left_x, right_x, left_x), x1 = c(left_x, right_x, right_x),
                  y0 = c(left_y, right_y, node_y), y1 = rep(node_y, 3L),
                  gp = grid::gpar(col = "#595959")
                )
              }
            }
          )
        }
        annotations$facet <- primary_strip
        if (show_subfacets) {
          secondary_labels <- panels$secondary
          previous <- seq_len(nrow(panels) - 1L)
          secondary_breaks <- previous[
            panels$primary[previous] == panels$primary[previous + 1L] &
              secondary_labels[previous] != secondary_labels[previous + 1L]
          ] + 1L
          secondary_strip <- ComplexHeatmap::AnnotationFunction(
            which = "column",
            n = ncol(assoc_mat),
            height = grid::unit(base_size * 1.6, "pt"),
            fun = function(index) {
              blocks <- rle(secondary_labels[index])
              ends <- cumsum(blocks$lengths)
              for (k in seq_along(ends)) {
                span <- blocks$lengths[k] / length(index)
                center <- ends[k] / length(index) - span / 2
                grid::pushViewport(grid::viewport(x = center, width = span, clip = "on"))
                grid::grid.rect(gp = grid::gpar(fill = "#FAFBFC", col = "#8A8A8A", lwd = 0.6))
                grid::grid.text(blocks$values[k], gp = grid::gpar(
                  fontfamily = font_family, fontsize = base_size * 0.72
                ))
                grid::popViewport()
              }
            }
          )
          annotations$subgroup <- secondary_strip
        }
        panel_titles <- do.call(ComplexHeatmap::HeatmapAnnotation,
                                c(annotations, list(show_annotation_name = FALSE)))
      }
      heatmap <- ComplexHeatmap::Heatmap(
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
        column_labels = if (by_sample || is.null(panels)) colnames(assoc_mat) else panels$secondary,
        show_row_dend = show_row_dend_opt,
        show_column_dend = show_col_dend_opt,
        row_dend_side = "left",
        row_names_side = if (show_row_dend_opt) "right" else "left",
        column_names_side = "bottom",
        column_names_rot = heatmap_x_label_angle(),
        row_names_gp = grid::gpar(fontfamily = font_family, fontsize = plot_axis_label_size(base_size, 10), col = "black"),
        column_names_gp = grid::gpar(fontfamily = font_family, fontsize = plot_axis_label_size(base_size, 10), col = "black"),
        top_annotation = panel_titles,
        row_names_max_width = grid::unit(min(60, canvas_width_mm * 0.35), "mm"),
        show_heatmap_legend = FALSE,
        height = body_height,
        column_title = NULL,
        layer_fun = if (show_values || show_marks || length(secondary_breaks)) function(j, i, x, y, width, height, fill) {
          boundaries <- which(!duplicated(j) & j %in% secondary_breaks)
          if (length(boundaries)) {
            edge <- x[boundaries] - width[boundaries] / 2
            grid::grid.segments(x0 = edge, x1 = edge, y0 = 0, y1 = 1,
                                gp = grid::gpar(col = "#595959", lwd = 1.2))
          }
          text <- labels[cbind(i, j)]
          drawn <- which(nzchar(text))
          if (!length(drawn)) return(invisible(NULL))
          rgb <- grDevices::col2rgb(fill[drawn])
          light <- as.vector(crossprod(c(0.299, 0.587, 0.114), rgb)) >= 150
          fontsize <- if (show_values) {
            # Scale labels down for narrow or short cells.
            pmin(
              max(6, base_size - 2),
              grid::convertWidth(width[drawn], "pt", valueOnly = TRUE) / (0.65 * nchar(text[drawn])),
              grid::convertHeight(height[drawn], "pt", valueOnly = TRUE) * 0.85
            )
          } else {
            sig_fontsize_pt
          }
          grid::grid.text(
            text[drawn], x = x[drawn],
            y = if (show_values) y[drawn] else y[drawn] - grid::unit(2.5, "pt"),
            gp = grid::gpar(
              fontfamily = font_family, fontsize = fontsize,
              fontface = if (show_values) "plain" else "bold",
              col = ifelse(light, "#1F2937", "#FFFFFF")
            )
          )
        } else NULL
      )
      list(heatmap = heatmap, colors = col_fun, legend = legend_params)
    }

    heatmap_legend_font_sizes <- reactive({
      base_size <- heatmap_base_size()
      c(title = max(6, base_size * 0.9), labels = max(5, base_size * 0.85))
    })

    plot_dims_px <- reactive({
      tree_space_px <- if (samples_view() && isTRUE(input$cluster_sample_columns)) {
        ceiling(sample_tree_height_pt(heatmap_base_size()) * 110 / 72) + 4L
      } else 0L
      if (identical(input$plot_size_mode, "manual")) {
        clamp_dimension <- function(value, fallback, minimum) {
          if (is.null(value) || !is.finite(value)) value <- fallback
          as.integer(max(minimum, min(2400, value)))
        }
        return(list(
          width_px = clamp_dimension(input$plot_width, 900L, 400L),
          height_px = clamp_dimension(input$plot_height, 700L, 300L) + tree_space_px
        ))
      }
      payload <- tryCatch(display_corr_payload(), error = function(e) NULL)
      if (is.null(payload)) {
        return(list(width_px = 900L, height_px = 700L))
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
      base_size <- heatmap_base_size()
      label_size <- plot_axis_label_size(base_size, 10)
      panels <- payload$column_panels
      by_sample <- isTRUE(payload$run_settings$samples)
      column_labels <- if (by_sample || is.null(panels)) colnames(assoc_mat) else panels$secondary

      # Measure label widths on a temporary device without creating a graphics file.
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off(), add = TRUE)
      text_width_px <- function(labels, size) {
        grid::pushViewport(grid::viewport(gp = grid::gpar(fontfamily = "sans", fontsize = size)))
        on.exit(grid::popViewport(), add = TRUE)
        max(grid::convertWidth(grid::stringWidth(labels), "inches", valueOnly = TRUE)) * 110
      }
      label_cap_px <- 60 / 25.4 * 110
      row_labels_px <- min(label_cap_px, text_width_px(rownames(assoc_mat), label_size))
      column_label_width <- text_width_px(column_labels, label_size)
      radians <- heatmap_x_label_angle() * pi / 180
      label_height <- label_size * 110 / 72
      column_labels_px <- min(label_cap_px, column_label_width * sin(radians) + label_height * cos(radians))
      legend_fonts <- heatmap_legend_font_sizes()
      legend_px <- max(
        text_width_px(heatmap_legend_title(payload), legend_fonts[["title"]]),
        text_width_px("-100.00", legend_fonts[["labels"]]) + 36
      )
      horizontal_space <- row_labels_px + legend_px + if (isTRUE(input$show_row_dendrogram)) 170 else 100
      vertical_space <- max(230, column_labels_px + if (isTRUE(input$show_col_dendrogram)) 120 else 80)
      if (by_sample && !is.null(payload$column_group)) {
        vertical_space <- vertical_space + ceiling(base_size * 1.6 * 110 / 72)
      }
      panel_gap_px <- if (is.null(panels)) 0 else 8 * (length(unique(panels$primary)) - 1L)
      label_span_px <- column_label_width * cos(radians) + label_height * sin(radians)
      target_cell_width <- max(if (by_sample) 24 else 36, ceiling(base_size * 2.8), min(120, label_span_px + 8))
      width_px <- as.integer(ceiling(max(500, n_cols * target_cell_width + horizontal_space + panel_gap_px)))
      height_px <- as.integer(ceiling(max(700, n_rows * auto_cell_height(n_rows) + vertical_space))) + tree_space_px

      list(width_px = width_px, height_px = height_px)
    })

    display_plot_width <- reactive({
      width <- plot_dims_px()$width_px
      if (samples_view() && !identical(input$plot_size_mode, "manual")) {
        # Keep the render size independent of the output container's layout.
        # The browser fits the resulting image to the available width.
        width <- min(width, 2400L)
      }
      width
    })

    draw_heatmap <- function() {
      grid::grid.newpage()
      canvas_width <- grid::convertWidth(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
      canvas_height <- grid::convertHeight(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
      spec <- heatmap_plot_spec(canvas_width)
      title_lines <- strsplit(spec$legend$title, "\n", fixed = TRUE)[[1L]]
      title_width <- max(vapply(title_lines, function(line) {
        title <- grid::textGrob(line, gp = spec$legend$title_gp)
        grid::convertWidth(grid::grobWidth(title), "mm", valueOnly = TRUE)
      }, numeric(1)))
      legend <- do.call(ComplexHeatmap::Legend, c(list(col_fun = spec$colors), spec$legend))
      vertical_width <- grid::convertWidth(grid::grobWidth(legend@grob), "mm", valueOnly = TRUE)
      vertical_height <- grid::convertHeight(grid::grobHeight(legend@grob), "mm", valueOnly = TRUE)
      legend_width <- grid::unit(max(vertical_width, title_width) + 12, "mm")
      place_right <- canvas_width - grid::convertWidth(legend_width, "mm", valueOnly = TRUE) >= 72 &&
        vertical_height + 10 <= canvas_height
      if (place_right) {
        layout <- grid::grid.layout(1, 2, widths = grid::unit.c(grid::unit(1, "null"), legend_width))
      } else {
        legend <- do.call(ComplexHeatmap::Legend, c(
          list(col_fun = spec$colors, direction = "horizontal",
               legend_width = grid::unit(min(60, canvas_width - 16), "mm")),
          spec$legend
        ))
        legend_height <- grid::grobHeight(legend@grob) + grid::unit(12, "mm")
        layout <- grid::grid.layout(2, 1, heights = grid::unit.c(grid::unit(1, "null"), legend_height))
      }

      grid::pushViewport(grid::viewport(layout = layout, name = "association_heatmap_layout"))
      grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1,
                                       clip = "on", name = "association_heatmap_panel"))
      left_pad_mm <- if (isTRUE(input$show_row_dendrogram)) 16 else 10
      ComplexHeatmap::draw(spec$heatmap, newpage = FALSE, show_heatmap_legend = FALSE,
                           padding = grid::unit(c(10, 4, 6, left_pad_mm), "mm"))
      grid::seekViewport("association_heatmap_layout")
      grid::pushViewport(grid::viewport(layout.pos.row = if (place_right) 1 else 2,
                                       layout.pos.col = if (place_right) 2 else 1,
                                       clip = "on", name = "association_heatmap_legend"))
      ComplexHeatmap::draw(legend, x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"),
                           just = c("center", "center"))
      grid::upViewport(0)
      invisible(NULL)
    }

    output$corr_heatmap_plot <- renderPlot(
      draw_heatmap(),
      width = function() display_plot_width(),
      height = function() plot_dims_px()$height_px,
      res = 110
    )
    # Draw while the result panel is still hidden so the "Preparing the heatmap"
    # prompt stays visible until the finished image is ready.
    outputOptions(output, "corr_heatmap_plot", suspendWhenHidden = FALSE)

    output$heatmap_status <- renderText({
      req(display_corr_payload())
      payload <- display_corr_payload()
      by_sample <- isTRUE(payload$run_settings$samples)
      sample_detection <- by_sample && identical(input$association_mode, "prevalence")
      primary_group_label <- if (is.null(input$group_var) || !nzchar(input$group_var)) "(Not selected)" else input$group_var
      current_secondary_var <- secondary_var_resolved()
      secondary_group_label <- if (is.null(current_secondary_var)) "(None)" else current_secondary_var
      analysis_scope_label <- if (by_sample) {
        if (is.null(payload$column_group)) {
          "Each column is one sample. Primary groups define panels."
        } else {
          "Each column is one sample. Primary groups define panels with labeled secondary subpanels."
        }
      } else if (is.null(current_secondary_var)) {
        "Primary groups are used directly."
      } else if (!is.null(payload$column_panels)) {
        "Secondary groups are analyzed separately within each displayed primary group."
      } else {
        "All primary levels are pooled for the continuous secondary variable."
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
          paste0("Column view: ", if (by_sample) "Individual Samples" else "Group Summary"),
          paste0("Association type: ", payload$assoc_type),
          paste0("Primary group variable: ", primary_group_label),
          paste0("Secondary group variable: ", secondary_group_label),
          paste0("Current analysis scope: ", analysis_scope_label),
          if (!by_sample) paste0("Analyzed group levels/axis: ", analyzed_group_levels) else NULL,
          paste0("Samples used: ", payload$n_samples),
          paste0("Taxa shown: ", n_taxa_display),
          paste0("Association mode: ", if (identical(input$association_mode, "prevalence")) "Prevalence-based (binary)" else "Abundance-based"),
          paste0("Transform: ", if (identical(input$association_mode, "prevalence")) "Not applied (binary prevalence)" else if (identical(input$transform_method, "clr")) "CLR" else "TSS"),
          if (!by_sample) paste0("Method: ", if (identical(input$association_mode, "prevalence")) "Pearson/phi for correlation context" else tools::toTitleCase(input$corr_method)) else NULL,
          paste0(
            "Heatmap scale: ",
            if (sample_detection) {
              "Binary detection (0 = absent, 1 = present)"
            } else if (identical(input$association_mode, "prevalence")) {
              "Raw (delta prevalence for categorical groups)"
            } else if (length(tss_value_scale_label(payload))) {
              tss_value_scale_label(payload)
            } else if (identical(input$value_scale, "zscore")) {
              "Z-score (by taxa)"
            } else {
              "Raw"
            }
          ),
          paste0("Color palette: ", heatmap_palette()$label),
          if (sample_detection) {
            "Color scale: absent to present"
          } else if (cell_values_are_percent(payload)) {
            if (by_sample) "Color scale: zero to maximum relative abundance" else
              "Color scale: zero to maximum mean relative abundance"
          } else if (length(tss_value_scale_label(payload))) {
            "Color scale: zero to maximum transformed value"
          } else {
            paste0("Color scale abs limit: ", if (is.finite(as.numeric(input$abs_limit)) && as.numeric(input$abs_limit) > 0) format(as.numeric(input$abs_limit), digits = 3) else "Auto (98th percentile)")
          },
          paste0("Plot size mode: ", if (identical(input$plot_size_mode, "manual")) "Manual" else "Auto"),
          if (identical(input$plot_size_mode, "manual")) paste0("Plot width (px): ", plot_dims_px()$width_px) else NULL,
          if (identical(input$plot_size_mode, "manual")) paste0("Plot height (px): ", plot_dims_px()$height_px) else NULL,
          if (!by_sample) paste0("Significance marks: ", if (isTRUE(input$show_significance_marks)) "Shown" else "Hidden") else NULL,
          paste0("Cell values: ", if (isTRUE(input$show_cell_values)) "Shown" else "Hidden"),
          if (!by_sample) paste0("Significance filter: ", if (isTRUE(input$apply_sig_filter)) "On" else "Off") else NULL,
          if (payload$sig_tested) paste0("Significance test method: ", test_method_label) else NULL,
          if (!by_sample) paste0("FDR cutoff (q): ", format(as.numeric(input$q_cutoff), digits = 3)) else NULL,
          if (payload$sig_tested) paste0("Significance note: ", payload$sig_note) else
            if (by_sample) "Significance tests do not apply to individual sample values." else "Significance tests not run.",
          if (payload$sig_tested) paste0("Significance tests: ", payload$n_tests) else NULL,
          if (payload$sig_tested) paste0("Significant cells: ", payload$n_sig) else NULL,
          paste0("Taxonomic level: ", input$tax_level),
          if (!by_sample && identical(analysis_target(), "taxa_group")) paste0("Group type: ", tools::toTitleCase(input$group_type)) else NULL,
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
      box_width <- display_plot_width()
      if (is.null(box_width) || !is.finite(box_width)) {
        box_width <- 600
      }
      tags$div(
        class = "app-result-card app-status-card",
        style = paste0("width: ", box_width, "px; max-width: 100%;"),
        verbatimTextOutput(session$ns("heatmap_status"))
      )
    })

    output$heatmap_legend_box <- renderUI({
      req_plot_legend_tab(session, "Association Heatmap")
      payload <- display_corr_payload()
      box_width <- display_plot_width()
      if (is.null(box_width) || !is.finite(box_width)) {
        box_width <- 600
      }
      tags$div(
        style = paste(sprintf("width: %spx;", box_width), "min-width: 0; height: auto; overflow: visible;"),
        tags$div(
          style = "color: #1f2937; font-size: 12.5px; line-height: 1.55; white-space: normal; overflow-wrap: anywhere;",
          heatmap_figure_legend(payload)
        )
      )
    })
    outputOptions(output, "heatmap_legend_box", suspendWhenHidden = FALSE, priority = -10)

    output$heatmap_status_separator <- renderUI({
      box_width <- display_plot_width()
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

    heatmap_figure_legend <- function(payload) {
      settings <- payload$run_settings
      by_sample <- isTRUE(settings$samples)
      assoc_type_label <- if (!is.null(payload) && "assoc_type" %in% names(payload)) payload$assoc_type else "N/A"
      significance_text <- if (by_sample) NULL else
        "Significance tests are calculated when significance marks or FDR masking is enabled."
      if (payload$sig_tested) {
        significance_text <- payload$sig_note
        if (isTRUE(input$show_significance_marks)) {
          significance_text <- paste0(
            "An asterisk (*) marks a cell with a Benjamini-Hochberg-adjusted p-value (q) at or below ",
            format(payload$q_cutoff), ". ", significance_text
          )
        }
      }
      if (payload$sig_tested && isTRUE(input$apply_sig_filter)) {
        significance_text <- paste0(
          significance_text,
          " Cells above the FDR cutoff (q = ", format(payload$q_cutoff),
          ") and cells without a valid q-value are masked."
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
      transformed_tss <- length(tss_value_scale_label(payload)) > 0L
      tss_description <- if (transformed_tss) {
        quantity <- if (by_sample) "sample relative abundance" else "mean relative abundance within each group"
        if (identical(payload$run_settings$value_scale, "hellinger")) {
          paste0("Colors show the square root of ", quantity, " expressed as a proportion.")
        } else {
          paste0("Colors show log10(x + 1), where x is ", quantity, " expressed as a percentage.")
        }
      }
      body_text <- if (by_sample) {
        measurement <- if (identical(settings$association_mode, "prevalence")) {
          "Cell color shows whether the taxon was detected in that sample."
        } else if (identical(settings$value_scale, "zscore")) {
          "Cell color shows the abundance Z-score for each taxon across the included samples."
        } else if (cell_values_are_percent(payload)) {
          "Cell color shows relative abundance in that sample."
        } else if (transformed_tss) {
          tss_description
        } else {
          "Cell color shows CLR-transformed abundance in that sample."
        }
        paste(
          "Each column is one sample, labeled by SampleID.",
          paste0("Panels show ", settings$group_var, " groups."),
          if (!is.null(payload$column_group)) paste0("Labeled subpanels show ", payload$column_group, " groups within each panel.") else "",
          measurement,
          "Panels share the same color scale. Gray marks unavailable values."
        )
      } else if (!is.null(payload$column_panels)) {
        cell_summary <- if (identical(input$association_mode, "prevalence")) {
          "Each cell shows the difference in prevalence between the indicated secondary group and the other secondary groups in that panel."
        } else if (transformed_tss) {
          tss_description
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
      } else if (cell_values_are_percent(payload)) {
        paste0("Cells show mean relative abundance by ", group_label,
               ". The color scale starts at zero and ends at the highest displayed mean.")
      } else if (transformed_tss) {
        tss_description
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
          if (isTRUE(input$show_cell_values)) {
            tags$div(
              if (cell_values_are_percent(payload)) {
                if (by_sample) "Cell labels show sample relative abundance as a percentage, rounded to one decimal place." else
                  "Cell labels show mean relative abundance as a percentage, rounded to one decimal place."
              } else {
                "Cell labels show the plotted values rounded to one decimal place."
              }
            )
          },
          if (!is.null(significance_text)) tags$div(style = "margin-top: 4px;", significance_text)
        )
      )
    }

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
        draw_heatmap()
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
