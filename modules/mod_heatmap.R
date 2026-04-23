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
    ")),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        h4(icon("table-cells"), "Corr. Heatmap"),
        hr(),
        selectInput(
          ns("group_var"),
          "1. Group variable",
          choices = character(0),
          selected = NULL
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
          choices = c("ASV", "Genus", "Species"),
          selected = "Genus"
        ),
        selectInput(
          ns("transform_method"),
          "4. Data transform",
          choices = c("TSS" = "tss", "CLR" = "clr"),
          selected = "clr"
        ),
        selectInput(
          ns("corr_method"),
          "5. Correlation method",
          choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
          selected = "pearson"
        ),
        selectInput(
          ns("value_scale"),
          "6. Heatmap value scale",
          choices = c("Raw" = "raw", "Z-score (by taxa)" = "zscore"),
          selected = "raw"
        ),
        tags$details(
          tags$summary("Advanced options"),
          numericInput(
            ns("prevalence_filter_pct"),
            "7. Prevalence filter (%)",
            value = 10,
            min = 0,
            max = 100,
            step = 1
          ),
          numericInput(
            ns("max_taxa"),
            "8. Max taxa",
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
          numericInput(ns("abs_limit"), "Color scale abs limit (0 = auto)", value = 0, min = 0, max = 100, step = 0.1),
          checkboxInput(ns("apply_sig_filter"), "Mask non-significant cells (FDR)", value = FALSE),
          numericInput(ns("q_cutoff"), "FDR cutoff (q)", value = 0.05, min = 0.0001, max = 1, step = 0.01)
        ),
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Cell Size"),
        numericInput(ns("cell_width"), "Cell width (px)", value = 30, min = 6, max = 80, step = 1),
        numericInput(ns("cell_height"), "Cell height (px)", value = 10, min = 6, max = 80, step = 1),
        numericInput(ns("base_size"), "Base Font Size:", value = 11, min = 6, max = 30, step = 1),
        actionButton(ns("run_heatmap"), "Run Heatmap", class = "btn-danger", style = "font-size: 12px;"),
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
             if (msg.label) btn.textContent = msg.label;
           });"
        ))
      ),
      mainPanel(
        h4("Correlation Heatmap"),
        plotOutput(ns("corr_heatmap_plot"), height = "auto"),
        uiOutput(ns("heatmap_legend_box")),
        uiOutput(ns("heatmap_status_separator")),
        h5(icon("circle-info"), "Heatmap Status"),
        uiOutput(ns("heatmap_status_box"))
      )
    )
  )
}

mod_heatmap_server <- function(id, ps_obj, meta_vars = NULL) {
  moduleServer(id, function(input, output, session) {
    analysis_target <- reactive("taxa_group")

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
    })

    taxa_choices_preview <- reactive({
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

    corr_payload <- eventReactive(input$run_heatmap, {
      heatmap_running(TRUE)
      session$sendCustomMessage(
        "toggle-heatmap-run-btn",
        list(id = session$ns("run_heatmap"), disabled = TRUE, label = "Running...")
      )
      on.exit(
        {
          heatmap_running(FALSE)
          session$sendCustomMessage(
            "toggle-heatmap-run-btn",
            list(id = session$ns("run_heatmap"), disabled = FALSE, label = "Run Heatmap")
          )
        },
        add = TRUE
      )
      req(ps_obj(), input$tax_level, input$transform_method, input$corr_method, input$value_scale)
      ps_data <- ps_obj()
      validate(
        need(phyloseq::nsamples(ps_data) >= 3, "At least 3 samples are required."),
        need(phyloseq::ntaxa(ps_data) >= 3, "At least 3 taxa are required.")
      )

      if (!identical(input$tax_level, "ASV")) {
        validate(
          need(!is.null(phyloseq::tax_table(ps_data, errorIfNULL = FALSE)), "Taxonomy table is required for taxonomic aggregation.")
        )
        tax_cols <- colnames(phyloseq::tax_table(ps_data))
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

      taxa_var <- apply(transformed, 1, stats::var, na.rm = TRUE)
      keep_var <- is.finite(taxa_var) & taxa_var > 0
      transformed <- transformed[keep_var, , drop = FALSE]
      validate(
        need(nrow(transformed) >= 3, "Too few taxa remain after variance filtering.")
      )

      meta_df <- metadata_df()
      sample_ids <- colnames(transformed)
      if (nrow(meta_df) > 0) {
        meta_df <- meta_df[sample_ids, , drop = FALSE]
      }

      p_mat <- NULL
      q_mat <- NULL
      sig_supported <- FALSE
      sig_note <- "Not computed."
      row_hclust <- NULL
      col_hclust <- NULL

      if (identical(analysis_target(), "taxa_taxa")) {
        assoc_mat <- stats::cor(t(transformed), method = input$corr_method, use = "pairwise.complete.obs")
        assoc_type <- "Taxa vs Taxa correlation"
        sig_supported <- TRUE
        sig_note <- "Computed from pairwise correlation tests (BH-adjusted)."

        n_tax <- nrow(transformed)
        p_mat <- matrix(NA_real_, nrow = n_tax, ncol = n_tax, dimnames = list(rownames(transformed), rownames(transformed)))
        diag(p_mat) <- 0
        if (n_tax >= 2) {
          for (i in seq_len(n_tax - 1)) {
            xi <- transformed[i, ]
            for (j in (i + 1):n_tax) {
              xj <- transformed[j, ]
              p_val <- tryCatch(
                suppressWarnings(stats::cor.test(xi, xj, method = input$corr_method, exact = FALSE)$p.value),
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

        group_vec <- meta_df[[input$group_var]]
        names(group_vec) <- rownames(meta_df)
        group_vec <- group_vec[sample_ids]

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
          assoc_vals <- apply(transformed, 1, function(taxa_x) {
            suppressWarnings(stats::cor(taxa_x, group_num, method = input$corr_method, use = "pairwise.complete.obs"))
          })
          assoc_mat <- matrix(as.numeric(assoc_vals), ncol = 1)
          rownames(assoc_mat) <- names(assoc_vals)
          colnames(assoc_mat) <- input$group_var
          assoc_type <- "Taxa vs continuous group correlation"
          sig_supported <- TRUE
          sig_note <- "Computed from taxa-group correlation tests (BH-adjusted)."

          p_vals <- apply(transformed, 1, function(taxa_x) {
            tryCatch(
              suppressWarnings(stats::cor.test(taxa_x, group_num, method = input$corr_method, exact = FALSE)$p.value),
              error = function(e) NA_real_
            )
          })
          q_vals <- rep(NA_real_, length(p_vals))
          ok <- is.finite(p_vals)
          if (any(ok)) {
            q_vals[ok] <- stats::p.adjust(p_vals[ok], method = "BH")
          }
          p_mat <- matrix(as.numeric(p_vals), ncol = 1, dimnames = list(names(p_vals), input$group_var))
          q_mat <- matrix(as.numeric(q_vals), ncol = 1, dimnames = list(names(q_vals), input$group_var))
        } else {
          group_chr <- as.character(group_vec)
          group_chr[is.na(group_chr) | !nzchar(group_chr)] <- "(Missing)"
          levels_keep <- unique(group_chr)
          validate(
            need(length(levels_keep) >= 2, "Categorical group needs at least 2 levels.")
          )
          assoc_mat <- sapply(levels_keep, function(lv) {
            idx <- which(group_chr == lv)
            rowMeans(transformed[, idx, drop = FALSE], na.rm = TRUE)
          })
          if (!is.matrix(assoc_mat)) {
            assoc_mat <- matrix(assoc_mat, ncol = length(levels_keep))
          }
          rownames(assoc_mat) <- rownames(transformed)
          colnames(assoc_mat) <- levels_keep
          assoc_type <- "Taxa vs categorical group level mean"
          sig_supported <- TRUE
          sig_note <- "Computed per cell by one-vs-rest Wilcoxon rank-sum test (taxa x group level), BH-adjusted."

          group_fac <- factor(group_chr, levels = levels_keep)
          p_mat <- matrix(
            NA_real_,
            nrow = nrow(assoc_mat),
            ncol = ncol(assoc_mat),
            dimnames = list(rownames(assoc_mat), colnames(assoc_mat))
          )
          for (j in seq_along(levels_keep)) {
            lv <- levels_keep[j]
            in_group <- which(group_fac == lv)
            out_group <- which(group_fac != lv)
            if (length(in_group) < 2 || length(out_group) < 2) {
              next
            }
            p_vec <- apply(transformed, 1, function(taxa_x) {
              tryCatch(
                stats::wilcox.test(
                  x = as.numeric(taxa_x[in_group]),
                  y = as.numeric(taxa_x[out_group]),
                  exact = FALSE
                )$p.value,
                error = function(e) NA_real_
              )
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

      assoc_mat[is.na(assoc_mat)] <- 0
      assoc_mat[!is.finite(assoc_mat)] <- 0
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

      if (identical(input$value_scale, "zscore")) {
        z_by_row <- t(apply(assoc_mat, 1, function(row_x) {
          row_sd <- stats::sd(row_x, na.rm = TRUE)
          if (!is.finite(row_sd) || row_sd <= 0) {
            return(rep(0, length(row_x)))
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

      if (isTRUE(input$show_col_dendrogram) && ncol(assoc_mat) >= 3) {
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
        heat_df = heat_df,
        n_samples = ncol(transformed),
        n_taxa = nrow(transformed),
        assoc_type = assoc_type,
        sig_supported = sig_supported,
        sig_note = sig_note,
        row_hclust = row_hclust,
        col_hclust = col_hclust,
        n_tests = sum(is.finite(heat_df$QValue)),
        n_sig = sum(is.finite(heat_df$QValue) & heat_df$QValue <= q_cutoff, na.rm = TRUE)
      )
    }, ignoreNULL = TRUE)

    observeEvent(corr_payload(), {
      taxa_choices <- rownames(corr_payload()$assoc_mat)
      current_selected <- isolate(input$selected_taxa)
      selected_keep <- if (length(current_selected) > 0) intersect(current_selected, taxa_choices) else character(0)
      updateSelectizeInput(session, "selected_taxa", choices = taxa_choices, selected = selected_keep, server = TRUE)
    }, ignoreInit = FALSE)

    cell_dims_px <- reactive({
      cell_width <- as.numeric(input$cell_width)
      if (!is.finite(cell_width)) cell_width <- 30
      cell_width <- max(6, min(80, cell_width))

      cell_height <- as.numeric(input$cell_height)
      if (!is.finite(cell_height)) cell_height <- 10
      cell_height <- max(6, min(80, cell_height))

      list(cell_width = cell_width, cell_height = cell_height)
    })

    heatmap_plot_obj <- reactive({
      req(corr_payload())
      heat_df <- corr_payload()$heat_df
      assoc_mat <- corr_payload()$assoc_mat
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

      fill_limit <- max(abs(heat_df$Association), na.rm = TRUE)
      if (!is.finite(fill_limit) || fill_limit <= 0) fill_limit <- 1
      abs_limit_input <- as.numeric(input$abs_limit)
      if (is.finite(abs_limit_input) && abs_limit_input > 0) {
        fill_limit <- abs_limit_input
      }

      base_size <- as.numeric(input$base_size)
      if (!is.finite(base_size)) base_size <- 11
      base_size <- max(6, min(30, base_size))
      sig_fontsize_pt <- max(8, base_size * 0.95)

      dims_cell <- cell_dims_px()
      cell_width_pt <- dims_cell$cell_width * 0.75
      cell_height_pt <- dims_cell$cell_height * 0.75

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
        c("#2166AC", "#F7F7F7", "#B2182B")
      )

      show_row_dend_opt <- isTRUE(input$show_row_dendrogram)
      show_col_dend_opt <- isTRUE(input$show_col_dendrogram)

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

      ComplexHeatmap::Heatmap(
        assoc_mat,
        name = "Value",
        col = col_fun,
        na_col = "#D9D9D9",
        cluster_rows = cluster_rows_opt,
        cluster_columns = cluster_cols_opt,
        show_row_dend = show_row_dend_opt,
        show_column_dend = show_col_dend_opt,
        row_dend_side = "left",
        row_names_side = if (show_row_dend_opt) "right" else "left",
        column_names_side = "bottom",
        column_names_rot = 90,
        row_names_gp = grid::gpar(fontsize = max(6, base_size - 2)),
        column_names_gp = grid::gpar(fontsize = max(6, base_size - 2)),
        heatmap_legend_param = list(
          title = "Value",
          title_gp = grid::gpar(fontsize = base_size),
          labels_gp = grid::gpar(fontsize = max(6, base_size - 1))
        ),
        width = grid::unit(ncol(assoc_mat) * cell_width_pt, "pt"),
        height = grid::unit(nrow(assoc_mat) * cell_height_pt, "pt"),
        column_title = NULL,
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (identical(sig_mat[i, j], "*")) {
            grid::grid.text("*", x = x, y = y - grid::unit(2.5, "pt"), gp = grid::gpar(fontsize = sig_fontsize_pt))
          }
        }
      )
    })

    plot_dims_px <- reactive({
      payload <- tryCatch(corr_payload(), error = function(e) NULL)
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
      width_px <- as.integer((n_cols * cell_width) + 340 + left_extra_px)
      height_px <- as.integer((n_rows * cell_height) + 230)

      list(width_px = width_px, height_px = height_px)
    })

    output$corr_heatmap_plot <- renderPlot(
      {
        if (is.null(input$run_heatmap) || input$run_heatmap < 1) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Click 'Run Heatmap' to start analysis.")
          return(invisible(NULL))
        }
        if (isTRUE(heatmap_running())) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Heatmap analysis is running. Please wait...")
          return(invisible(NULL))
        }
        ht <- heatmap_plot_obj()
        grid::grid.newpage()
        left_pad_mm <- if (isTRUE(input$show_row_dendrogram)) 16 else 10
        top_pad_mm <- 10
        ComplexHeatmap::draw(
          ht,
          heatmap_legend_side = "right",
          padding = grid::unit(c(top_pad_mm, 6, 6, left_pad_mm), "mm")
        )
      },
      width = function() plot_dims_px()$width_px,
      height = function() plot_dims_px()$height_px,
      res = 110
    )

    output$heatmap_status <- renderText({
      req(corr_payload())
      payload <- corr_payload()
      test_method_label <- if (identical(input$group_type, "continuous")) {
        paste0("Correlation test (", tools::toTitleCase(input$corr_method), ") + BH correction")
      } else if (identical(input$group_type, "categorical")) {
        "One-vs-rest Wilcoxon rank-sum test (taxa x group) + BH correction"
      } else {
        "Auto: continuous uses correlation test, categorical uses one-vs-rest Wilcoxon + BH correction"
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
          paste0("Samples used: ", payload$n_samples),
          paste0("Taxa shown: ", n_taxa_display),
          paste0("Transform: ", if (identical(input$transform_method, "clr")) "CLR" else "TSS"),
          paste0("Method: ", tools::toTitleCase(input$corr_method)),
          paste0("Heatmap scale: ", if (identical(input$value_scale, "zscore")) "Z-score (by taxa)" else "Raw"),
          paste0("Color scale abs limit: ", if (is.finite(as.numeric(input$abs_limit)) && as.numeric(input$abs_limit) > 0) format(as.numeric(input$abs_limit), digits = 3) else "Auto"),
          paste0("Cell width (px): ", format(as.numeric(input$cell_width), digits = 3)),
          paste0("Cell height (px): ", format(as.numeric(input$cell_height), digits = 3)),
          paste0("Auto plot width (px): ", plot_dims_px()$width_px),
          paste0("Auto plot height (px): ", plot_dims_px()$height_px),
          paste0("Significance mark size: linked to base font size (", format(as.numeric(input$base_size), digits = 3), ")"),
          paste0("Significance filter: ", if (isTRUE(input$apply_sig_filter)) "On" else "Off"),
          paste0("Significance test method: ", test_method_label),
          paste0("FDR cutoff (q): ", format(as.numeric(input$q_cutoff), digits = 3)),
          paste0("Significance note: ", payload$sig_note),
          if (isTRUE(payload$sig_supported)) paste0("Significance tests: ", payload$n_tests) else NULL,
          if (isTRUE(payload$sig_supported)) paste0("Significant cells (*): ", payload$n_sig) else NULL,
          paste0("Taxonomic level: ", input$tax_level),
          if (identical(analysis_target(), "taxa_group")) paste0("Group type: ", tools::toTitleCase(input$group_type)) else NULL,
          if (identical(analysis_target(), "taxa_group")) paste0("Group variable: ", input$group_var) else NULL
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
      req(input$tax_level, input$transform_method, input$corr_method, input$value_scale)
      group_label <- if (is.null(input$group_var) || !nzchar(input$group_var)) "selected group variable" else input$group_var
      title_text <- "Association heatmap"
      body_text <- paste0(
        "Cells represent taxa-to-group association values at ",
        tolower(input$tax_level),
        " level. The analysis uses ",
        toupper(input$transform_method),
        " transformed abundance and ",
        tools::toTitleCase(input$corr_method),
        " statistics. Color indicates association direction and magnitude, with value scale set to ",
        if (identical(input$value_scale, "zscore")) "row-wise z-score" else "raw association value",
        ". Group variable is ",
        group_label,
        if (isTRUE(input$apply_sig_filter)) {
          paste0(", and non-significant cells are masked at FDR q < ", input$q_cutoff, ".")
        } else {
          "."
        }
      )
      tags$div(
        tags$div(style = "font-weight: 600; margin-bottom: 4px;", title_text),
        tags$div(body_text)
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
        req(corr_payload())
        dims <- plot_dims_px()
        export_scale <- 2
        grDevices::png(
          filename = file,
          width = as.integer(dims$width_px * export_scale),
          height = as.integer(dims$height_px * export_scale),
          units = "px",
          res = as.integer(110 * export_scale),
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
          padding = grid::unit(c(top_pad_mm, 6, 6, left_pad_mm), "mm")
        )
      }
    )

    output$download_corr_tsv <- downloadHandler(
      filename = function() paste0("association_matrix_", analysis_target(), "_", Sys.Date(), ".tsv"),
      content = function(file) {
        req(corr_payload())
        mat_df <- as.data.frame(corr_payload()$assoc_mat, stringsAsFactors = FALSE)
        mat_df <- cbind(Taxon = rownames(mat_df), mat_df, stringsAsFactors = FALSE)
        utils::write.table(mat_df, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    )
  })
}
