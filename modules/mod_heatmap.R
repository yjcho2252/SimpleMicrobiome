library(shiny)
library(phyloseq)
library(ggplot2)

mod_heatmap_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        width = 2,
        h4(icon("table-cells"), "Correlation Heatmap"),
        hr(),
        selectInput(
          ns("analysis_target"),
          "1. Analysis target",
          choices = c("Taxa vs Taxa" = "taxa_taxa", "Taxa vs Group" = "taxa_group"),
          selected = "taxa_taxa"
        ),
        uiOutput(ns("group_var_ui")),
        uiOutput(ns("group_type_ui")),
        selectInput(
          ns("tax_level"),
          "2. Taxonomic level",
          choices = c("ASV", "Genus", "Species"),
          selected = "Genus"
        ),
        selectInput(
          ns("transform_method"),
          "3. Data transform",
          choices = c("TSS" = "tss", "CLR" = "clr"),
          selected = "clr"
        ),
        selectInput(
          ns("corr_method"),
          "4. Correlation method",
          choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
          selected = "spearman"
        ),
        selectInput(
          ns("value_scale"),
          "5. Heatmap value scale",
          choices = c("Raw" = "raw", "Z-score (by taxa)" = "zscore"),
          selected = "zscore"
        ),
        numericInput(
          ns("prevalence_filter_pct"),
          "6. Prevalence filter cutoff (%)",
          value = 10,
          min = 0,
          max = 100,
          step = 1
        ),
        numericInput(
          ns("max_taxa"),
          "7. Max taxa",
          value = 80,
          min = 10,
          max = 300,
          step = 5
        ),
        checkboxInput(ns("show_diagonal"), "Show diagonal (Taxa vs Taxa only)", value = FALSE),
        checkboxInput(ns("cluster_order"), "Cluster row/column order", value = TRUE),
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Plot Dimensions"),
        numericInput(ns("plot_width"), "Plot width (px)", value = 900, min = 400, max = 2400, step = 50),
        numericInput(ns("plot_height"), "Plot height (px)", value = 780, min = 300, max = 2400, step = 50),
        actionButton(ns("run_heatmap"), "Run Heatmap", class = "btn-danger", style = "font-size: 12px;")
      ),
      mainPanel(
        h4("Correlation Heatmap"),
        tags$div(
          style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
          downloadButton(ns("download_heatmap_png"), "Download Heatmap (PNG)"),
          downloadButton(ns("download_corr_tsv"), "Download Matrix (TSV)")
        ),
        plotOutput(ns("corr_heatmap_plot"), height = "auto"),
        verbatimTextOutput(ns("heatmap_status"))
      )
    )
  )
}

mod_heatmap_server <- function(id, ps_obj, meta_vars = NULL) {
  moduleServer(id, function(input, output, session) {
    sanitize_taxa_names <- function(ps_data, rank_name) {
      tt <- phyloseq::tax_table(ps_data, errorIfNULL = FALSE)
      if (is.null(tt) || !rank_name %in% colnames(tt)) {
        return(ps_data)
      }
      tax_df <- as.data.frame(tt, stringsAsFactors = FALSE)
      rank_values <- as.character(tax_df[[rank_name]])
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

    output$group_var_ui <- renderUI({
      req(input$analysis_target)
      if (!identical(input$analysis_target, "taxa_group")) {
        return(NULL)
      }
      meta_df <- metadata_df()
      candidate_cols <- colnames(meta_df)
      candidate_cols <- setdiff(candidate_cols, "SampleID")
      selectInput(
        session$ns("group_var"),
        "Group variable (for Taxa vs Group)",
        choices = candidate_cols,
        selected = if (length(candidate_cols) > 0) candidate_cols[1] else NULL
      )
    })

    output$group_type_ui <- renderUI({
      req(input$analysis_target)
      if (!identical(input$analysis_target, "taxa_group")) {
        return(NULL)
      }
      selectInput(
        session$ns("group_type"),
        "Group type",
        choices = c("Auto" = "auto", "Continuous" = "continuous", "Categorical" = "categorical"),
        selected = "auto"
      )
    })

    corr_payload <- eventReactive(input$run_heatmap, {
      req(ps_obj(), input$tax_level, input$transform_method, input$corr_method, input$analysis_target, input$value_scale)
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
      if (is.na(max_taxa) || max_taxa < 3) max_taxa <- 80
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

      if (identical(input$analysis_target, "taxa_taxa")) {
        assoc_mat <- stats::cor(t(transformed), method = input$corr_method, use = "pairwise.complete.obs")
        assoc_type <- "Taxa vs Taxa correlation"
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
        is_numeric_group <- if (identical(group_type_selected, "continuous")) {
          TRUE
        } else if (identical(group_type_selected, "categorical")) {
          FALSE
        } else {
          is.numeric(group_vec) || is.integer(group_vec)
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
        }
      }

      assoc_mat[is.na(assoc_mat)] <- 0
      assoc_mat[!is.finite(assoc_mat)] <- 0

      if (identical(input$analysis_target, "taxa_taxa")) {
        assoc_mat <- pmax(-1, pmin(1, assoc_mat))
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

      if (isTRUE(input$cluster_order) && nrow(assoc_mat) >= 3) {
        row_dist <- stats::dist(assoc_mat, method = "euclidean")
        row_ord <- tryCatch(
          stats::hclust(row_dist, method = "average")$order,
          error = function(e) seq_len(nrow(assoc_mat))
        )
        assoc_mat <- assoc_mat[row_ord, , drop = FALSE]

        if (ncol(assoc_mat) >= 3) {
          col_dist <- stats::dist(t(assoc_mat), method = "euclidean")
          col_ord <- tryCatch(
            stats::hclust(col_dist, method = "average")$order,
            error = function(e) seq_len(ncol(assoc_mat))
          )
          assoc_mat <- assoc_mat[, col_ord, drop = FALSE]
        }
      }

      heat_df <- as.data.frame(as.table(assoc_mat), stringsAsFactors = FALSE)
      colnames(heat_df) <- c("AxisX", "AxisY", "Association")
      if (identical(input$analysis_target, "taxa_taxa") && !isTRUE(input$show_diagonal)) {
        heat_df <- heat_df[heat_df$AxisX != heat_df$AxisY, , drop = FALSE]
      }

      list(
        assoc_mat = assoc_mat,
        heat_df = heat_df,
        n_samples = ncol(transformed),
        n_taxa = nrow(transformed),
        assoc_type = assoc_type
      )
    }, ignoreNULL = FALSE)

    heatmap_plot <- reactive({
      req(corr_payload())
      heat_df <- corr_payload()$heat_df
      assoc_mat <- corr_payload()$assoc_mat
      validate(
        need(nrow(heat_df) > 0, "No cells to plot with current options.")
      )

      heat_df$AxisX <- factor(heat_df$AxisX, levels = rownames(assoc_mat))
      heat_df$AxisY <- factor(heat_df$AxisY, levels = colnames(assoc_mat))
      fill_limit <- max(abs(heat_df$Association), na.rm = TRUE)
      if (!is.finite(fill_limit) || fill_limit <= 0) fill_limit <- 1

      if (identical(input$analysis_target, "taxa_group")) {
        heat_df$x_plot <- factor(as.character(heat_df$AxisY), levels = colnames(assoc_mat))
        heat_df$y_plot <- factor(as.character(heat_df$AxisX), levels = rownames(assoc_mat))
      } else {
        heat_df$x_plot <- heat_df$AxisX
        heat_df$y_plot <- heat_df$AxisY
      }

      ggplot2::ggplot(heat_df, ggplot2::aes(x = x_plot, y = y_plot, fill = Association)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
          low = "#2166AC",
          mid = "#F7F7F7",
          high = "#B2182B",
          midpoint = 0,
          limits = c(-fill_limit, fill_limit)
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y = ggplot2::element_text(size = 7),
          panel.grid = ggplot2::element_blank()
        ) +
        ggplot2::labs(
          title = paste0(
            "Heatmap (",
            if (identical(input$transform_method, "clr")) "CLR" else "TSS",
            " + ",
            tools::toTitleCase(input$corr_method), ", ",
            if (identical(input$value_scale, "zscore")) "Z-score" else "Raw",
            ")"
          ),
          x = if (identical(input$analysis_target, "taxa_taxa")) "Taxa" else "Group / Level",
          y = "Taxa",
          fill = "Value"
        ) +
        ggplot2::coord_fixed()
    })

    output$corr_heatmap_plot <- renderPlot(
      {
        heatmap_plot()
      },
      width = function() input$plot_width,
      height = function() input$plot_height,
      res = 110
    )

    output$heatmap_status <- renderText({
      req(corr_payload())
      payload <- corr_payload()
      paste(
        c(
          paste0("Analysis target: ", if (identical(input$analysis_target, "taxa_taxa")) "Taxa vs Taxa" else "Taxa vs Group"),
          paste0("Association type: ", payload$assoc_type),
          paste0("Samples used: ", payload$n_samples),
          paste0("Taxa shown: ", payload$n_taxa),
          paste0("Transform: ", if (identical(input$transform_method, "clr")) "CLR" else "TSS"),
          paste0("Method: ", tools::toTitleCase(input$corr_method)),
          paste0("Heatmap scale: ", if (identical(input$value_scale, "zscore")) "Z-score (by taxa)" else "Raw"),
          paste0("Taxonomic level: ", input$tax_level),
          if (identical(input$analysis_target, "taxa_group")) paste0("Group type: ", tools::toTitleCase(input$group_type)) else NULL,
          if (identical(input$analysis_target, "taxa_group")) paste0("Group variable: ", input$group_var) else NULL
        ),
        collapse = "\n"
      )
    })

    output$download_heatmap_png <- downloadHandler(
      filename = function() {
        paste0(
          "association_heatmap_",
          input$analysis_target,
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
        dpi_val <- 300
        width_in <- input$plot_width / dpi_val
        height_in <- input$plot_height / dpi_val
        ggplot2::ggsave(
          filename = file,
          plot = heatmap_plot(),
          device = "png",
          width = width_in,
          height = height_in,
          units = "in",
          dpi = dpi_val
        )
      }
    )

    output$download_corr_tsv <- downloadHandler(
      filename = function() paste0("association_matrix_", input$analysis_target, "_", Sys.Date(), ".tsv"),
      content = function(file) {
        req(corr_payload())
        mat_df <- as.data.frame(corr_payload()$assoc_mat, stringsAsFactors = FALSE)
        mat_df <- cbind(Taxon = rownames(mat_df), mat_df, stringsAsFactors = FALSE)
        utils::write.table(mat_df, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    )
  })
}
