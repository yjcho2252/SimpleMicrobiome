library(shiny)
library(phyloseq)
library(ggplot2)
library(dplyr)

mod_biplot_ui <- function(id) {
  ns <- NS(id)
  tagList(
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
          ns("distance_metric"),
          "3. Distance metric",
          choices = c("Bray-Curtis" = "bray", "Aitchison" = "aitchison"),
          selected = "bray"
        ),
        numericInput(
          ns("prevalence_filter_pct"),
          "4. Prevalence filter cutoff (%)",
          value = 10,
          min = 0,
          max = 100,
          step = 1
        ),
        numericInput(
          ns("max_taxa"),
          "5. Max taxa",
          value = 100,
          min = 10,
          max = 500,
          step = 10
        ),
        numericInput(
          ns("top_taxa_vectors"),
          "6. Top taxa vectors",
          value = 15,
          min = 0,
          max = 100,
          step = 1
        ),
        checkboxInput(ns("show_taxa_vectors"), "Show taxa vectors", value = TRUE),
        checkboxInput(ns("show_group_vectors"), "Show group vectors", value = TRUE),
        checkboxInput(ns("show_group_centroid"), "Show group centroids", value = TRUE),
        checkboxInput(ns("show_sample_names"), "Show sample names", value = FALSE),
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Plot Dimensions"),
        numericInput(ns("plot_width"), "Plot width (px)", value = 900, min = 400, max = 2400, step = 50),
        numericInput(ns("plot_height"), "Plot height (px)", value = 700, min = 300, max = 2400, step = 50),
        actionButton(ns("run_biplot"), "Run Biplot", class = "btn-danger", style = "font-size: 12px;")
      ),
      mainPanel(
        h4("Association Biplot"),
        tags$div(
          style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
          downloadButton(ns("download_biplot_png"), "Download Biplot (PNG)"),
          downloadButton(ns("download_biplot_scores"), "Download Scores (TSV)")
        ),
        plotOutput(ns("biplot_plot"), height = "auto"),
        verbatimTextOutput(ns("biplot_status"))
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

    biplot_payload <- eventReactive(input$run_biplot, {
      req(ps_obj(), input$group_var, input$tax_level, input$distance_metric)
      validate(
        need(requireNamespace("vegan", quietly = TRUE), "Package 'vegan' is required for dbRDA.")
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

      if (identical(input$distance_metric, "aitchison")) {
        log_mat <- log(otu + 1)
        otu_use <- apply(log_mat, 2, function(x) x - mean(x, na.rm = TRUE))
        if (!is.matrix(otu_use)) otu_use <- matrix(otu_use, nrow = nrow(log_mat))
        rownames(otu_use) <- rownames(log_mat)
        colnames(otu_use) <- colnames(log_mat)
      } else {
        totals <- colSums(otu, na.rm = TRUE)
        totals[totals <= 0 | is.na(totals)] <- 1
        otu_use <- sweep(otu, 2, totals, "/")
      }

      taxa_var <- apply(otu_use, 1, stats::var, na.rm = TRUE)
      otu_use <- otu_use[is.finite(taxa_var) & taxa_var > 0, , drop = FALSE]
      validate(need(nrow(otu_use) >= 3, "Too few taxa remain after variance filtering."))

      sd <- phyloseq::sample_data(ps_data, errorIfNULL = FALSE)
      validate(need(!is.null(sd), "Metadata is required for biplot."))
      meta_df <- as.data.frame(sd, stringsAsFactors = FALSE)
      validate(need(input$group_var %in% colnames(meta_df), "Selected group variable is not available in metadata."))

      sample_ids <- colnames(otu_use)
      meta_df <- meta_df[sample_ids, , drop = FALSE]
      group_raw <- meta_df[[input$group_var]]
      group_chr <- as.character(group_raw)
      keep_samples <- !is.na(group_raw) & nzchar(group_chr)
      validate(need(sum(keep_samples) >= 4, "At least 4 samples with valid group values are required."))
      otu_use <- otu_use[, keep_samples, drop = FALSE]
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
      comm_df <- as.data.frame(t(otu_use), stringsAsFactors = FALSE)
      validate(
        need(nrow(model_df) == nrow(comm_df), "Model data and community matrix row counts do not match."),
        need(all(rownames(model_df) == rownames(comm_df)), "Sample IDs between model data and community matrix do not match.")
      )

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

      anova_all <- tryCatch(vegan::anova.cca(model, permutations = 999), error = function(e) NULL)
      p_all <- if (!is.null(anova_all) && nrow(as.data.frame(anova_all)) >= 1) as.data.frame(anova_all)[1, "Pr(>F)"] else NA_real_

      site_scores <- vegan::scores(model, display = "sites", choices = 1:2, scaling = 2)
      species_scores <- tryCatch(vegan::scores(model, display = "species", choices = 1:2, scaling = 2), error = function(e) NULL)
      biplot_scores <- tryCatch(vegan::scores(model, display = "bp", choices = 1:2, scaling = 2), error = function(e) NULL)

      site_df <- as.data.frame(site_scores, stringsAsFactors = FALSE)
      validate(
        need(ncol(site_df) >= 2, "dbRDA did not return 2D site scores.")
      )
      colnames(site_df)[1:2] <- c("Axis1", "Axis2")
      site_df$SampleID <- rownames(site_df)
      site_df$Group <- as.character(model_df$group_model)

      species_df <- data.frame()
      if (!is.null(species_scores)) {
        species_df <- as.data.frame(species_scores, stringsAsFactors = FALSE)
        if (ncol(species_df) >= 2) {
          colnames(species_df)[1:2] <- c("Axis1", "Axis2")
          species_df$Taxon <- rownames(species_df)
          species_df$Length <- sqrt(species_df$Axis1^2 + species_df$Axis2^2)
          species_df <- species_df[order(species_df$Length, decreasing = TRUE), , drop = FALSE]
          n_top <- as.integer(input$top_taxa_vectors)
          if (is.na(n_top) || n_top < 0) n_top <- 15
          if (n_top == 0) {
            species_df <- species_df[0, , drop = FALSE]
          } else if (nrow(species_df) > n_top) {
            species_df <- species_df[seq_len(n_top), , drop = FALSE]
          }
        } else {
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
          dplyr::mutate(Variable = paste0("Group: ", Group))
      }

      eig <- vegan::eigenvals(model)
      axis1_var <- if (length(eig) >= 1) round(100 * eig[1] / sum(eig), 1) else NA_real_
      axis2_var <- if (length(eig) >= 2) round(100 * eig[2] / sum(eig), 1) else NA_real_

      list(
        site_df = site_df,
        species_df = species_df,
        bp_df = bp_df,
        group_arrow_df = group_arrow_df,
        metric_label = if (identical(input$distance_metric, "aitchison")) "Aitchison" else "Bray-Curtis",
        axis1_var = axis1_var,
        axis2_var = axis2_var,
        model_p = as.numeric(p_all),
        n_samples = nrow(site_df),
        n_taxa = nrow(otu_use)
      )
    }, ignoreNULL = FALSE)

    biplot_plot <- reactive({
      req(biplot_payload())
      payload <- biplot_payload()
      site_df <- payload$site_df
      taxa_df <- payload$species_df
      bp_df <- payload$bp_df
      group_arrow_df <- payload$group_arrow_df

      p <- ggplot2::ggplot(site_df, ggplot2::aes(x = Axis1, y = Axis2, color = Group)) +
        ggplot2::geom_point(size = 3, alpha = 0.9) +
        ggplot2::theme_bw() +
        ggplot2::labs(
          title = paste0("Association Biplot (dbRDA, ", payload$metric_label, ")"),
          x = paste0("CAP1", if (is.finite(payload$axis1_var)) paste0(" (", payload$axis1_var, "%)") else ""),
          y = paste0("CAP2", if (is.finite(payload$axis2_var)) paste0(" (", payload$axis2_var, "%)") else ""),
          color = input$group_var
        )

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
            size = 3,
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
            size = 3,
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
              size = 3.2,
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
              size = 3.2,
              vjust = -0.4
            )
        }
      }

      if (isTRUE(input$show_sample_names)) {
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_text_repel(
            data = site_df,
            ggplot2::aes(x = Axis1, y = Axis2, label = SampleID, color = Group),
            size = 2.8,
            show.legend = FALSE,
            box.padding = 0.2,
            point.padding = 0.15
          )
        } else {
          p <- p + ggplot2::geom_text(
            data = site_df,
            ggplot2::aes(x = Axis1, y = Axis2, label = SampleID, color = Group),
            size = 2.8,
            hjust = -0.1,
            show.legend = FALSE
          )
        }
      }

      p
    })

    output$biplot_plot <- renderPlot(
      {
        biplot_plot()
      },
      width = function() input$plot_width,
      height = function() input$plot_height,
      res = 110
    )

    output$biplot_status <- renderText({
      req(biplot_payload())
      payload <- biplot_payload()
      pval_txt <- if (is.finite(payload$model_p)) format(payload$model_p, digits = 3) else "NA"
      paste(
        c(
          paste0("Distance metric: ", payload$metric_label),
          paste0("Group variable: ", input$group_var),
          paste0("Samples used: ", payload$n_samples),
          paste0("Taxa used: ", payload$n_taxa),
          paste0("Top taxa vectors: ", nrow(payload$species_df)),
          paste0("dbRDA model p-value (perm): ", pval_txt)
        ),
        collapse = "\n"
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
