## UI
mod_alpha_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML("
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
        h4(icon("circle-nodes"), "Alpha Diversity"),
        hr(),
        uiOutput(ns("local_group_selector")),
        uiOutput(ns("secondary_group_selector")), 
        hr(),
        
        h4(icon("list-check"), "Index Selection"),
        tags$style(HTML(sprintf(
          "#%s > label.control-label { display: none; margin: 0; } #%s { margin-bottom: 4px; }",
          ns("alpha_methods"), ns("alpha_methods")
        ))),
        checkboxGroupInput(ns("alpha_methods"), "",
                           choices = c(
                             "Observed" = "Observed",
                             "Chao1" = "Chao1",
                             "Shannon" = "Shannon",
                             "Simpson (1-D)" = "Simpson"
                           ),
                           selected = c("Chao1", "Shannon")),
        hr(),
        h4(icon("sliders"), "Plot Settings"),
        selectInput(
          ns("plot_type"),
          "Plot Type:",
          choices = c("Box plot" = "boxplot", "Bar plot" = "barplot"),
          selected = "boxplot"
        ),
        selectInput(
          ns("color_palette"),
          "Color Palette:",
          choices = c(
            "Set2" = "set2",
            "Dark2" = "dark2",
            "Paired" = "paired",
            "Grayscale" = "gray"
          ),
          selected = "set2"
        ),
        checkboxInput(ns("use_ggpattern"), "Use ggpattern (Bar plot)", value = FALSE),
        
        hr(),
        h4(icon("sliders"), "P-value Options"),
        tags$style(HTML(sprintf(
          "#%s { margin-bottom: 2px; } #%s, #%s { margin-top: 0; margin-bottom: 4px; }",
          ns("show_p_val"), ns("show_sig_only"), ns("show_sig_as_marks")
        ))),
        checkboxInput(ns("show_p_val"), "Show p-value bars", value = TRUE),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("show_p_val")),
          checkboxInput(ns("show_sig_only"), "Show only p-value < 0.05", value = TRUE),
          checkboxInput(ns("show_sig_as_marks"), "Show as significant marks", value = TRUE)
        ),
        
        selectInput(ns("stat_method"), "Statistical Method:",
                    choices = c("Wilcoxon" = "wilcox.test", 
                                "T-test" = "t.test"),
                    selected = "wilcox.test"),
        selectInput(
          ns("p_adjust_method"),
          "Pairwise Multiple Testing Correction:",
          choices = c(
            "None" = "none",
            "Holm" = "holm",
            "BH" = "BH",
            "Bonferroni" = "bonferroni"
          ),
          selected = "BH"
        ),
        
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Plot Dimensions"),
        numericInput(ns("plot_width"), "Plot Width:", value = 600, min = 300, step = 50),
        numericInput(ns("plot_height"), "Plot Height:", value = 500, min = 300, step = 50),
        numericInput(ns("base_size"), "Base Font Size:", value = 11, min = 6, max = 30, step = 1),
        
        hr(),
        h5(icon("download"), "Download"),
        div(
          style = "display: flex; gap: 4px; flex-wrap: nowrap;",
          downloadButton(
            ns("download_alpha_plot"),
            "Download Plot (PNG)",
            style = "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
          ),
          downloadButton(
            ns("download_alpha_data"),
            "Download Data (TSV)",
            style = "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
          )
        )
      ),
      mainPanel(
        h4("Alpha Diversity"),
        plotOutput(ns("alpha_plot_out"), height = "auto"), 
        uiOutput(ns("alpha_legend_box")),
        br(),
        textOutput(ns("rarefy_size_text"))
      )
    )
  )
}

## Server
mod_alpha_server <- function(id, ps_obj, meta_cols, active_tab = NULL) { 
  moduleServer(id, function(input, output, session) {
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

    is_active_tab <- reactive({
      if (is.null(active_tab)) {
        return(TRUE)
      }
      current_tab <- tryCatch(active_tab(), error = function(e) NULL)
      identical(current_tab, "Alpha Diversity")
    })
    
    output$local_group_selector <- renderUI({
      req(meta_cols())
      group_choices <- setdiff(meta_cols(), "SampleID")
      selectInput(session$ns("group_var"), "Primary Group:", 
                  choices = group_choices, selected = group_choices[1]) 
    })
    
    output$secondary_group_selector <- renderUI({
      req(meta_cols())
      primary_input <- input$group_var
      if (is.null(primary_input) || !nzchar(primary_input)) {
        group_base <- setdiff(meta_cols(), "SampleID")
        primary_input <- if (length(group_base) > 0) group_base[1] else NULL
      }
      resolved_primary <- resolve_meta_colname(primary_input, meta_cols())
      group_choices <- setdiff(meta_cols(), c("SampleID", primary_input, resolved_primary))
      selectInput(session$ns("secondary_group_var"), "Secondary Group (Optional):",
                  choices = c("(None)" = "none", group_choices), selected = "none")
    })
    
    rarefaction_size <- reactive({
      req(isTRUE(is_active_tab()))
      req(ps_obj())
      min(phyloseq::sample_sums(ps_obj()))
    })
    
    alpha_long_reactive <- reactive({
      req(isTRUE(is_active_tab()))
      req(ps_obj(), input$group_var, input$alpha_methods)
      physeq_rarefied <- phyloseq::rarefy_even_depth(
        ps_obj(), sample.size = rarefaction_size(), rngseed = 42, replace = FALSE, verbose = FALSE
      )
      meta_data <- as(phyloseq::sample_data(physeq_rarefied), "data.frame")
      primary_col <- resolve_meta_colname(input$group_var, colnames(meta_data))
      secondary_col <- resolve_meta_colname(input$secondary_group_var, colnames(meta_data))
      is_secondary <- !is.null(secondary_col) && secondary_col != "none"
      validate(need(primary_col %in% colnames(meta_data), paste0("Primary grouping variable '", input$group_var, "' was not found in metadata.")))
      if (is_secondary) {
        validate(need(secondary_col %in% colnames(meta_data), paste0("Secondary grouping variable '", input$secondary_group_var, "' was not found in metadata.")))
      }
      valid_cols <- if (is_secondary) c(primary_col, secondary_col) else primary_col
      valid_samples <- complete.cases(meta_data[, valid_cols, drop = FALSE])
      physeq_filtered <- phyloseq::prune_samples(valid_samples, physeq_rarefied)
      meta_filtered <- meta_data[valid_samples, , drop = FALSE]
      alpha_df <- phyloseq::estimate_richness(physeq_filtered, measures = input$alpha_methods)
      alpha_df[[primary_col]] <- factor(meta_filtered[[primary_col]])
      if (is_secondary) alpha_df[[secondary_col]] <- factor(meta_filtered[[secondary_col]])
      alpha_df$SampleID <- rownames(alpha_df)
      out_long <- tidyr::pivot_longer(
        alpha_df,
        cols = tidyselect::any_of(input$alpha_methods),
        names_to = "Alpha_Index",
        values_to = "Value"
      )
      out_long$Alpha_Index <- ifelse(
        out_long$Alpha_Index == "Simpson",
        "Simpson (1-D)",
        out_long$Alpha_Index
      )
      out_long
    })
    
    alpha_plot_reactive <- reactive({
      req(isTRUE(is_active_tab()))
      req(alpha_long <- alpha_long_reactive())
      primary_col <- resolve_meta_colname(input$group_var, colnames(alpha_long))
      secondary_col <- resolve_meta_colname(input$secondary_group_var, colnames(alpha_long))
      is_secondary <- !is.null(secondary_col) && secondary_col != "none"
      x_axis_col <- if (is_secondary) secondary_col else primary_col
      group_levels <- levels(factor(alpha_long[[x_axis_col]]))
      num_groups <- length(group_levels)
      plot_type <- input$plot_type
      if (is.null(plot_type) || !plot_type %in% c("boxplot", "barplot")) {
        plot_type <- "boxplot"
      }
      palette_key <- input$color_palette
      if (is.null(palette_key) || !palette_key %in% c("set2", "dark2", "paired", "gray")) {
        palette_key <- "set2"
      }
      base_size <- input$base_size
      if (is.null(base_size) || !is.finite(base_size)) {
        base_size <- 11
      }
      n_groups <- max(1, length(group_levels))
      fill_values <- switch(
        palette_key,
        "dark2" = grDevices::hcl.colors(n_groups, palette = "Dark 3"),
        "paired" = grDevices::hcl.colors(n_groups, palette = "Set 3"),
        "gray" = grDevices::gray.colors(n_groups, start = 0.2, end = 0.85),
        grDevices::hcl.colors(n_groups, palette = "Set 2")
      )
      names(fill_values) <- group_levels
      use_pattern <- plot_type == "barplot" && isTRUE(input$use_ggpattern) && requireNamespace("ggpattern", quietly = TRUE)
      pattern_values <- rep(c("stripe", "stripe", "crosshatch", "none"), length.out = n_groups)
      names(pattern_values) <- group_levels
      pattern_angle_values <- rep(c(45, -45, 0, 0), length.out = n_groups)
      names(pattern_angle_values) <- group_levels
      
      p <- ggplot2::ggplot(alpha_long, ggplot2::aes_string(x = x_axis_col, y = "Value", fill = x_axis_col)) +
        ggplot2::theme_bw(base_size = base_size) +
        ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
        ggplot2::labs(
          title = "Alpha Diversity Comparison",
          y = "Alpha Diversity"
        ) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = base_size + 3, face = "bold"),
          axis.title.x = ggplot2::element_text(face = "bold", size = base_size + 1),
          axis.title.y = ggplot2::element_text(face = "bold", size = base_size + 1),
          plot.margin = ggplot2::margin(t = 26, r = 8, b = 8, l = 8)
        )

      if (plot_type == "barplot") {
        facet_group_cols <- c("Alpha_Index")
        if (is_secondary) {
          facet_group_cols <- c(facet_group_cols, primary_col)
        }
        # facet_grid(Alpha_Index ~ PrimaryGroup, scales = "free_y") shares y-scale by row.
        # Use row-level minima so bar baselines align with the effective y-domain.
        scale_group_cols <- c("Alpha_Index")
        summary_group_cols <- c(facet_group_cols, x_axis_col)
        summary_df <- alpha_long %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(summary_group_cols))) %>%
          dplyr::summarise(
            mean_val = mean(Value, na.rm = TRUE),
            se_val = stats::sd(Value, na.rm = TRUE) / sqrt(dplyr::n()),
            .groups = "drop"
          )
        scale_baseline_df <- alpha_long %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(scale_group_cols))) %>%
          dplyr::summarise(
            scale_min = {
              vals <- Value[is.finite(Value)]
              if (length(vals) == 0) 0 else min(vals, na.rm = TRUE)
            },
            .groups = "drop"
          )
        summary_df <- dplyr::left_join(summary_df, scale_baseline_df, by = scale_group_cols)
        summary_df$mean_val[!is.finite(summary_df$mean_val)] <- NA_real_
        summary_df$se_val[!is.finite(summary_df$se_val)] <- 0
        summary_df$scale_min[!is.finite(summary_df$scale_min)] <- 0
        summary_df <- summary_df[is.finite(summary_df$mean_val), , drop = FALSE]

        if (use_pattern) {
          if ("geom_crossbar_pattern" %in% getNamespaceExports("ggpattern")) {
            p <- p +
              ggpattern::geom_crossbar_pattern(
                data = summary_df,
                mapping = ggplot2::aes(
                  x = .data[[x_axis_col]],
                  y = mean_val,
                  ymin = scale_min,
                  ymax = mean_val,
                  fill = .data[[x_axis_col]],
                  pattern = .data[[x_axis_col]],
                  pattern_angle = .data[[x_axis_col]]
                ),
                inherit.aes = FALSE,
                color = "black",
                linewidth = 0.25,
                alpha = 0.85,
                width = 0.7,
                fatten = 0,
                pattern_fill = "white",
                pattern_colour = "#222222",
                pattern_density = 0.01,
                pattern_spacing = 0.03
              ) +
              ggpattern::scale_pattern_manual(values = pattern_values, drop = FALSE) +
              ggpattern::scale_pattern_angle_manual(values = pattern_angle_values, drop = FALSE) +
              ggplot2::guides(pattern = "none", pattern_angle = "none")
          } else {
            p <- p +
              ggplot2::geom_crossbar(
                data = summary_df,
                mapping = ggplot2::aes(
                  x = .data[[x_axis_col]],
                  y = mean_val,
                  ymin = scale_min,
                  ymax = mean_val,
                  fill = .data[[x_axis_col]]
                ),
                inherit.aes = FALSE,
                alpha = 0.8,
                width = 0.7,
                color = "#333333",
                linewidth = 0.25,
                fatten = 0
              )
          }
          p <- p +
            ggplot2::geom_errorbar(
              data = summary_df,
              mapping = ggplot2::aes(
                x = .data[[x_axis_col]],
                ymin = pmax(mean_val - se_val, scale_min),
                ymax = mean_val + se_val
              ),
              inherit.aes = FALSE,
              width = 0.2,
              linewidth = 0.4
            )
        } else {
          p <- p +
            ggplot2::geom_crossbar(
              data = summary_df,
              mapping = ggplot2::aes(
                x = .data[[x_axis_col]],
                y = mean_val,
                ymin = scale_min,
                ymax = mean_val,
                fill = .data[[x_axis_col]]
              ),
              inherit.aes = FALSE,
              alpha = 0.8,
              width = 0.7,
              color = "#333333",
              linewidth = 0.25,
              fatten = 0
            ) +
            ggplot2::geom_errorbar(
              data = summary_df,
              mapping = ggplot2::aes(
                x = .data[[x_axis_col]],
                ymin = pmax(mean_val - se_val, scale_min),
                ymax = mean_val + se_val
              ),
              inherit.aes = FALSE,
              width = 0.2,
              linewidth = 0.4
            )
        }
        p <- p + ggplot2::geom_point(
          position = ggplot2::position_jitter(width = 0.1, height = 0),
          alpha = 0.5,
          size = 1.5
        )
      } else {
        p <- p +
          ggplot2::geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.7) +
          ggplot2::geom_jitter(width = 0.1, alpha = 0.5, size = 1.5)
      }
      
      if (is_secondary) {
        p <- p + ggplot2::facet_grid(stats::as.formula(paste("Alpha_Index ~", primary_col)), scales = "free_y")
      } else {
        p <- p + ggplot2::facet_wrap(~ Alpha_Index, scales = "free_y", ncol = 2)
      }
      
      if (input$show_p_val && num_groups >= 2) {
        all_comps <- utils::combn(group_levels, 2, simplify = FALSE)
        sig_only <- isTRUE(input$show_sig_only)
        sig_as_marks <- isTRUE(input$show_sig_as_marks)
        p_adjust_method <- input$p_adjust_method
        if (is.null(p_adjust_method) || !p_adjust_method %in% c("none", "holm", "BH", "bonferroni")) {
          p_adjust_method <- "none"
        }

        format_p_value <- function(pv) {
          if (!is.finite(pv)) return(NA_character_)
          if (pv < 0.001) return("< 0.001")
          sprintf("%.3f", pv)
        }
        signif_mark <- function(pv) {
          if (!is.finite(pv)) return("ns")
          if (pv < 0.001) return("***")
          if (pv < 0.01) return("**")
          if (pv < 0.05) return("*")
          "ns"
        }

        facet_vars <- c("Alpha_Index")
        if (is_secondary) facet_vars <- c(facet_vars, primary_col)
        facet_split <- interaction(alpha_long[, facet_vars, drop = FALSE], drop = TRUE, lex.order = TRUE)
        facet_data <- split(alpha_long, facet_split)
        annotation_rows <- vector("list", length(facet_data))
        row_idx <- 0L

        for (sub_data in facet_data) {
          if (nrow(sub_data) == 0) next
          comp_pvals <- vapply(all_comps, function(this_comp) {
            d1 <- sub_data$Value[sub_data[[x_axis_col]] == this_comp[1]]
            d2 <- sub_data$Value[sub_data[[x_axis_col]] == this_comp[2]]
            if (length(d1) <= 1 || length(d2) <= 1) {
              return(NA_real_)
            }
            tryCatch(
              if (input$stat_method == "wilcox.test") stats::wilcox.test(d1, d2)$p.value else stats::t.test(d1, d2)$p.value,
              error = function(e) NA_real_
            )
          }, numeric(1))

          adj_pvals <- comp_pvals
          if (p_adjust_method != "none") {
            ok <- is.finite(adj_pvals)
            if (any(ok)) {
              adj_pvals[ok] <- stats::p.adjust(adj_pvals[ok], method = p_adjust_method)
            }
          }

          valid_idx <- which(is.finite(adj_pvals))
          if (sig_only) valid_idx <- valid_idx[adj_pvals[valid_idx] < 0.05]
          if (length(valid_idx) == 0) next

          y_max <- suppressWarnings(max(sub_data$Value, na.rm = TRUE))
          if (!is.finite(y_max) || y_max <= 0) y_max <- 1
          base_mult <- if (sig_as_marks) 1.12 else 1.20
          step_mult <- if (sig_as_marks) 0.08 else 0.12

          ann_df <- data.frame(
            group1 = vapply(all_comps[valid_idx], `[[`, character(1), 1),
            group2 = vapply(all_comps[valid_idx], `[[`, character(1), 2),
            p.adj = adj_pvals[valid_idx],
            stringsAsFactors = FALSE
          )
          ann_df$p.signif <- vapply(ann_df$p.adj, signif_mark, character(1))
          ann_df$p.format <- vapply(ann_df$p.adj, format_p_value, character(1))
          ann_df$y.position <- y_max * (base_mult + step_mult * (seq_len(nrow(ann_df)) - 1))
          for (fv in facet_vars) {
            ann_df[[fv]] <- sub_data[[fv]][1]
          }

          row_idx <- row_idx + 1L
          annotation_rows[[row_idx]] <- ann_df
        }

        annotation_rows <- annotation_rows[seq_len(row_idx)]
        if (length(annotation_rows) > 0) {
          annotation_df <- do.call(rbind, annotation_rows)
          p <- p + ggpubr::stat_pvalue_manual(
            annotation_df,
            label = if (sig_as_marks) "p.signif" else "p.format",
            xmin = "group1",
            xmax = "group2",
            y.position = "y.position",
            tip.length = 0.01,
            size = max(2, base_size / 3.2),
            bracket.size = 0.3,
            vjust = if (sig_as_marks) 0.25 else -0.15,
            hide.ns = FALSE,
            inherit.aes = FALSE
          )
        }
      }
      return(
        p + ggplot2::theme(
          legend.position = "none",
          strip.background = ggplot2::element_rect(fill = "grey85", color = "grey20")
        )
      )
    })
    
    output$alpha_plot_out <- renderPlot({
      ps_now <- ps_obj()
      if (is.null(ps_now)) {
        graphics::plot.new()
        graphics::text(
          0.5, 0.5,
          "Applying selected samples and preparing alpha diversity plot.\nPlease wait...",
          cex = 0.85
        )
        return(invisible(NULL))
      }
      if (phyloseq::nsamples(ps_now) == 0) {
        graphics::plot.new()
        graphics::text(
          0.5, 0.5,
          "No samples are currently selected.\nPlease select at least one sample in Preprocessing.",
          cex = 0.85
        )
        return(invisible(NULL))
      }
      tryCatch(
        alpha_plot_reactive(),
        error = function(e) {
          graphics::plot.new()
          graphics::text(
            0.5, 0.5,
            paste0("Alpha diversity plot is not ready yet. Please wait...\n", conditionMessage(e)),
            cex = 0.85
          )
        }
      )
    },
    height = function() { 
      input$plot_height 
    },
    width = function() { 
      input$plot_width 
    })
    
    output$rarefy_size_text <- renderText({
      req(rarefaction_size())
      paste("Rarefaction depth:", scales::comma(rarefaction_size()))
    })
    
    output$download_alpha_plot <- downloadHandler(
      filename = function() { 
        paste0("alpha_plot_", Sys.Date(), ".png") 
      },
      content = function(file) {
        ggplot2::ggsave(file, plot = alpha_plot_reactive(), device = "png", 
                        width = input$plot_width / 72, height = input$plot_height / 72, dpi = 300)
      }
    )
    
    output$download_alpha_data <- downloadHandler(
      filename = function() { 
        paste0("alpha_data_", Sys.Date(), ".tsv") 
      },
      content = function(file) {
        utils::write.table(alpha_long_reactive(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )

    output$alpha_legend_box <- renderUI({
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
          uiOutput(session$ns("alpha_figure_legend"))
        )
      )
    })

    output$alpha_figure_legend <- renderUI({
      req(input$alpha_methods, input$plot_type)
      selected_indices <- input$alpha_methods
      if (is.null(selected_indices) || length(selected_indices) == 0) {
        selected_indices <- c("Chao1", "Shannon")
      }
      meta_data <- tryCatch(
        as(phyloseq::sample_data(ps_obj()), "data.frame"),
        error = function(e) NULL
      )
      primary_col <- if (!is.null(meta_data)) resolve_meta_colname(input$group_var, colnames(meta_data)) else NULL
      secondary_col <- if (!is.null(meta_data)) resolve_meta_colname(input$secondary_group_var, colnames(meta_data)) else NULL
      is_secondary <- !is.null(secondary_col) && secondary_col != "none" && !is.null(meta_data) && secondary_col %in% colnames(meta_data)
      x_axis_col <- if (is_secondary) secondary_col else primary_col
      n_groups <- if (!is.null(meta_data) && !is.null(x_axis_col) && x_axis_col %in% colnames(meta_data)) {
        length(unique(stats::na.omit(meta_data[[x_axis_col]])))
      } else {
        0
      }
      stat_method_label <- if (identical(input$stat_method, "t.test")) "t-test" else "Wilcoxon rank-sum test"
      p_adjust_label <- switch(
        input$p_adjust_method,
        "holm" = "Holm",
        "BH" = "Benjamini-Hochberg (BH)",
        "bonferroni" = "Bonferroni",
        "none" = "no multiple-testing correction",
        "BH"
      )
      sig_sentence <- if (n_groups >= 3) {
        paste0(
          " Pairwise significance is annotated above each comparison using ",
          stat_method_label,
          " with ",
          p_adjust_label,
          "."
        )
      } else if (n_groups == 2) {
        paste0(
          " Pairwise significance is annotated using ",
          stat_method_label,
          " (single comparison)."
        )
      } else {
        ""
      }
      index_label <- paste(selected_indices, collapse = ", ")
      plot_type_label <- if (identical(input$plot_type, "barplot")) {
        "bar plots summarize group means with standard error bars"
      } else {
        "box plots summarize group distributions with overlaid sample points"
      }
      tags$div(
        tags$div(
          style = "font-weight: 600; margin-bottom: 4px;",
          "Alpha diversity comparison plot"
        ),
        tags$div(
          paste0(
            "This figure shows alpha diversity indices ",
            index_label,
            " across groups; ",
            plot_type_label,
            ".",
            sig_sentence
          )
        )
      )
    })
  })
}
