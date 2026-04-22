## UI
mod_alpha_ui <- function(id) {
  ns <- NS(id)
  tagList(
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
        selectInput(
          ns("plot_type"),
          "Plot Type:",
          choices = c("Boxplot" = "boxplot", "Barplot (mean ± SE)" = "barplot"),
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
        checkboxInput(ns("use_ggpattern"), "Use ggpattern (Barplot)", value = FALSE),
        
        hr(),
        h4(icon("sliders"), "P-value Options"),
        tags$style(HTML(sprintf(
          "#%s { margin-bottom: 2px; } #%s { margin-top: 0; margin-bottom: 4px; }",
          ns("show_p_val"), ns("only_sig")
        ))),
        checkboxInput(ns("show_p_val"), "Show P-value Comparison Bars", value = TRUE),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("show_p_val")),
          checkboxInput(ns("only_sig"), "Show Only Significant (p < 0.05)", value = TRUE)
        ),
        
        selectInput(ns("stat_method"), "Statistical Method:",
                    choices = c("Wilcoxon (Non-parametric)" = "wilcox.test", 
                                "T-test (Parametric)" = "t.test"),
                    selected = "wilcox.test"),
        selectInput(
          ns("p_adjust_method"),
          "Pairwise Multiple Testing Correction:",
          choices = c(
            "None" = "none",
            "Holm" = "holm",
            "Benjamini-Hochberg (FDR)" = "BH",
            "Bonferroni" = "bonferroni"
          ),
          selected = "BH"
        ),
        
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Plot Dimensions"),
        numericInput(ns("plot_width"), "Plot Width:", value = 800, min = 300, step = 50),
        numericInput(ns("plot_height"), "Plot Height:", value = 600, min = 300, step = 50),
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
        plotOutput(ns("alpha_plot_out"), height = "auto"), 
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
      selectInput(session$ns("group_var"), "Primary Grouping Variable:", 
                  choices = group_choices, selected = group_choices[1]) 
    })
    
    output$secondary_group_selector <- renderUI({
      req(input$group_var)
      resolved_primary <- resolve_meta_colname(input$group_var, meta_cols())
      group_choices <- setdiff(meta_cols(), c("SampleID", input$group_var, resolved_primary))
      selectInput(session$ns("secondary_group_var"), "Secondary Grouping Variable (Optional):",
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
      valid_cols <- if(is_secondary) c(primary_col, secondary_col) else primary_col
      valid_samples <- complete.cases(meta_data[, valid_cols, drop=FALSE])
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
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(plot.margin = ggplot2::margin(t = 26, r = 8, b = 8, l = 8))

      if (plot_type == "barplot") {
        if (use_pattern) {
          summary_group_cols <- c("Alpha_Index", x_axis_col)
          if (is_secondary) {
            summary_group_cols <- c(summary_group_cols, primary_col)
          }
          summary_df <- alpha_long %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(summary_group_cols))) %>%
            dplyr::summarise(
              mean_val = mean(Value, na.rm = TRUE),
              se_val = stats::sd(Value, na.rm = TRUE) / sqrt(dplyr::n()),
              .groups = "drop"
            )
          summary_df$se_val[!is.finite(summary_df$se_val)] <- 0

          p <- p +
            ggpattern::geom_col_pattern(
              data = summary_df,
              mapping = ggplot2::aes(
                x = .data[[x_axis_col]],
                y = mean_val,
                fill = .data[[x_axis_col]],
                pattern = .data[[x_axis_col]],
                pattern_angle = .data[[x_axis_col]]
              ),
              inherit.aes = FALSE,
              color = "black",
              linewidth = 0.25,
              alpha = 0.85,
              width = 0.7,
              pattern_fill = "white",
              pattern_colour = "#222222",
              pattern_density = 0.01,
              pattern_spacing = 0.03
            ) +
            ggplot2::geom_errorbar(
              data = summary_df,
              mapping = ggplot2::aes(
                x = .data[[x_axis_col]],
                ymin = pmax(mean_val - se_val, 0),
                ymax = mean_val + se_val
              ),
              inherit.aes = FALSE,
              width = 0.2,
              linewidth = 0.4
            ) +
            ggpattern::scale_pattern_manual(values = pattern_values, drop = FALSE) +
            ggpattern::scale_pattern_angle_manual(values = pattern_angle_values, drop = FALSE) +
            ggplot2::guides(pattern = "none", pattern_angle = "none")
        } else {
          p <- p +
            ggplot2::stat_summary(fun = mean, geom = "col", alpha = 0.8, width = 0.7, color = "black", linewidth = 0.25) +
            ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 0.2, linewidth = 0.4)
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
        current_label <- if (input$only_sig) "p.signif" else "p.format"
        p_adjust_method <- input$p_adjust_method
        if (is.null(p_adjust_method) || !p_adjust_method %in% c("none", "holm", "BH", "bonferroni")) {
          p_adjust_method <- "none"
        }
        
        if (input$only_sig) {
          final_comps <- list()
          for (comp in all_comps) {
            is_sig_anywhere <- FALSE
            for (idx in unique(alpha_long$Alpha_Index)) {
              sub_data <- alpha_long[alpha_long$Alpha_Index == idx, , drop = FALSE]
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
              if (p_adjust_method != "none") {
                ok <- is.finite(comp_pvals)
                if (any(ok)) {
                  comp_pvals[ok] <- stats::p.adjust(comp_pvals[ok], method = p_adjust_method)
                }
              }
              comp_idx <- which(vapply(all_comps, function(this_comp) identical(this_comp, comp), logical(1)))
              if (length(comp_idx) == 1) {
                p_val <- comp_pvals[comp_idx]
                if (!is.na(p_val) && p_val < 0.05) {
                  is_sig_anywhere <- TRUE
                  break
                }
              }
            }
            if (is_sig_anywhere) final_comps[[length(final_comps) + 1]] <- comp
          }
        } else {
          final_comps <- all_comps
        }

        if (length(final_comps) > 0) {
          stat_args <- list(
            comparisons = final_comps,
            method = input$stat_method,
            label = current_label,
            step.increase = 0.14,
            label.y.npc = 1.22,
            hide.ns = input$only_sig,
            digits = 3,
            size = 3.2
          )
          if (p_adjust_method != "none") {
            stat_args$p.adjust.method <- p_adjust_method
          }
          p <- p + do.call(ggpubr::stat_compare_means, stat_args)
        }
      }
      return(
        p + ggplot2::theme(
          legend.position = "none",
          strip.background = ggplot2::element_rect(fill = "white")
        )
      )
    })
    
    output$alpha_plot_out <- renderPlot({ alpha_plot_reactive() },
                                        height = function() { input$plot_height },
                                        width = function() { input$plot_width })
    
    output$rarefy_size_text <- renderText({
      req(rarefaction_size())
      paste("Rarefaction depth:", scales::comma(rarefaction_size()))
    })
    
    output$download_alpha_plot <- downloadHandler(
      filename = function() { paste0("alpha_plot_", Sys.Date(), ".png") },
      content = function(file) {
        ggplot2::ggsave(file, plot = alpha_plot_reactive(), device = "png", 
                        width = input$plot_width/72, height = input$plot_height/72, dpi = 300)
      }
    )
    
    output$download_alpha_data <- downloadHandler(
      filename = function() { paste0("alpha_data_", Sys.Date(), ".tsv") },
      content = function(file) {
        utils::write.table(alpha_long_reactive(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
  })
}
