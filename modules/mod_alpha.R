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
        
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Plot Dimensions"),
        numericInput(ns("plot_width"), "Plot Width:", value = 800, min = 300, step = 50),
        numericInput(ns("plot_height"), "Plot Height:", value = 600, min = 300, step = 50),
        
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
mod_alpha_server <- function(id, ps_obj, meta_cols) { 
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
      req(ps_obj())
      min(phyloseq::sample_sums(ps_obj()))
    })
    
    alpha_long_reactive <- reactive({
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
      req(alpha_long <- alpha_long_reactive())
      primary_col <- resolve_meta_colname(input$group_var, colnames(alpha_long))
      secondary_col <- resolve_meta_colname(input$secondary_group_var, colnames(alpha_long))
      is_secondary <- !is.null(secondary_col) && secondary_col != "none"
      x_axis_col <- if (is_secondary) secondary_col else primary_col
      group_levels <- levels(factor(alpha_long[[x_axis_col]]))
      num_groups <- length(group_levels)
      
      p <- ggplot2::ggplot(alpha_long, ggplot2::aes_string(x = x_axis_col, y = "Value", fill = x_axis_col)) +
        ggplot2::geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.7) +
        ggplot2::geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
        ggplot2::theme_bw() +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.45)))
      
      if (is_secondary) {
        p <- p + ggplot2::facet_grid(stats::as.formula(paste("Alpha_Index ~", primary_col)), scales = "free_y")
      } else {
        p <- p + ggplot2::facet_wrap(~ Alpha_Index, scales = "free_y", ncol = 2)
      }
      
      if (input$show_p_val && num_groups >= 2) {
        all_comps <- utils::combn(group_levels, 2, simplify = FALSE)
        current_label <- if (input$only_sig) "p.signif" else "p.format"
        
        if (input$only_sig) {
          final_comps <- list()
          for (comp in all_comps) {
            is_sig_anywhere <- FALSE
            for (idx in unique(alpha_long$Alpha_Index)) {
              sub_data <- alpha_long[alpha_long$Alpha_Index == idx, ]
              d1 <- sub_data$Value[sub_data[[x_axis_col]] == comp[1]]
              d2 <- sub_data$Value[sub_data[[x_axis_col]] == comp[2]]
              if (length(d1) > 1 && length(d2) > 1) {
                p_val <- if (input$stat_method == "wilcox.test") stats::wilcox.test(d1, d2)$p.value else stats::t.test(d1, d2)$p.value
                if (!is.na(p_val) && p_val < 0.05) { is_sig_anywhere <- TRUE; break }
              }
            }
            if (is_sig_anywhere) final_comps[[length(final_comps) + 1]] <- comp
          }
        } else {
          final_comps <- all_comps
        }
        
        if (length(final_comps) > 0) {
          p <- p + ggpubr::stat_compare_means(
            comparisons = final_comps,
            method = input$stat_method,
            label = current_label,
            step.increase = 0.08,
            label.y.npc = 0.96,
            hide.ns = input$only_sig,
            digits = 3
          )
        }
        
        p <- p + ggpubr::stat_compare_means(
          method = if (num_groups > 2) "kruskal.test" else input$stat_method,
          label.x.npc = "center",
          label.y.npc = 1.03,
          size = 3
        )
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
