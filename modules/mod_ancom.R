mod_ancom_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        
        selectInput(ns("group_var"), "1. Metadata grouping variable", choices = NULL),
        
        selectizeInput(ns("control_groups"), "2. Control Group Levels", choices = NULL, multiple = TRUE, 
                       options = list(placeholder = 'Select metadata levels for Control Group', 
                                      plugins = list('remove_button'))), 
        
        selectizeInput(ns("comparison_groups"), "3. Comparison Group Levels", choices = NULL, multiple = TRUE,
                       options = list(placeholder = 'Select metadata levels for Comparison Group',
                                      plugins = list('remove_button'))), 
        uiOutput(ns("group_sample_counts")),
        hr(),
        
        selectInput(ns("tax_level"), "4. Taxonomic level",
                    choices = c("ASV", "Genus", "Species"), selected = "Genus"),
        
        selectInput(ns("volcano_y_axis"), "5. Volcano Plot Y-axis Metric",
                    choices = c("FDR-adjusted p-value (q-value)" = "q_val", 
                                "Raw p-value (p-value)" = "p_val"),
                    selected = "q_val"),
        
        numericInput(ns("plot_height"), "Plot height (px)", value = 600, min = 300, max = 2000, step = 50),
        numericInput(ns("plot_width"), "Plot width (px)", value = 900, min = 400, max = 2000, step = 50),
        actionButton(ns("run_ancom_btn"), "Run ANCOM-BC2", class = "btn-danger")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Table",
                   downloadButton(ns("download_ancom_table"), "Download Table (TSV)"),
                   DTOutput(ns("ancom_table"))),
          tabPanel("Volcano Plot",
                   downloadButton(ns("download_volcano"), "Download Plot (PNG)"),
                   plotOutput(ns("ancom_plot"), height = "auto")),
          tabPanel("Bar Plot",
                   downloadButton(ns("download_barplot"), "Download Plot (PNG)"),
                   plotOutput(ns("ancom_barplot"), height = "auto"))
        )
      )
    )
  )
}

mod_ancom_server <- function(id, ps_obj) {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(ps_obj(), {
      req(ps_obj())
      meta_df <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
      meta_cols <- colnames(meta_df)
      group_choices <- setdiff(meta_cols, "SampleID")
      selected_group <- if (length(group_choices) > 0) group_choices[1] else NULL
      
      updateSelectInput(session, "group_var",
                        choices = group_choices,
                        selected = selected_group)
    }, ignoreNULL = FALSE)
    
    observeEvent(list(ps_obj(), input$group_var), {
      req(ps_obj(), input$group_var)
      meta_df <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
      validate(
        need(input$group_var %in% colnames(meta_df), paste0("ERROR: Metadata variable '", input$group_var, "' not found."))
      )
      
      level_choices <- sort(unique(as.character(meta_df[[input$group_var]])))
      level_choices <- level_choices[!is.na(level_choices) & nzchar(level_choices)]
      
      updateSelectizeInput(session, "control_groups",
                           choices = level_choices,
                           selected = NULL,
                           server = TRUE)
      updateSelectizeInput(session, "comparison_groups",
                           choices = level_choices,
                           selected = NULL,
                           server = TRUE)
    }, ignoreNULL = FALSE)
    
    selected_group_labels <- reactive({
      req(input$control_groups, input$comparison_groups)
      list(
        control = paste(input$control_groups, collapse = " + "),
        comparison = paste(input$comparison_groups, collapse = " + ")
      )
    })
    
    group_selection_info <- reactive({
      req(ps_obj(), input$group_var)
      ps <- ps_obj()
      meta_df <- as.data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE)
      validate(
        need(input$group_var %in% colnames(meta_df), paste0("ERROR: Metadata variable '", input$group_var, "' not found in sample_data."))
      )
      
      control_levels <- if (is.null(input$control_groups)) character(0) else input$control_groups
      comparison_levels <- if (is.null(input$comparison_groups)) character(0) else input$comparison_groups
      overlapping_levels <- intersect(control_levels, comparison_levels)
      
      group_values <- as.character(meta_df[[input$group_var]])
      sample_ids <- rownames(meta_df)
      control_ids <- sample_ids[group_values %in% control_levels]
      comparison_ids <- sample_ids[group_values %in% comparison_levels]
      
      list(
        control_levels = control_levels,
        comparison_levels = comparison_levels,
        overlapping_levels = overlapping_levels,
        control_ids = control_ids,
        comparison_ids = comparison_ids
      )
    })
    
    output$group_sample_counts <- renderUI({
      req(ps_obj(), input$group_var)
      info <- group_selection_info()
      
      tagList(
        h5("Selected Sample Counts"),
        tags$div(paste0("Control: ", length(info$control_ids), " sample(s)")),
        tags$div(paste0("Comparison: ", length(info$comparison_ids), " sample(s)")),
        if (length(info$overlapping_levels) > 0) {
          tags$div(
            style = "color:#b30000; font-weight:600;",
            paste0("Overlapping levels: ", paste(info$overlapping_levels, collapse = ", "))
          )
        }
      )
    })
    
    ps_filtered <- reactive({
      req(ps_obj(), input$group_var, input$control_groups, input$comparison_groups, input$tax_level)
      
      info <- group_selection_info()
      control_levels <- info$control_levels
      comparison_levels <- info$comparison_levels
      
      validate(
        need(length(control_levels) > 0, "ERROR: Select at least one metadata level for the Control Group."),
        need(length(comparison_levels) > 0, "ERROR: Select at least one metadata level for the Comparison Group.")
      )
      
      validate(
        need(length(info$overlapping_levels) == 0, 
             paste("ERROR: Levels overlap. The following levels are in both groups:", paste(info$overlapping_levels, collapse = ", ")))
      )
      
      ps <- ps_obj()
      control_ids <- info$control_ids
      comparison_ids <- info$comparison_ids
      validate(
        need(length(control_ids) > 0, "ERROR: No samples found for selected Control group levels."),
        need(length(comparison_ids) > 0, "ERROR: No samples found for selected Comparison group levels.")
      )
      target_ids <- c(control_ids, comparison_ids)
      
      ps_sub <- phyloseq::prune_samples(target_ids, ps)
      
      current_group_var <- "UserGroup"
      
      sample_groups <- character(length(target_ids))
      sample_groups[phyloseq::sample_names(ps_sub) %in% control_ids] <- "Control"
      sample_groups[phyloseq::sample_names(ps_sub) %in% comparison_ids] <- "Comparison"
      
      group_factor <- factor(sample_groups, levels = c("Control", "Comparison"))
      
      phyloseq::sample_data(ps_sub)[[current_group_var]] <- group_factor
      
      if (input$tax_level != "ASV") {
        tax_rank <- input$tax_level
        validate(
          need(!is.null(phyloseq::tax_table(ps_sub)), "ERROR: Taxonomic table missing for aggregation.")
        )
        tax_cols <- colnames(phyloseq::tax_table(ps_sub))
        validate(
          need(tax_rank %in% tax_cols, paste("ERROR: Taxonomic rank", tax_rank, "not found in taxonomy table."))
        )
        ps_sub <- phyloseq::tax_glom(ps_sub, taxrank = tax_rank, NArm = TRUE)
        
        validate(
          need(phyloseq::ntaxa(ps_sub) > 0, paste("ERROR: After aggregating to", tax_rank, ", no taxa remain."))
        )
        
        tax_table_sub <- phyloseq::tax_table(ps_sub)
        if (!is.null(tax_table_sub) && tax_rank %in% colnames(tax_table_sub)) {
          tax_df <- as.data.frame(tax_table_sub)
          rank_values <- as.character(tax_df[[tax_rank]])
          rank_values[is.na(rank_values) | rank_values == ""] <- phyloseq::taxa_names(ps_sub)[is.na(rank_values) | rank_values == ""]
          phyloseq::taxa_names(ps_sub) <- make.unique(rank_values)
        }
      }
      
      ps_sub
    })
    
    ancom_res <- eventReactive(input$run_ancom_btn, {
      req(ps_filtered())
      
      # ⭐ 수정: current_group_var를 명확히 'UserGroup'으로 정의 ⭐
      current_group_var <- "UserGroup" 
      
      result <- tryCatch({
        withProgress(message = 'Running ANCOM-BC2...', value = 0, {
          
          ps_current <- ps_filtered()
          tax_level_arg <- if (input$tax_level == "ASV") NULL else input$tax_level
          
          out <- ANCOMBC::ancombc2(data = ps_current,
                                   tax_level = tax_level_arg,
                                   # ⭐ 수정: fix_formula에 그룹 변수를 문자열로 지정 ⭐
                                   fix_formula = current_group_var,
                                   # ⭐ 수정: group 인자에 그룹 변수를 문자열로 지정 (에러 해결) ⭐
                                   group = current_group_var,
                                   struc_zero = TRUE, neg_lb = TRUE,
                                   alpha = 0.05, n_cl = 1)
          
          res <- out$res
          res$feature_id <- res$taxon
          
          tax_table_current <- phyloseq::tax_table(ps_current)
          if (!is.null(tax_table_current)) {
            tax_df <- as.data.frame(tax_table_current)
            tax_df$feature_id <- rownames(tax_df)
            res <- dplyr::left_join(res, tax_df, by = "feature_id")
          }
          
          preferred_order <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
          res$taxa_label <- NA_character_
          for (col in preferred_order) {
            if (col %in% colnames(res)) {
              idx <- is.na(res$taxa_label) & !is.na(res[[col]]) & res[[col]] != ""
              res$taxa_label[idx] <- as.character(res[[col]][idx])
            }
          }
          missing_idx <- is.na(res$taxa_label) | res$taxa_label == ""
          res$taxa_label[missing_idx] <- res$feature_id[missing_idx]
          
          res$taxon <- res$taxa_label
          
          rownames(res) <- res$feature_id
          return(res)
        })
      }, error = function(e) {
        showNotification(paste("ANCOM-BC2 Execution Failed:", e$message), type = "error", duration = NULL)
        NULL 
      })
      
      validate(
        need(!is.null(result), "ANCOM-BC2 analysis failed. Check the error message above.")
      )
      
      return(result)
    })
    
    get_metric_column <- function(res, prefixes, suffixes, fallbacks = character()) {
      candidates <- unlist(lapply(prefixes, function(pref) paste0(pref, suffixes)))
      candidates <- unique(c(candidates, fallbacks))
      match <- candidates[candidates %in% colnames(res)]
      if (length(match) == 0) return(NA_character_)
      match[1]
    }
    
    ancom_processed <- reactive({
      req(ancom_res())
      res <- ancom_res()
      
      current_group_var <- "UserGroup" 
      comparison_group_key <- "Comparison"
      comparison_group_label <- selected_group_labels()$comparison
      
      suffixes <- unique(c(
        paste0(current_group_var, comparison_group_key),
        paste0(current_group_var, "_", comparison_group_key),
        paste0(current_group_var, comparison_group_key, "V"),
        paste0(current_group_var, "_", comparison_group_key, "V"),
        comparison_group_key,
        paste0(comparison_group_key, "V")
      ))
      
      lfc_col <- get_metric_column(res, c("lfc_"), suffixes, fallbacks = c("lfc"))
      se_col <- get_metric_column(res, c("se_"), suffixes)
      W_col <- get_metric_column(res, c("W_"), suffixes)
      p_col <- get_metric_column(res, c("p_"), suffixes, fallbacks = c("p_val", "pvalue"))
      q_col <- get_metric_column(res, c("q_", "q_val_"), suffixes, fallbacks = c("q_val", "padj"))
      diff_col <- get_metric_column(res, c("diff_"), suffixes, fallbacks = c("diff", "diff_abn"))
      
      validate(
        need(!is.na(lfc_col), "LFC column not found in ANCOM results."),
        need(!is.na(p_col) || !is.na(q_col), "P-value or q-value columns not found in ANCOM results.")
      )
      
      to_numeric <- function(values) suppressWarnings(as.numeric(values))
      
      lfc_values <- to_numeric(res[[lfc_col]])
      se_values <- if (!is.na(se_col)) to_numeric(res[[se_col]]) else rep(NA_real_, nrow(res))
      W_values <- if (!is.na(W_col)) to_numeric(res[[W_col]]) else rep(NA_real_, nrow(res))
      p_values <- if (!is.na(p_col)) to_numeric(res[[p_col]]) else rep(NA_real_, nrow(res))
      q_values <- if (!is.na(q_col)) to_numeric(res[[q_col]]) else p.adjust(p_values, method = "holm")
      
      diff_values <- if (!is.na(diff_col) && diff_col %in% colnames(res)) res[[diff_col]] else res[["diff_abn"]]
      if (is.null(diff_values)) {
        diff_values <- rep(FALSE, nrow(res))
      }
      diff_values <- as.logical(diff_values)
      diff_values[is.na(diff_values)] <- FALSE
      
      taxonomy_cols <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(res))
      taxonomy_df <- if (length(taxonomy_cols) > 0) res[, taxonomy_cols, drop = FALSE] else NULL
      
      base_df <- data.frame(
        taxon = if ("taxon" %in% colnames(res)) res$taxon else res$taxa_label,
        taxa_label = res$taxa_label,
        feature_id = res$feature_id,
        lfc = lfc_values,
        se = se_values,
        W = W_values,
        p_val = p_values,
        q_val = q_values,
        diff = diff_values,
        stringsAsFactors = FALSE
      )
      
      if (!is.null(taxonomy_df)) {
        base_df <- cbind(taxonomy_df, base_df)
      }
      
      base_df$direction <- ifelse(base_df$lfc > 0,
                                  paste("Increase in", comparison_group_label),
                                  paste("Decrease in", comparison_group_label))
      base_df
    })

    output$ancom_table <- renderDT({
      req(ancom_processed())
      res <- ancom_processed()
      leading_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "taxon", "taxa_label", "feature_id")
      existing_leading <- intersect(leading_cols, colnames(res))
      remaining_cols <- setdiff(colnames(res), existing_leading)
      res <- res[, c(existing_leading, remaining_cols), drop = FALSE]
      datatable(res, options = list(scrollX = TRUE))
    })
    
    output$download_ancom_table <- downloadHandler(
      filename = function() {
        control_tag <- gsub("[^A-Za-z0-9_]+", "_", selected_group_labels()$control)
        comparison_tag <- gsub("[^A-Za-z0-9_]+", "_", selected_group_labels()$comparison)
        paste0("ancombc2_table_", control_tag, "_vs_", comparison_tag, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        req(ancom_processed())
        readr::write_tsv(ancom_processed(), file)
      }
    )
    
    volcano_plot_reactive <- reactive({
      req(ancom_processed(), input$volcano_y_axis)
      res <- ancom_processed()
      
      y_col_name <- input$volcano_y_axis 
      y_axis_label <- if (y_col_name == "q_val") "-log10(FDR-adjusted p-value)" else "-log10(Raw p-value)"
      y_data <- res[[y_col_name]]
      
      # plot_df는 lfc와 선택된 p-value/q-value가 유효한 행만 사용
      plot_df <- res[!is.na(res$lfc) & !is.na(y_data) & y_data > 0, , drop = FALSE]
      
      validate(
        need(nrow(plot_df) > 0, "No valid rows to plot. Check ANCOM results.")
      )
      
      plot_df$y_metric <- -log10(y_data[!is.na(res$lfc) & !is.na(y_data) & y_data > 0])
      
      top_labels <- plot_df[order(plot_df[[y_col_name]]), , drop = FALSE]
      
      if (nrow(top_labels) > 10) {
        top_labels <- top_labels[seq_len(10), , drop = FALSE]
      }

      comparison_group_name <- selected_group_labels()$comparison
      control_group_name <- selected_group_labels()$control
      
      p <- ggplot(plot_df, aes(x = lfc, y = y_metric, color = diff)) + 
        geom_point(alpha = 0.8, size = 2.5) +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50")) +
        geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey") +
        labs(title = paste("ANCOM-BC2 Volcano Plot:", control_group_name, "vs", comparison_group_name),
             x = paste("Log Fold Change (", control_group_name, "vs", comparison_group_name, ")", sep = ""),
             y = y_axis_label, 
             color = "Significant",
             subtitle = paste("Positive LFC = increased in", comparison_group_name)) +
        theme_minimal()
      
      if (nrow(top_labels) > 0) {
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_text_repel(
            data = top_labels,
            aes(label = taxa_label),
            size = 3,
            show.legend = FALSE
          )
        } else {
          p <- p + geom_text(
            data = top_labels,
            aes(label = taxa_label),
            size = 3,
            vjust = -0.8,
            show.legend = FALSE,
            check_overlap = TRUE
          )
        }
      }
      
      p
    })
    
    bar_plot_reactive <- reactive({
      req(ancom_processed())
      res <- ancom_processed()
      plot_df <- res[!res$diff & !is.na(res$lfc), , drop = FALSE]
      if (nrow(plot_df) == 0) {
        plot_df <- res[!is.na(res$lfc), , drop = FALSE]
      }
      validate(
        need(nrow(plot_df) > 0, "No valid rows to plot for bar chart.")
      )
      
      plot_df <- plot_df[order(abs(plot_df$lfc), decreasing = TRUE), , drop = FALSE]
      top_n <- min(10, nrow(plot_df))
      res_top <- plot_df[seq_len(top_n), , drop = FALSE]
      res_top$direction_flag <- res_top$lfc > 0
      
      comparison_group_name <- selected_group_labels()$comparison
      control_group_name <- selected_group_labels()$control
      
      ggplot(res_top, aes(x = reorder(taxa_label, lfc), y = lfc, fill = direction_flag)) +
        geom_col(color = "black") +
        coord_flip() +
        scale_fill_manual(
          values = c("TRUE" = "red", "FALSE" = "blue"),
          labels = c(
            "TRUE" = paste("Increase in", comparison_group_name),
            "FALSE" = paste("Decrease in", comparison_group_name)
          )
        ) +
        labs(title = paste("Top", top_n, "Differential Taxa by LFC:", control_group_name, "vs", comparison_group_name),
             x = input$tax_level,
             y = paste("Log Fold Change (", control_group_name, "vs", comparison_group_name, ")", sep = ""),
             fill = "Direction",
             subtitle = paste("Positive LFC = increased in", comparison_group_name)) +
        theme_bw()
        
    })
    
    output$ancom_plot <- renderPlot(
      { volcano_plot_reactive() },
      height = function() { req(input$plot_height); input$plot_height },
      width = function() { req(input$plot_width); input$plot_width }
    )
    output$ancom_barplot <- renderPlot(
      { bar_plot_reactive() },
      height = function() { req(input$plot_height); input$plot_height },
      width = function() { req(input$plot_width); input$plot_width }
    )
    
    output$download_volcano <- downloadHandler(
      filename = function() {
        control_tag <- gsub("[^A-Za-z0-9_]+", "_", selected_group_labels()$control)
        comparison_tag <- gsub("[^A-Za-z0-9_]+", "_", selected_group_labels()$comparison)
        paste0("ancombc2_volcano_", control_tag, "_vs_", comparison_tag, ".png")
      },
      content = function(file) {
        req(input$plot_height, input$plot_width)
        ggplot2::ggsave(file, plot = volcano_plot_reactive(), device = "png",
                        width = input$plot_width / 100, height = input$plot_height / 100)
      }
    )
    
    output$download_barplot <- downloadHandler(
      filename = function() {
        control_tag <- gsub("[^A-Za-z0-9_]+", "_", selected_group_labels()$control)
        comparison_tag <- gsub("[^A-Za-z0-9_]+", "_", selected_group_labels()$comparison)
        paste0("ancombc2_barplot_", control_tag, "_vs_", comparison_tag, ".png")
      },
      content = function(file) {
        req(input$plot_height, input$plot_width)
        ggplot2::ggsave(file, plot = bar_plot_reactive(), device = "png",
                        width = input$plot_width / 100, height = input$plot_height / 100)
      }
    )
  })
}
