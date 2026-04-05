library(shiny)
library(ggplot2)
library(dplyr)
library(phyloseq)

## UI
mod_taxa_comparison_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4(icon("square-poll-vertical"), "Taxa Comparison"),
        hr(),
        uiOutput(ns("group_selector")),
        selectInput(
          ns("tax_level"),
          "Taxonomic Level:",
          choices = c("ASV", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
          selected = "Genus"
        ),
        selectizeInput(
          ns("taxa_selected"),
          "Taxa to Compare:",
          choices = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Select one or more taxa",
            plugins = list("remove_button")
          )
        ),
        selectInput(
          ns("abundance_mode"),
          "Abundance Scale:",
          choices = c(
            "CLR Abundance" = "clr",
            "Relative Abundance (%)" = "relative",
            "Log TSS (log10(% + 1))" = "log_tss"
          ),
          selected = "clr"
        ),
        hr(),
        numericInput(ns("plot_width"), "Plot Width (px):", value = 400, min = 300, step = 50),
        numericInput(ns("plot_height"), "Plot Height (px):", value = 500, min = 300, step = 50),
        hr(),
        h5(icon("download"), "Download Plot"),
        downloadButton(ns("download_taxa_plot"), "Download Plot (PNG)", style = "font-size: 12px;"),
        downloadButton(ns("download_taxa_data"), "Download Data (TSV)", style = "font-size: 12px;"),
        hr(),
        h5(icon("table"), "Download Full Matrix"),
        downloadButton(ns("download_raw_matrix"), "Raw Counts Matrix (TSV)", style = "font-size: 12px;"),
        downloadButton(ns("download_rel_matrix"), "Relative Abundance Matrix (TSV)", style = "font-size: 12px;")
      ),
      mainPanel(
        width = 9,
        plotOutput(ns("taxa_comparison_plot"), height = "auto")
      )
    )
  )
}

## Server
mod_taxa_comparison_server <- function(id, ps_obj, meta_cols) {
  moduleServer(id, function(input, output, session) {
    
    output$group_selector <- renderUI({
      req(meta_cols())
      group_choices <- setdiff(meta_cols(), "SampleID")
      selected_col <- if (length(group_choices) > 0) group_choices[1] else NULL
      selectInput(
        session$ns("group_var"),
        "Grouping Variable:",
        choices = group_choices,
        selected = selected_col
      )
    })
    
    taxa_long_data <- reactive({
      req(ps_obj(), input$group_var, input$tax_level)
      
      ps <- ps_obj()
      tax_level <- input$tax_level
      
      if (tax_level != "ASV") {
        validate(need(!is.null(phyloseq::tax_table(ps)), "Taxonomy table is missing."))
        tax_ranks <- phyloseq::rank_names(ps)
        validate(need(tax_level %in% tax_ranks, paste0("Level '", tax_level, "' not found.")))
        ps <- phyloseq::tax_glom(ps, taxrank = tax_level, NArm = FALSE)
      }
      
      abundance_mode <- input$abundance_mode
      if (is.null(abundance_mode) || !abundance_mode %in% c("clr", "relative", "log_tss")) {
        abundance_mode <- "clr"
      }

      if (abundance_mode %in% c("relative", "log_tss")) {
        ps <- phyloseq::transform_sample_counts(ps, function(x) {
          s <- sum(x)
          if (s == 0) x else x / s
        })
      } else {
        otu_mat <- as(phyloseq::otu_table(ps), "matrix")
        if (!phyloseq::taxa_are_rows(ps)) {
          otu_mat <- t(otu_mat)
        }
        
        otu_clr <- apply(otu_mat, 2, function(x) {
          x <- as.numeric(x) + 1
          log_x <- log(x)
          log_x - mean(log_x, na.rm = TRUE)
        })
        
        otu_clr <- as.matrix(otu_clr)
        rownames(otu_clr) <- rownames(otu_mat)
        colnames(otu_clr) <- colnames(otu_mat)
        
        otu_obj <- phyloseq::otu_table(otu_clr, taxa_are_rows = TRUE)
        sdata_obj <- phyloseq::sample_data(ps)
        tax_obj <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
        ps <- if (is.null(tax_obj)) {
          phyloseq::phyloseq(otu_obj, sdata_obj)
        } else {
          phyloseq::phyloseq(otu_obj, tax_obj, sdata_obj)
        }
      }
      
      df <- phyloseq::psmelt(ps)
      
      meta_df <- as.data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE)
      for (meta_col in colnames(meta_df)) {
        syntactic_col <- make.names(meta_col)
        if (!meta_col %in% names(df) && syntactic_col %in% names(df)) {
          names(df)[names(df) == syntactic_col] <- meta_col
        }
      }
      
      if (tax_level == "ASV") {
        df$Taxa <- as.character(df$OTU)
      } else {
        taxa_values <- as.character(df[[tax_level]])
        taxa_values[is.na(taxa_values) | taxa_values == ""] <- as.character(df$OTU[is.na(taxa_values) | taxa_values == ""])
        df$Taxa <- taxa_values
      }
      
      df$Group <- as.factor(df[[input$group_var]])
      df <- df[!is.na(df$Group), , drop = FALSE]
      
      if (abundance_mode == "relative") {
        df$AbundancePlot <- df$Abundance * 100
      } else if (abundance_mode == "log_tss") {
        df$AbundancePlot <- log10((df$Abundance * 100) + 1)
      } else {
        df$AbundancePlot <- df$Abundance
      }
      df <- df[is.finite(df$AbundancePlot), , drop = FALSE]
      
      df
    })
    
    tax_level_ps_raw <- reactive({
      req(ps_obj(), input$tax_level)
      ps <- ps_obj()
      tax_level <- input$tax_level
      
      if (tax_level != "ASV") {
        validate(
          need(!is.null(phyloseq::tax_table(ps)), "Taxonomy table is missing.")
        )
        tax_ranks <- phyloseq::rank_names(ps)
        validate(
          need(tax_level %in% tax_ranks, paste0("Taxonomic level '", tax_level, "' is not available in this dataset."))
        )
        ps <- phyloseq::tax_glom(ps, taxrank = tax_level, NArm = FALSE)
      }
      ps
    })
    
    tax_matrices <- reactive({
      req(tax_level_ps_raw(), input$tax_level)
      ps <- tax_level_ps_raw()
      tax_level <- input$tax_level
      
      otu_mat <- as(phyloseq::otu_table(ps), "matrix")
      if (!phyloseq::taxa_are_rows(ps)) {
        otu_mat <- t(otu_mat)
      }
      
      if (tax_level == "ASV") {
        taxa_labels <- rownames(otu_mat)
      } else {
        tax_tab <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE)
        taxa_labels <- as.character(tax_tab[[tax_level]])
        missing_idx <- is.na(taxa_labels) | taxa_labels == ""
        taxa_labels[missing_idx] <- rownames(otu_mat)[missing_idx]
      }
      taxa_labels <- make.unique(as.character(taxa_labels))
      rownames(otu_mat) <- taxa_labels
      
      raw_df <- as.data.frame(t(otu_mat), stringsAsFactors = FALSE, check.names = FALSE)
      raw_df <- tibble::rownames_to_column(raw_df, var = "SampleID")
      
      rel_mat <- apply(otu_mat, 2, function(x) {
        s <- sum(x, na.rm = TRUE)
        if (s == 0) rep(0, length(x)) else (x / s) * 100
      })
      rel_mat <- as.matrix(rel_mat)
      rownames(rel_mat) <- rownames(otu_mat)
      colnames(rel_mat) <- colnames(otu_mat)
      
      rel_df <- as.data.frame(t(rel_mat), stringsAsFactors = FALSE, check.names = FALSE)
      rel_df <- tibble::rownames_to_column(rel_df, var = "SampleID")
      
      list(raw = raw_df, rel = rel_df)
    })
    
    observeEvent(list(taxa_long_data(), input$tax_level), {
      df <- taxa_long_data()
      taxa_summary <- df %>%
        dplyr::group_by(Taxa) %>%
        dplyr::summarise(mean_abundance = mean(AbundancePlot, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(mean_abundance))
      
      taxa_choices <- as.character(taxa_summary$Taxa)
      taxa_choices <- taxa_choices[!is.na(taxa_choices) & nzchar(taxa_choices)]
      
      current_selected <- input$taxa_selected
      default_selected <- if (length(current_selected) > 0) {
        intersect(current_selected, taxa_choices)
      } else {
        head(taxa_choices, 1)
      }
      if (length(default_selected) == 0 && length(taxa_choices) > 0) {
        default_selected <- head(taxa_choices, 1)
      }
      
      updateSelectizeInput(
        session,
        "taxa_selected",
        choices = taxa_choices,
        selected = default_selected,
        server = TRUE
      )
    })
    
    taxa_plot_data <- reactive({
      req(taxa_long_data(), input$group_var)
      df <- taxa_long_data()
      selected_taxa <- if (is.null(input$taxa_selected) || length(input$taxa_selected) == 0) {
        character(0)
      } else {
        input$taxa_selected
      }
      
      validate(need(length(selected_taxa) > 0, "Select at least one taxa in 'Taxa to Compare'."))
      
      df_sub <- df[df$Taxa %in% selected_taxa, , drop = FALSE]
      
      taxa_order <- df_sub %>%
        dplyr::group_by(Taxa) %>%
        dplyr::summarise(med = stats::median(AbundancePlot, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(med)) %>%
        dplyr::pull(Taxa)
      df_sub$Taxa <- factor(df_sub$Taxa, levels = taxa_order)
      
      df_sub
    })
    
    taxa_plot_reactive <- reactive({
      req(taxa_plot_data())
      df <- taxa_plot_data()
      abundance_mode <- input$abundance_mode
      if (is.null(abundance_mode) || !abundance_mode %in% c("clr", "relative", "log_tss")) {
        abundance_mode <- "clr"
      }
      y_label <- switch(
        abundance_mode,
        "relative" = "Relative Abundance (%)",
        "log_tss" = "Log TSS (log10(% + 1))",
        "CLR Abundance"
      )
      integer_breaks <- function(x) {
        if (length(x) == 0 || all(!is.finite(x))) {
          return(NULL)
        }
        rng <- range(x, na.rm = TRUE)
        lo <- 0
        hi <- ceiling(rng[2])
        if (!is.finite(lo) || !is.finite(hi)) {
          return(NULL)
        }
        if (hi < 0) {
          hi <- 0
        }
        if (lo == hi) {
          return(lo)
        }
        step <- max(1, ceiling((hi - lo) / 6))
        seq(lo, hi, by = step)
      }
      
      pval_df <- df %>%
        dplyr::group_by(Taxa) %>%
        dplyr::group_modify(~{
          n_groups <- dplyr::n_distinct(.x$Group)
          p_val <- tryCatch({
            if (n_groups == 2) {
              stats::wilcox.test(AbundancePlot ~ Group, data = .x)$p.value
            } else if (n_groups > 2) {
              stats::kruskal.test(AbundancePlot ~ Group, data = .x)$p.value
            } else {
              NA_real_
            }
          }, error = function(e) NA_real_)
          dplyr::tibble(
            p_label = if (is.na(p_val)) "p-value: NA" else paste0("p-value = ", format(p_val, digits = 3))
          )
        }) %>%
        dplyr::ungroup()
      
      df <- dplyr::left_join(df, pval_df, by = "Taxa")
      df$Taxa_with_p <- paste0(as.character(df$Taxa), "\n", df$p_label)
      
      ggplot2::ggplot(df, ggplot2::aes(x = Group, y = AbundancePlot, fill = Group)) +
        ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.65) +
        ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
        ggplot2::scale_y_continuous(breaks = integer_breaks, limits = c(0, NA)) +
        ggplot2::facet_wrap(~ Taxa_with_p, scales = "free_y") +
        ggplot2::theme_bw() +
        ggplot2::labs(title = paste0("Taxa Comparison: ", input$tax_level),
                      x = input$group_var, y = y_label) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
    })
    
    output$taxa_comparison_plot <- renderPlot({
      taxa_plot_reactive()
    }, height = function() { req(input$plot_height); input$plot_height },
       width = function() { req(input$plot_width); input$plot_width })
    
    output$download_taxa_plot <- downloadHandler(
      filename = function() { paste0("taxa_plot_", Sys.Date(), ".png") },
      content = function(file) {
        ggplot2::ggsave(file, plot = taxa_plot_reactive(), device = "png",
                        width = input$plot_width/100, height = input$plot_height/100)
      }
    )
    
    output$download_taxa_data <- downloadHandler(
      filename = function() { paste0("taxa_data_", Sys.Date(), ".tsv") },
      content = function(file) {
        out_df <- taxa_plot_data()[, c("Sample", "Group", "Taxa", "AbundancePlot")]
        colnames(out_df) <- c(
          "SampleID",
          input$group_var,
          "Taxa",
          if (input$abundance_mode == "relative") {
            "Relative_Abundance_percent"
          } else if (input$abundance_mode == "log_tss") {
            "Log_TSS_log10_percent_plus1"
          } else {
            "CLR_Abundance"
          }
        )
        readr::write_tsv(out_df, file)
      }
    )
    
    output$download_raw_matrix <- downloadHandler(
      filename = function() {
        paste0("taxa_raw_matrix_", input$tax_level, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        req(tax_matrices())
        readr::write_tsv(tax_matrices()$raw, file)
      }
    )
    
    output$download_rel_matrix <- downloadHandler(
      filename = function() {
        paste0("taxa_relative_abundance_matrix_", input$tax_level, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        req(tax_matrices())
        readr::write_tsv(tax_matrices()$rel, file)
      }
    )
  })
}
