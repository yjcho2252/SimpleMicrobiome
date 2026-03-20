library(shiny)
library(ggplot2)
library(dplyr)
library(phyloseq)

## UI Module ------------------------------------------------------------------
mod_taxa_comparison_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        width = 3,
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
          ns("plot_mode"),
          "Plot Mode:",
          choices = c(
            "Multiple Taxa (Facet)" = "multi",
            "Single Taxa Detail (p-value)" = "single"
          ),
          selected = "single"
        ),
        uiOutput(ns("single_taxa_selector")),
        checkboxInput(ns("use_relative"), "Use Relative Abundance (%)", value = FALSE),
        hr(),
        numericInput(ns("plot_width"), "Plot Width (px):", value = 1100, min = 400, step = 50),
        numericInput(ns("plot_height"), "Plot Height (px):", value = 700, min = 300, step = 50),
        hr(),
        downloadButton(ns("download_taxa_plot"), "Download Plot (PNG)"),
        downloadButton(ns("download_taxa_data"), "Download Data (TSV)"),
        hr(),
        h5("Download Full Matrix"),
        downloadButton(ns("download_raw_matrix"), "Raw Counts Matrix (TSV)"),
        downloadButton(ns("download_rel_matrix"), "Relative Abundance Matrix (TSV)")
      ),
      mainPanel(
        width = 9,
        plotOutput(ns("taxa_comparison_plot"), height = "auto")
      )
    )
  )
}

## Server Module --------------------------------------------------------------
mod_taxa_comparison_server <- function(id, ps_obj, meta_cols) {
  moduleServer(id, function(input, output, session) {
    
    # 1. Grouping Variable Selector UI 생성
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
    
    # 2. Phyloseq 객체를 긴 형태(Long format)의 데이터프레임으로 변환
    taxa_long_data <- reactive({
      req(ps_obj(), input$group_var, input$tax_level)
      
      ps <- ps_obj()
      tax_level <- input$tax_level
      
      # Taxonomy Glom (ASV가 아닌 경우)
      if (tax_level != "ASV") {
        validate(need(!is.null(phyloseq::tax_table(ps)), "Taxonomy table is missing."))
        tax_ranks <- phyloseq::rank_names(ps)
        validate(need(tax_level %in% tax_ranks, paste0("Level '", tax_level, "' not found.")))
        ps <- phyloseq::tax_glom(ps, taxrank = tax_level, NArm = FALSE)
      }
      
      # Relative Abundance 변환 (필요시)
      if (isTRUE(input$use_relative)) {
        ps <- phyloseq::transform_sample_counts(ps, function(x) {
          s <- sum(x)
          if (s == 0) x else x / s
        })
      } else {
        # raw count 대신 CLR 변환값 사용 (pseudo-count +1)
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
      
      # 메타데이터 컬럼명 복원 (make.names 방지)
      meta_df <- as.data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE)
      for (meta_col in colnames(meta_df)) {
        syntactic_col <- make.names(meta_col)
        if (!meta_col %in% names(df) && syntactic_col %in% names(df)) {
          names(df)[names(df) == syntactic_col] <- meta_col
        }
      }
      
      # Taxa 이름 정리
      if (tax_level == "ASV") {
        df$Taxa <- as.character(df$OTU)
      } else {
        taxa_values <- as.character(df[[tax_level]])
        taxa_values[is.na(taxa_values) | taxa_values == ""] <- as.character(df$OTU[is.na(taxa_values) | taxa_values == ""])
        df$Taxa <- taxa_values
      }
      
      # 그룹 설정 및 NA 필터링
      df$Group <- as.factor(df[[input$group_var]])
      df <- df[!is.na(df$Group), , drop = FALSE]
      
      # Plot용 Abundance 값 계산
      df$AbundancePlot <- if (isTRUE(input$use_relative)) df$Abundance * 100 else df$Abundance
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
    
    # 3. Taxa 선택 UI 업데이트 (상위 8개 자동 선택)
    observeEvent(list(taxa_long_data(), input$tax_level), {
      df <- taxa_long_data()
      taxa_summary <- df %>%
        dplyr::group_by(Taxa) %>%
        dplyr::summarise(mean_abundance = mean(AbundancePlot, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(mean_abundance))
      
      taxa_choices <- as.character(taxa_summary$Taxa)
      taxa_choices <- taxa_choices[!is.na(taxa_choices) & nzchar(taxa_choices)]
      
      # 초기값이 없으면 상위 8개를 기본값으로 세팅
      current_selected <- input$taxa_selected
      default_selected <- if (length(current_selected) > 0) {
        intersect(current_selected, taxa_choices)
      } else {
        head(taxa_choices, min(8, length(taxa_choices)))
      }
      
      updateSelectizeInput(
        session,
        "taxa_selected",
        choices = taxa_choices,
        selected = default_selected,
        server = TRUE
      )
    })
    
    # 4. Single Taxa Selector UI 생성
    output$single_taxa_selector <- renderUI({
      req(input$plot_mode == "single", taxa_long_data())
      df <- taxa_long_data()
      taxa_choices <- df %>%
        dplyr::group_by(Taxa) %>%
        dplyr::summarise(mean_abundance = mean(AbundancePlot, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(mean_abundance)) %>%
        dplyr::pull(Taxa) %>%
        as.character()
      
      selectInput(
        session$ns("single_taxa_selected"),
        "Single Taxa:",
        choices = taxa_choices,
        selected = head(taxa_choices, 1)
      )
    })
    
    # 5. 최종 플롯용 데이터 필터링 (초기 로딩 시 Fallback 포함)
    taxa_plot_data <- reactive({
      req(taxa_long_data(), input$group_var, input$plot_mode)
      df <- taxa_long_data()
      
      # UI가 로딩되기 전이라도 데이터를 표시하기 위한 Fallback 로직
      if (input$plot_mode == "single") {
        selected_taxa <- if (is.null(input$single_taxa_selected) || input$single_taxa_selected == "") {
          head(unique(df$Taxa), 1)
        } else {
          input$single_taxa_selected
        }
      } else {
        selected_taxa <- if (is.null(input$taxa_selected) || length(input$taxa_selected) == 0) {
          head(unique(df$Taxa), 8)
        } else {
          input$taxa_selected
        }
      }
      
      validate(need(length(selected_taxa) > 0, "Wait for taxa selection..."))
      
      df_sub <- df[df$Taxa %in% selected_taxa, , drop = FALSE]
      
      # 중앙값 기준 정렬
      taxa_order <- df_sub %>%
        dplyr::group_by(Taxa) %>%
        dplyr::summarise(med = stats::median(AbundancePlot, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(med)) %>%
        dplyr::pull(Taxa)
      df_sub$Taxa <- factor(df_sub$Taxa, levels = taxa_order)
      
      df_sub
    })
    
    # 6. 플롯 생성 로직
    taxa_plot_reactive <- reactive({
      req(taxa_plot_data())
      df <- taxa_plot_data()
      y_label <- if (isTRUE(input$use_relative)) "Relative Abundance (%)" else "CLR Abundance"
      
      if (input$plot_mode == "single") {
        n_groups <- length(unique(df$Group))
        p_val <- tryCatch({
          if (n_groups == 2) {
            stats::wilcox.test(AbundancePlot ~ Group, data = df)$p.value
          } else if (n_groups > 2) {
            stats::kruskal.test(AbundancePlot ~ Group, data = df)$p.value
          } else { NA_real_ }
        }, error = function(e) NA_real_)
        
        p_label <- if (is.na(p_val)) "p-value: NA" else paste0("p-value = ", format(p_val, digits = 3))
        y_annot <- max(df$AbundancePlot, na.rm = TRUE) * 1.1
        
        ggplot2::ggplot(df, ggplot2::aes(x = Group, y = AbundancePlot, fill = Group)) +
          ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.65) +
          ggplot2::geom_jitter(width = 0.18, alpha = 0.7, size = 1.8) +
          ggplot2::theme_bw() +
          ggplot2::labs(title = paste0("Taxa Detail: ", unique(df$Taxa)[1]),
                        subtitle = p_label, x = input$group_var, y = y_label) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
      } else {
        ggplot2::ggplot(df, ggplot2::aes(x = Group, y = AbundancePlot, fill = Group)) +
          ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.65) +
          ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
          ggplot2::facet_wrap(~ Taxa, scales = "free_y") +
          ggplot2::theme_bw() +
          ggplot2::labs(title = paste0("Taxa Comparison: ", input$tax_level),
                        x = input$group_var, y = y_label) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
      }
    })
    
    # 7. Outputs (Plot & Download)
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
          if (isTRUE(input$use_relative)) "Relative_Abundance_percent" else "CLR_Abundance"
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
