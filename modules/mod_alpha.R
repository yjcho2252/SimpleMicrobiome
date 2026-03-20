mod_alpha_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        # ⭐ 1. 주 그룹 변수 선택 UI (Primary) ⭐
        uiOutput(ns("local_group_selector")),
        
        # ⭐ 2. 보조 그룹 변수 선택 UI (Secondary) ⭐
        uiOutput(ns("secondary_group_selector")), 
        hr(),
        # -------------------------------
        
        h4("Index Selection"),
        checkboxGroupInput(ns("alpha_methods"), "Indices:",
                           choices = c("Observed", "Chao1", "Shannon", "Simpson"),
                           selected = c("Chao1", "Shannon")),
        
        # ⭐ P-값 조정 방법 선택 UI ⭐
        hr(),
        h4("Statistical Comparison"),
        selectInput(ns("p_adj_method"), "P-value Correction Method:",
                    choices = c(
                      "FDR (Default)" = "fdr",
                      "Bonferroni" = "bonferroni",
                      "Holm" = "holm",
                      "None (No Correction)" = "none"
                    ),
                    selected = "fdr"),
        
        # ⭐ Plot 크기 조절 Numeric Input ⭐
        hr(),
        h4("Plot Dimensions (px)"),
        numericInput(ns("plot_width"), "Plot Width:", 
                     value = 800, min = 300, step = 50),
        numericInput(ns("plot_height"), "Plot Height:", 
                     value = 700, min = 300, step = 50),
        
        # 다운로드 버튼 추가 (ID 변경 및 데이터 다운로드 추가)
        hr(),
        h4("Download Results"),
        downloadButton(ns("download_alpha_plot"), "Download Plot (PNG)"), # ID 변경
        downloadButton(ns("download_alpha_data"), "Download Data (TSV)") # 데이터 다운로드 추가
      ),
      mainPanel(
        # 플롯 출력
        plotOutput(ns("alpha_plot_out"), height = "auto"), 
        br(),
        # 희귀화 깊이 출력
        textOutput(ns("rarefy_size_text"))
      )
    )
  )
}
mod_alpha_server <- function(id, ps_obj, meta_cols) { 
  moduleServer(id, function(input, output, session) {
    
    # 필수 패키지 확인
    if (!requireNamespace("ggpubr", quietly = TRUE)) {
      showNotification("Package 'ggpubr' is required for statistical testing on plot.", type = "warning")
    }
    
    # ----------------------------------------------------
    # 모듈 내 주/보조 그룹 변수 선택 로직
    # ----------------------------------------------------
    
    output$local_group_selector <- renderUI({
      req(meta_cols())
      group_choices <- setdiff(meta_cols(), "SampleID")
      selected_col <- if (length(group_choices) > 0) group_choices[1] else meta_cols()[1]
      
      selectInput(session$ns("group_var"), "Select Primary Grouping Variable:", 
                  choices = group_choices, selected = selected_col) 
    })
    
    group_var <- reactive({
      req(input$group_var)
      input$group_var
    })
    
    output$secondary_group_selector <- renderUI({
      req(meta_cols())
      group_choices <- setdiff(meta_cols(), c("SampleID", input$group_var))
      select_options <- c("(None)" = "none", group_choices)
      
      selectInput(session$ns("secondary_group_var"), "Select Secondary Grouping Variable:",
                  choices = select_options,
                  selected = "none")
    })
    
    secondary_group_var <- reactive({
      input$secondary_group_var
    })
    
    # ----------------------------------------------------
    # 데이터 준비: 희귀화 (Rarefaction)
    # ----------------------------------------------------
    rarefaction_size <- reactive({
      req(ps_obj())
      # 최소 시퀀싱 깊이를 희귀화 깊이로 설정
      min(phyloseq::sample_sums(ps_obj()))
    })
    
    alpha_data_reactive <- reactive({
      req(ps_obj())
      min_size <- rarefaction_size()
      
      tryCatch({
        # Rarefaction 수행
        physeq_rarefied <- phyloseq::rarefy_even_depth(
          ps_obj(),
          sample.size = min_size,
          rngseed = 42,
          replace = FALSE,
          verbose = FALSE
        )
        
        # Rarefaction 후 0인 샘플 제거
        na_samples <- which(is.na(phyloseq::sample_sums(physeq_rarefied)) | is.infinite(phyloseq::sample_sums(physeq_rarefied)) | phyloseq::sample_sums(physeq_rarefied) == 0)
        if (length(na_samples) > 0) {
          physeq_rarefied <- phyloseq::prune_samples(phyloseq::sample_names(physeq_rarefied)[-na_samples], physeq_rarefied)
        }
        
        return(physeq_rarefied)
      }, error = function(e) {
        showNotification(paste("Rarefaction failed:", e$message), type = "error")
        return(NULL)
      })
    })
    
    # ----------------------------------------------------
    # ⭐ 데이터 준비: 알파 다양성 값 계산 및 Tidy Data (Long Format) 생성 ⭐
    # (Plot 및 Download에 사용)
    # ----------------------------------------------------
    alpha_long_reactive <- reactive({
      req(physeq_rarefied <- alpha_data_reactive(), group_var(), input$alpha_methods)
      
      primary_col <- group_var()
      secondary_col <- secondary_group_var() 
      is_secondary <- secondary_col != "none"
      
      meta_data <- as(phyloseq::sample_data(physeq_rarefied), "data.frame")
      
      # 1. 유효성 검사 및 필터링
      valid_cols <- c(primary_col)
      if (is_secondary) valid_cols <- c(valid_cols, secondary_col)
      
      # 필수 변수에 NA가 없는 샘플만 사용
      valid_samples <- apply(meta_data[, valid_cols, drop=FALSE], 1, function(x) !any(is.na(x)))
      physeq_filtered <- phyloseq::prune_samples(valid_samples, physeq_rarefied)
      meta_filtered <- meta_data[valid_samples, , drop = FALSE]
      
      validate(need(nrow(meta_filtered) > 0, "No samples available after filtering by the selected group(s)."))
      
      # 그룹 변수를 무조건 factor로 강제 변환
      meta_filtered[[primary_col]] <- factor(meta_filtered[[primary_col]])
      if (is_secondary) {
        meta_filtered[[secondary_col]] <- factor(meta_filtered[[secondary_col]])
      }
      
      measures_req <- unique(input$alpha_methods)
      alpha_df <- phyloseq::estimate_richness(physeq_filtered, measures = measures_req)
      
      # 2. 그룹 변수 추가 및 타입 안정화
      alpha_df[[primary_col]] <- meta_filtered[[primary_col]]
      if (is_secondary) {
        alpha_df[[secondary_col]] <- meta_filtered[[secondary_col]]
      }
      alpha_df$SampleID <- rownames(alpha_df) # SampleID 열 이름 명확화
      
      measures_to_plot <- intersect(measures_req, colnames(alpha_df))
      validate(need(length(measures_to_plot) > 0, "Selected alpha diversity measures are unavailable."))
      
      # 3. ggplot 전용 tidy 데이터 생성 (Long Format)
      alpha_long_list <- lapply(measures_to_plot, function(measure) {
        df <- data.frame(
          Alpha_Index = rep(measure, nrow(alpha_df)), # Index 이름을 명확히 변경
          Value = alpha_df[[measure]], # Value 이름을 명확히 변경
          SampleID = alpha_df$SampleID,
          stringsAsFactors = FALSE
        )
        df[[primary_col]] <- alpha_df[[primary_col]]
        if (is_secondary) {
          df[[secondary_col]] <- alpha_df[[secondary_col]]
        }
        return(df)
      })
      alpha_long <- do.call(rbind, alpha_long_list)
      
      alpha_long <- alpha_long[!is.na(alpha_long$Value), , drop = FALSE]
      validate(need(nrow(alpha_long) > 0, "No valid alpha diversity values available."))
      
      return(alpha_long) # 최종 롱 포맷 데이터 반환
    })
    
    # ----------------------------------------------------
    # 그래프를 그리는 reactive expression
    # ----------------------------------------------------
    alpha_plot_reactive <- eventReactive(list(alpha_long_reactive(), group_var(), secondary_group_var(), input$p_adj_method), {
      req(alpha_long <- alpha_long_reactive(), group_var(), secondary_group_var(), input$p_adj_method)
      
      primary_col <- group_var()
      secondary_col <- secondary_group_var() 
      is_secondary <- secondary_col != "none"
      
      # Plot 설정 변수 계산
      if (is_secondary) {
        x_axis_col <- secondary_col
        fill_col <- secondary_col
        facet_col <- primary_col
        stat_compare_col <- secondary_col
      } else {
        x_axis_col <- primary_col
        fill_col <- primary_col
        facet_col <- "NULL" 
        stat_compare_col <- primary_col
      }
      
      # 통계 비교를 위한 그룹 레벨 확인 (비교 대상 그룹)
      group_levels <- sort(unique(as.character(alpha_long[[stat_compare_col]])))
      num_groups <- length(group_levels)
      
      validate(need(num_groups >= 2, paste0("Selected comparison group '", stat_compare_col, "' has only one unique value. Cannot perform comparison.")))
      
      stat_method <- if (num_groups == 2) "t.test" else "kruskal.test"
      
      # 5. ggplot 생성
      p <- ggplot2::ggplot(alpha_long, ggplot2::aes_string(x = x_axis_col, y = "Value", fill = fill_col)) + # Value 사용
        ggplot2::geom_boxplot(color = "black", width = 0.6, outlier.shape = NA) +
        ggplot2::geom_jitter(color = "black", width = 0.2, alpha = 0.6, size = 1.5, show.legend = FALSE) +
        ggplot2::theme_bw() +
        ggplot2::labs(
          title = paste("Alpha Diversity Indices by", primary_col, if(is_secondary) paste0("(", secondary_col, " comparisons)") else "", "(Rarefied)"),
          y = "Diversity Index",
          x = x_axis_col
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold"),
          axis.title = ggplot2::element_text(size = 12),
          strip.text = ggplot2::element_text(size = 11, face = "bold"),
          legend.position = "bottom"
        ) +
        ggplot2::guides(color = "none") +
        ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.1))
      
      # Facet (Alpha_Index ~ Primary_Group)
      if (is_secondary) {
        p <- p + ggplot2::facet_grid(stats::as.formula(paste("Alpha_Index", "~", facet_col)), scales = "free_y")
      } else {
        p <- p + ggplot2::facet_wrap(~ Alpha_Index, scales = "free_y", ncol = 2)
      }
      
      # P-값 표시
      if (requireNamespace("ggpubr", quietly = TRUE) && num_groups >= 2) {
        
        comparison_levels <- unique(as.character(alpha_long[[stat_compare_col]]))
        comparison_levels <- comparison_levels[!is.na(comparison_levels) & nzchar(comparison_levels)]
        
        # stat_signif 내부 factor 연산 오류 방지를 위해 문자열 레벨로 고정
        alpha_long[[stat_compare_col]] <- factor(as.character(alpha_long[[stat_compare_col]]), levels = comparison_levels)
        
        final_comparisons <- utils::combn(comparison_levels, 2, simplify = FALSE)
        adj_method <- input$p_adj_method
        group_by_col <- if (is_secondary) primary_col else NULL 
        if (length(final_comparisons) > 0) {
          p <- p + ggpubr::stat_compare_means(
            comparisons = final_comparisons,
            method = stat_method,
            label = "p.format",
            label.x.tip = TRUE,
            tip.length = 0.01,
            size = 4,
            p.adjust.method = adj_method,
            group.by = group_by_col,
            label.y.npc = 0.95,
            vjust = 0
          )
        }
      }
      
      return(p)
    })
    
    # ----------------------------------------------------
    # 출력 및 다운로드 핸들러
    # ----------------------------------------------------
    
    # 1. 플롯 출력
    output$alpha_plot_out <- renderPlot({ alpha_plot_reactive() },
                                        height = function() { input$plot_height },
                                        width = function() { input$plot_width }
    )
    
    # 2. 희귀화 깊이 텍스트 출력
    output$rarefy_size_text <- renderText({
      req(rarefaction_size())
      paste("Rarefaction depth (minimum sequencing depth):", scales::comma(rarefaction_size()))
    })
    
    # 3. 플롯 다운로드 핸들러 (ID 변경됨)
    output$download_alpha_plot <- downloadHandler(
      filename = function() { paste0("alpha_diversity_plot_", Sys.Date(), ".png") },
      content = function(file) {
        plot_width_in <- input$plot_width / 72
        plot_height_in <- input$plot_height / 72
        
        ggplot2::ggsave(file, 
                        plot = alpha_plot_reactive(), 
                        device = "png", 
                        width = plot_width_in, 
                        height = plot_height_in,
                        dpi = 300 
        )
      }
    )
    
    # ⭐ 4. 데이터 테이블 다운로드 핸들러 (CSV) ⭐
    # ⭐ 4. 데이터 테이블 다운로드 핸들러 (TSV로 수정) ⭐
    output$download_alpha_data <- downloadHandler(
      filename = function() {
        # 파일 확장자를 .tsv로 변경
        paste0("alpha_diversity_data_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        data_to_write <- alpha_long_reactive() 
        
        # write.table을 사용하여 탭(\t)으로 구분된 파일(TSV) 생성
        utils::write.table(data_to_write, 
                           file, 
                           sep = "\t", 
                           row.names = FALSE, 
                           quote = FALSE)
      }
    )
  })
}
