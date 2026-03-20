### UI 함수
mod_barplot_ui <- function(id) {
  ns <- NS(id)
  tax_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  palette_choices <- c(
    "Paired (Default, Max 12)" = "Paired",
    "Custom (HCL 30+ Colors)" = "Custom_HCL",
    "Set1 (Max 9)" = "Set1",
    "Set2 (Max 8)" = "Set2",
    "Set3 (Max 12)" = "Set3",
    "Dark2 (Max 8)" = "Dark2",
    "Accent (Max 8)" = "Accent",
    "Pastel1 (Max 9)" = "Pastel1",
    "Pastel2 (Max 8)" = "Pastel2"
  )
  
  tagList(
    sidebarLayout(
      sidebarPanel(
        width = 3,
        
        #h4("Plot Settings"),
        # Primary Group (기존의 Facet 변수)
        uiOutput(ns("local_group_selector")), 
        
        # ⭐ 추가: Secondary Group 선택 UI ⭐
        uiOutput(ns("secondary_group_selector")),
        hr(),
        
        # Plot Mode 선택 UI: 기본값을 "Sample Level"로 설정
        selectInput(ns("plot_mode"), "Plot Display Mode:",
                    choices = c("Sample Level" = "Sample",
                                "Group Mean" = "Group_Mean"),
                    selected = "Sample"), # 기본값: 샘플별 표시
        hr(),
        
        selectInput(ns("tax_level"), "Taxonomic Level:",
                    choices = tax_ranks,
                    selected = "Genus"),
        
        numericInput(ns("top_n_taxa"), "Number of Top Taxa to Display:",
                     value = 10, min = 1, max = 50, step = 1),
        
        selectInput(ns("color_palette"), "Color Palette Selection:",
                    choices = palette_choices,
                    selected = "Paired"),
        
        selectInput(ns("name_display_mode"), "Taxa Name Display:",
                    choices = c("Full Hierarchy" = "full",
                                "Current Rank Only" = "current"),
                    selected = "full"),
        
        hr(),
        h4("Sample Sorting"),
        
        # ⭐ 조건부: Plot Mode가 Sample Level일 때만 샘플 정렬 UI 표시 ⭐
        conditionalPanel(
          condition = paste0("input['", ns("plot_mode"), "'] == 'Sample'"),
          
          selectInput(ns("sort_method"), "Primary Sort Criterion (within Group):",
                      choices = c(
                        "Top Taxa Abundance (Default)" = "taxa_top1",
                        "Sample Name (Alphabetical)" = "sample_name",
                        "Metadata Variable" = "metadata_var"
                      ),
                      selected = "taxa_top1"),
          
          conditionalPanel(
            condition = paste0("input['", ns("sort_method"), "'] == 'metadata_var'"),
            hr(),
            uiOutput(ns("sort_by_metadata_ui")),
            selectInput(ns("secondary_sort"), "Secondary Sort (after Metadata):",
                        choices = c(
                          "None (Metadata only)" = "none",
                          "Top Taxa Abundance (Descending)" = "taxa_top1",
                          "Sample Name (Alphabetical)" = "sample_name_asc"
                        ),
                        selected = "taxa_top1")
          )
        ),
        
        hr(),
        numericInput(ns("plot_width"), "Plot Width (px):",
                     value = 1000, min = 300, step = 50),
        numericInput(ns("plot_height"), "Plot Height (px):",
                     value = 600, min = 300, step = 50),
        
        downloadButton(ns("download_barplot"), "Download Plot (PNG)")
      ),
      mainPanel(
        width = 9,
        plotOutput(ns("barplot_out"), height = "auto")
      )
    )
  )
}
### Server 함수
mod_barplot_server <- function(id, ps_obj, meta_cols) {
  moduleServer(id, function(input, output, session) {
    
    # ------------------------------------
    # UI Render 및 Reactive 변수 정의
    # ------------------------------------
    output$local_group_selector <- renderUI({
      req(meta_cols())
      group_choices <- setdiff(meta_cols(), "SampleID")
      selected_col <- if (length(group_choices) > 0) group_choices[1] else meta_cols()[1]
      selectInput(session$ns("group_var"), "1. Select Primary Grouping Variable (Facet):",
                  choices = group_choices, selected = selected_col)
    })
    
    # ⭐ 추가: Secondary Group Selector ⭐
    output$secondary_group_selector <- renderUI({
      req(meta_cols(), input$group_var)
      # Primary Group으로 선택된 변수는 제외
      group_choices <- setdiff(meta_cols(), c("SampleID", input$group_var))
      # 기본값은 첫 번째 사용 가능한 변수 또는 NULL
      selected_col <- if (length(group_choices) > 0) group_choices[1] else NULL
      selectInput(session$ns("secondary_var"), "2. Select Secondary Grouping Variable (Sub-Group):",
                  choices = c("None" = "None", group_choices), 
                  selected = "None")
    })
    
    group_var <- reactive({
      req(input$group_var)
      input$group_var
    })
    
    secondary_var <- reactive({
      req(input$secondary_var)
      if (input$secondary_var == "None") {
        return(NULL)
      } else {
        return(input$secondary_var)
      }
    })
    
    output$sort_by_metadata_ui <- renderUI({
      req(meta_cols(), group_var())
      selectInput(session$ns("sort_metadata_var"), "Metadata Variable for Primary Sort:",
                  choices = meta_cols(),
                  selected = group_var())
    })
    
    # Custom HCL 색상 생성 함수 (이전과 동일)
    generate_taxa_colors <- function(n) {
      if (n <= 1) return(character(0))
      hues <- seq(15, 375, length.out = n + 1)[1:n]
      colors <- grDevices::hcl(h = hues, c = 100, l = 65)
      return(colors)
    }
    
    # ------------------------------------
    # 그래프를 그리는 reactive expression
    # ------------------------------------
    barplot_reactive <- reactive({
      req(ps_obj(), group_var(), input$tax_level, input$name_display_mode, 
          input$plot_mode, input$top_n_taxa, input$color_palette)
      
      # 필요한 변수 설정
      current_rank <- input$tax_level
      primary_var <- group_var()
      secondary_var_val <- secondary_var()
      
      # Facet 변수는 항상 Primary 변수
      facet_var <- primary_var 
      topN <- input$top_n_taxa
      plot_mode <- input$plot_mode
      
      # 1. 상대 풍부도로 변환 및 선택된 분류 단위로 통합 (이전과 동일)
      ps_rel <- phyloseq::transform_sample_counts(ps_obj(), function(x) x / sum(x))
      ps_glom <- phyloseq::tax_glom(ps_rel, taxrank = current_rank)
      
      # Full Hierarchy Taxa Name 컬럼 생성 (이전과 동일)
      tax_df <- as.data.frame(phyloseq::tax_table(ps_glom))
      tax_ranks <- phyloseq::rank_names(ps_glom)
      # ... (Full_Taxa_Name 생성 로직은 생략. 이전과 동일함) ...
      tax_df$Full_Taxa_Name <- apply(tax_df, 1, function(taxa) {
        valid_ranks <- taxa[tax_ranks]
        full_name_parts <- c()
        for (i in 1:length(valid_ranks)) {
          part <- valid_ranks[i]
          if (is.na(part) || part == "") {
            part <- paste0(names(valid_ranks)[i], "__Unassigned")
          }
          full_name_parts <- c(full_name_parts, part)
        }
        return(paste(full_name_parts, collapse = ";"))
      })
      phyloseq::tax_table(ps_glom) <- phyloseq::tax_table(as.matrix(tax_df))
      
      # 2. 상위 N개 분류군 추출 및 OTU ID 결정 (이전과 동일)
      sorted_taxa <- names(sort(phyloseq::taxa_sums(ps_glom), decreasing = TRUE))
      top_otu_ids <- sorted_taxa[1:min(topN, length(sorted_taxa))]
      
      # 3. psmelt를 사용하여 데이터프레임으로 변환 (샘플별 원본 데이터)
      df_sample <- phyloseq::psmelt(ps_glom)
      meta_df_base <- as.data.frame(phyloseq::sample_data(ps_rel), stringsAsFactors = FALSE)
      meta_cols <- colnames(meta_df_base)
      
      # psmelt가 비표준 컬럼명(예: patient-id)을 make.names 형태(patient.id)로 바꾼 경우 복원
      for (meta_col in meta_cols) {
        syntactic_col <- make.names(meta_col)
        if (!meta_col %in% names(df_sample) && syntactic_col %in% names(df_sample)) {
          names(df_sample)[names(df_sample) == syntactic_col] <- meta_col
        }
      }
      
      if (!primary_var %in% names(df_sample)) {
        showNotification(
          paste0("Error: Primary grouping variable '", primary_var, "' was not found in melted data."),
          type = "error"
        )
        return(NULL)
      }
      df_sample[[primary_var]] <- factor(df_sample[[primary_var]])
      
      if (!is.null(secondary_var_val)) {
        if (!secondary_var_val %in% names(df_sample)) {
          showNotification(
            paste0("Error: Secondary grouping variable '", secondary_var_val, "' was not found in melted data."),
            type = "error"
          )
          return(NULL)
        }
        df_sample[[secondary_var_val]] <- factor(df_sample[[secondary_var_val]])
      }
      
      # 4 & 5. Taxa Name 및 Group 결정 (이전과 동일)
      if (input$name_display_mode == "full") {
        df_sample$Taxa_Name_Display <- df_sample$Full_Taxa_Name
        fill_label <- paste(current_rank, "Taxa (Full Hierarchy)")
      } else {
        df_sample$Taxa_Name_Display <- df_sample[[current_rank]]
        fill_label <- current_rank
      }
      
      df_sample$Taxa_Group <- ifelse(df_sample$OTU %in% top_otu_ids,
                                     as.character(df_sample$Taxa_Name_Display),
                                     "Others")
      
      # 6. Taxa Group의 레벨 설정 (이전과 동일)
      taxa_mean_abundance <- df_sample %>%
        dplyr::group_by(Taxa_Group) %>%
        dplyr::summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = 'drop') %>%
        dplyr::arrange(desc(Mean_Abundance))
      
      ordered_taxa_levels <- taxa_mean_abundance$Taxa_Group
      if ("Others" %in% ordered_taxa_levels) {
        ordered_taxa_levels <- c(setdiff(ordered_taxa_levels, "Others"), "Others")
      }
      df_sample$Taxa_Group <- factor(df_sample$Taxa_Group, levels = ordered_taxa_levels)
      
      # ⭐ Combined Group 변수 생성 ⭐
      # Group Mean 모드에서 X축이 될 변수
      if (!is.null(secondary_var_val)) {
        
        # 그룹 변수들을 Factor로 변환하여 순서를 올바르게 설정
        meta_df_base[[primary_var]] <- factor(meta_df_base[[primary_var]])
        meta_df_base[[secondary_var_val]] <- factor(meta_df_base[[secondary_var_val]])
        
        # df_sample의 그룹 변수들도 factor로 변환
        df_sample[[primary_var]] <- factor(df_sample[[primary_var]])
        df_sample[[secondary_var_val]] <- factor(df_sample[[secondary_var_val]])
        
        
        # 1. Combined Group Level 순서 결정 (meta_df_base 사용)
        # meta_df_base를 tibble로 변환하여 dplyr 연산 안정화
        meta_df_tibble <- tibble::as_tibble(meta_df_base)
        
        unique_groups <- meta_df_tibble %>%
          dplyr::arrange(.data[[primary_var]], .data[[secondary_var_val]]) %>%
          dplyr::mutate(Combined_Group = paste(.data[[primary_var]], .data[[secondary_var_val]], sep = "_")) %>%
          # ⭐ 수정: pull() 대신 select()를 사용하여 벡터 추출 및 unique 적용을 명확히 함 ⭐
          dplyr::select(Combined_Group) %>% 
          unique() %>%
          dplyr::pull(Combined_Group)
        
        # 2. df_sample에 Combined_Group 변수 생성 및 레벨 적용
        df_sample$Combined_Group <- paste(df_sample[[primary_var]], df_sample[[secondary_var_val]], sep = "_")
        df_sample$Combined_Group <- factor(df_sample$Combined_Group, levels = unique_groups)
        
        combined_title_suffix <- paste0(" by ", primary_var, " & ", secondary_var_val)    
      } else {
        # Secondary 변수가 없으면 Primary 변수가 Combined Group이 됨
        df_sample$Combined_Group <- df_sample[[primary_var]]
        combined_title_suffix <- paste0(" by ", primary_var)
      }      
      
      # ----------------------------------------------------
      # ⭐ Plot Mode에 따른 데이터프레임 및 축 설정 ⭐
      # ----------------------------------------------------
      if (plot_mode == "Group_Mean") {
        # 그룹 평균 계산 (Combined_Group 및 Primary Variable 기준으로 집계)
        df_plot <- df_sample %>%
          # ⭐ 수정: primary_var를 그룹 변수에 추가 ⭐
          dplyr::group_by(Combined_Group, .data[[primary_var]], Taxa_Group) %>% 
          dplyr::summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = 'drop')
        
        # X축은 결합된 그룹 변수가 됨
        x_var <- "Combined_Group"
        x_label <- paste(primary_var, if (!is.null(secondary_var_val)) secondary_var_val else "")
        plot_title <- paste("Mean Relative Abundance of Top", topN, current_rank, "Taxa", combined_title_suffix)
        
        # df_plot$Combined_Group은 이미 factor이지만, 확실하게 한 번 더 지정합니다.
        df_plot[[x_var]] <- factor(df_plot[[x_var]])
        
      } else { # Sample Level (Default)
        df_plot <- df_sample
        x_var <- "Sample"
        x_label <- ""
        plot_title <- paste("Top", topN, current_rank, "Taxa + Others", combined_title_suffix)
        
        # 7. 샘플 정렬 로직 (Sample Mode에서만 실행)
        req(input$sort_method)
        sort_method <- input$sort_method
        ordered_sample_levels <- c()
        meta_df_base <- as.data.frame(phyloseq::sample_data(ps_rel), stringsAsFactors = FALSE)
        if (!"Sample" %in% names(meta_df_base)) { meta_df_base$Sample <- rownames(meta_df_base) }
        
        # 샘플 정렬 시에도 Secondary Group 변수를 고려하여 정렬 기준에 포함
        sorting_vars <- c("Sample", primary_var)
        if (!is.null(secondary_var_val)) { sorting_vars <- c(sorting_vars, secondary_var_val) }
        
        # ... (샘플 정렬 로직은 Primary/Secondary 변수를 모두 고려하도록 수정되어야 함) ...
        # 여기서는 편의상 Primary/Secondary 변수 순서로 그룹화하여 정렬하도록 단순화합니다.
        
        required_vars <- c("Sample", primary_var)
        if (!is.null(secondary_var_val)) { required_vars <- c(required_vars, secondary_var_val) }
        
        if (sort_method == "metadata_var") {
          req(input$sort_metadata_var)
          if (input$sort_metadata_var %in% names(meta_df_base)) {
            required_vars <- c(required_vars, input$sort_metadata_var)
          } else {
            showNotification(paste0("Error: Metadata variable '", input$sort_metadata_var, "' not found in sample data."), type = "error")
            return(NULL)
          }
        }
        
        valid_vars <- required_vars[required_vars %in% names(meta_df_base)]
        meta_df_temp <- meta_df_base[, valid_vars, drop = FALSE]
        meta_df <- unique(meta_df_temp)
        meta_df <- tibble::as_tibble(meta_df)
        req(nrow(meta_df) > 0)
        
        # 정렬 인자 생성 (Primary > Secondary > Sort_Criterion)
        sort_args <- list(dplyr::sym(primary_var))
        if (!is.null(secondary_var_val)) {
          sort_args <- c(sort_args, dplyr::sym(secondary_var_val))
        }
        
        if (sort_method == "taxa_top1") {
          main_top_taxa_name <- as.character(ordered_taxa_levels[1])
          df_sort <- df_plot %>% 
            dplyr::filter(Taxa_Group == main_top_taxa_name) %>% 
            dplyr::group_by_at(c("Sample", primary_var, if(!is.null(secondary_var_val)) secondary_var_val)) %>% 
            dplyr::summarise(Sort_Abundance = sum(Abundance, na.rm = TRUE), .groups = 'drop') 
          
          sort_args_taxa <- c(sort_args, dplyr::quo(dplyr::desc(!!dplyr::sym("Sort_Abundance"))), dplyr::sym("Sample"))
          df_sort_ordered <- df_sort %>% dplyr::arrange(!!!sort_args_taxa)
          ordered_sample_levels <- unique(df_sort_ordered$Sample)
          
        } else if (sort_method == "sample_name") {
          sort_args_name <- c(sort_args, dplyr::sym("Sample"))
          ordered_sample_levels <- meta_df %>% dplyr::arrange(!!!sort_args_name) %>% dplyr::pull(Sample)
          
        } else if (sort_method == "metadata_var") {
          req(input$secondary_sort)
          sort_meta_var <- input$sort_metadata_var; secondary_sort_crit <- input$secondary_sort; meta_df_sortable <- meta_df
          
          if (secondary_sort_crit == "taxa_top1") {
            main_top_taxa_name <- as.character(ordered_taxa_levels[1])
            taxa_abun_data <- df_plot %>% dplyr::filter(Taxa_Group == main_top_taxa_name) %>% dplyr::group_by(Sample) %>% dplyr::summarise(Sort_Abundance = sum(Abundance, na.rm = TRUE), .groups = 'drop')
            meta_df_sortable <- meta_df_sortable %>% dplyr::left_join(taxa_abun_data, by = "Sample")
          }
          
          sort_args_meta <- c(sort_args, dplyr::sym(sort_meta_var))
          if (secondary_sort_crit == "taxa_top1") { 
            sort_args_meta <- c(sort_args_meta, dplyr::quo(dplyr::desc(!!dplyr::sym("Sort_Abundance"))))
          } else if (secondary_sort_crit == "sample_name_asc") { 
            sort_args_meta <- c(sort_args_meta, dplyr::sym("Sample")) 
          }
          
          meta_df_sorted <- meta_df_sortable %>% dplyr::arrange(!!!sort_args_meta)
          ordered_sample_levels <- meta_df_sorted %>% dplyr::pull(Sample) %>% unique()
          
        } else { # Fallback
          sort_args_fallback <- c(sort_args, dplyr::sym("Sample"))
          ordered_sample_levels <- meta_df %>% dplyr::arrange(!!!sort_args_fallback) %>% dplyr::pull(Sample) 
        }
        
        # 샘플 레벨 적용
        df_plot$Sample <- factor(df_plot$Sample, levels = ordered_sample_levels)
      }
      # ----------------------------------------------------
      
      # 색상 로직 (이전과 동일)
      palette_choice <- input$color_palette
      num_top_taxa <- length(setdiff(ordered_taxa_levels, "Others"))
      top_taxa_names <- setdiff(ordered_taxa_levels, "Others")
      # ... (색상 계산 및 할당 로직은 생략. 이전과 동일함) ...
      
      if (palette_choice == "Custom_HCL") {
        taxa_base_colors <- generate_taxa_colors(num_top_taxa)
      } else {
        tryCatch({
          pal_info <- RColorBrewer::brewer.pal.info[palette_choice, ]
          max_colors <- pal_info$maxcolors
          n_colors_to_use <- min(num_top_taxa, max_colors)
          brewer_colors <- RColorBrewer::brewer.pal(n_colors_to_use, palette_choice)
          taxa_base_colors <- rep(brewer_colors, length.out = num_top_taxa)
        }, error = function(e) {
          showNotification(paste0("Warning: RColorBrewer function failed (", e$message, "). Using HCL fallback."), type = "warning", duration = 3)
          taxa_base_colors <- generate_taxa_colors(num_top_taxa)
        })
      }
      
      names(taxa_base_colors) <- top_taxa_names
      taxa_colors <- c(taxa_base_colors, "Others" = "gray80")
      missing_levels <- setdiff(ordered_taxa_levels, names(taxa_colors))
      if (length(missing_levels) > 0) {
        extra_colors <- rep("gray50", length(missing_levels))
        names(extra_colors) <- missing_levels
        taxa_colors <- c(taxa_colors, extra_colors)
      }
      
      # 9. ggplot 생성
      fill_var <- "Taxa_Group"
      
      p <- ggplot2::ggplot(df_plot, ggplot2::aes_string(x = x_var, y = "Abundance", fill = fill_var)) +
        ggplot2::geom_bar(stat = "identity", position = "fill") +
        ggplot2::scale_fill_manual(values = taxa_colors) + 
        
        # Facet 설정 (df_plot에 primary_var가 포함되었으므로 이제 정상 작동)
        if (plot_mode == "Sample" || !is.null(secondary_var_val)) {
          # Sample Mode이거나 Secondary 변수가 있는 경우 Primary 변수로 Facet
          ggplot2::facet_wrap(as.formula(paste0("~`", primary_var, "`")), scales = "free_x")
        } else {
          # Secondary 변수가 없는 Group Mean (Primary Group만 있는 경우)는 Facet을 사용하지 않음
          NULL 
        }
      
      p <- p +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold"),
          # Group Mean 모드일 경우 X축 레이블을 45도 회전
          axis.text.x = ggplot2::element_text(angle = if(plot_mode == "Group_Mean") 45 else 8, hjust = 1, size = if(plot_mode == "Group_Mean") 10 else 8),
          panel.grid.major = ggplot2::element_line(color = "white", linewidth = 0.1),
          panel.grid.minor = ggplot2::element_line(color = "white", linewidth = 0.1),
          axis.title.x = ggplot2::element_blank(), # x축 제목은 비워둠 (Combine Group Label로 대체)
          legend.title = ggplot2::element_text(size = 10)
        ) +
        ggplot2::labs(
          title = plot_title,
          y = "Relative Abundance",
          fill = fill_label
        )
      
      return(p)
    })
    
    # ------------------------------------
    # 출력 및 다운로드 로직 (이전과 동일)
    # ------------------------------------
    output$barplot_out <- renderPlot({ barplot_reactive() },
                                     height = function() { input$plot_height },
                                     width = function() { input$plot_width }
    )
    
    output$download_barplot <- downloadHandler(
      filename = function() { paste0("taxa_barplot_", input$tax_level, "_", input$plot_mode, ".png") },
      content = function(file) {
        plot_width_in <- input$plot_width / 72
        plot_height_in <- input$plot_height / 72
        
        ggplot2::ggsave(file,
                        plot = barplot_reactive(),
                        device = "png",
                        width = plot_width_in,
                        height = plot_height_in,
                        dpi = 300
        )
      }
    )
  })
}
