library(shiny)
library(phyloseq)
library(dplyr)
library(DT)
library(microbiome) # CLR 변환 최적화를 위해 권장

### UI 함수 -------------------------------------------------------------------
mod_preprocessing_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        width = 4,
        wellPanel(
          h4(icon("filter"), "Sample Selection Overview"),
          p("All samples are selected by default."),
          p("Unselect rows in the table to **EXCLUDE** samples."),
          hr(),
          # 정보 표시창
          div(style = "background: #f8f9fa; padding: 10px; border-radius: 5px;",
              verbatimTextOutput(ns("sample_count_info"))
          ),
          hr(),
          # h4(icon("vial"), "Normalization"),
          # selectInput(ns("norm_method"), "Select Transformation Method:",
          #             choices = c(
          #               "TSS (Total Sum Scaling / %)" = "tss",
          #               "CLR (Centered Log-Ratio)" = "clr",
          #               "None (Raw Counts / Filtered)" = "none"
          #             ),
          #             selected = "tss"
          # ),
          # uiOutput(ns("norm_description")),
          # hr(),
          actionButton(ns("reset_selection"), "Reset: Select All Samples", 
                       icon = icon("sync"), class = "btn-secondary btn-sm", width = "100%")
        )
      ),
      column(
        width = 8,
        wellPanel(
          h4(icon("list-check"), "Individual Sample Selection (Click to Toggle)"),
          # 샘플 목록 표
          DTOutput(ns("sample_table"))
        )
      )
    )
  )
}

### Server 함수 ---------------------------------------------------------------
mod_preprocessing_server <- function(id, ps_obj_initial, active_tab) {
  moduleServer(id, function(input, output, session) {
    
    # 내부 반응형 값 저장소
    ps_filtered <- reactiveVal(NULL) 
    ps_normalized <- reactiveVal(NULL) 
    
    # 1. 메타데이터 가공 (테이블용)
    meta_df_for_dt <- reactive({
      req(ps_obj_initial())
      ps <- ps_obj_initial()
      
      meta_df <- data.frame(phyloseq::sample_data(ps))
      
      # SampleID 중복 제거 및 정리
      if ("SampleID" %in% colnames(meta_df)) {
        meta_df <- meta_df %>% dplyr::select(-SampleID)
      }
      
      read_counts <- phyloseq::sample_sums(ps)
      
      df_for_dt <- meta_df %>%
        tibble::rownames_to_column("SampleID") %>% 
        dplyr::mutate(`Read Count` = as.numeric(read_counts[SampleID])) %>%
        dplyr::rename_with(~gsub("_", " ", .x)) %>%
        dplyr::select(`SampleID`, `Read Count`, everything())
      
      return(df_for_dt)
    })
    
    # 2. 샘플 테이블 렌더링 (초기 전체 선택 핵심)
    output$sample_table <- renderDT({
      df <- meta_df_for_dt()
      req(df)
      
      datatable(
        df,
        rownames = FALSE, 
        selection = list(
          mode = 'multiple', 
          selected = seq_len(nrow(df)), # <--- 모든 행을 초기에 강제 선택
          target = 'row'
        ), 
        options = list(
          dom = 'ftip', 
          pageLength = 20,
          columnDefs = list(list(className = 'dt-center', targets = "_all"))
        )
      ) 
    }, server = TRUE)
    
    # 3. 리셋 버튼 (전체 다시 선택)
    proxy <- DT::dataTableProxy('sample_table')
    observeEvent(input$reset_selection, {
      req(meta_df_for_dt())
      DT::selectRows(proxy, seq_len(nrow(meta_df_for_dt())))
    })
    
    # 4. 필터링 로직 (사용자 선택 반영)
    observe({
      # 탭이 활성화되거나 선택이 바뀔 때 실행
      req(ps_obj_initial(), meta_df_for_dt())
      
      # 만약 아무것도 선택되지 않았다면 (초기 로딩 포함) 초기 데이터 사용 시도
      selected_indices <- input$sample_table_rows_selected
      
      if (is.null(selected_indices)) {
        # 초기 상태: 모든 샘플 포함
        ps_new <- ps_obj_initial()
      } else if (length(selected_indices) == 0) {
        # 사용자가 수동으로 다 해제했을 때
        ps_new <- phyloseq::prune_samples(character(0), ps_obj_initial())
      } else {
        # 선택된 샘플만 추출
        selected_samples <- meta_df_for_dt()$SampleID[selected_indices]
        ps_new <- phyloseq::prune_samples(selected_samples, ps_obj_initial())
      }
      
      # Abundance가 0인 Taxa 제거 (Data Clean-up)
      if (phyloseq::nsamples(ps_new) > 0) {
        ps_new <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_new) > 0, ps_new)
      }
      
      ps_filtered(ps_new)
    })
    
    # 5. 정규화 로직
    observe({
      ps_obj <- ps_filtered()
      req(ps_obj)
      
      method <- input$norm_method
      
      if (phyloseq::nsamples(ps_obj) == 0) {
        ps_normalized(ps_obj)
        return()
      }
      
      # 정규화 계산
      ps_res <- tryCatch({
        if (method == "clr") {
          # microbiome 패키지가 있다면 최우선 사용 (안전함)
          if (requireNamespace("microbiome", quietly = TRUE)) {
            microbiome::transform(ps_obj, "clr")
          } else {
            # ALDEx2 방식 (기존 코드 유지)
            otu_mat <- as.matrix(phyloseq::otu_table(ps_obj))
            # ALDEx2는 sample이 열이어야 함
            if (!phyloseq::taxa_are_rows(ps_obj)) otu_mat <- t(otu_mat)
            clr_res <- ALDEx2::aldex.clr(otu_mat, mc.samples = 1, verbose = FALSE)
            # 간략화된 CLR 추출 로직
            clr_mat <- sapply(clr_res@analysisData, function(x) x[,1])
            ps_t <- ps_obj
            phyloseq::otu_table(ps_t) <- phyloseq::otu_table(clr_mat, taxa_are_rows = TRUE)
            ps_t
          }
        } else if (method == "tss") {
          phyloseq::transform_sample_counts(ps_obj, function(x) (x / sum(x)) * 100)
        } else {
          ps_obj
        }
      }, error = function(e) {
        showNotification(paste("Transformation Error:", e$message), type = "error")
        ps_obj
      })
      
      ps_normalized(ps_res)
    })
    
    # 6. UI 보조 출력
    output$norm_description <- renderUI({
      desc <- switch(input$norm_method,
                     "tss" = "Relative abundance (%). Best for Barplots.",
                     "clr" = "Log-ratio transform. Best for PCA/ANCOM.",
                     "none" = "Filtered raw counts."
      )
      p(em(desc), style = "color: gray; font-size: 0.9em;")
    })
    
    output$sample_count_info <- renderText({
      req(ps_obj_initial())
      ps_now <- ps_normalized()
      
      res <- sprintf("Initial Total: %d samples\n", phyloseq::nsamples(ps_obj_initial()))
      if (!is.null(ps_now)) {
        res <- paste0(res, sprintf("Currently Selected: %d samples\n", phyloseq::nsamples(ps_now)))
        res <- paste0(res, sprintf("Active Taxa (ASVs): %d", phyloseq::ntaxa(ps_now)))
      } else {
        res <- paste0(res, "Initializing selection...")
      }
      res
    })
    
    # 7. 모듈 외부로 데이터 전달
    return(list(
      ps_filtered_raw = ps_filtered,
      ps_normalized = ps_normalized
    ))
  })
}