library(shiny)
library(phyloseq)
library(tidyverse)
library(ANCOMBC)
library(DT)
library(bslib)
library(ggplot2)
library(ALDEx2)
library(microbiome)
library(ggpubr)
library(shinyWidgets) # 추가적인 UI 위젯 활용

# 모듈 로드
mod_path <- file.path(getwd(), "modules")
if (dir.exists(mod_path)) {
  list.files(mod_path, pattern = "\\.R$", full.names = TRUE) %>% walk(source)
}

## 1. UI 정의 -----------------------------------------------------------------
ui <- page_navbar(
  title = span(
    tagList(icon("microscope"), " SimpleMicrobiome"),
    style = "font-weight: bold; letter-spacing: 0.5px;"
  ),
  id = "tab_panel_main",
  theme = bs_theme(
    version = 5, # Bootstrap 5 활용
    bootswatch = "flatly", # 전문적이고 깨끗한 테마 (또는 "lux")
    primary = "#2c3e50"
  ),
  
  # 상단 탭 구성
  nav_panel(
    title = "Data Loading", 
    icon = icon("upload"),
    mod_fileload_ui("mod_fileload")
  ),
  nav_panel(
    title = "Preprocessing", 
    icon = icon("sliders"),
    mod_preprocessing_ui("mod_preprocessing")
  ),
  
  # 분석 메뉴 그룹화 (드롭다운으로 묶어 전문성 향상)
  nav_menu(
    title = "Visualization",
    icon = icon("chart-line"),
    nav_panel("Taxa Barplot", icon = icon("chart-bar"), mod_barplot_ui("mod_barplot")),
    nav_panel("Taxa Comparison", icon = icon("square-poll-vertical"), mod_taxa_comparison_ui("mod_taxa_comparison"))
  ),
  
  nav_menu(
    title = "Diversity Analysis",
    icon = icon("dna"),
    nav_panel("Alpha Diversity", icon = icon("circle-nodes"), mod_alpha_ui("mod_alpha")),
    nav_panel("Beta Diversity", icon = icon("project-diagram"), mod_beta_ui("mod_beta"))
  ),
  
  nav_panel(
    title = "Differential Abundance", 
    icon = icon("vial-circle-check"),
    mod_ancom_ui("mod_ancom")
  )
  
  # 우측 상단 정보 (선택 사항)
  #nav_spacer(),
  #nav_item(
  #  tags$a(icon("github"), "Source", href = "#", target = "_blank", class = "nav-link")
  #)
)

## 2. Server 정의 -------------------------------------------------------------
server <- function(input, output, session) {
  
  # 데이터 로딩 모듈
  file_data <- mod_fileload_server("mod_fileload")
  ps_obj_initial <- file_data$ps_initial
  meta_vars <- file_data$meta_vars
  
  # 현재 활성화된 탭 추적
  active_tab <- reactive({ input$tab_panel_main })
  
  # 데이터가 로드되면 자동으로 전처리 탭으로 이동 (사용자 경험 개선)
  observeEvent(ps_obj_initial(), {
    req(ps_obj_initial())
    # SweetAlert 등으로 로딩 완료 알림을 주면 더 전문적입니다.
    # shinyWidgets::sendSweetAlert(session, title = "Data Loaded", type = "success")
    # nav_select(id = "tab_panel_main", selected = "Preprocessing")
    updateTabsetPanel(session, inputId = "tab_panel_main", selected = "Data Preprocessing")
  }, ignoreInit = TRUE)
  
  # 전처리 모듈
  preprocessing_data <- mod_preprocessing_server("mod_preprocessing", ps_obj_initial, active_tab)
  
  ps_obj_filtered_raw <- preprocessing_data$ps_filtered_raw
  ps_obj_normalized <- preprocessing_data$ps_normalized
  
  # 분석 모듈들 호출
  mod_barplot_server("mod_barplot", ps_obj_normalized, meta_vars)
  mod_taxa_comparison_server("mod_taxa_comparison", ps_obj_filtered_raw, meta_vars)
  mod_alpha_server("mod_alpha", ps_obj_filtered_raw, meta_vars)
  mod_beta_server("mod_beta", ps_obj_normalized, meta_vars)
  mod_ancom_server("mod_ancom", ps_obj_filtered_raw)
  
}

shinyApp(ui, server)