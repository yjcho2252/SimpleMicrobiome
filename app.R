library(shiny)
library(phyloseq)
library(tidyverse)
library(ANCOMBC)
library(DT)
library(bslib)
library(ggplot2)
library(ALDEx2)
library(ggpubr)
library(shinyWidgets)
library(Maaslin2)

## Module Loading
mod_path <- file.path(getwd(), "modules")
if (dir.exists(mod_path)) {
  list.files(mod_path, pattern = "\\.R$", full.names = TRUE) %>% walk(source)
}

## App UI
ui <- page_navbar(
  title = actionLink(
    "go_home",
    label = tagList(icon("microscope"), " SimpleMicrobiome"),
    style = "font-weight: 700; letter-spacing: 0.3px; font-size: 16px; color: inherit; text-decoration: none; padding: 0;"
  ),
  id = "tab_panel_main",
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2c3e50"
  ),
  header = tags$head(
    tags$style(HTML("
      body { font-size: 14px; }
      .navbar { min-height: 40px; }
      .navbar-brand { font-size: 15px; padding-top: 6px; padding-bottom: 6px; line-height: 1.1; }
      .navbar-nav .nav-link { font-size: 13px; padding-top: 6px; padding-bottom: 6px; line-height: 1.1; }
      .dropdown-menu { font-size: 13px; }
      .btn, .form-control, .form-select, .form-check-label { font-size: 13px; }
      h4 { font-size: 1.15rem; }
      h5 { font-size: 1.02rem; }
    "))
  ),
  
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
  
  nav_menu(
    title = "Differential Abundance",
    icon = icon("vial-circle-check"),
    nav_panel("ANCOM-BC2", icon = icon("vial-circle-check"), mod_ancom_ui("mod_ancom")),
    nav_panel("MaAsLin2", icon = icon("flask"), mod_maaslin2_ui("mod_maaslin2"))
  ),
  
  nav_menu(
    title = "Network Analysis",
    icon = icon("diagram-project"),
    nav_panel("SpiecEasi", icon = icon("diagram-project"), mod_spieceasi_ui("mod_spieceasi"))
  )
  
)

## App Server
server <- function(input, output, session) {
  observeEvent(input$go_home, {
    updateTabsetPanel(session, inputId = "tab_panel_main", selected = "Data Loading")
  })
  
  file_data <- mod_fileload_server("mod_fileload")
  ps_obj_initial <- file_data$ps_initial
  meta_vars <- file_data$meta_vars
  
  active_tab <- reactive({ input$tab_panel_main })
  
  observeEvent(ps_obj_initial(), {
    req(ps_obj_initial())
    updateTabsetPanel(session, inputId = "tab_panel_main", selected = "Preprocessing")
  }, ignoreInit = TRUE)
  
  preprocessing_data <- mod_preprocessing_server("mod_preprocessing", ps_obj_initial, active_tab)
  
  ps_obj_filtered_raw <- preprocessing_data$ps_filtered_raw
  ps_obj_normalized <- preprocessing_data$ps_normalized
  
  mod_barplot_server("mod_barplot", ps_obj_normalized, meta_vars)
  mod_taxa_comparison_server("mod_taxa_comparison", ps_obj_filtered_raw, meta_vars)
  mod_alpha_server("mod_alpha", ps_obj_filtered_raw, meta_vars)
  mod_beta_server("mod_beta", ps_obj_normalized, meta_vars)
  mod_ancom_server("mod_ancom", ps_obj_filtered_raw)
  mod_maaslin2_server("mod_maaslin2", ps_obj_filtered_raw)
  mod_spieceasi_server("mod_spieceasi", ps_obj_filtered_raw)
  
  observeEvent(input[["mod_fileload-reset_all_app"]], {
    session$reload()
  })
  
}

shinyApp(ui, server)

