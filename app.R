library(shiny)
library(phyloseq)
library(tidyverse)
library(ANCOMBC)
library(DT)
library(bslib)
library(ggplot2)
library(ggpubr)
library(shinyWidgets)
library(Maaslin2)
library(microbiome)
library(cluster)
library(vegan)
library(dplyr)
library(randomForest)
library(ComplexHeatmap)
library(NetCoMi)
library(igraph)
library(ggraph)
library(ggrepel)
library(shapr)

app_script_path <- tryCatch({
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE)
}, error = function(e) {
  ""
})
app_root_dir <- if (nzchar(app_script_path)) dirname(app_script_path) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
options(simplemicrobiome_app_dir = app_root_dir)
options(shiny.error = function(...) {
  args <- list(...)
  msg <- if (length(args) >= 1 && inherits(args[[1]], "condition")) {
    conditionMessage(args[[1]])
  } else if (length(args) >= 1) {
    as.character(args[[1]])
  } else {
    geterrmessage()
  }
  cat(
    sprintf(
      "[%s] [SHINY-ERROR] %s\n",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      msg
    ),
    file = stderr()
  )
})

## Module Loading
mod_path <- file.path(app_root_dir, "modules")
if (dir.exists(mod_path)) {
  list.files(mod_path, pattern = "\\.R$", full.names = TRUE) %>% walk(source)
}

## App UI
ui <- page_navbar(
  title = actionLink(
    "go_home",
    label = tags$img(
      src = "brand-wordmark.svg",
      alt = "SimpleMicrobiome",
      style = "height: 30px; width: auto; display: block;"
    ),
    style = "display: inline-flex; align-items: center; text-decoration: none; padding: 0;"
  ),
  window_title = "SimpleMicrobiome",
  id = "tab_panel_main",
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2c3e50"
  ),
  header = tags$head(
    tags$title("SimpleMicrobiome"),
    tags$link(rel = "icon", type = "image/svg+xml", href = "favicon.svg"),
    tags$style(HTML("
      body { font-size: 14px; }
      .navbar { min-height: 40px; }
      .navbar-brand { font-size: 15px; padding-top: 6px; padding-bottom: 6px; line-height: 1.1; }
      .navbar .container-fluid { gap: 0.2rem; }
      .navbar-brand { margin-right: 0.2rem !important; }
      .navbar-nav .nav-link { font-size: 13px; padding-top: 6px; padding-bottom: 6px; line-height: 1.1; }
      .navbar-nav .nav-link[data-value='Top'] { display: none !important; }
      .dropdown-menu { font-size: 13px; }
      .btn, .form-control, .form-select, .form-check-label { font-size: 13px; }
      h4 { font-size: 1.15rem; }
      h5 { font-size: 1.02rem; }
      .sidebar h4 { font-size: 0.98rem; margin-bottom: 0.45rem; }
      .sidebar h5 { font-size: 0.90rem; }
      .sidebar .control-label,
      .sidebar .form-label,
      .sidebar .form-check-label,
      .sidebar .help-block,
      .sidebar details summary,
      .sidebar .shiny-input-container label {
        font-size: 12px;
      }
      .sidebar .form-control,
      .sidebar .form-select,
      .sidebar .btn {
        font-size: 12px;
      }
      .selectize-input,
      .selectize-dropdown,
      .selectize-dropdown-content,
      .selectize-dropdown .option,
      .selectize-dropdown .optgroup-header,
      select.form-select,
      .irs-grid-text {
        font-size: 12px !important;
      }
      .sidebar .shiny-input-container {
        margin-bottom: 5px;
      }
      .sidebar .control-label,
      .sidebar .form-label {
        margin-bottom: 2px;
      }
      .sidebar .form-group {
        margin-bottom: 5px;
      }
      .sidebar hr {
        margin-top: 6px;
        margin-bottom: 6px;
      }
      .sidebar details {
        margin-top: 2px !important;
        margin-bottom: 2px !important;
      }
      .sidebar .help-block {
        margin-top: 1px;
        margin-bottom: 2px;
      }
    "))
  ),
  
  nav_panel(
    title = "Top",
    icon = icon("house"),
    mod_top_ui("mod_top")
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
    nav_panel("MaAsLin2", icon = icon("flask"), mod_maaslin2_ui("mod_maaslin2")),
    nav_panel("Random Forest", icon = icon("tree"), mod_randomforest_ui("mod_randomforest"))
  ),
  
  nav_menu(
    title = "Network Analysis",
    icon = icon("diagram-project"),
    nav_panel("SparCC", icon = icon("share-nodes"), mod_sparcc_ui("mod_sparcc")),
    nav_panel("SpiecEasi", icon = icon("diagram-project"), mod_spieceasi_ui("mod_spieceasi"))
  ),
  nav_panel(
    title = "Citation",
    icon = icon("book"),
    mod_citation_ui("mod_citation")
  )
  
)

## App Server
server <- function(input, output, session) {
  observeEvent(input$go_home, {
    updateTabsetPanel(session, inputId = "tab_panel_main", selected = "Top")
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
  
  mod_barplot_server("mod_barplot", ps_obj_filtered_raw, meta_vars)
  mod_taxa_comparison_server("mod_taxa_comparison", ps_obj_filtered_raw, meta_vars)
  mod_alpha_server("mod_alpha", ps_obj_filtered_raw, meta_vars)
  mod_beta_server("mod_beta", ps_obj_filtered_raw, meta_vars)
  mod_ancom_server("mod_ancom", ps_obj_filtered_raw)
  mod_maaslin2_server("mod_maaslin2", ps_obj_filtered_raw)
  mod_randomforest_server("mod_randomforest", ps_obj_filtered_raw)
  mod_spieceasi_server("mod_spieceasi", ps_obj_filtered_raw)
  mod_sparcc_server("mod_sparcc", ps_obj_filtered_raw)
  mod_citation_server("mod_citation")
  mod_top_server("mod_top")
  
  observeEvent(input[["mod_fileload-reset_all_app"]], {
    session$reload()
  })
  
}

shinyApp(ui, server)

