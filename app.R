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
library(maaslin3)
library(microbiome)
library(cluster)
library(vegan)
library(randomForest)
library(ComplexHeatmap)
library(NetCoMi)
library(igraph)
library(ggraph)
library(ggrepel)
library(shapr)
library(ggpattern)
library(httr2)
library(readr)
library(dplyr)

app_script_path <- tryCatch({
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE)
}, error = function(e) {
  ""
})
app_root_dir <- if (nzchar(app_script_path)) dirname(app_script_path) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
options(simplemicrobiome_app_dir = app_root_dir)
options(shiny.maxRequestSize = 100 * 1024^2) # 100MB upload limit
addResourcePath("docs", file.path(app_root_dir, "docs"))
addResourcePath("manual", file.path(app_root_dir, "docs", "manual_html"))
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
shared_path <- file.path(mod_path, "shared")
list.files(shared_path, pattern = "\\.R$", full.names = TRUE) %>% walk(source)
list.files(mod_path, pattern = "^mod_.*\\.R$", full.names = TRUE) %>% walk(source)

## App UI
ui <- page_navbar(
  title = actionLink(
    "go_home",
    label = tagList(
      tags$img(
        src = "brand-wordmark.svg",
        alt = "SimpleMicrobiome",
        style = "height: 30px; width: auto; display: block;"
      ),
      tags$span("v1.2.1", class = "navbar-app-version")
    ),
    style = "display: inline-flex; align-items: center; gap: 0.4rem; text-decoration: none; padding: 0;"
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
    tags$link(rel = "icon", type = "image/png", sizes = "64x64", href = "260617_SimpleMicrobiome_icon_64.png"),
    tags$link(rel = "shortcut icon", type = "image/png", href = "260617_SimpleMicrobiome_icon_64.png"),
    tags$script(src = "plot-scroll-layout.js", defer = NA),
    tags$script(src = "plot-legend-readiness.js", defer = NA),
    tags$script(src = "profile-group-selectors.js", defer = NA),
    tags$style(HTML("
      body { font-size: 14px; }
      .app-plot-legend:not([data-legend-ready]) { display: none; }
      .app-result-card,
      .app-plot-legend > div {
        max-width: 100%;
        margin-top: 12px;
        padding: 12px 14px;
        border: 1px solid #e5e7eb;
        border-left: 4px solid #6b7280;
        border-radius: 8px;
        background: linear-gradient(180deg, #fcfcfd 0%, #f7f8fa 100%);
        box-shadow: 0 1px 2px rgba(0, 0, 0, 0.04);
        box-sizing: border-box;
        color: #1f2937;
        font-size: 12.5px;
        line-height: 1.55;
        overflow-wrap: anywhere;
      }
      .app-status-card {
        width: 600px;
        height: 130px;
        overflow: auto;
        margin-bottom: 10px;
      }
      .app-status-card pre {
        margin: 0;
        padding: 0;
        border: 0;
        background: transparent;
        color: inherit;
        font-family: inherit;
        font-size: inherit;
        line-height: inherit;
        white-space: pre-wrap;
        overflow-wrap: anywhere;
      }
      .app-status-card p { margin-bottom: 0.4rem; }
      .profile-group-selectors[disabled] > .shiny-input-container { pointer-events: none; }
      :root {
        --taxa-plot-loading-width: 900px;
        --taxa-plot-loading-height: 700px;
      }
      .taxa-comparison-result > .taxa-plot-stage:has(:is(.taxa-plot-container, .shiny-plot-output):is(:empty, .recalculating, .shiny-output-error)) + div,
      :is(.taxa-plot-stage, .taxa-barplot-view):has(.plotly.html-widget-output:is(.recalculating, .shiny-output-error, :not([data-taxa-ready]))) + .taxa-plot-details,
      :is(.taxa-plot-stage, .taxa-barplot-view):has(.taxa-plot-container:is(:empty, .recalculating, .shiny-output-error)) + .taxa-plot-details {
        display: none;
      }
      .taxa-plot-stage { position: relative; }
      .taxa-plot-stage .shiny-output-error-taxa-plot-wait {
        display: flex;
        align-items: center;
        justify-content: center;
        width: min(var(--taxa-plot-loading-width), 100%);
        min-height: var(--taxa-plot-loading-height);
        text-align: center;
        color: #334155;
      }
      .taxa-plot-stage:has(.plotly.html-widget-output:not(.shiny-output-error):not([data-taxa-ready])),
      .taxa-plot-stage:has(.shiny-plot-output:empty:not(.shiny-output-error)),
      .taxa-plot-stage:has(.taxa-plot-container:empty:not(.shiny-output-error)) {
        min-height: var(--taxa-plot-loading-height);
      }
      .taxa-plot-stage:has(.taxa-plot-container:empty:not(.shiny-output-error)) {
        width: min(var(--taxa-plot-loading-width), 100%);
      }
      .taxa-plot-loading {
        position: absolute;
        inset: 0;
        max-width: var(--taxa-plot-loading-width);
        display: none;
        align-items: center;
        justify-content: center;
        pointer-events: none;
        z-index: 1;
      }
      .taxa-plot-loading span {
        padding: 8px 12px;
        border-radius: 6px;
        background: rgba(255, 255, 255, 0.94);
        color: #334155;
      }
      .taxa-plot-stage:has(.plotly.html-widget-output:not(.shiny-output-error):not([data-taxa-ready])) > .taxa-plot-loading,
      .taxa-plot-stage:has(.plotly.html-widget-output.recalculating:not(.shiny-output-error)) > .taxa-plot-loading,
      .taxa-plot-stage:has(.shiny-plot-output:is(:empty, .recalculating):not(.shiny-output-error)) > .taxa-plot-loading,
      .taxa-plot-stage:has(.taxa-plot-container.recalculating:not(.shiny-output-error)) > .taxa-plot-loading,
      .taxa-plot-stage:has(.taxa-plot-container:empty:not(.shiny-output-error)) > .taxa-plot-loading {
        display: flex;
      }
      @media (min-width: 768px) and (min-height: 480px) {
        .plot-scroll-layout > div {
          max-height: var(--plot-scroll-height, calc(100vh - 5rem));
          overflow: auto;
          overscroll-behavior: contain;
        }
      }
      .navbar { min-height: 40px; }
      .navbar-brand { font-size: 15px; padding-top: 6px; padding-bottom: 6px; line-height: 1.1; }
      .navbar-app-version { color: #fff; font-size: 11px; white-space: nowrap; }
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
      .selectize-input > input,
      .selectize-input > input::placeholder,
      .selectize-dropdown,
      .selectize-dropdown-content,
      .selectize-dropdown .option,
      .selectize-dropdown .optgroup-header,
      select.form-select,
      select.form-control,
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
    title = "Taxa Profiles",
    icon = icon("chart-line"),
    nav_panel("Taxa Bar Plot", icon = icon("chart-bar"), mod_barplot_ui("mod_barplot")),
    nav_panel("Taxa Bubble Plot", icon = icon("circle"), mod_taxa_summary_ui("mod_taxa_summary")),
    nav_panel("Taxa Comparison", icon = icon("square-poll-vertical"), mod_taxa_comparison_ui("mod_taxa_comparison"))
  ),
  
  nav_menu(
    title = "Diversity",
    icon = icon("dna"),
    nav_panel("Alpha Diversity", icon = icon("circle-nodes"), mod_alpha_ui("mod_alpha")),
    nav_panel("Beta Diversity", icon = icon("project-diagram"), mod_beta_ui("mod_beta")),
    nav_panel("Constrained Ordination", icon = icon("arrows-up-down-left-right"), mod_biplot_ui("mod_biplot"))
  ),
  
  nav_menu(
    title = "Differential Abundance",
    icon = icon("vial-circle-check"),
    nav_panel("ANCOM-BC2", icon = icon("vial-circle-check"), mod_ancom_ui("mod_ancom")),
    nav_panel("MaAsLin2", icon = icon("flask"), mod_maaslin2_ui("mod_maaslin2")),
    nav_panel("MaAsLin3", icon = icon("flask"), mod_maaslin3_ui("mod_maaslin3")),
    nav_panel("Random Forest", icon = icon("tree"), mod_randomforest_ui("mod_randomforest"))
  ),
  
  nav_menu(
    title = "Network",
    icon = icon("diagram-project"),
    nav_panel("SparCC", icon = icon("share-nodes"), mod_sparcc_ui("mod_sparcc")),
    nav_panel("SPIEC-EASI", icon = icon("diagram-project"), mod_spieceasi_ui("mod_spieceasi"))
  ),
  nav_menu(
    title = "Association",
    icon = icon("diagram-next"),
    nav_panel("Association Heatmap", icon = icon("table-cells"), mod_heatmap_ui("mod_heatmap")),
    nav_panel("Taxa Correlogram", icon = icon("braille"), mod_taxa_correlogram_ui("mod_taxa_correlogram")),
    nav_panel("Upset Plot", icon = icon("layer-group"), mod_upset_ui("mod_upset"))
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
  session$userData$active_plot_tab <- active_tab
  
  preprocessing_data <- mod_preprocessing_server("mod_preprocessing", ps_obj_initial, active_tab)
  
  ps_obj_filtered_raw <- preprocessing_data$ps_filtered_raw
  
  plot_order <- preprocessing_data$plot_order
  mod_barplot_server("mod_barplot", preprocessing_data$ps_selected_raw, meta_vars, plot_order)
  mod_taxa_summary_server("mod_taxa_summary", preprocessing_data$ps_selected_raw, plot_order)
  mod_taxa_comparison_server("mod_taxa_comparison", ps_obj_filtered_raw, meta_vars, active_tab, plot_order)
  mod_alpha_server("mod_alpha", ps_obj_filtered_raw, meta_vars, active_tab, plot_order)
  mod_beta_server("mod_beta", ps_obj_filtered_raw, meta_vars, plot_order)
  mod_ancom_server("mod_ancom", ps_obj_filtered_raw, plot_order)
  mod_maaslin2_server("mod_maaslin2", ps_obj_filtered_raw, plot_order)
  mod_maaslin3_server("mod_maaslin3", ps_obj_filtered_raw, plot_order)
  mod_randomforest_server("mod_randomforest", ps_obj_filtered_raw)
  mod_spieceasi_server("mod_spieceasi", ps_obj_filtered_raw)
  mod_sparcc_server("mod_sparcc", ps_obj_filtered_raw)
  mod_heatmap_server("mod_heatmap", ps_obj_filtered_raw, meta_vars, plot_order,
                     active = reactive(identical(active_tab(), "Association Heatmap")))
  mod_taxa_correlogram_server("mod_taxa_correlogram", ps_obj_filtered_raw,
                              active = reactive(identical(active_tab(), "Taxa Correlogram")))
  mod_upset_server("mod_upset", ps_obj_filtered_raw, plot_order,
                   active = reactive(identical(active_tab(), "Upset Plot")))
  mod_biplot_server("mod_biplot", ps_obj_filtered_raw, meta_vars)
  mod_citation_server("mod_citation")
  mod_top_server("mod_top")
  
  observeEvent(input[["mod_fileload-reset_all_app"]], {
    session$reload()
  })
  
}

shinyApp(ui, server)
