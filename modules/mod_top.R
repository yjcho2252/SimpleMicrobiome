## UI
mod_top_ui <- function(id) {
  ns <- NS(id)
  converter_url <- Sys.getenv("CONVERTER_APP_URL", "https://simplemicrobiome.mglab.org/convert/ui")

  tagList(
    tags$style(HTML(
      paste0(
        "#", ns("top_wrap"), " .top-hero {",
        "  background: linear-gradient(135deg, #f1f6fd 0%, #e5eefb 100%);",
        "  border: 1px solid #c8d8ef;",
        "  border-radius: 16px;",
        "  padding: 22px 24px;",
        "  margin-bottom: 16px;",
        "  box-shadow: 0 4px 14px rgba(20, 56, 99, 0.08);",
        "}",
        "#", ns("top_wrap"), " .top-title {",
        "  font-size: 28px;",
        "  font-weight: 700;",
        "  color: #0f3568;",
        "  margin-bottom: 6px;",
        "}",
        "#", ns("top_wrap"), " .top-subtitle {",
        "  font-size: 14px;",
        "  color: #355074;",
        "  margin-bottom: 6px;",
        "}",
        "#", ns("top_wrap"), " .top-grid { margin: 0 -6px; }",
        "#", ns("top_wrap"), " .top-col { padding: 0 6px; }",
        "#", ns("top_wrap"), " .top-card {",
        "  border: 1px solid #d0ddef;",
        "  border-radius: 14px;",
        "  background: #ffffff;",
        "  padding: 14px 16px 12px 16px;",
        "  margin-bottom: 12px;",
        "  box-shadow: 0 2px 8px rgba(20, 56, 99, 0.05);",
        "}",
        "#", ns("top_wrap"), " .top-card h4 {",
        "  margin-top: 0;",
        "  margin-bottom: 8px;",
        "  font-size: 16px;",
        "  color: #123b70;",
        "}",
        "#", ns("top_wrap"), " .top-list { margin: 0; padding-left: 18px; }",
        "#", ns("top_wrap"), " .top-start {",
        "  background: #e8f1ff;",
        "  border: 1px solid #b8d1f2;",
        "  border-radius: 10px;",
        "  padding: 10px 12px;",
        "  color: #0f3568;",
        "  margin-top: 10px;",
        "  font-size: 13px;",
        "}",
        "#", ns("top_wrap"), " .top-note {",
        "  background: #fff8ea;",
        "  border: 1px solid #f0d9a2;",
        "  border-radius: 10px;",
        "  padding: 10px 12px;",
        "  color: #6b4f00;",
        "}",
        "#", ns("top_wrap"), " details summary { cursor: pointer; color: #0f3568; }",
        "#", ns("top_wrap"), " a { color: #0f4c92; font-weight: 600; text-decoration: none; }",
        "#", ns("top_wrap"), " a:hover { color: #0a3a70; text-decoration: underline; }"
      )
    )),
    div(
      id = ns("top_wrap"),
      div(
        class = "top-hero",
        div(class = "top-title", tagList(icon("microscope"), " SimpleMicrobiome")),
        p(
          class = "top-subtitle",
          "SimpleMicrobiome supports microbiome data loading, preprocessing, visualization, diversity analysis, differential abundance, and network analysis."
        ),
        p(
          class = "top-subtitle",
          "Developed by the Microbial Genomics Laboratory, Kangwon National University."
        ),
        div(
          class = "top-start",
          tagList(
            icon("play"),
            strong(" Start here: "),
            "Go to the ",
            tags$a(
              href = "#",
              onclick = "var el=document.querySelector('.navbar-nav .nav-link[data-value=\"Data Loading\"]'); if(el){el.click();} return false;",
              "Data Loading"
            ),
            " menu first."
          )
        ),
        div(
          class = "top-start",
          tagList(
            icon("right-to-bracket"),
            strong(" QZA Converter (External Service): "),
            tags$a(
              href = converter_url,
              target = "_blank",
              rel = "noopener noreferrer",
              "Open Converter Tool"
            )
          )
        )
      ),
      fluidRow(
        class = "top-grid",
        column(
          width = 7,
          class = "top-col",
          div(
            class = "top-card",
            h4(icon("route"), " Quick Workflow"),
            tags$ol(
              class = "top-list",
              tags$li("Data Loading: Upload feature table, taxonomy, and metadata."),
              tags$li("Preprocessing: Filter data by read counts and user-selected options."),
              tags$li("Analysis: Run visualization/diversity/DA/network modules as needed.")
            )
          )
        ),
        column(
          width = 5,
          class = "top-col",
          div(
            class = "top-card",
            h4(icon("shield-halved"), " Notes"),
            div(
              class = "top-note",
              tags$ul(
                class = "top-list",
                tags$li("Use enough samples per group for stable statistics."),
                tags$li("Network outputs are sensitive to filtering settings.")
              )
            )
          )
        )
      ),
      div(
        class = "top-card",
        h4(icon("clock-rotate-left"), " Version History"),
        tags$details(
          tags$summary("Show updates"),
          tags$ul(
            class = "top-list",
            tags$li("260408: Network Analysis modules added."),
            tags$li("260407: Random Forest module added."),
            tags$li("260406: MaAsLin2 module added."),
            tags$li("260321: Taxa Comparison feature added."),
            tags$li("260310: Preprocessing functionality improved."),
            tags$li("251210: First public version released.")
          )
        )
      )
    )
  )
}

## Server
mod_top_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    invisible(NULL)
  })
}
