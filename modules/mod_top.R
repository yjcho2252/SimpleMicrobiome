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
        "  max-width: 920px;",
        "  margin-left: 0;",
        "  margin-right: 0;",
        "  box-shadow: 0 4px 14px rgba(20, 56, 99, 0.08);",
        "}",
        "#", ns("top_wrap"), " .top-title {",
        "  display: inline-flex;",
        "  align-items: center;",
        "  margin-bottom: 14px;",
        "}",
        "#", ns("top_wrap"), " .top-title img {",
        "  height: 42px;",
        "  width: auto;",
        "  display: block;",
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
        "#", ns("top_wrap"), " .top-card-readable {",
        "  max-width: 920px;",
        "  margin-left: 0;",
        "  margin-right: 0;",
        "}",
        "#", ns("top_wrap"), " .top-card h4 {",
        "  margin-top: 0;",
        "  margin-bottom: 8px;",
        "  font-size: 16px;",
        "  color: #123b70;",
        "}",
        "#", ns("top_wrap"), " .top-workflow-img {",
        "  width: 100%;",
        "  height: auto;",
        "  display: block;",
        "  border: 1px solid #d7e4f5;",
        "  border-radius: 12px;",
        "  margin-bottom: 10px;",
        "  background: #f8fbff;",
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
        "#", ns("top_wrap"), " a:hover { color: #0a3a70; text-decoration: underline; }",
        "#", ns("top_wrap"), " .top-footer {",
        "  display: flex;",
        "  align-items: center;",
        "  justify-content: center;",
        "  gap: 8px;",
        "  margin-top: 4px;",
        "  max-width: 920px;",
        "  width: 100%;",
        "  color: #355074;",
        "  font-size: 13px;",
        "}",
        "#", ns("top_wrap"), " .top-footer img {",
        "  height: 18px;",
        "  width: auto;",
        "  display: block;",
        "}"
      )
    )),
    div(
      id = ns("top_wrap"),
      div(
        class = "top-hero",
        div(
          class = "top-title",
          tags$img(src = "brand-wordmark-dark.svg", alt = "SimpleMicrobiome")
        ),
        p(
          class = "top-subtitle",
          "SimpleMicrobiome enables key microbiome analyses and visualizations."
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
        
      ),
      fluidRow(
        class = "top-grid",
        column(
          width = 12,
          class = "top-col",
          div(
            class = "top-card top-card-readable",
            h4(icon("route"), " Quick Workflow"),
            tags$ol(
              class = "top-list",
              tags$li("Data Loading: Upload feature table, taxonomy, and metadata."),
              tags$li("Preprocessing: Filter data by read counts and user-selected options."),
              tags$li("Analysis: Run Taxa Profiles/Diversity/DA/Network/Association modules as needed.")
            )
          )
        )
      ),
      div(
        class = "top-card top-card-readable",
        h4(icon("clock-rotate-left"), " Version History"),
        tags$details(
          tags$summary("Show updates"),
          tags$ul(
            class = "top-list",
            tags$li("260423: UI improvements and bug fixes."),
            tags$li("260421: Bar plot in Taxa Profiles and Alpha Diversity added."),
            tags$li("260420: Association analysis modules added."),
            tags$li("260417: K-fold CV in RF added."),
            tags$li("260415: QZA Converter added."),
            tags$li("260414: Network Comparison function added."),
            tags$li("260408: Network Analysis modules added."),
            tags$li("260407: Random Forest module added."),
            tags$li("260406: MaAsLin2 module added."),
            tags$li("260321: Taxa Comparison feature added."),
            tags$li("260310: Preprocessing functionality improved."),
            tags$li("251210: First public version released.")
          )
        )
      ),
      div(
        class = "top-footer",
        tags$img(
          src = "https://www.kangwon.ac.kr/assets/ko/images/sub/symbol1.webp",
          alt = "Kangwon National University"
        ),
        tags$span("\u00a9 2025-2026 SimpleMicrobiome | Microbial Genomics Lab"),
        tags$br(),
        tags$span("Contact: yongjoon (at) kangwon.ac.kr")
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
