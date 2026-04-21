library(shiny)
library(DT)
library(readr)
library(httr2)

ui <- fluidPage(
  titlePanel("QIIME2 qza to TSV Converter (API Client)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("tax_qza", "Upload Taxonomy qza", accept = ".qza"),
      fileInput("feat_qza", "Upload Feature Table qza", accept = ".qza"),
      actionButton("convert_btn", "Convert via Service"),
      tags$div(
        style = "margin-top: 8px;",
        textOutput("progress_hint")
      ),
      br(), br(),
      downloadButton("download_tax", "Download Taxonomy TSV"),
      br(), br(),
      downloadButton("download_feat", "Download Feature Table TSV"),
      br(), br(),
      verbatimTextOutput("status")
    ,
      tags$script(HTML(
        "Shiny.addCustomMessageHandler('toggle-convert-btn', function(msg) {
           var btn = document.getElementById(msg.id);
           if (!btn) return;
           btn.disabled = !!msg.disabled;
           if (msg.label) btn.textContent = msg.label;
         });"
      ))
    ),
    mainPanel(
      h4("Taxonomy Preview"),
      DTOutput("tax_preview"),
      br(),
      h4("Feature Table Preview"),
      DTOutput("feat_preview")
    )
  )
)

server <- function(input, output, session) {
  tax_data <- reactiveVal(NULL)
  feat_data <- reactiveVal(NULL)
  tax_file_path <- reactiveVal(NULL)
  feat_file_path <- reactiveVal(NULL)
  status_text <- reactiveVal("Ready.")
  is_running <- reactiveVal(FALSE)

  output$status <- renderText(status_text())
  output$progress_hint <- renderText({
    if (is_running()) "Converting files... please wait." else ""
  })

  normalize_base_url <- function(x) {
    out <- trimws(as.character(x))
    out <- sub("/+$", "", out)
    out
  }
  service_url <- normalize_base_url(Sys.getenv("QIIME_CONVERTER_URL", "https://simplemicrobiome.mglab.org/convert"))

  post_file_and_save <- function(base_url, endpoint, field_name, source_file, out_path) {
    req <- request(paste0(base_url, endpoint)) |>
      req_method("POST") |>
      req_body_multipart(!!field_name := upload_file(source_file)) |>
      req_timeout(300)

    resp <- req_perform(req)
    if (resp_status(resp) >= 400) {
      stop(resp_body_string(resp))
    }

    raw <- resp_body_raw(resp)
    writeBin(raw, out_path)
    out_path
  }

  observeEvent(input$convert_btn, {
    req(input$tax_qza, input$feat_qza)

    base_url <- service_url
    validate(need(nzchar(base_url), "Converter service URL is not configured."))

    is_running(TRUE)
    session$sendCustomMessage(
      "toggle-convert-btn",
      list(id = session$ns("convert_btn"), disabled = TRUE, label = "Converting...")
    )
    on.exit({
      is_running(FALSE)
      session$sendCustomMessage(
        "toggle-convert-btn",
        list(id = session$ns("convert_btn"), disabled = FALSE, label = "Convert via Service")
      )
    }, add = TRUE)

    status_text("Checking converter service...")
    tax_data(NULL)
    feat_data(NULL)
    tax_file_path(NULL)
    feat_file_path(NULL)

    # Create per-run temp outputs for download handlers.
    tmp_dir <- tempfile("qza_client_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    tax_out <- file.path(tmp_dir, "taxonomy.tsv")
    feat_out <- file.path(tmp_dir, "feature_table.tsv")

    health_ok <- tryCatch({
      resp <- request(paste0(base_url, "/health")) |>
        req_method("GET") |>
        req_timeout(15) |>
        req_perform()
      resp_status(resp) < 400
    }, error = function(e) {
      FALSE
    })

    if (!health_ok) {
      status_text("Converter service is unreachable. Check service URL and server status.")
      return()
    }

    status_text("Converting taxonomy via service...")
    tax_res <- tryCatch({
      post_file_and_save(
        base_url = base_url,
        endpoint = "/convert/taxonomy",
        field_name = "taxonomy_qza",
        source_file = input$tax_qza$datapath,
        out_path = tax_out
      )
    }, error = function(e) {
      NULL
    })

    if (is.null(tax_res) || !file.exists(tax_out)) {
      status_text("Taxonomy conversion failed.")
      return()
    }

    tax_df <- tryCatch(
      read_tsv(tax_out, show_col_types = FALSE),
      error = function(e) NULL
    )
    if (is.null(tax_df)) {
      status_text("Taxonomy conversion succeeded, but preview parsing failed.")
      tax_file_path(tax_out)
      return()
    }
    tax_data(tax_df)
    tax_file_path(tax_out)

    status_text("Converting feature table via service...")
    feat_res <- tryCatch({
      post_file_and_save(
        base_url = base_url,
        endpoint = "/convert/feature-table",
        field_name = "feature_table_qza",
        source_file = input$feat_qza$datapath,
        out_path = feat_out
      )
    }, error = function(e) {
      NULL
    })

    if (is.null(feat_res) || !file.exists(feat_out)) {
      status_text("Feature table conversion failed.")
      return()
    }

    feat_df <- tryCatch(
      read.delim(feat_out, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(feat_df)) {
      status_text("Feature table conversion succeeded, but preview parsing failed.")
      feat_file_path(feat_out)
      return()
    }

    feat_data(feat_df)
    feat_file_path(feat_out)
    status_text("Conversion completed.")
  })

  output$tax_preview <- renderDT({
    req(tax_data())
    datatable(tax_data(), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$feat_preview <- renderDT({
    req(feat_data())
    datatable(feat_data(), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$download_tax <- downloadHandler(
    filename = function() "taxonomy.tsv",
    content = function(file) {
      req(tax_file_path())
      file.copy(tax_file_path(), file)
    }
  )

  output$download_feat <- downloadHandler(
    filename = function() "feature_table.tsv",
    content = function(file) {
      req(feat_file_path())
      file.copy(feat_file_path(), file)
    }
  )
}

shinyApp(ui, server)
