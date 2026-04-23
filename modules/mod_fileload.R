## UI
mod_fileload_ui <- function(id) {
  ns <- NS(id)
  converter_url <- Sys.getenv("CONVERTER_APP_URL", "https://simplemicrobiome.mglab.org/convert/ui")
  tagList(
    tags$style(HTML("
      .file-input-row {
        display: flex;
        align-items: flex-start;
        gap: 10px;
        margin-bottom: 6px;
        justify-content: flex-start;
      }
      .file-input-label {
        flex: 0 0 200px;
        width: 200px;
        height: 38px;
        display: flex;
        align-items: center;
        margin: 0;
        line-height: 1.25;
        font-size: 12px;
        transform: translateY(-3px);
      }
      .file-input-control {
        flex: 0 0 360px;
        width: 360px;
      }
      .file-input-control .shiny-input-container,
      .file-input-control .form-group,
      .file-input-control .input-group {
        margin-bottom: 0 !important;
      }
      .file-input-control .shiny-file-input-progress {
        display: none !important;
      }
      .shiny-file-input-progress {
        display: none !important;
      }
      .file-input-control .btn-file {
        font-size: 11px;
      }
    ")),
    h4(icon("upload"), "Data Input"),
    div(
      class = "file-input-row",
      tags$p("1. ASV/OTU Abundance Matrix", class = "file-input-label"),
      div(
        class = "file-input-control",
        fileInput(ns("otu_file"), NULL, accept = c(".csv", ".tsv", ".txt"))
      )
    ),
    div(
      class = "file-input-row",
      tags$p("2. Taxonomy Table", class = "file-input-label"),
      div(
        class = "file-input-control",
        fileInput(ns("tax_file"), NULL, accept = c(".csv", ".tsv", ".txt"))
      )
    ),
    div(
      class = "file-input-row",
      tags$p("3. Metadata File", class = "file-input-label"),
      div(
        class = "file-input-control",
        fileInput(ns("meta_file"), NULL, accept = c(".csv", ".tsv", ".txt"))
      )
    ),
    div(
      style = "display: flex; align-items: center; gap: 6px; flex-wrap: nowrap; overflow-x: auto; margin-bottom: 4px;",
      actionButton(
        ns("load_example"),
        "Load Example",
        icon = icon("flask"),
        class = "btn btn-info btn-sm",
        style = "font-size: 12px; padding: 3px 8px; width: 110px; white-space: nowrap; flex: 0 0 auto;"
      ),
      downloadButton(
        ns("download_example"),
        "Download Example",
        class = "btn btn-secondary btn-sm",
        style = "font-size: 12px; padding: 3px 8px; width: 130px; white-space: nowrap; flex: 0 0 auto;"
      ),
      actionButton(
        ns("reset_all_app"),
        "Reset All",
        icon = icon("rotate-left"),
        class = "btn btn-warning btn-sm",
        style = "font-size: 12px; padding: 3px 8px; width: 110px; white-space: nowrap; flex: 0 0 auto;"
      )
    ),
    tags$div(
      style = "font-size: 12px; color: #555; margin-top: 2px; margin-bottom: 8px;",
      tags$p("This tool requires 3 files:", style = "margin: 0 0 4px 0;"),
      tags$ul(
        style = "margin: 0 0 0 18px; padding: 0;",
        tags$li("ASV/OTU Abundance Matrix from table.qza of QIIME2"),
        tags$li("Taxonomy table from taxonomy.qza of QIIME2"),
        tags$li("Metadata table prepared manually by the user (see the example file)")
      ),
      tags$p(
        tagList(
          "To convert QZA files to input files, use ",
          tags$a(
            href = converter_url,
            target = "_blank",
            rel = "noopener noreferrer",
            "QZA Converter"
          ),
          ""
        ),
        style = "margin: 6px 0 0 0;"
      )
    ),
    uiOutput(ns("load_status"))
  )
}

## Server
mod_fileload_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    load_completed <- reactiveVal(FALSE) 
    example_files <- reactiveVal(NULL)

    resolve_sample_dir <- function() {
      file.path(getwd(), "sample")
    }

    observeEvent(input$load_example, {
      sample_dir <- resolve_sample_dir()
      if (!dir.exists(sample_dir)) {
        showNotification(
          "Example load failed: sample folder was not found.",
          type = "error",
          duration = 8
        )
        return(NULL)
      }
      otu_path <- file.path(sample_dir, "1_ASV_table.txt")
      tax_path <- file.path(sample_dir, "2_taxonomy_table.txt")
      meta_path <- file.path(sample_dir, "3_metadata.txt")

      required_paths <- c(otu_path, tax_path, meta_path)
      if (!all(file.exists(required_paths))) {
        showNotification(
          "Example load failed: required files are missing.",
          type = "error",
          duration = 8
        )
        return(NULL)
      }

      example_files(list(
        otu_file = list(datapath = otu_path, name = basename(otu_path)),
        tax_file = list(datapath = tax_path, name = basename(tax_path)),
        meta_file = list(datapath = meta_path, name = basename(meta_path))
      ))

      showNotification(
        "Example files loaded!",
        type = "message",
        duration = 3
      )
    })

    output$download_example <- downloadHandler(
      filename = function() {
        "example_data.zip"
      },
      content = function(file) {
        sample_dir <- resolve_sample_dir()
        if (!dir.exists(sample_dir)) {
          stop("Example download failed: sample folder was not found.")
        }
        required_files <- c(
          "1_ASV_table.txt",
          "2_taxonomy_table.txt",
          "3_metadata.txt"
        )
        required_paths <- file.path(sample_dir, required_files)

        if (!all(file.exists(required_paths))) {
          stop("Example download failed: required files are missing in the sample folder.")
        }

        temp_export_dir <- file.path(tempdir(), paste0("example_data_", as.integer(Sys.time())))
        dir.create(temp_export_dir, recursive = TRUE, showWarnings = FALSE)
        file.copy(required_paths, temp_export_dir, overwrite = TRUE)

        old_wd <- getwd()
        on.exit(setwd(old_wd), add = TRUE)
        setwd(temp_export_dir)

        utils::zip(
          zipfile = file,
          files = required_files,
          flags = "-r9Xq"
        )
      }
    )
    
    read_microbiome_file <- function(file_input, is_meta = FALSE) {
      if (is.null(file_input)) return(NULL)
      
      file_path <- file_input$datapath
      data_content <- readLines(file_path, n = 5)
      
      if (length(data_content) < 2) {
        return(NULL)  
      }
      
      if (grepl("\t", data_content[2])) {
        sep_val <- "\t"  
      } else {
        sep_val <- ","  
      }
      
      rn_val <- if (is_meta) NULL else 1
      
      data <- tryCatch({
        read.table(file_path, header = TRUE, row.names = rn_val, sep = sep_val, 
                   comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
      }, error = function(e) {
        message(paste("File read failed:", file_input$name, "Error:", conditionMessage(e)))
        return(NULL)
      })
      return(data)
    }
    
    ps_obj_initial <- reactive({ 
      selected_files <- example_files()

      otu_input <- if (!is.null(input$otu_file)) input$otu_file else if (!is.null(selected_files)) selected_files$otu_file else NULL
      tax_input <- if (!is.null(input$tax_file)) input$tax_file else if (!is.null(selected_files)) selected_files$tax_file else NULL
      meta_input <- if (!is.null(input$meta_file)) input$meta_file else if (!is.null(selected_files)) selected_files$meta_file else NULL

      req(otu_input, tax_input, meta_input)
      
      otu_df <- read_microbiome_file(otu_input)
      tax_df <- read_microbiome_file(tax_input)
      meta_df <- read_microbiome_file(meta_input, is_meta = TRUE) 
      
      if(is.null(otu_df) | is.null(tax_df) | is.null(meta_df)) {
        load_completed(FALSE)
        
        showNotification("File read error: Check file format (csv/tsv), delimiter, and data structure. Please re-upload all files.",
                         type = "error", duration = 10)
        return(NULL)
      }
      
      ps <- tryCatch({
        colnames(meta_df)[1] <- "SampleID"
        rownames(meta_df) <- meta_df$SampleID
        
        otu <- otu_table(as.matrix(otu_df), taxa_are_rows = TRUE) 
        tax <- tax_table(as.matrix(tax_df))
        meta <- sample_data(meta_df) 
        
        ps <- phyloseq(otu, tax, meta)
        
        sample_keep <- sample_names(otu) %in% sample_names(meta)
        if (length(sample_keep) != nsamples(ps)) {
          stop(sprintf(
            paste0(
              "Length mismatch before prune_samples: ",
              "length(sample_names(otu) %%in%% sample_names(meta)) = %d, ",
              "nsamples(ps) = %d, ",
              "length(sample_names(otu)) = %d, ",
              "length(sample_names(meta)) = %d"
            ),
            length(sample_keep),
            nsamples(ps),
            length(sample_names(otu)),
            length(sample_names(meta))
          ))
        }
        ps <- prune_samples(sample_keep, ps)

        taxa_keep <- taxa_names(otu) %in% taxa_names(tax)
        if (length(taxa_keep) != ntaxa(ps)) {
          stop(sprintf(
            paste0(
              "Length mismatch before prune_taxa: ",
              "length(taxa_names(otu) %%in%% taxa_names(tax)) = %d, ",
              "ntaxa(ps) = %d, ",
              "length(taxa_names(otu)) = %d, ",
              "length(taxa_names(tax)) = %d"
            ),
            length(taxa_keep),
            ntaxa(ps),
            length(taxa_names(otu)),
            length(taxa_names(tax))
          ))
        }
        ps <- prune_taxa(taxa_keep, ps)
        ps <- prune_taxa(taxa_sums(ps) > 0, ps)
        ps <- prune_samples(sample_sums(ps) > 0, ps)
        
        if (ntaxa(ps) == 0) {
          stop("No valid taxa found. Check ASV/Taxonomy ID matching. (0 taxa)")
        }
        if (nsamples(ps) == 0) {
          stop("No valid samples found. Check Sample ID matching. (0 samples)")
        }
        
        load_completed(TRUE) 
        
        return(ps)
        
      }, error = function(e) {
        load_completed(FALSE)
        showNotification(paste0("Data construction error: ", conditionMessage(e),
                                "\nPlease verify file contents and try again."),
                         type = "error", duration = 10)
        return(NULL)
      })
    }) |> bindEvent(
      input$otu_file,
      input$tax_file,
      input$meta_file,
      example_files(),
      ignoreInit = TRUE
    )
    
    meta_vars <- reactive({
      req(ps_obj_initial())
      colnames(sample_data(ps_obj_initial()))
    })
    
    output$load_status <- renderUI({
      if (load_completed()) {
        h5(
          span("✅ Go to Preprocessing Tab!", style = "color: green; font-weight: bold;"),
          style = "margin-top: 2px; margin-bottom: 0;"
        )
      } else {
        h5(
          span("Waiting for 3 data files...", style = "color: orange;"),
          style = "margin-top: 2px; margin-bottom: 0;"
        )
      }
    })
    
    return(list(
      ps_initial = ps_obj_initial,
      meta_vars = meta_vars,
      load_completed = load_completed
    ))
  })
}
