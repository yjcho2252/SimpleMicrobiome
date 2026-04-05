## UI
mod_fileload_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML(sprintf(
      "#%s .btn-file, #%s .btn-file, #%s .btn-file { font-size: 11px; }",
      ns("otu_file"), ns("tax_file"), ns("meta_file")
    ))),
    h4(icon("upload"), "Data Input"),
    fileInput(ns("otu_file"),
              "1. ASV/OTU Abundance Matrix",
              accept = c(".csv", ".tsv", ".txt")),
    fileInput(ns("tax_file"),
              "2. Taxonomy Table",
              accept = c(".csv", ".tsv", ".txt")),
    fileInput(ns("meta_file"),
              "3. Metadata File", 
              accept = c(".csv", ".tsv", ".txt")),
    actionButton(
      ns("reset_all_app"),
      "Reset All",
      icon = icon("rotate-left"),
      class = "btn btn-warning btn-sm",
      style = "font-size: 12px; padding: 3px 8px; width: 110px; display: inline-block; white-space: nowrap;"
    ),
    hr(),
    uiOutput(ns("load_status"))
  )
}

## Server
mod_fileload_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    load_completed <- reactiveVal(FALSE) 
    
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
      req(input$otu_file, input$tax_file, input$meta_file)
      
      otu_df <- read_microbiome_file(input$otu_file)
      tax_df <- read_microbiome_file(input$tax_file)
      meta_df <- read_microbiome_file(input$meta_file, is_meta = TRUE) 
      
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
    })
    
    meta_vars <- reactive({
      req(ps_obj_initial())
      colnames(sample_data(ps_obj_initial()))
    })
    
    output$load_status <- renderUI({
      if (load_completed()) {
        h5(span("✅ Go to Preprocessing Tab!", style = "color: green; font-weight: bold;"))
      } else {
        h5(span("Waiting for 3 data files...", style = "color: orange;"))
      }
    })
    
    return(list(
      ps_initial = ps_obj_initial,
      meta_vars = meta_vars,
      load_completed = load_completed
    ))
  })
}
