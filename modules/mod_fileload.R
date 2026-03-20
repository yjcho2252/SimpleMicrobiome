# --- UI 함수 ---
mod_fileload_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Data Input (Raw Counts Required)"),
    fileInput(ns("otu_file"),
              "1. ASV/OTU Abundance Matrix (.csv/.tsv)",
              accept = c(".csv", ".tsv", ".txt")),
    fileInput(ns("tax_file"),
              "2. Taxonomy Table (.csv/.tsv)",
              accept = c(".csv", ".tsv", ".txt")),
    fileInput(ns("meta_file"),
              "3. Metadata File (.csv/.tsv)", 
              accept = c(".csv", ".tsv", ".txt")),
    hr(),
    uiOutput(ns("load_status"))
  )
}

# --- Server 함수 ---
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
        
        showNotification("🚨 파일 읽기 오류: 파일 형식(csv/tsv), 구분자 또는 데이터 구조를 확인하세요. 모든 파일을 다시 업로드해주세요.",
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
        
        ps <- prune_samples(sample_names(otu) %in% sample_names(meta), ps)
        ps <- prune_taxa(taxa_names(otu) %in% taxa_names(tax), ps)
        ps <- prune_taxa(taxa_sums(ps) > 0, ps)
        ps <- prune_samples(sample_sums(ps) > 0, ps)
        
        if (ntaxa(ps) == 0) {
          stop("Taxa가 유효하지 않습니다. ASV/Taxonomy ID 일치를 확인하세요. (0 taxa)")
        }
        if (nsamples(ps) == 0) {
          stop("샘플이 유효하지 않습니다. Sample ID 일치를 확인하세요. (0 samples)")
        }
        
        load_completed(TRUE) 
        
        return(ps)
        
      }, error = function(e) {
        load_completed(FALSE)
        showNotification(paste0("❌ 데이터 생성 오류: ", conditionMessage(e),
                                "\n파일 내용을 확인하고 다시 시도하세요."),
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
        h5(span("✅ Go to Data Preprocessing Tab!", style = "color: green; font-weight: bold;"))
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