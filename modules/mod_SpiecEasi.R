library(shiny)
library(phyloseq)
library(ggplot2)
library(DT)

## UI
mod_spieceasi_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4(icon("diagram-project"), "SpiecEasi Network"),
        hr(),
        selectInput(ns("tax_level"), "1. Taxonomic level", choices = c("ASV", "Genus", "Species"), selected = "Genus"),
        numericInput(ns("prevalence_filter_pct"), "2. Prevalence filter cutoff (%)", value = 10, min = 0, max = 100, step = 1),
        numericInput(ns("max_taxa"), "3. Max taxa for network", value = 200, min = 20, max = 2000, step = 10),
        hr(),
        selectInput(ns("method"), "4. Inference method", choices = c("mb", "glasso"), selected = "mb"),
        numericInput(ns("nlambda"), "5. Number of lambdas", value = 20, min = 5, max = 100, step = 1),
        numericInput(ns("lambda_min_ratio"), "6. Lambda min ratio", value = 1e-2, min = 1e-4, max = 0.5, step = 1e-3),
        numericInput(ns("pulsar_rep_num"), "7. Pulsar rep.num", value = 20, min = 5, max = 200, step = 1),
        numericInput(ns("stars_thresh"), "8. StARS threshold", value = 0.05, min = 0.01, max = 0.2, step = 0.01),
        numericInput(ns("seed"), "9. Seed", value = 1001, min = 1, step = 1),
        numericInput(ns("plot_width"), "Plot width (px)", value = 900, min = 400, max = 2400, step = 50),
        numericInput(ns("plot_height"), "Plot height (px)", value = 700, min = 300, max = 2400, step = 50),
        actionButton(ns("run_network_btn"), "Run Network Analysis", class = "btn-danger", style = "font-size: 12px;"),
        tags$script(HTML(
          "Shiny.addCustomMessageHandler('toggle-network-run-btn', function(msg) {
             var btn = document.getElementById(msg.id);
             if (!btn) return;
             btn.disabled = !!msg.disabled;
             if (msg.label) btn.textContent = msg.label;
           });"
        ))
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Network Plot",
            downloadButton(ns("download_network_plot"), "Download Plot (PNG)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
            plotOutput(ns("network_plot"), height = "auto")
          ),
          tabPanel(
            "Edge Table",
            downloadButton(ns("download_edge_table"), "Download Table (TSV)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
            DTOutput(ns("edge_table"))
          ),
          tabPanel(
            "Summary",
            verbatimTextOutput(ns("network_summary"))
          )
        )
      )
    )
  )
}

## Server
mod_spieceasi_server <- function(id, ps_obj) {
  moduleServer(id, function(input, output, session) {

    build_network_inputs <- reactive({
      req(ps_obj())
      ps <- ps_obj()
      validate(
        need(phyloseq::nsamples(ps) > 2, "At least 3 samples are required for network analysis."),
        need(phyloseq::ntaxa(ps) > 2, "At least 3 taxa are required for network analysis.")
      )

      ps_use <- ps

      if (input$tax_level != "ASV") {
        validate(
          need(!is.null(phyloseq::tax_table(ps_use)), "Taxonomy table is required to aggregate by taxonomic level.")
        )
        tax_cols <- colnames(phyloseq::tax_table(ps_use))
        validate(
          need(input$tax_level %in% tax_cols, paste("Taxonomic rank", input$tax_level, "not found in taxonomy table."))
        )

        ps_use <- phyloseq::tax_glom(ps_use, taxrank = input$tax_level, NArm = TRUE)

        tax_df <- as.data.frame(phyloseq::tax_table(ps_use))
        rank_values <- as.character(tax_df[[input$tax_level]])
        rank_values[is.na(rank_values) | rank_values == ""] <- phyloseq::taxa_names(ps_use)[is.na(rank_values) | rank_values == ""]
        phyloseq::taxa_names(ps_use) <- make.unique(rank_values)
      }

      otu_mat <- as.matrix(phyloseq::otu_table(ps_use))
      if (!phyloseq::taxa_are_rows(phyloseq::otu_table(ps_use))) {
        otu_mat <- t(otu_mat)
      }

      prevalence_cutoff <- as.numeric(input$prevalence_filter_pct)
      if (is.na(prevalence_cutoff)) prevalence_cutoff <- 10
      prevalence_cutoff <- max(0, min(100, prevalence_cutoff))
      taxa_prevalence <- rowMeans(otu_mat > 0) * 100
      keep_prev <- taxa_prevalence >= prevalence_cutoff
      otu_mat <- otu_mat[keep_prev, , drop = FALSE]

      validate(
        need(nrow(otu_mat) >= 3, "Too few taxa remain after prevalence filtering. Lower the cutoff.")
      )

      max_taxa <- as.integer(input$max_taxa)
      if (is.na(max_taxa) || max_taxa < 3) max_taxa <- 200
      if (nrow(otu_mat) > max_taxa) {
        taxa_order <- order(rowSums(otu_mat), decreasing = TRUE)
        otu_mat <- otu_mat[taxa_order[seq_len(max_taxa)], , drop = FALSE]
      }

      list(otu_mat = otu_mat)
    })

    network_result <- eventReactive(input$run_network_btn, {
      req(build_network_inputs())

      validate(
        need(requireNamespace("SpiecEasi", quietly = TRUE), "SpiecEasi package is not installed. Please install 'SpiecEasi'."),
        need(requireNamespace("igraph", quietly = TRUE), "igraph package is not installed. Please install 'igraph'.")
      )

      session$sendCustomMessage("toggle-network-run-btn", list(
        id = session$ns("run_network_btn"),
        disabled = TRUE,
        label = "Running..."
      ))
      on.exit({
        session$sendCustomMessage("toggle-network-run-btn", list(
          id = session$ns("run_network_btn"),
          disabled = FALSE,
          label = "Run Network Analysis"
        ))
      }, add = TRUE)

      set.seed(as.integer(input$seed))
      otu_mat <- build_network_inputs()$otu_mat

      result <- tryCatch({
        withProgress(message = "Running SpiecEasi...", value = 0, {
          fit <- SpiecEasi::spiec.easi(
            otu_mat,
            method = input$method,
            nlambda = as.integer(input$nlambda),
            lambda.min.ratio = as.numeric(input$lambda_min_ratio),
            pulsar.params = list(
              rep.num = as.integer(input$pulsar_rep_num),
              thresh = as.numeric(input$stars_thresh)
            ),
            verbose = FALSE
          )

          se_adj <- SpiecEasi::getRefit(fit)
          graph <- SpiecEasi::adj2igraph(se_adj)
          edge_df <- igraph::as_data_frame(graph, what = "edges")

          if (nrow(edge_df) > 0 && "weight" %in% colnames(edge_df)) {
            edge_df$weight <- as.numeric(edge_df$weight)
            edge_df$sign <- ifelse(edge_df$weight >= 0, "Positive", "Negative")
          } else {
            edge_df$weight <- NA_real_
            edge_df$sign <- NA_character_
          }

          list(
            graph = graph,
            edge_table = edge_df,
            n_taxa = nrow(otu_mat),
            n_samples = ncol(otu_mat)
          )
        })
      }, error = function(e) {
        showNotification(paste("SpiecEasi execution failed:", e$message), type = "error", duration = NULL)
        NULL
      })

      validate(
        need(!is.null(result), "Network analysis failed. Check the error message above.")
      )

      result
    })

    output$network_summary <- renderText({
      req(network_result())
      res <- network_result()
      g <- res$graph

      paste(
        c(
          "SpiecEasi Network Summary",
          paste0("Samples used: ", res$n_samples),
          paste0("Taxa used: ", res$n_taxa),
          paste0("Nodes: ", igraph::vcount(g)),
          paste0("Edges: ", igraph::ecount(g)),
          paste0("Connected components: ", igraph::components(g)$no),
          paste0("Method: ", input$method)
        ),
        collapse = "\n"
      )
    })

    network_plot_reactive <- reactive({
      req(network_result())
      g <- network_result()$graph

      validate(
        need(igraph::vcount(g) > 0, "No nodes available for plotting."),
        need(igraph::ecount(g) > 0, "No edges were inferred. Try changing parameters or filtering settings.")
      )

      layout_mat <- igraph::layout_with_fr(g)
      edge_weight <- igraph::E(g)$weight
      edge_col <- ifelse(edge_weight >= 0, "#D73027", "#4575B4")
      edge_col[is.na(edge_col)] <- "#999999"

      plot(
        g,
        layout = layout_mat,
        vertex.size = 5,
        vertex.label.cex = 0.7,
        vertex.label.color = "#222222",
        vertex.color = "#F1C40F",
        edge.width = 1.2,
        edge.color = edge_col,
        main = "SpiecEasi Inferred Association Network"
      )
    })

    output$network_plot <- renderPlot(
      { network_plot_reactive() },
      height = function() { req(input$plot_height); input$plot_height },
      width = function() { req(input$plot_width); input$plot_width }
    )

    output$edge_table <- renderDT({
      req(network_result())
      datatable(network_result()$edge_table, options = list(scrollX = TRUE))
    })

    output$download_network_plot <- downloadHandler(
      filename = function() {
        paste0("spieceasi_network_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(input$plot_height, input$plot_width)
        grDevices::png(file, width = input$plot_width, height = input$plot_height, res = 120)
        network_plot_reactive()
        grDevices::dev.off()
      }
    )

    output$download_edge_table <- downloadHandler(
      filename = function() {
        paste0("spieceasi_edges_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        req(network_result())
        readr::write_tsv(network_result()$edge_table, file)
      }
    )
  })
}
