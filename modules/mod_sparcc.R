## UI
mod_sparcc_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        width = 2,
        h4(icon("diagram-project"), "SparCC Network"),
        hr(),
        selectInput(ns("tax_level"), "1. Taxonomic level", choices = c("ASV", "Genus", "Species"), selected = "Genus"),
        numericInput(ns("prevalence_filter_pct"), "2. Prevalence filter cutoff (%)", value = 10, min = 0, max = 100, step = 1),
        numericInput(ns("max_taxa"), "3. Max taxa for network", value = 200, min = 20, max = 2000, step = 10),
        numericInput(ns("seed"), "4. Seed", value = 1001, min = 1, step = 1),
        numericInput(ns("plot_width"), "Plot width (px)", value = 900, min = 400, max = 2400, step = 50),
        numericInput(ns("plot_height"), "Plot height (px)", value = 700, min = 300, max = 2400, step = 50),
        actionButton(ns("run_network_btn"), "Run SparCC", class = "btn-danger", style = "font-size: 12px;"),
        tags$script(HTML(
          "Shiny.addCustomMessageHandler('toggle-sparcc-run-btn', function(msg) {
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
          tabPanel("Summary", verbatimTextOutput(ns("network_summary")))
        )
      )
    )
  )
}

## Server
mod_sparcc_server <- function(id, ps_obj) {
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
      if (!phyloseq::taxa_are_rows(phyloseq::otu_table(ps_use))) otu_mat <- t(otu_mat)

      prevalence_cutoff <- as.numeric(input$prevalence_filter_pct)
      if (is.na(prevalence_cutoff)) prevalence_cutoff <- 10
      prevalence_cutoff <- max(0, min(100, prevalence_cutoff))
      keep_prev <- rowMeans(otu_mat > 0) * 100 >= prevalence_cutoff
      otu_mat <- otu_mat[keep_prev, , drop = FALSE]
      validate(need(nrow(otu_mat) >= 3, "Too few taxa remain after prevalence filtering. Lower the cutoff."))

      max_taxa <- as.integer(input$max_taxa)
      if (is.na(max_taxa) || max_taxa < 3) max_taxa <- 200
      if (nrow(otu_mat) > max_taxa) {
        otu_mat <- otu_mat[order(rowSums(otu_mat), decreasing = TRUE)[seq_len(max_taxa)], , drop = FALSE]
      }

      list(otu_mat = otu_mat)
    })

    network_result <- eventReactive(input$run_network_btn, {
      req(build_network_inputs())
      validate(
        need(requireNamespace("NetCoMi", quietly = TRUE), "NetCoMi package is not installed. Please install 'NetCoMi'."),
        need(requireNamespace("igraph", quietly = TRUE), "igraph package is not installed. Please install 'igraph'.")
      )

      session$sendCustomMessage("toggle-sparcc-run-btn", list(id = session$ns("run_network_btn"), disabled = TRUE, label = "Running..."))
      on.exit({
        session$sendCustomMessage("toggle-sparcc-run-btn", list(id = session$ns("run_network_btn"), disabled = FALSE, label = "Run SparCC"))
      }, add = TRUE)

      set.seed(as.integer(input$seed))
      otu_mat <- build_network_inputs()$otu_mat

      result <- tryCatch({
        withProgress(message = "Running NetCoMi SparCC...", value = 0, {
          net_obj <- NetCoMi::netConstruct(data = t(otu_mat), dataType = "counts", measure = "sparcc", verbose = 0)
          adja <- net_obj$adjaMat1
          if (is.null(adja)) adja <- net_obj$adjaMat
          validate(need(!is.null(adja), "NetCoMi returned no adjacency matrix for SparCC."))
          graph <- igraph::graph_from_adjacency_matrix(as.matrix(adja), mode = "undirected", weighted = TRUE, diag = FALSE)
          edge_df <- igraph::as_data_frame(graph, what = "edges")
          if (nrow(edge_df) > 0 && "weight" %in% colnames(edge_df)) {
            edge_df$weight <- as.numeric(edge_df$weight)
            edge_df$sign <- ifelse(edge_df$weight >= 0, "Positive", "Negative")
          } else {
            edge_df$weight <- NA_real_
            edge_df$sign <- NA_character_
          }
          list(graph = graph, edge_table = edge_df, n_taxa = nrow(otu_mat), n_samples = ncol(otu_mat))
        })
      }, error = function(e) {
        showNotification(paste("SparCC execution failed:", e$message), type = "error", duration = NULL)
        NULL
      })

      validate(need(!is.null(result), "Network analysis failed. Check the error message above."))
      result
    })

    output$network_summary <- renderText({
      req(network_result())
      g <- network_result()$graph
      paste(c(
        "SparCC Network Summary",
        paste0("Samples used: ", network_result()$n_samples),
        paste0("Taxa used: ", network_result()$n_taxa),
        paste0("Nodes: ", igraph::vcount(g)),
        paste0("Edges: ", igraph::ecount(g)),
        paste0("Connected components: ", igraph::components(g)$no)
      ), collapse = "\n")
    })

    network_plot_reactive <- reactive({
      req(network_result())
      g <- network_result()$graph
      validate(
        need(igraph::vcount(g) > 0, "No nodes available for plotting."),
        need(igraph::ecount(g) > 0, "No edges were inferred. Try changing parameters or filtering settings.")
      )
      edge_weight <- igraph::E(g)$weight
      edge_col <- ifelse(edge_weight >= 0, "#4575B4", "#D73027")
      edge_col[is.na(edge_col)] <- "#999999"
      edge_abs <- abs(edge_weight)
      edge_abs[!is.finite(edge_abs)] <- 0
      ew_rng <- range(edge_abs, na.rm = TRUE)
      if (!all(is.finite(ew_rng)) || diff(ew_rng) == 0) {
        edge_width <- rep(2.2, length(edge_abs))
      } else {
        edge_width <- 1 + 5 * (edge_abs - ew_rng[1]) / (ew_rng[2] - ew_rng[1])
      }
      plot(g, layout = igraph::layout_with_fr(g), vertex.size = 5, vertex.label.cex = 0.7, vertex.label.color = "#222222", vertex.color = "#F1C40F", edge.width = edge_width, edge.color = edge_col, main = "NetCoMi SparCC Network")
      graphics::legend(
        "topleft",
        legend = c("Positive", "Negative"),
        col = c("#4575B4", "#D73027"),
        lty = 1,
        lwd = 3,
        bty = "n",
        cex = 0.9,
        title = "Edge sign"
      )
      width_levels <- c(ew_rng[1], mean(ew_rng), ew_rng[2])
      if (!all(is.finite(width_levels))) width_levels <- c(0, 0, 0)
      if (!all(is.finite(ew_rng)) || diff(ew_rng) == 0) {
        width_legend <- rep(2.2, 3)
      } else {
        width_legend <- 1 + 5 * (width_levels - ew_rng[1]) / (ew_rng[2] - ew_rng[1])
      }
      graphics::legend(
        "bottomleft",
        legend = paste0(c("Low", "Mid", "High"), " |w|=", sprintf("%.2f", width_levels)),
        col = "#4D4D4D",
        lty = 1,
        lwd = width_legend,
        bty = "n",
        cex = 0.85,
        title = "Edge strength"
      )
    })

    output$network_plot <- renderPlot({ network_plot_reactive() }, height = function() input$plot_height, width = function() input$plot_width)
    output$edge_table <- renderDT({ datatable(network_result()$edge_table, options = list(scrollX = TRUE)) })

    output$download_network_plot <- downloadHandler(
      filename = function() paste0("sparcc_network_", Sys.Date(), ".png"),
      content = function(file) {
        grDevices::png(file, width = input$plot_width, height = input$plot_height, res = 120)
        network_plot_reactive()
        grDevices::dev.off()
      }
    )

    output$download_edge_table <- downloadHandler(
      filename = function() paste0("sparcc_edges_", Sys.Date(), ".tsv"),
      content = function(file) readr::write_tsv(network_result()$edge_table, file)
    )
  })
}
