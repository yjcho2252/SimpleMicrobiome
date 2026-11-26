mod_upset_ui <- function(id) {
  ns <- NS(id)
  plot_scroll_layout(
    sidebarPanel(
      width = 2,
      h4(icon("layer-group"), "Upset Plot"),
      hr(),
      selectInput(ns("group"), "Grouping Variable", choices = c("Choose a variable" = "")),
      selectizeInput(ns("levels"), "Groups to Include", choices = NULL, multiple = TRUE),
      hr(),
      h4(icon("sliders"), "Analysis & Plot Settings"),
      selectInput(ns("rank"), "Taxonomic Level", choices = c("Genus", "ASV"), selected = "Genus"),
      tags$details(
        tags$summary("Advanced Options"),
        numericInput(ns("min_prevalence"), "Minimum Group Prevalence (%)", 10, min = 0, max = 100, step = 1)
      ),
      hr(),
      h4(icon("sliders"), "Display Options"),
      selectInput(ns("plot_size_mode"), "Plot Size", choices = c("Auto" = "auto", "Manual" = "manual")),
      conditionalPanel(
        "input.plot_size_mode === 'manual'", ns = ns,
        numericInput(ns("plot_width"), "Plot Width (px)", 900, min = 400, max = 3000, step = 50),
        numericInput(ns("plot_height"), "Plot Height (px)", 550, min = 300, max = 2400, step = 50)
      ),
      numericInput(ns("base_size"), "Base Font Size", 11, min = 6, max = 30),
      numericInput(ns("max_intersections"), "Max Intersections to Plot", 25, min = 1, max = 100),
      hr(),
      h5(icon("download"), "Download"),
      div(
        style = "display: flex; flex-wrap: wrap; gap: 4px;",
        downloadButton(ns("download_plot"), "Plot (PNG)"),
        downloadButton(ns("download_table"), "Intersections (TSV)")
      )
    ),
    mainPanel(
      width = 9,
      h4("Upset Plot"),
      automatic_analysis_view(ns,
        tabsetPanel(
          tabPanel("Upset Plot",
            div(
              class = "taxa-plot-stage",
              uiOutput(ns("plot_container"), class = "taxa-plot-container", style = "overflow-x: auto;"),
              taxa_plot_loading("Upset Plot")
            ),
            plot_legend_output(ns, "figure_legend", c("plot"))
          ),
          tabPanel("Intersection Table", DT::DTOutput(ns("intersections")))
        ),
        uiOutput(ns("status_box"))
      )
    )
  )
}

mod_upset_server <- function(id, ps_obj, plot_order = function() NULL, active = function() TRUE) {
  moduleServer(id, function(input, output, session) {
    metadata <- reactive({
      ps <- ps_obj()
      req(!is.null(ps))
      as.data.frame(phyloseq::sample_data(ps))
    })

    observeEvent(ps_obj(), {
      columns <- names(metadata())
      columns <- columns[gsub("[^a-z0-9]", "", tolower(columns)) != "sampleid"]
      group <- isolate(input$group)
      default_group <- if (length(columns)) columns[1L] else ""
      updateSelectInput(session, "group", choices = c("Choose a variable" = "", stats::setNames(columns, columns)),
                        selected = if (isTruthy(group) && group %in% columns) group else default_group)
      taxonomy <- phyloseq::tax_table(ps_obj(), errorIfNULL = FALSE)
      ranks <- unique(c("ASV", colnames(taxonomy)))
      rank <- isolate(input$rank)
      default_rank <- if ("Genus" %in% ranks) "Genus" else "ASV"
      updateSelectInput(session, "rank", choices = ranks,
                        selected = if (isTruthy(rank) && rank %in% ranks) rank else default_rank)
    })

    group_values <- reactive({
      meta <- metadata()
      if (!isTruthy(input$group) || !input$group %in% names(meta)) return(character())
      values <- trimws(as.character(meta[[input$group]]))
      stats::setNames(values, rownames(meta))
    })

    observeEvent(group_values(), {
      values <- group_values()
      choices <- sort(unique(values[!is.na(values) & nzchar(values)]))
      updateSelectizeInput(session, "levels", choices = choices,
                          selected = if (length(choices) <= 12L) choices else character())
    }, ignoreNULL = FALSE)

    settings <- reactive({
      validate(need(isTruthy(input$group), "Choose a grouping variable."))
      groups <- unique(input$levels)
      validate(need(
        length(groups) >= 2L && length(groups) <= 12L,
        "Select between 2 and 12 groups."
      ))
      list(group = input$group, groups = groups, rank = req(input$rank), cutoff = req(input$min_prevalence))
    })

    analysis_inputs <- reactive({
      ps <- ps_obj()
      validate(need(!is.null(ps), "Load and preprocess data to view the Upset Plot."))
      list(ps = ps, config = settings())
    })

    displayed_result <- automatic_analysis(active, analysis_inputs, function(request) {
      ps <- request$ps
      config <- request$config
      ranked <- taxa_rank_counts(ps, config$rank)
      counts <- ranked$counts
      samples <- colnames(counts)
      values <- unname(group_values()[samples])
      totals <- unname(ranked$sample_totals[samples])
      keep <- is.finite(totals) & totals > 0 & values %in% config$groups
      sample_counts <- table(factor(values[keep], levels = config$groups))
      validate(need(all(sample_counts > 0L), "Each selected group needs at least one sample with positive read counts."))

      membership <- vapply(config$groups, function(group) {
        detected <- rowSums(counts[, keep & values == group, drop = FALSE] > 0)
        as.integer(detected > 0L & 100 * detected / as.integer(sample_counts[group]) >= config$cutoff)
      }, integer(nrow(counts)))
      membership <- matrix(membership, nrow = nrow(counts),
                           dimnames = list(ranked$labels, config$groups))
      membership <- membership[rowSums(membership) > 0L, , drop = FALSE]
      validate(need(nrow(membership) > 0L, "No taxa passed the prevalence threshold in the selected groups."))
      list(membership = membership, sample_counts = sample_counts, config = config)
    })

    combinations <- reactive({
      data <- displayed_result()
      groups <- colnames(data$membership)
      ordered <- levels(plot_order_factor(groups, data$config$group, plot_order(), groups))
      membership <- data$membership[, ordered, drop = FALSE]
      ComplexHeatmap::make_comb_mat(membership, mode = "distinct")
    })

    plot_combinations <- reactive({
      combinations <- combinations()
      sizes <- ComplexHeatmap::comb_size(combinations)
      limit <- if (isTruthy(input$max_intersections)) input$max_intersections else 25L
      shown <- head(order(sizes, decreasing = TRUE), max(1L, min(100L, as.integer(limit))))
      combinations[shown]
    })

    intersection_table <- reactive({
      combinations <- combinations()
      sizes <- ComplexHeatmap::comb_size(combinations)
      codes <- ComplexHeatmap::comb_name(combinations)
      groups <- ComplexHeatmap::set_name(combinations)
      cases <- lapply(order(sizes, decreasing = TRUE), function(i) {
        included <- strsplit(codes[i], "", fixed = TRUE)[[1L]] == "1"
        taxa <- sort(ComplexHeatmap::extract_comb(combinations, codes[i]))
        data.frame(
          Taxon = taxa, Groups = paste(groups[included], collapse = " + "),
          GroupCount = sum(included), IntersectionSize = sizes[i]
        )
      })
      do.call(rbind, cases)
    })

    dimensions <- reactive({
      if (identical(input$plot_size_mode, "manual")) {
        return(c(width = req(input$plot_width), height = req(input$plot_height)))
      }
      m <- plot_combinations()
      c(width = max(650, 320 + 22 * length(ComplexHeatmap::comb_size(m))),
        height = max(430, 270 + 24 * length(ComplexHeatmap::set_size(m))))
    })

    draw_upset <- function() {
      m <- plot_combinations()
      font <- if (isTruthy(input$base_size)) input$base_size else 11
      sizes <- ComplexHeatmap::comb_size(m)
      size <- dimensions()
      show_counts <- length(sizes) <= 12L && size[["width"]] >= 260 + 22 * length(sizes)
      top_height <- grid::unit(max(100, size[["height"]] * 0.35) / 96, "inches")
      intersection_color <- "#246A73"
      annotation_font <- grid::gpar(fontfamily = "sans", fontsize = font, col = "#334155")
      plot <- ComplexHeatmap::UpSet(
        m, set_order = seq_along(ComplexHeatmap::set_size(m)),
        comb_order = seq_along(sizes),
        row_names_gp = grid::gpar(fontfamily = "sans", fontsize = font, fontface = "bold", col = "#243446"),
        comb_col = intersection_color, pt_size = grid::unit(2.5, "mm"), lwd = 1.7,
        bg_col = c("#F3F6F8", "#FFFFFF"), bg_pt_col = "#CFD8DF",
        top_annotation = ComplexHeatmap::upset_top_annotation(
          m, height = top_height, gp = grid::gpar(fill = intersection_color), annotation_name_gp = annotation_font,
          add_numbers = show_counts, numbers_rot = 0,
          numbers_gp = grid::gpar(fontfamily = "sans", fontsize = max(7, font - 2), col = "#243446"),
          numbers_offset = grid::unit(1.5, "mm"),
          ylim = if (show_counts) c(0, max(sizes) * 1.2) else NULL
        ),
        left_annotation = ComplexHeatmap::upset_left_annotation(
          m, gp = grid::gpar(fill = "#A8B9C6"), annotation_name_gp = annotation_font
        ),
        right_annotation = NULL
      )
      ComplexHeatmap::draw(plot)
    }

    output$plot_container <- renderUI({
      size <- dimensions()
      plotOutput(session$ns("plot"), width = size[["width"]], height = size[["height"]])
    })
    register_automatic_analysis_status(output, active, displayed_result,
                                       "Choose data and at least two groups to view the Upset Plot.",
                                       updating_message = "Preparing the Upset Plot. Please wait...")
    output$plot <- renderPlot(draw_upset(), res = 96)
    output$figure_legend <- renderUI({
      req_plot_legend_tab(session, "Upset Plot")
      config <- displayed_result()$config
      shown <- length(ComplexHeatmap::comb_size(plot_combinations()))
      total <- length(ComplexHeatmap::comb_size(combinations()))
      plot_legend_box(
        dimensions()[["width"]], "Upset Plot",
        paste(
          paste0("Each row is a level of ", config$group, ". A taxon at the ", config$rank,
                 " level belongs to a group if it was detected in at least one sample and reached ",
                 format(config$cutoff), "% prevalence in that group."),
          "Side bars show the number of included taxa per group. Each top bar counts taxa that belong to exactly the groups marked by filled dots among the selected groups.",
          if (shown < total) {
            paste0("The plot shows ", shown, " of ", total,
                   " intersections, ordered by size. The Intersection Table lists all intersections.")
          } else {
            "The plot and Intersection Table include all observed intersections."
          }
        )
      )
    })
    outputOptions(output, "figure_legend", suspendWhenHidden = FALSE, priority = -10)
    output$intersections <- DT::renderDT({
      DT::datatable(intersection_table(), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
    })
    output$status_box <- renderUI({
      displayed_result()
      tags$div(
        style = sprintf("width: %dpx; max-width: 100%%;", as.integer(dimensions()[["width"]])),
        hr(),
        h5(icon("circle-info"), "Upset Plot Status"),
        tags$div(class = "app-result-card app-status-card", style = "width: 100%;",
                 verbatimTextOutput(session$ns("status")))
      )
    })
    output$status <- renderText({
      data <- displayed_result()
      sample_counts <- data$sample_counts
      paste(
        paste("Grouping variable:", data$config$group),
        paste("Taxonomic level:", data$config$rank),
        paste0("Minimum group prevalence: ", format(data$config$cutoff), "%"),
        paste("Taxa included:", nrow(data$membership)),
        paste("Samples per group:", paste(paste0(names(sample_counts), " (", as.integer(sample_counts), ")"), collapse = ", ")),
        paste("Intersections plotted:", length(ComplexHeatmap::comb_size(plot_combinations())), "of",
              length(ComplexHeatmap::comb_size(combinations()))),
        sep = "\n"
      )
    })

    output$download_plot <- downloadHandler(
      filename = function() "upset_plot.png",
      content = function(file) {
        size <- dimensions()
        grDevices::png(file, width = size[["width"]] / 96, height = size[["height"]] / 96,
                       units = "in", res = 300, bg = "white")
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_upset()
      }
    )
    output$download_table <- downloadHandler(
      filename = function() "upset_plot_intersections.tsv",
      content = function(file) readr::write_tsv(intersection_table(), file)
    )
  })
}
