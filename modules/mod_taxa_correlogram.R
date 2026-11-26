mod_taxa_correlogram_ui <- function(id) {
  ns <- NS(id)
  plot_scroll_layout(
    sidebarPanel(
      width = 2,
      h4(icon("braille"), "Taxa Correlogram"),
      hr(),
      selectInput(ns("group_var"), "Primary Grouping Variable", choices = c("None" = "")),
      conditionalPanel(
        "input.group_var && input.group_var !== ''", ns = ns,
        selectInput(ns("group_level"), "Group Level to Include", choices = c("All samples" = ""))
      ),
      hr(),
      h4(icon("sliders"), "Analysis & Plot Settings"),
      selectInput(ns("rank"), "Taxonomic Level", choices = c("Genus", "ASV"), selected = "Genus"),
      selectInput(ns("mode"), "Association Mode", choices = c(
        "Abundance-based" = "abundance", "Prevalence-based (binary)" = "prevalence"
      )),
      tags$details(
        tags$summary("Advanced Options"),
        conditionalPanel(
          "input.mode === 'abundance'", ns = ns,
          selectInput(ns("transform"), "Data Transform", choices = c("CLR" = "clr", "TSS" = "tss")),
          selectInput(ns("method"), "Correlation Method", choices = c("Pearson" = "pearson", "Spearman" = "spearman"))
        ),
        numericInput(ns("min_prevalence"), "Prevalence Filter (%)", 10, min = 0, max = 100),
        numericInput(ns("max_taxa"), "Max Taxa", 20, min = 2, max = 100),
        checkboxInput(ns("mask"), "Mask Non-Significant Pairs (FDR)", FALSE),
        numericInput(ns("q_cutoff"), "FDR Cutoff (q)", 0.05, min = 0.0001, max = 1, step = 0.01),
        selectizeInput(
          ns("selected_taxa"), "Selected Taxa (Optional)",
          choices = NULL, multiple = TRUE,
          options = list(placeholder = "If selected, only these taxa are shown")
        )
      ),
      hr(),
      h4(icon("sliders"), "Display Options"),
      selectInput(ns("plot_size_mode"), "Plot Size", choices = c("Auto" = "auto", "Manual" = "manual")),
      conditionalPanel(
        "input.plot_size_mode === 'manual'", ns = ns,
        numericInput(ns("plot_width"), "Plot Width (px)", 850, min = 400, max = 3000, step = 50),
        numericInput(ns("plot_height"), "Plot Height (px)", 750, min = 400, max = 3000, step = 50)
      ),
      numericInput(ns("base_size"), "Base Font Size", 11, min = 6, max = 30),
      selectInput(ns("x_label_angle"), "X-Axis Label Angle", choices = c(
        "0° Horizontal" = "0", "45° Diagonal" = "45", "90° Vertical" = "90"
      ), selected = "45"),
      selectInput(ns("color_palette"), "Color Palette", diverging_plot_palettes, selected = "blue_red"),
      checkboxInput(ns("cluster"), "Cluster Taxa", TRUE),
      checkboxInput(ns("show_significance_marks"), "Show Significance Marks (*)", FALSE),
      numericInput(ns("max_circle_size"), "Maximum Circle Size", 8, min = 1, max = 20, step = 0.5),
      hr(),
      h5(icon("download"), "Download"),
      div(
        style = "display: flex; flex-wrap: wrap; gap: 4px;",
        downloadButton(ns("download_plot"), "Plot (PNG)"),
        downloadButton(ns("download_matrix"), "Matrix (TSV)")
      )
    ),
    mainPanel(
      width = 9,
      h4("Taxa Correlogram"),
      automatic_analysis_view(ns,
        div(
          class = "taxa-plot-stage",
          uiOutput(ns("plot_container"), class = "taxa-plot-container", style = "overflow-x: auto;"),
          taxa_plot_loading("correlogram")
        ),
        plot_legend_output(ns, "figure_legend", c("plot")),
        uiOutput(ns("status_box"))
      )
    )
  )
}

mod_taxa_correlogram_server <- function(id, ps_obj, active = function() TRUE) {
  moduleServer(id, function(input, output, session) {
    metadata <- reactive({
      ps <- ps_obj()
      req(ps)
      as.data.frame(phyloseq::sample_data(ps, errorIfNULL = FALSE))
    })

    observeEvent(ps_obj(), {
      taxonomy <- phyloseq::tax_table(ps_obj(), errorIfNULL = FALSE)
      ranks <- unique(c("ASV", colnames(taxonomy)))
      selected <- isolate(input$rank)
      default_rank <- if ("Genus" %in% ranks) "Genus" else "ASV"
      updateSelectInput(session, "rank", choices = ranks,
                        selected = if (isTruthy(selected) && selected %in% ranks) selected else default_rank)
      columns <- names(metadata())
      columns <- columns[gsub("[^a-z0-9]", "", tolower(columns)) != "sampleid"]
      group <- isolate(input$group_var)
      default_group <- if (length(columns)) columns[1L] else ""
      updateSelectInput(session, "group_var", choices = c(stats::setNames(columns, columns), "None" = ""),
                        selected = if (isTruthy(group) && group %in% columns) group else default_group)
    })

    group_values <- reactive({
      meta <- metadata()
      if (!isTruthy(input$group_var) || !input$group_var %in% names(meta)) return(character())
      stats::setNames(trimws(as.character(meta[[input$group_var]])), rownames(meta))
    })

    observeEvent(group_values(), {
      values <- group_values()
      choices <- sort(unique(values[!is.na(values) & nzchar(values)]))
      default_level <- if (length(choices)) choices[1L] else ""
      updateSelectInput(session, "group_level",
                        choices = c(stats::setNames(choices, choices), "All samples" = ""),
                        selected = default_level)
    }, ignoreNULL = FALSE)

    settings <- reactive({
      req(input$rank, input$mode)
      binary <- identical(input$mode, "prevalence")
      group_level <- if (isTruthy(input$group_var) && isTruthy(input$group_level)) input$group_level else NULL
      list(
        rank = input$rank, binary = binary,
        group_var = if (!is.null(group_level)) input$group_var else NULL,
        group_level = group_level,
        transform = if (binary) NULL else req(input$transform),
        method = if (binary) "pearson" else req(input$method),
        min_prevalence = req(input$min_prevalence),
        max_taxa = req(input$max_taxa)
      )
    })

    analysis_inputs <- reactive({
      ps <- ps_obj()
      validate(need(!is.null(ps), "Load and preprocess data to view the correlogram."))
      list(ps = ps, config = settings())
    })

    displayed_result <- automatic_analysis(active, analysis_inputs, function(request) {
      ps <- request$ps
      config <- request$config
      ranked <- taxa_rank_counts(ps, config$rank)
      counts <- ranked$counts
      totals <- unname(ranked$sample_totals[colnames(counts)])
      keep_samples <- is.finite(totals) & totals > 0
      if (!is.null(config$group_level)) {
        values <- unname(group_values()[colnames(counts)])
        keep_samples <- keep_samples & values %in% config$group_level
      }
      validate(need(sum(keep_samples) >= 3L, "At least three samples with positive read counts are required in the selected sample set."))
      counts <- counts[, keep_samples, drop = FALSE]
      rownames(counts) <- ranked$labels
      relative <- sweep(counts, 2, totals[keep_samples], "/")
      values <- if (config$binary) {
        1 * (counts > 0)
      } else if (identical(config$transform, "clr")) {
        logged <- log1p(counts)
        sweep(logged, 2, colMeans(logged), "-")
      } else {
        relative
      }
      variable <- apply(values, 1, stats::var)
      keep_taxa <- which(rowMeans(counts > 0) * 100 >= config$min_prevalence &
                           is.finite(variable) & variable > 0)
      validate(need(length(keep_taxa) >= 2L, "Fewer than two variable taxa passed the prevalence filter."))
      ranked_taxa <- keep_taxa[order(rowMeans(relative)[keep_taxa], decreasing = TRUE)]
      values <- values[head(ranked_taxa, max(2L, min(100L, as.integer(config$max_taxa)))), , drop = FALSE]
      correlation <- stats::cor(t(values), method = config$method)
      pairs <- which(upper.tri(correlation), arr.ind = TRUE)
      p <- vapply(seq_len(nrow(pairs)), function(i) {
        suppressWarnings(stats::cor.test(
          values[pairs[i, 1L], ], values[pairs[i, 2L], ], method = config$method, exact = FALSE
        )$p.value)
      }, numeric(1))
      q <- matrix(NA_real_, nrow(correlation), ncol(correlation), dimnames = dimnames(correlation))
      q[pairs] <- stats::p.adjust(p, method = "BH")
      q[cbind(pairs[, 2L], pairs[, 1L])] <- q[pairs]
      list(correlation = correlation, q = q, samples = ncol(values), config = config)
    })

    observeEvent(displayed_result(), {
      choices <- rownames(displayed_result()$correlation)
      selected <- intersect(input$selected_taxa, choices)
      updateSelectizeInput(session, "selected_taxa", choices = choices, selected = selected, server = TRUE)
    })

    selected_correlation <- reactive({
      matrix <- displayed_result()$correlation
      selected <- intersect(rownames(matrix), input$selected_taxa)
      if (!length(selected)) return(matrix)
      matrix[selected, selected, drop = FALSE]
    })

    plot_data <- reactive({
      data <- displayed_result()
      matrix <- selected_correlation()
      if (isTRUE(input$cluster) && nrow(matrix) > 1L) {
        distance <- 1 - matrix
        distance[distance < 0] <- 0
        order <- stats::hclust(stats::as.dist(distance), method = "average")$order
        matrix <- matrix[order, order, drop = FALSE]
      }
      labels <- rownames(matrix)
      pairs <- which(lower.tri(matrix, diag = length(labels) == 1L), arr.ind = TRUE)
      q <- data$q[labels, labels, drop = FALSE][pairs]
      cutoff <- req(input$q_cutoff)
      cells <- data.frame(
        x = factor(labels[pairs[, 2L]], levels = labels),
        y = factor(labels[pairs[, 1L]], levels = rev(labels)),
        correlation = matrix[pairs],
        significant = is.finite(q) & q <= cutoff
      )
      list(matrix = matrix, cells = cells, labels = labels)
    })

    correlogram <- reactive({
      view <- plot_data()
      font <- if (isTruthy(input$base_size)) input$base_size else 11
      circle_size <- if (isTruthy(input$max_circle_size)) input$max_circle_size else 8
      colors <- diverging_plot_colors(input$color_palette)
      angle <- if (isTruthy(input$x_label_angle)) as.numeric(input$x_label_angle) else 45
      visible <- if (isTRUE(input$mask) && length(view$labels) > 1L) {
        view$cells[view$cells$significant, , drop = FALSE]
      } else {
        view$cells
      }
      plot <- ggplot2::ggplot(view$cells, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_tile(width = 1, height = 1, fill = "white", colour = "#D8DEE4", linewidth = 0.25) +
        ggplot2::geom_point(data = visible, ggplot2::aes(size = abs(correlation), colour = correlation))
      if (isTRUE(input$show_significance_marks)) {
        plot <- plot + ggplot2::geom_text(
          data = visible[visible$significant, , drop = FALSE],
          label = "*", fontface = "bold", vjust = 0.75, size = font / 3
        )
      }
      plot +
        ggplot2::scale_size_area(
          limits = c(0, 1), max_size = circle_size, name = "|Correlation|",
          breaks = c(0.25, 0.50, 0.75, 1), labels = c("0.25", "0.50", "0.75", "1.00")
        ) +
        ggplot2::scale_colour_gradient2(
          low = colors[1L], mid = colors[2L], high = colors[3L],
          limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), name = "Correlation",
          guide = ggplot2::guide_colourbar(theme = ggplot2::theme(
            legend.key.width = grid::unit(3, "mm"), legend.key.height = grid::unit(26, "mm")
          ))
        ) +
        ggplot2::scale_x_discrete(limits = view$labels) +
        ggplot2::scale_y_discrete(limits = rev(view$labels)) +
        ggplot2::coord_fixed() +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = font) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            angle = angle, hjust = if (angle == 0) 0.5 else 1,
            vjust = if (angle == 90) 0.5 else 1
          ),
          legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
          legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
          legend.spacing.y = grid::unit(2, "mm"),
          panel.grid = ggplot2::element_blank()
        )
    })

    dimensions <- reactive({
      if (identical(input$plot_size_mode, "manual")) {
        return(c(width = req(input$plot_width), height = req(input$plot_height)))
      }
      n <- nrow(selected_correlation())
      c(width = max(620, 310 + 21 * n), height = max(480, 180 + 21 * n))
    })

    output$plot_container <- renderUI({
      size <- dimensions()
      plotOutput(session$ns("plot"), width = size[["width"]], height = size[["height"]])
    })
    register_automatic_analysis_status(output, active, displayed_result,
                                       "Choose data and analysis settings to view the correlogram.",
                                       updating_message = "Preparing the correlogram. Please wait...")
    output$plot <- renderPlot(correlogram(), res = 96)
    output$figure_legend <- renderUI({
      req_plot_legend_tab(session, "Taxa Correlogram")
      config <- displayed_result()$config
      method <- if (config$binary) "phi correlation on presence/absence" else
        paste(tools::toTitleCase(config$method), "correlation on", toupper(config$transform), "abundance")
      scope <- if (is.null(config$group_level)) "all included samples" else
        paste0("samples with ", config$group_var, " = ", config$group_level)
      cutoff <- format(req(input$q_cutoff))
      plot_legend_box(
        dimensions()[["width"]], "Taxa Correlogram",
        paste(
          paste0("Each cell shows ", method, " for a taxon pair across ", scope, "."),
          paste("Circle area shows the absolute correlation.",
                "The color scale shows negative to positive correlations, with near-white at zero."),
          if (isTRUE(input$cluster)) "Taxa are ordered by hierarchical clustering.",
          if (isTRUE(input$show_significance_marks)) {
            paste0("An asterisk (*) marks a pair with a Benjamini-Hochberg-adjusted p-value (q) at or below ", cutoff, ".")
          },
          if (isTRUE(input$mask)) {
            paste0("Pairs above q = ", cutoff, " or without a valid q-value are shown as empty cells.")
          },
          "Correlation alone does not establish a direct interaction."
        )
      )
    })
    outputOptions(output, "figure_legend", suspendWhenHidden = FALSE, priority = -10)
    output$status_box <- renderUI({
      displayed_result()
      tags$div(
        style = sprintf("width: %dpx; max-width: 100%%;", as.integer(dimensions()[["width"]])),
        hr(),
        h5(icon("circle-info"), "Taxa Correlogram Status"),
        tags$div(class = "app-result-card app-status-card", style = "width: 100%;",
                 verbatimTextOutput(session$ns("status")))
      )
    })
    output$status <- renderText({
      data <- displayed_result()
      config <- data$config
      method <- if (config$binary) "Phi correlation on presence/absence" else
        paste(tools::toTitleCase(config$method), "correlation on", toupper(config$transform), "abundance")
      scope <- if (is.null(config$group_level)) "All included samples" else
        paste(config$group_var, "=", config$group_level)
      q <- data$q[upper.tri(data$q)]
      cutoff <- req(input$q_cutoff)
      paste(
        paste("Analysis scope:", scope),
        paste("Method:", method),
        paste("Taxonomic level:", config$rank),
        paste("Samples used:", data$samples),
        paste("Taxa shown:", nrow(selected_correlation())),
        paste("Unique pairs tested (BH-adjusted):", sum(is.finite(q))),
        paste0("Pairs meeting the FDR cutoff (q <= ", format(cutoff), "): ", sum(is.finite(q) & q <= cutoff)),
        if (isTRUE(input$show_significance_marks)) "An asterisk (*) marks a pair that meets the cutoff." else "Significance marks are hidden.",
        if (isTRUE(input$mask)) "Pairs that do not meet the cutoff are hidden." else "All pairs are shown.",
        sep = "\n"
      )
    })

    output$download_plot <- downloadHandler(
      filename = function() "taxa_correlogram.png",
      content = function(file) {
        size <- dimensions()
        ggplot2::ggsave(file, correlogram(), width = size[["width"]] / 96,
                        height = size[["height"]] / 96, units = "in", dpi = 300,
                        device = "png", bg = "white")
      }
    )
    output$download_matrix <- downloadHandler(
      filename = function() "taxa_correlation_matrix.tsv",
      content = function(file) {
        matrix <- plot_data()$matrix
        table <- data.frame(Taxon = rownames(matrix), matrix, check.names = FALSE)
        readr::write_tsv(table, file)
      }
    )
  })
}
