mod_taxa_summary_ui <- function(id) {
  ns <- NS(id)
  plot_scroll_layout(
    sidebarPanel(
      width = 2,
      h4(icon("chart-line"), "Taxa Bubble Plot"),
      hr(),
      profile_group_selectors(ns, "primary", "Primary Group (Facet)",
                              "secondary", "Secondary Group (Optional)"),
      hr(),
      h4(icon("sliders"), "Analysis & Plot Settings"),
      selectInput(ns("rank"), "Taxonomic Level", choices = c("Genus", "ASV"), selected = "Genus"),
      numericInput(ns("top_n"), "Top Taxa by Mean Relative Abundance", 20, min = 1, max = 100, step = 1),
      hr(),
      h4(icon("sliders"), "Display Options"),
      selectInput(ns("plot_size_mode"), "Plot Size", choices = c("Auto" = "auto", "Manual" = "manual")),
      conditionalPanel(
        "input.plot_size_mode === 'manual'", ns = ns,
        numericInput(ns("plot_width"), "Plot Width (px)", 500, min = 300, max = 4000, step = 50),
        numericInput(ns("plot_height"), "Plot Height (px)", 700, min = 300, max = 4000, step = 50)
      ),
      numericInput(ns("base_size"), "Base Font Size", 12, min = 6, max = 30, step = 1),
      selectInput(ns("x_label_angle"), "X-Axis Label Angle", choices = c(
        "0° Horizontal" = "0", "45° Diagonal" = "45", "90° Vertical" = "90"
      ), selected = "45"),
      selectInput(ns("color_palette"), "Color Palette", sequential_plot_palettes, selected = "plasma"),
      selectInput(ns("bubble_size"), "Bubble Size", choices = c(
        "Prevalence" = "prevalence", "Mean Relative Abundance" = "abundance",
        "Fixed" = "fixed"
      )),
      selectInput(ns("bubble_color"), "Bubble Color", choices = c(
        "Auto" = "auto", "Mean Relative Abundance" = "abundance", "Prevalence" = "prevalence"
      )),
      checkboxInput(ns("show_bubble_outline"), "Show Bubble Outline", value = FALSE),
      checkboxInput(ns("hide_bubble_grid"), "Hide Vertical Grid Lines", value = FALSE),
      conditionalPanel(
        "input.bubble_size === 'fixed'", ns = ns,
        numericInput(ns("fixed_bubble_size"), "Fixed Bubble Size", 7, min = 1, max = 20, step = 0.5)
      ),
      conditionalPanel(
        "input.bubble_size !== 'fixed'", ns = ns,
        numericInput(ns("max_bubble_size"), "Maximum Bubble Size", 7, min = 1, max = 20, step = 0.5)
      ),
      hr(),
      h5(icon("download"), "Download"),
      div(
        style = "display: flex; flex-wrap: wrap; gap: 4px;",
        downloadButton(ns("download_plot"), "Plot (PNG)"),
        downloadButton(ns("download_summary"), "Summary (TSV)")
      )
    ),
    mainPanel(
      width = 9,
      h4("Taxa Bubble Plot"),
      div(
        class = "taxa-plot-stage",
        uiOutput(ns("bubble_container"), class = "taxa-plot-container"),
        taxa_plot_loading("taxa bubble plot")
      ),
      plot_legend_output(ns, "figure_legend", c("bubble"), class = "taxa-plot-details")
    )
  )
}

mod_taxa_summary_server <- function(id, ps_obj, plot_order = function() NULL) {
  moduleServer(id, function(input, output, session) {
    metadata <- reactive({
      ps <- ps_obj()
      req(!is.null(ps))
      as.data.frame(phyloseq::sample_data(ps))
    })

    grouping_columns <- reactive({
      if (is.null(ps_obj())) return(character())
      columns <- names(metadata())
      columns[gsub("[^a-z0-9]", "", tolower(columns)) != "sampleid"]
    })
    register_profile_group_selectors(
      input, output, session, grouping_columns, "primary", "secondary",
      primary_extra = c("All samples" = ""), secondary_none = c("None" = "")
    )

    observeEvent(ps_obj(), {
      taxonomy <- phyloseq::tax_table(ps_obj(), errorIfNULL = FALSE)
      ranks <- unique(c("ASV", colnames(taxonomy)))
      selected <- isolate(input$rank)
      default_rank <- if ("Genus" %in% ranks) "Genus" else "ASV"
      updateSelectInput(session, "rank", choices = ranks,
                        selected = if (isTruthy(selected) && selected %in% ranks) selected else default_rank)
    })

    rank_counts <- reactive({
      metadata()
      taxa_rank_counts(ps_obj(), req(input$rank))
    })

    taxa_summary <- reactive({
      rank <- rank_counts()
      sample_ids <- colnames(rank$counts)
      meta <- metadata()[sample_ids, , drop = FALSE]
      primary <- input$primary
      secondary <- input$secondary
      has_primary <- isTruthy(primary) && primary %in% names(meta)
      has_secondary <- has_primary && isTruthy(secondary) && secondary %in% setdiff(names(meta), primary)
      primary_values <- if (has_primary) trimws(as.character(meta[[primary]])) else rep("All samples", nrow(meta))
      secondary_values <- if (has_secondary) trimws(as.character(meta[[secondary]])) else primary_values
      totals <- unname(rank$sample_totals[sample_ids])
      keep <- is.finite(totals) & totals > 0 &
        !is.na(primary_values) & nzchar(primary_values) &
        !is.na(secondary_values) & nzchar(secondary_values)
      validate(need(any(keep), "No samples have positive read counts and complete group labels."))

      order_config <- plot_order()
      primary_groups <- plot_order_factor(primary_values[keep], primary, order_config)
      x_groups <- if (has_secondary) {
        plot_order_factor(secondary_values[keep], secondary, order_config)
      } else {
        primary_groups
      }
      groups <- if (has_secondary) {
        expand.grid(x = levels(x_groups), panel = levels(primary_groups), stringsAsFactors = FALSE)
      } else {
        data.frame(x = levels(x_groups), panel = "All samples")
      }

      counts <- rank$counts[, keep, drop = FALSE]
      relative <- sweep(counts, 2, totals[keep], "/") * 100
      top_n <- input$top_n
      if (!isTruthy(top_n) || !is.finite(top_n)) top_n <- 20L
      selected <- head(order(rowMeans(relative), decreasing = TRUE), max(1L, min(100L, as.integer(top_n))))
      taxa <- rank$labels[selected]
      counts <- counts[selected, , drop = FALSE]
      relative <- relative[selected, , drop = FALSE]

      rows <- lapply(seq_len(nrow(groups)), function(j) {
        in_group <- x_groups == groups$x[j]
        if (has_secondary) in_group <- in_group & primary_groups == groups$panel[j]
        n <- sum(in_group)
        detected <- if (n > 0L) rowSums(counts[, in_group, drop = FALSE] > 0) else NA_real_
        data.frame(
          taxon = taxa, x = groups$x[j], panel = groups$panel[j], samples = n,
          detected = detected,
          prevalence = if (n > 0L) 100 * detected / n else NA_real_,
          abundance = if (n > 0L) rowMeans(relative[, in_group, drop = FALSE]) else NA_real_
        )
      })
      result <- do.call(rbind, rows)
      result$taxon <- factor(result$taxon, levels = rev(taxa))
      result$x <- factor(result$x, levels = levels(x_groups))
      result$panel <- factor(result$panel, levels = unique(groups$panel))
      result$detail <- ifelse(result$samples > 0L,
        sprintf("Mean relative abundance: %.3g%%<br>Prevalence: %.3g%% (%d/%d samples)",
                result$abundance, result$prevalence, result$detected, result$samples),
        "No samples in this group combination"
      )
      group_label <- if (has_secondary) {
        paste(htmltools::htmlEscape(result$panel), htmltools::htmlEscape(result$x), sep = " / ")
      } else {
        htmltools::htmlEscape(result$x)
      }
      result$detail <- paste(htmltools::htmlEscape(result$taxon), group_label, result$detail, sep = "<br>")
      list(cells = result, faceted = has_secondary, x_label = if (has_secondary) secondary else if (has_primary) primary else NULL)
    })

    bubble_metrics <- reactive({
      size <- switch(
        if (isTruthy(input$bubble_size)) input$bubble_size else "prevalence",
        fixed = NULL, abundance = "abundance", "prevalence"
      )
      color <- input$bubble_color
      if (!isTruthy(color) || !color %in% c("abundance", "prevalence")) {
        color <- if (identical(size, "abundance")) "prevalence" else "abundance"
      }
      list(size = size, color = color)
    })

    metric_label <- function(metric) {
      if (identical(metric, "prevalence")) "Prevalence (%)" else "Mean relative abundance (%)"
    }

    metric_guide_label <- function(metric) {
      if (identical(metric, "abundance")) "Mean relative\nabundance (%)" else metric_label(metric)
    }

    bubble_legend_data <- reactiveVal(NULL)

    output$figure_legend <- renderUI({
      req_plot_legend_tab(session, "Taxa Bubble Plot")
      info <- req(bubble_legend_data())
      summary <- info$summary
      dimensions <- info$dimensions
      width <- if (is.numeric(dimensions$width)) dimensions$width else dimensions$max_width
      grouping <- if (summary$faceted) {
        paste0("Panels show ", info$primary, " groups. Columns show ", summary$x_label, " levels.")
      } else if (isTruthy(summary$x_label)) {
        paste0("Columns show ", summary$x_label, " groups.")
      } else {
        "The column summarizes all included samples."
      }
      metrics <- info$metrics
      size_text <- if (is.null(metrics$size)) "All bubbles have the same area." else
        paste0("Bubble area shows ", tolower(metric_label(metrics$size)), ".")
      interpretation <- paste(
        size_text, paste0("Color shows ", tolower(metric_label(metrics$color)), "."),
        "No bubble is drawn when a taxon was not detected. A gray cross marks a group combination with no samples."
      )
      tags$div(
        style = sprintf("width: %dpx;", as.integer(width)),
        tags$div(style = "font-weight: 600; margin-bottom: 4px;",
                 "Taxa Bubble Plot"),
        tags$div(paste0(summary$n_taxa, " taxa at the ", info$rank, " level. ",
                        grouping, " ", interpretation)),
        tags$div(
          style = "margin-top: 4px;",
          "Prevalence is the percentage of samples in which a taxon was detected. Mean relative abundance includes zeros from samples where the taxon was not detected."
        )
      )
    })
    outputOptions(output, "figure_legend", suspendWhenHidden = FALSE, priority = -10)

    label_angle <- reactive(if (isTruthy(input$x_label_angle)) as.numeric(input$x_label_angle) else 45)
    base_font <- reactive(if (isTruthy(input$base_size)) input$base_size else 12)
    bubble_marker_size <- reactive({
      fixed <- is.null(bubble_metrics()$size)
      size <- if (fixed) input$fixed_bubble_size else input$max_bubble_size
      if (isTruthy(size)) size else 7
    })

    vertical_title_offset <- reactive({
      summary <- taxa_summary()
      if (label_angle() != 90 || !isTruthy(summary$x_label)) return(0)
      label_width <- max(nchar(levels(summary$cells$x), type = "width"))
      ceiling(label_width * base_font() * 0.8) + 16
    })

    summary_plot <- function() {
      summary <- taxa_summary()
      cells <- summary$cells
      observed <- cells[cells$samples > 0L, , drop = FALSE]
      missing <- cells[cells$samples == 0L, , drop = FALSE]
      palette <- if (isTruthy(input$color_palette)) input$color_palette else "plasma"
      metric_limits <- function(metric) {
        if (identical(metric, "prevalence")) return(c(0, 100))
        upper <- max(c(0, cells$abundance), na.rm = TRUE)
        c(0, if (upper > 0) upper else 1)
      }
      metric_color_scale <- function(metric, aesthetic) {
        ggplot2::scale_colour_gradientn(
          name = metric_guide_label(metric), colours = sequential_plot_colors(palette),
          limits = metric_limits(metric), na.value = "grey85", aesthetics = aesthetic
        )
      }
      angle <- label_angle()
      radius_px <- bubble_marker_size() * 96 / 25.4 / 2 + 5
      panel_height <- plot_dimensions()$height - vertical_title_offset() - 145
      y_padding <- max(0.85, (nlevels(cells$taxon) - 1L) * radius_px /
                         max(75, panel_height - 2 * radius_px))
      p <- ggplot2::ggplot(cells, ggplot2::aes(x = x, y = taxon, text = detail)) +
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::scale_y_discrete(drop = FALSE, expand = ggplot2::expansion(add = y_padding)) +
        ggplot2::labs(x = summary$x_label, y = NULL) +
        ggplot2::theme_minimal(base_size = base_font(), base_family = "sans") +
        ggplot2::theme(
          legend.title = ggplot2::element_text(size = ggplot2::rel(0.78)),
          legend.text = ggplot2::element_text(size = ggplot2::rel(0.75)),
          legend.key.size = grid::unit(4, "mm"),
          legend.spacing.y = grid::unit(2, "mm"),
          strip.background = ggplot2::element_rect(fill = "#F1F3F5", colour = "#595959", linewidth = 0.8),
          strip.text.x = ggplot2::element_text(face = "bold", vjust = 0.5),
          axis.text.x = ggplot2::element_text(
            angle = angle, hjust = if (angle == 0) 0.5 else 1, vjust = if (angle == 90) 0.5 else 1
          ),
          axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = if (angle == 90) 10 else 5))
        )
      if (summary$faceted) p <- p + ggplot2::facet_grid(cols = ggplot2::vars(panel))

      metrics <- bubble_metrics()
      size <- bubble_marker_size()
      points <- observed[observed$detected > 0, , drop = FALSE]
      points$bubble_color <- points[[metrics$color]]
      outline <- isTRUE(input$show_bubble_outline)
      stroke <- if (outline) 0.35 else 0
      p <- p + ggplot2::theme(
        panel.border = ggplot2::element_rect(colour = "#333333", fill = NA, linewidth = 0.9),
        plot.margin = ggplot2::margin(t = 14, r = 8, b = 8, l = 8)
      ) +
        ggplot2::geom_blank() +
        ggplot2::geom_point(data = missing, shape = 4, colour = "grey60", size = 2)
      if (isTRUE(input$hide_bubble_grid)) {
        p <- p + ggplot2::theme(
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank()
        )
      }
      if (is.null(metrics$size)) {
        p <- p + ggplot2::geom_point(
          data = points, ggplot2::aes(fill = bubble_color),
          shape = 21, colour = "#333333", stroke = stroke, size = size
        )
      } else {
        points$bubble_area <- points[[metrics$size]]
        p <- p + ggplot2::geom_point(
          data = points, ggplot2::aes(size = bubble_area, fill = bubble_color),
          shape = 21, colour = "#333333", stroke = stroke
        ) +
          ggplot2::scale_size_area(
            name = metric_guide_label(metrics$size), limits = metric_limits(metrics$size), max_size = size
          )
      }
      p + metric_color_scale(metrics$color, "fill")
    }

    plot_dimensions <- reactive({
      if (identical(input$plot_size_mode, "manual")) {
        req(isTruthy(input$plot_width), isTruthy(input$plot_height))
        return(list(width = input$plot_width, height = input$plot_height))
      }
      cells <- taxa_summary()$cells
      n_columns <- nlevels(cells$x) * nlevels(cells$panel)
      per_column <- min(50, max(28, 110 / sqrt(n_columns)))
      auto_width <- ceiling(min(2200, max(500, 350 + n_columns * per_column)))
      list(width = "100%", max_width = auto_width,
           height = max(380, 700 + 18 * (nlevels(cells$taxon) - 20)) + vertical_title_offset())
    })

    plot_container <- function(id) {
      dimensions <- plot_dimensions()
      div(
        style = if (is.null(dimensions$max_width)) "overflow-x: auto;" else
          sprintf("overflow-x: auto; max-width: %dpx;", dimensions$max_width),
        plotly::plotlyOutput(session$ns(id), width = dimensions$width, height = dimensions$height)
      )
    }

    output$bubble_container <- renderUI({
      validate_taxa_plot_data(ps_obj(), "taxa bubble plot")
      plot_container("bubble")
    })
    interactive_plot <- function() {
      dimensions <- plot_dimensions()
      manual <- is.numeric(dimensions$width)
      width <- if (manual) dimensions$width else NULL
      plot <- summary_plot()
      offset <- vertical_title_offset()
      if (offset > 0) plot <- plot + ggplot2::labs(x = NULL)
      widget <- plotly::ggplotly(
        plot, width = width, height = dimensions$height, tooltip = "text"
      )
      font_family <- "sans-serif"
      compact_colorbar <- function(bar) {
        bar$lenmode <- "pixels"
        bar$len <- 135
        bar$thickness <- 12
        bar$tickfont$size <- base_font() * 0.85
        bar$tickfont$family <- font_family
        if (is.list(bar$title)) {
          bar$title$text <- gsub("\n", "<br>", bar$title$text, fixed = TRUE)
          bar$title$font$size <- base_font() * 0.9
          bar$title$font$family <- font_family
        } else if (is.character(bar$title)) {
          bar$title <- gsub("\n", "<br>", bar$title, fixed = TRUE)
          bar$titlefont$size <- base_font() * 0.9
          bar$titlefont$family <- font_family
        }
        bar
      }
      widget$x$data <- lapply(widget$x$data, function(trace) {
        if (!is.null(trace$colorbar)) trace$colorbar <- compact_colorbar(trace$colorbar)
        if (!is.null(trace$marker$colorbar)) trace$marker$colorbar <- compact_colorbar(trace$marker$colorbar)
        trace
      })
      axes <- grep("^[xy]axis[0-9]*$", names(widget$x$layout), value = TRUE)
      for (axis in axes) {
        widget$x$layout[[axis]]$tickfont$family <- font_family
        if (startsWith(axis, "yaxis")) {
          # Reserve space for the labels without drawing visible tick marks.
          widget$x$layout[[axis]]$ticks <- "outside"
          widget$x$layout[[axis]]$ticklen <- 10
          widget$x$layout[[axis]]$tickcolor <- "rgba(0,0,0,0)"
        }
        if (is.list(widget$x$layout[[axis]]$title)) {
          widget$x$layout[[axis]]$title$font$family <- font_family
        }
      }
      widget$x$layout$shapes <- lapply(widget$x$layout$shapes, function(shape) {
        # ggplotly can place panel borders under the data and drop their line width.
        if (identical(shape$type, "rect") && identical(shape$fillcolor, "transparent") &&
            is.null(shape$ysizemode)) {
          shape$layer <- "above"
          shape$line$width <- 2
        } else if (identical(shape$type, "rect") && identical(shape$ysizemode, "pixel")) {
          shape$line$width <- 1
        }
        shape
      })
      strips <- Filter(function(shape) identical(shape$ysizemode, "pixel"), widget$x$layout$shapes)
      annotations <- lapply(widget$x$layout$annotations, function(annotation) {
        annotation$font$family <- font_family
        # Plotly converts facet labels to annotations anchored at the bottom of the strip.
        if (is.null(annotation$annotationType) && identical(annotation$yanchor, "bottom")) {
          for (strip in strips) {
            if (isTRUE(abs(annotation$x - mean(c(strip$x0, strip$x1))) < 1e-6) &&
                isTRUE(abs(annotation$y - strip$yanchor) < 1e-6)) {
              annotation$yanchor <- "middle"
              annotation$yshift <- mean(c(strip$y0, strip$y1))
              annotation$text <- paste0("<b>", htmltools::htmlEscape(annotation$text), "</b>")
              break
            }
          }
        }
        annotation
      })
      if (offset > 0) {
        # Position the shared title below the vertical tick labels in single and faceted plots.
        axis_title <- list(
          text = htmltools::htmlEscape(taxa_summary()$x_label),
          xref = "paper", yref = "paper", x = 0.5, y = 0, yshift = -offset,
          yanchor = "top", showarrow = FALSE,
          font = list(family = font_family, size = base_font() * 96 / 72)
        )
        margin <- widget$x$layout$margin
        margin$b <- max(margin$b, offset + 2 * base_font() + 12)
        widget <- plotly::layout(widget, margin = margin)
        annotations <- c(annotations, list(axis_title))
      }
      widget$x$layout$annotations <- annotations
      widget <- plotly::layout(widget, autosize = !manual,
                               font = list(family = font_family),
                               legend = list(font = list(family = font_family, size = base_font() * 0.85)))
      widget <- plotly::config(widget, responsive = !manual)
      summary <- taxa_summary()
      bubble_legend_data(list(
        summary = list(faceted = summary$faceted, x_label = summary$x_label,
                       n_taxa = nlevels(summary$cells$taxon)),
        dimensions = dimensions, primary = input$primary, rank = input$rank,
        metrics = bubble_metrics()
      ))
      htmlwidgets::onRender(widget, "function(el) {
        el.setAttribute('data-taxa-ready', 'true');
        el.dispatchEvent(new Event('plotlegend:widgetready', {bubbles: true}));
      }")
    }

    output$bubble <- plotly::renderPlotly(interactive_plot())

    output$download_plot <- downloadHandler(
      filename = function() "taxa_summary_bubble.png",
      content = function(file) {
        dimensions <- plot_dimensions()
        width <- dimensions$width
        if (!is.numeric(width)) {
          width <- session$clientData[[paste0("output_", session$ns("bubble"), "_width")]]
          if (length(width) != 1L || !is.finite(width) || width <= 0) width <- dimensions$max_width
          width <- min(width, dimensions$max_width)
        }
        ggplot2::ggsave(
          file, plot = summary_plot(), device = "png",
          width = width / 96, height = dimensions$height / 96, units = "in",
          dpi = 300, bg = "white"
        )
      }
    )
    output$download_summary <- downloadHandler(
      filename = function() "taxa_summary.tsv",
      content = function(file) {
        summary <- taxa_summary()
        cells <- summary$cells
        export <- data.frame(
          Taxon = cells$taxon, PrimaryGroup = if (summary$faceted) cells$panel else cells$x,
          SecondaryGroup = if (summary$faceted) cells$x else NA_character_,
          SampleCount = cells$samples, DetectedCount = cells$detected,
          MeanRelativeAbundancePct = cells$abundance, PrevalencePct = cells$prevalence
        )
        readr::write_tsv(export, file, na = "NA")
      }
    )
  })
}
