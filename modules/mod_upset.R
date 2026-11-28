mod_upset_ui <- function(id) {
  ns <- NS(id)
  table_style <- "max-width: 100%; overflow-x: auto; font-size: 0.8em;"
  download_row <- function(plot_id, table_id, label) {
    button_style <- "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
    div(
      style = "display: flex; gap: 4px; align-items: center; flex-wrap: nowrap;",
      downloadButton(ns(plot_id), paste(label, "PNG"), style = button_style),
      downloadButton(ns(table_id), paste(label, "TSV"), style = button_style)
    )
  }
  plot_scroll_layout(
    sidebarPanel(
      width = 2,
      h4(icon("layer-group"), "Core & Indicator Taxa"),
      hr(),
      selectInput(ns("group"), "Grouping Variable", choices = c("Choose a variable" = "")),
      selectizeInput(ns("levels"), "Groups to Include", choices = NULL, multiple = TRUE),
      hr(),
      h4(icon("sliders"), "Analysis & Plot Settings"),
      selectInput(ns("analysis_mode"), "Analysis Mode",
                  choices = c("Core Taxa" = "core", "Indicator Taxa (IndVal)" = "indval")),
      selectInput(ns("rank"), "Taxonomic Level", choices = c("Genus", "ASV"), selected = "Genus"),
      conditionalPanel("input.analysis_mode === 'indval'", ns = ns,
        selectInput(ns("indval_max_pattern_size"), "Maximum Groups per Pattern",
                    choices = c("1", "2", "3"), selected = "3"),
        numericInput(ns("q_cutoff"), "Maximum BH-adjusted q-value", 0.05, min = 0, max = 1, step = 0.01)
      ),
      conditionalPanel("input.analysis_mode === 'core'", ns = ns,
        numericInput(ns("min_relative_abundance"), "Minimum Relative Abundance (%)", 0.01,
                     min = 0, max = 100, step = 0.01),
        numericInput(ns("min_prevalence"), "Minimum Group Prevalence (%)", 50, min = 0, max = 100, step = 1)
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
      selectInput(ns("pattern_size"), "Groups per Pattern", choices = c("All" = "all")),
      selectizeInput(ns("pattern_members"), "Groups in Pattern", choices = NULL,
                     multiple = TRUE, options = list(placeholder = "Any group")),
      selectInput(ns("pattern_limit"), "Patterns to Show",
                  choices = c("Top 15" = "15", "Top 25" = "25", "All" = "all"), selected = "25"),
      hr(),
      h5(icon("download"), "Download"),
      conditionalPanel("input.analysis_mode === 'core'", ns = ns,
        download_row("download_plot", "download_core_taxa", "Core")
      ),
      conditionalPanel("input.analysis_mode === 'indval'", ns = ns,
        download_row("download_indval_plot", "download_indval_results", "IndVal")
      )
    ),
    mainPanel(
      width = 9,
      h4("Core & Indicator Taxa"),
      automatic_analysis_view(ns,
        conditionalPanel("input.analysis_mode === 'core'", ns = ns,
          tabsetPanel(
            tabPanel("Plot",
              div(
                class = "taxa-plot-stage",
                uiOutput(ns("plot_container"), class = "taxa-plot-container", style = "overflow-x: auto;"),
                taxa_plot_loading("Core Taxa UpSet")
              ),
              plot_legend_output(ns, "figure_legend", c("plot"))
            ),
            tabPanel("Results Table",
              div(id = ns("core_taxa_frame"), style = table_style,
                  DT::DTOutput(ns("core_taxa"))))
          ),
          uiOutput(ns("core_table_width")),
          uiOutput(ns("status_box"))
        ),
        conditionalPanel("input.analysis_mode === 'indval'", ns = ns,
          tabsetPanel(
            id = ns("indval_result_tab"),
            tabPanel("Plot",
              div(class = "taxa-plot-stage",
                uiOutput(ns("indval_plot_container"), class = "taxa-plot-container", style = "overflow-x: auto;"),
                taxa_plot_loading("IndVal Plot")
              ),
              plot_legend_output(ns, "indval_figure_legend", c("indval_plot"))
            ),
            tabPanel("Results Table",
              div(id = ns("indval_table_frame"), style = table_style,
                  DT::DTOutput(ns("indval_table"))))
          ),
          uiOutput(ns("indval_table_width")),
          uiOutput(ns("indval_status_box"))
        )
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

    indval_pattern_cap <- reactive({
      max(1L, min(3L, length(unique(input$levels)) - 1L))
    })
    indval_pattern_size <- reactive({
      min(as.integer(req(input$indval_max_pattern_size)), indval_pattern_cap())
    })
    observeEvent(input$levels, {
      req(length(unique(input$levels)) >= 2L)
      updateSelectInput(session, "indval_max_pattern_size",
                        choices = as.character(seq_len(indval_pattern_cap())),
                        selected = as.character(indval_pattern_cap()))
    })

    settings <- function(mode) {
      validate(need(isTruthy(input$group), "Choose a grouping variable."))
      groups <- unique(input$levels)
      validate(need(
        length(groups) >= 2L && length(groups) <= 12L,
        "Select between 2 and 12 groups."
      ))
      config <- list(group = input$group, groups = groups, rank = req(input$rank),
                     mode = mode)
      if (identical(config$mode, "core")) {
        config$detection <- req(input$min_relative_abundance)
        config$cutoff <- req(input$min_prevalence)
      } else {
        config$max_pattern_size <- indval_pattern_size()
      }
      config
    }

    analysis_inputs <- function(mode) {
      ps <- ps_obj()
      validate(need(!is.null(ps), "Load and preprocess data to analyze core and indicator taxa."))
      list(ps = ps, config = settings(mode))
    }
    core_inputs <- reactive(analysis_inputs("core"))
    indval_inputs <- reactive(analysis_inputs("indval"))

    compute_result <- function(request) {
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

      if (identical(config$mode, "indval")) {
        validate(need(requireNamespace("indicspecies", quietly = TRUE),
                      "Install the indicspecies package to run Indicator Taxa."))
        validate(need(all(sample_counts >= 2L), "Indicator Taxa requires at least two samples per selected group."))
        detected <- rowSums(counts[, keep, drop = FALSE] > 0) > 0L
        validate(need(any(detected), "No taxa were detected in the selected groups."))
        relative <- sweep(counts[detected, keep, drop = FALSE], 2L, totals[keep], "/")
        rownames(relative) <- ranked$labels[detected]
        groups <- factor(values[keep], levels = config$groups)
        fit <- indicspecies::multipatt(t(relative), groups, func = "IndVal.g",
                                     max.order = config$max_pattern_size,
                                     control = permute::how(nperm = 999))
        group_rows <- match(config$groups, as.character(groups))
        pattern_groups <- lapply(seq_len(ncol(fit$comb)), function(i) {
          config$groups[as.logical(fit$comb[group_rows, i])]
        })
        pattern_labels <- make.unique(vapply(pattern_groups, paste, character(1), collapse = " + "))
        best <- fit$sign$index
        taxa <- seq_along(best)
        index <- cbind(taxa, best)
        results <- data.frame(
          Taxon = rownames(fit$sign),
          Group = pattern_labels[best],
          IndVal = fit$sign$stat,
          Specificity = fit$A[index],
          Fidelity = fit$B[index],
          PValue = fit$sign$p.value,
          QValue = stats::p.adjust(fit$sign$p.value, method = "BH")
        )
        for (group in config$groups) {
          results[[paste0("Prevalence (%) - ", group)]] <-
            100 * rowMeans(relative[, groups == group, drop = FALSE] > 0)
        }
        return(list(indicators = results, config = config, sample_counts = sample_counts,
                    patterns = stats::setNames(pattern_groups, pattern_labels)))
      }

      detected <- counts[, keep, drop = FALSE] > 0
      if (config$detection > 0) {
        relative <- 100 * sweep(counts[, keep, drop = FALSE], 2L, totals[keep], "/")
        detected <- relative >= config$detection
      }
      detected_counts <- vapply(config$groups, function(group) {
        rowSums(detected[, values[keep] == group, drop = FALSE])
      }, numeric(nrow(counts)))
      detected_counts <- matrix(detected_counts, nrow = nrow(counts),
                                dimnames = list(ranked$labels, config$groups))
      prevalence <- 100 * sweep(detected_counts, 2L, as.integer(sample_counts), "/")
      membership <- 1L * (detected_counts > 0L & prevalence >= config$cutoff)
      core <- rowSums(membership) > 0L
      validate(need(any(core), "No taxa met the detection and prevalence criteria in the selected groups."))
      relative <- 100 * sweep(counts[core, keep, drop = FALSE], 2L, totals[keep], "/")
      mean_relative_abundance <- vapply(config$groups, function(group) {
        rowMeans(relative[, values[keep] == group, drop = FALSE])
      }, numeric(sum(core)))
      mean_relative_abundance <- matrix(mean_relative_abundance, nrow = sum(core),
                                        dimnames = list(ranked$labels[core], config$groups))
      list(membership = membership[core, , drop = FALSE],
           detected_counts = detected_counts[core, , drop = FALSE],
           prevalence = prevalence[core, , drop = FALSE],
           mean_relative_abundance = mean_relative_abundance,
           sample_counts = sample_counts, config = config)
    }
    core_result <- automatic_analysis(
      function() isTRUE(active()) && identical(input$analysis_mode, "core"), core_inputs, compute_result
    )
    indval_result <- automatic_analysis(
      function() isTRUE(active()) && identical(input$analysis_mode, "indval"), indval_inputs, compute_result
    )
    displayed_result <- reactive({
      if (identical(req(input$analysis_mode), "indval")) indval_result() else core_result()
    })

    combinations <- reactive({
      data <- displayed_result()
      req(identical(data$config$mode, "core"))
      groups <- colnames(data$membership)
      ordered <- levels(plot_order_factor(groups, data$config$group, plot_order(), groups))
      membership <- data$membership[, ordered, drop = FALSE]
      ComplexHeatmap::make_comb_mat(membership, mode = "distinct")
    })

    observeEvent(list(input$levels, input$analysis_mode, input$indval_max_pattern_size), {
      groups <- input$levels
      max_size <- if (identical(input$analysis_mode, "indval")) {
        min(length(groups), indval_pattern_size())
      } else length(groups)
      sizes <- c("all", as.character(seq_len(max_size)))
      current <- input$pattern_size
      updateSelectInput(session, "pattern_size", choices = stats::setNames(sizes, c("All", sizes[-1L])),
                        selected = if (isTruthy(current) && current %in% sizes) current else "all")
      updateSelectizeInput(session, "pattern_members", choices = groups,
                           selected = intersect(input$pattern_members, groups))
    }, ignoreNULL = FALSE)

    filter_patterns <- function(members) {
      size <- req(input$pattern_size)
      required <- intersect(input$pattern_members, input$levels)
      keep <- if (identical(size, "all")) rep(TRUE, length(members)) else lengths(members) == as.integer(size)
      keep & vapply(members, function(pattern) all(required %in% pattern), logical(1))
    }
    limit_patterns <- function(indices) {
      limit <- req(input$pattern_limit)
      if (identical(limit, "all")) indices else utils::head(indices, as.integer(limit))
    }
    core_patterns <- reactive({
      m <- combinations()
      codes <- ComplexHeatmap::comb_name(m)
      groups <- ComplexHeatmap::set_name(m)
      members <- lapply(strsplit(codes, "", fixed = TRUE), function(code) groups[code == "1"])
      patterns <- data.frame(Code = codes,
                             Groups = make.unique(vapply(members, paste, character(1), collapse = " + ")),
                             GroupCount = lengths(members), IntersectionSize = as.integer(ComplexHeatmap::comb_size(m)))
      patterns <- patterns[filter_patterns(members), , drop = FALSE]
      patterns[order(-patterns$IntersectionSize, patterns$GroupCount, patterns$Groups), , drop = FALSE]
    })
    plotted_core_patterns <- reactive({
      patterns <- core_patterns()
      patterns[limit_patterns(seq_len(nrow(patterns))), , drop = FALSE]
    })
    plot_combinations <- reactive({
      codes <- plotted_core_patterns()$Code
      validate(need(length(codes) > 0L, "No intersections match the current filters."))
      m <- combinations()
      m[match(codes, ComplexHeatmap::comb_name(m))]
    })

    core_taxa_table <- reactive({
      data <- displayed_result()
      req(identical(data$config$mode, "core"))
      taxa <- rownames(data$membership)
      groups <- colnames(data$membership)
      rows <- data.frame(
        Taxon = rep(taxa, times = length(groups)),
        Group = rep(groups, each = length(taxa)),
        DetectedSamples = as.vector(data$detected_counts),
        TotalSamples = rep(as.integer(data$sample_counts), each = length(taxa)),
        PrevalencePct = as.vector(data$prevalence),
        MeanRelativeAbundancePct = as.vector(data$mean_relative_abundance),
        Core = as.vector(data$membership) == 1L
      )
      rows[order(rows$Taxon, match(rows$Group, groups)), , drop = FALSE]
    })

    indicator_result <- reactive({
      data <- displayed_result()
      req(identical(data$config$mode, "indval"))
      data
    })
    significant_indicators <- reactive({
      results <- indicator_result()$indicators
      results <- results[which(results$QValue <= req(input$q_cutoff)), , drop = FALSE]
      results[order(-results$IndVal, results$QValue, results$Taxon), , drop = FALSE]
    })
    indicator_group_order <- reactive({
      config <- indicator_result()$config
      levels(plot_order_factor(config$groups, config$group, plot_order(), config$groups))
    })
    matching_patterns <- reactive({
      result <- indicator_result()
      counts <- table(significant_indicators()$Group)
      groups <- as.character(names(counts))
      members <- result$patterns[groups]
      keep <- filter_patterns(members)
      patterns <- data.frame(Group = groups[keep], Count = as.integer(counts[keep]), Size = lengths(members)[keep])
      patterns[order(-patterns$Count, patterns$Size, patterns$Group), , drop = FALSE]
    })
    indicator_groups <- reactive({
      patterns <- matching_patterns()
      patterns[limit_patterns(seq_len(nrow(patterns))), , drop = FALSE]
    })
    selected_group <- reactiveVal(NULL)
    selected_core <- reactiveVal(NULL)
    indval_hit_map <- NULL
    core_hit_map <- NULL
    observeEvent(indval_inputs(), selected_group(NULL), ignoreInit = TRUE)
    observeEvent(core_inputs(), selected_core(NULL), ignoreInit = TRUE)
    observeEvent(matching_patterns(), {
      if (!is.null(selected_group()) && !selected_group() %in% matching_patterns()$Group) {
        selected_group(NULL)
      }
    })
    observeEvent(core_patterns(), {
      if (!is.null(selected_core()) && !selected_core() %in% core_patterns()$Code) selected_core(NULL)
    })
    clicked_pattern <- function(hit, click, labels) {
      req(!is.null(hit), identical(hit$labels, labels),
          is.finite(click$x), is.finite(click$y),
          click$x >= hit$left, click$x < hit$right, click$y <= hit$bottom)
      column <- floor((click$x - hit$left) / (hit$right - hit$left) * length(hit$labels)) + 1L
      hit$labels[column]
    }
    observeEvent(input$indval_click, {
      group <- clicked_pattern(indval_hit_map, input$indval_click$coords_css, indicator_groups()$Group)
      selected_group(if (identical(group, selected_group())) NULL else group)
    })
    observeEvent(input$core_click, {
      code <- clicked_pattern(core_hit_map, input$core_click$coords_css, plotted_core_patterns()$Code)
      selected_core(if (identical(code, selected_core())) NULL else code)
    })
    indval_width <- reactive({
      if (identical(input$plot_size_mode, "manual")) return(req(input$plot_width))
      max(650L, 320L + 22L * nrow(indicator_groups()))
    })
    indval_height <- reactive({
      if (identical(input$plot_size_mode, "manual")) return(req(input$plot_height))
      max(430L, 270L + 24L * length(indicator_group_order()))
    })
    output$indval_plot_container <- renderUI({
      req(identical(input$analysis_mode, "indval"))
      plotOutput(session$ns("indval_plot"), width = indval_width(), height = "auto",
                 click = session$ns("indval_click"))
    })
    outputOptions(output, "indval_plot_container", suspendWhenHidden = FALSE)
    output$indval_status_box <- renderUI({
      indicator_result()
      div(
        style = sprintf("width: %dpx; max-width: 100%%;", as.integer(indval_width())),
        hr(),
        h5(icon("circle-info"), "IndVal Analysis Status"),
        div(class = "app-result-card app-status-card", style = "width: 100%;",
            verbatimTextOutput(session$ns("indval_status")))
      )
    })
    output$indval_figure_legend <- renderUI({
      req_plot_legend_tab(session, "mod_upset")
      shown <- nrow(indicator_groups())
      total <- nrow(matching_patterns())
      group <- selected_group()
      plot_legend_box(
        indval_width(), "Indicator Taxa by Group Pattern",
        p("IndVal scores are based on relative abundance. Significance is assessed with 999 permutations and assumes independent samples."),
        p("Each top bar counts significant indicator taxa assigned to that pattern. Filled dots show the groups in the pattern. Left bars show significant indicator taxa per group across all patterns at the current q-value cutoff. A taxon assigned to a multi-group pattern contributes to each member group. Use Groups in Pattern to require all selected groups. Click a plot column to filter the Results Table."),
        sprintf("Showing %d of %d matching patterns with significant indicator taxa.", shown, total),
        if (!is.null(group)) div(paste("Selected pattern:", group))
      )
    })
    outputOptions(output, "indval_figure_legend", suspendWhenHidden = FALSE, priority = -10)
    indval_plot <- reactive({
      groups <- indicator_groups()
      validate(need(nrow(groups) > 0L, "No significant patterns match the current filters."))
      members <- indicator_group_order()
      patterns <- indicator_result()$patterns
      significant <- significant_indicators()
      assignments <- patterns[significant$Group]
      membership <- vapply(members, function(member) {
        vapply(assignments, function(pattern) member %in% pattern, logical(1))
      }, logical(nrow(significant)))
      membership <- matrix(as.integer(membership), nrow = nrow(significant),
                           dimnames = list(NULL, members))
      combinations <- ComplexHeatmap::make_comb_mat(membership, mode = "distinct")
      codes <- vapply(patterns[groups$Group], function(pattern) {
        paste0(as.integer(members %in% pattern), collapse = "")
      }, character(1))
      list(matrix = combinations[match(codes, ComplexHeatmap::comb_name(combinations))],
           labels = groups$Group)
    })
    draw_indval <- function(size) {
      plot <- indval_plot()
      selected <- selected_group()
      colors <- if (is.null(selected)) "#246A73" else
        ifelse(plot$labels == selected, "#246A73", "#A8B9C6")
      draw_upset_matrix(plot$matrix, size, "indval_patterns", colors, force_counts = TRUE)
      invisible(plot$labels)
    }
    output$indval_plot <- renderPlot({
      indval_hit_map <<- NULL
      size <- c(width = indval_width(), height = indval_height())
      labels <- draw_indval(size)
      indval_hit_map <<- pattern_hit_map("indval_patterns", size, labels)
    }, width = function() indval_width(), height = function() indval_height(), res = 96)
    outputOptions(output, "indval_plot", suspendWhenHidden = FALSE)
    displayed_indicators <- reactive({
      results <- significant_indicators()
      group <- selected_group()
      labels <- if (is.null(group)) matching_patterns()$Group else group
      results <- results[results$Group %in% labels, , drop = FALSE]
      results
    })
    output$indval_status <- renderText({
      data <- indicator_result()
      config <- data$config
      sample_counts <- data$sample_counts
      group <- selected_group()
      paste(
        paste("Grouping variable:", config$group),
        paste("Taxonomic level:", config$rank),
        paste("Samples per group:", paste(sprintf("%s (%d)", names(sample_counts),
                                                   as.integer(sample_counts)), collapse = ", ")),
        paste("Maximum groups per pattern:", config$max_pattern_size),
        paste("BH-adjusted q-value cutoff:", format(req(input$q_cutoff))),
        paste("Taxa tested:", nrow(data$indicators)),
        paste("Significant indicator taxa:", nrow(significant_indicators())),
        paste("Matching patterns plotted:", nrow(indicator_groups()), "of", nrow(matching_patterns())),
        paste("Result table:", if (is.null(group)) "All matching patterns" else group),
        sep = "\n"
      )
    })
    output$indval_table_width <- renderUI({
      tags$style(sprintf("#%s { width: %dpx; }", session$ns("indval_table_frame"),
                         as.integer(indval_width())))
    })
    outputOptions(output, "indval_table_width", suspendWhenHidden = FALSE)
    output$core_table_width <- renderUI({
      tags$style(sprintf("#%s { width: %dpx; }",
                         session$ns("core_taxa_frame"), as.integer(dimensions()[["width"]])))
    })
    outputOptions(output, "core_table_width", suspendWhenHidden = FALSE)
    output$indval_table <- DT::renderDT({
      table <- DT::datatable(displayed_indicators(), rownames = FALSE,
                             options = list(pageLength = 10, scrollX = TRUE, order = list()))
      table <- DT::formatRound(table, "IndVal", digits = 3)
      table <- DT::formatPercentage(table, c("Specificity", "Fidelity"), digits = 1)
      table <- DT::formatSignif(table, c("PValue", "QValue"), digits = 3)
      DT::formatRound(table, grep("^Prevalence", names(displayed_indicators()), value = TRUE), digits = 1)
    })
    outputOptions(output, "indval_table", suspendWhenHidden = FALSE)
    output$download_indval_results <- downloadHandler(
      filename = function() "indicator_taxa_current_results.tsv",
      content = function(file) readr::write_tsv(displayed_indicators(), file)
    )

    dimensions <- reactive({
      if (identical(input$plot_size_mode, "manual")) {
        return(c(width = req(input$plot_width), height = req(input$plot_height)))
      }
      c(width = max(650, 320 + 22 * nrow(plotted_core_patterns())),
        height = max(430, 270 + 24 * length(ComplexHeatmap::set_size(combinations()))))
    })

    pattern_hit_map <- function(name, size, labels) {
      bounds <- NULL
      ComplexHeatmap::decorate_heatmap_body(name, {
        left <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"), valueOnly = TRUE)
        right <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(0, "npc"), valueOnly = TRUE)
        device <- grDevices::dev.size("in")
        bounds <- list(left = left$x / device[1L] * size[["width"]],
                       right = right$x / device[1L] * size[["width"]],
                       bottom = (1 - left$y / device[2L]) * size[["height"]], labels = labels)
      }, envir = environment())
      bounds
    }

    draw_upset_matrix <- function(m, size, name, intersection_color = "#246A73", force_counts = FALSE) {
      font <- if (isTruthy(input$base_size)) input$base_size else 11
      sizes <- ComplexHeatmap::comb_size(m)
      show_counts <- force_counts || (length(sizes) <= 12L && size[["width"]] >= 260 + 22 * length(sizes))
      top_height <- grid::unit(max(100, size[["height"]] * 0.35) / 96, "inches")
      annotation_font <- grid::gpar(fontfamily = "sans", fontsize = font, col = "#334155")
      axis_font <- grid::gpar(fontfamily = "sans", fontsize = plot_axis_label_size(font, 11), col = "black")
      plot <- ComplexHeatmap::UpSet(
        m, name = name, set_order = seq_along(ComplexHeatmap::set_size(m)),
        comb_order = seq_along(sizes),
        row_names_gp = grid::gpar(fontfamily = "sans", fontsize = plot_axis_label_size(font, 11),
                                  fontface = "bold", col = "black"),
        comb_col = intersection_color, pt_size = grid::unit(2.5, "mm"), lwd = 1.7,
        bg_col = c("#F3F6F8", "#FFFFFF"), bg_pt_col = "#CFD8DF",
        top_annotation = ComplexHeatmap::upset_top_annotation(
          m, height = top_height, gp = grid::gpar(fill = intersection_color), annotation_name_gp = annotation_font,
          axis_param = list(gp = axis_font),
          add_numbers = show_counts, numbers_rot = 0,
          numbers_gp = grid::gpar(fontfamily = "sans", fontsize = max(7, font - 2), col = "#243446"),
          numbers_offset = grid::unit(1.5, "mm"),
          ylim = if (show_counts) c(0, max(sizes) * 1.2) else NULL
        ),
        left_annotation = ComplexHeatmap::upset_left_annotation(
          m, gp = grid::gpar(fill = "#A8B9C6"), annotation_name_gp = annotation_font,
          axis_param = list(direction = "reverse", gp = axis_font)
        ),
        right_annotation = NULL
      )
      ComplexHeatmap::draw(plot)
    }
    draw_upset <- function() {
      codes <- plotted_core_patterns()$Code
      selected <- selected_core()
      colors <- if (is.null(selected) || !selected %in% codes) "#246A73" else
        ifelse(codes == selected, "#246A73", "#A8B9C6")
      draw_upset_matrix(plot_combinations(), dimensions(), "core_patterns", colors)
    }

    output$plot_container <- renderUI({
      size <- dimensions()
      plotOutput(session$ns("plot"), width = size[["width"]], height = size[["height"]],
                 click = session$ns("core_click"))
    })
    outputOptions(output, "plot_container", suspendWhenHidden = FALSE)
    register_automatic_analysis_status(output, active, displayed_result,
                                       "Choose data and at least two groups to analyze core or indicator taxa.",
                                       updating_message = "Analyzing taxa. Please wait...")
    output$plot <- renderPlot({
      core_hit_map <<- NULL
      draw_upset()
      core_hit_map <<- pattern_hit_map("core_patterns", dimensions(), plotted_core_patterns()$Code)
    }, width = function() dimensions()[["width"]], height = function() dimensions()[["height"]], res = 96)
    outputOptions(output, "plot", suspendWhenHidden = FALSE)
    output$figure_legend <- renderUI({
      req_plot_legend_tab(session, "mod_upset")
      req(identical(input$analysis_mode, "core"))
      config <- displayed_result()$config
      shown <- nrow(plotted_core_patterns())
      total <- nrow(core_patterns())
      selected <- core_patterns()$Groups[match(selected_core(), core_patterns()$Code)]
      plot_legend_box(
        dimensions()[["width"]], "Core Taxa UpSet",
        paste(
          paste0("Each row is a level of ", config$group, ". A taxon at the ", config$rank,
                 " level belongs to a group if detected in at least one sample and its prevalence reaches ",
                 format(config$cutoff), "% in that group."),
          "Relative abundance is calculated using each sample's total read count after preprocessing filters.",
          if (config$detection > 0) {
            paste0("A taxon is detected in a sample when its relative abundance is at least ", format(config$detection), "%.")
          } else {
            "At a 0% threshold, any positive read count counts as detection."
          },
          "Side bars show the number of included taxa per group. Each top bar counts taxa that belong to exactly the groups marked by filled dots among the selected groups.",
          sprintf("Showing %d of %d matching patterns.", shown, total),
          "Click a plot column to see its taxa in the Results Table."
        ),
        if (length(selected)) div(paste("Selected pattern:", selected))
      )
    })
    outputOptions(output, "figure_legend", suspendWhenHidden = FALSE, priority = -10)
    displayed_core_taxa <- reactive({
      patterns <- core_patterns()
      if (!is.null(selected_core())) {
        patterns <- patterns[patterns$Code == selected_core(), , drop = FALSE]
      }
      m <- combinations()
      taxa <- unlist(lapply(patterns$Code, function(code) ComplexHeatmap::extract_comb(m, code)),
                     use.names = FALSE)
      rows <- core_taxa_table()
      rows[rows$Taxon %in% taxa & rows$Core, names(rows) != "Core", drop = FALSE]
    })
    output$core_taxa <- DT::renderDT({
      table <- DT::datatable(displayed_core_taxa(), rownames = FALSE,
                             options = list(pageLength = 10, scrollX = TRUE))
      table <- DT::formatRound(table, "PrevalencePct", digits = 2)
      DT::formatRound(table, "MeanRelativeAbundancePct", digits = 4)
    })
    output$status_box <- renderUI({
      displayed_result()
      tags$div(
        style = sprintf("width: %dpx; max-width: 100%%;", as.integer(dimensions()[["width"]])),
        hr(),
        h5(icon("circle-info"), "Core Taxa UpSet Plot Status"),
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
        paste0("Minimum relative abundance: ", format(data$config$detection), "%"),
        paste0("Minimum group prevalence: ", format(data$config$cutoff), "%"),
        paste("Taxa included:", nrow(data$membership)),
        paste("Samples per group:", paste(paste0(names(sample_counts), " (", as.integer(sample_counts), ")"), collapse = ", ")),
        paste("Matching patterns plotted:", nrow(plotted_core_patterns()), "of", nrow(core_patterns())),
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
    output$download_core_taxa <- downloadHandler(
      filename = function() "upset_core_taxa.tsv",
      content = function(file) readr::write_tsv(displayed_core_taxa(), file)
    )
    output$download_indval_plot <- downloadHandler(
      filename = function() "indicator_taxa_plot.png",
      content = function(file) {
        size <- c(width = indval_width(), height = indval_height())
        grDevices::png(file, width = size[["width"]] / 96, height = size[["height"]] / 96,
                       units = "in", res = 300, bg = "white")
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_indval(size)
      }
    )
  })
}
