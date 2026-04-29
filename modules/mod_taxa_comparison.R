library(shiny)
library(ggplot2)
library(dplyr)
library(phyloseq)

## UI
mod_taxa_comparison_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML("
      .taxa-result-card {
        height: 130px;
        overflow-y: auto;
        font-size: 12px;
        line-height: 1.45;
        background: #f8fafc;
        border: 1px solid #d9e2ec;
        border-radius: 10px;
        padding: 10px 12px;
        box-shadow: 0 1px 2px rgba(16, 24, 40, 0.04);
        margin-bottom: 10px;
      }
      .taxa-result-card pre {
        margin: 0;
        padding: 0;
        border: 0;
        background: transparent;
        font-size: 12px;
        line-height: 1.45;
        white-space: pre-wrap;
      }
      .well h4 { font-size: 16px; }
      .well h5 { font-size: 13px; }
      .well .control-label { font-size: 12px; }
      .well .checkbox label { font-size: 12px; }
      .well .form-control { font-size: 12px; }
      .well .btn { font-size: 11px; }
    ")),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        h4(icon("square-poll-vertical"), "Taxa Comparison"),
        hr(),
        uiOutput(ns("group_selector")),
        uiOutput(ns("secondary_group_selector")),
        hr(),
        div(
          style = "display: flex; align-items: center; gap: 6px; margin-bottom: 6px;",
          tags$input(
            id = ns("enable_longitudinal"),
            type = "checkbox",
            style = "margin: 0;"
          ),
          tags$span(
            style = "display: inline-flex; align-items: center; gap: 4px;",
            icon("timeline"),
            "Within-Subject Pairing"
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("enable_longitudinal")),
          selectInput(
            ns("subject_id_var"),
            "Subject Identifier:",
            choices = character(0),
            selected = NULL
          )
        ),
        div(
          style = "display: flex; align-items: center; gap: 6px; margin-bottom: 6px;",
          tags$input(
            id = ns("show_trend_line"),
            type = "checkbox",
            style = "margin: 0;"
          ),
          tags$span(
            style = "display: inline-flex; align-items: center; gap: 4px;",
            icon("chart-line"),
            "Show trend line"
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("show_trend_line")),
          selectInput(
            ns("trend_line_method"),
            "Trend Line Method:",
            choices = c(
              "Spearman" = "spearman",
              "Pearson" = "pearson",
              "Linear regression (lm)" = "lm"
            ),
            selected = "spearman"
          )
        ),
        hr(),
        h4(icon("sliders"), "Plot Settings"),
        selectInput(
          ns("tax_level"),
          "Taxonomic Level:",
          choices = c("ASV", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
          selected = "Genus"
        ),
        selectizeInput(
          ns("taxa_selected"),
          "Taxa to Compare:",
          choices = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Select one or more taxa",
            plugins = list("remove_button")
          )
        ),
        checkboxInput(
          ns("show_p_lt_0_05_only"),
          "Show only taxa with p-value < 0.05",
          value = FALSE
        ),
        selectInput(
          ns("abundance_mode"),
          "Transformation:",
          choices = c(
            "CLR Abundance" = "clr",
            "Relative Abundance (%)" = "relative",
            "Log TSS (log10(% + 1))" = "log_tss"
          ),
          selected = "clr"
        ),
        selectInput(
          ns("plot_type"),
          "Plot Type:",
          choices = c("Box plot" = "boxplot", "Bar plot" = "barplot", "Scatter plot" = "scatter"),
          selected = "boxplot"
        ),
        selectInput(
          ns("color_palette"),
          "Color Palette:",
          choices = c(
            "Set2" = "set2",
            "Dark2" = "dark2",
            "Paired" = "paired",
            "Grayscale" = "gray"
          ),
          selected = "set2"
        ),
        checkboxInput(ns("use_ggpattern"), "Use ggpattern (Bar plot)", value = FALSE),
        hr(),
        h4(icon("sliders"), "P-value Options"),
        checkboxInput(ns("show_p_val"), "Show p-value bars", value = TRUE),
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("show_p_val")),
          checkboxInput(ns("show_sig_only"), "Show only p-value < 0.05", value = TRUE),
          checkboxInput(ns("show_sig_as_marks"), "Show as significant marks", value = TRUE)
        ),
        selectInput(
          ns("stat_method"),
          "Statistical Method:",
          choices = c(
            "Wilcoxon" = "wilcox.test",
            "T-test" = "t.test"
          ),
          selected = "t.test"
        ),
        selectInput(
          ns("p_adjust_method"),
          "Pairwise Multiple Testing Correction:",
          choices = c(
            "None" = "none",
            "Holm" = "holm",
            "BH (FDR)" = "BH",
            "Bonferroni" = "bonferroni"
          ),
          selected = "BH"
        ),
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Plot Dimensions"),
        checkboxInput(ns("manual_facet_layout"), "Manual facet rows/cols", value = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("manual_facet_layout")),
          numericInput(ns("facet_ncol"), "Facet Columns (ncol):", value = 0, min = 0, step = 1),
          numericInput(ns("facet_nrow"), "Facet Rows (nrow):", value = 0, min = 0, step = 1),
          tags$small("Use 0 to keep automatic value for each field.")
        ),
        numericInput(ns("plot_width"), "Plot Width (px):", value = 400, min = 300, step = 50),
        numericInput(ns("plot_height"), "Plot Height (px):", value = 500, min = 300, step = 50),
        numericInput(ns("base_size"), "Base Font Size:", value = 11, min = 6, max = 30, step = 1),
        hr(),
        h5(icon("download"), "Download Plot"),
        div(
          style = "display: flex; gap: 4px; flex-wrap: nowrap;",
          downloadButton(
            ns("download_taxa_plot"),
            "Plot (PNG)",
            style = "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
          ),
          downloadButton(
            ns("download_taxa_data"),
            "Data (TSV)",
            style = "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
          )
        ),
        hr(),
        h5(icon("table"), "Download Full Matrix"),
        div(
          style = "display: flex; gap: 4px; flex-wrap: nowrap;",
          downloadButton(
            ns("download_raw_matrix"),
            "Raw Counts (TSV)",
            style = "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
          ),
          downloadButton(
            ns("download_rel_matrix"),
            "Relative Abundance (TSV)",
            style = "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
          )
        )
      ),
      mainPanel(
        width = 9,
        h4("Taxa Comparison"),
        plotOutput(ns("taxa_comparison_plot"), height = "auto"),
        uiOutput(ns("taxa_legend_box")),
        uiOutput(ns("taxa_status_separator")),
        h5(icon("circle-info"), "Comparison Status"),
        uiOutput(ns("comparison_status_box"))
      )
    )
  )
}

## Server
mod_taxa_comparison_server <- function(id, ps_obj, meta_cols, active_tab = NULL) {
  moduleServer(id, function(input, output, session) {
    resolve_meta_colname <- function(requested, available) {
      if (is.null(requested) || is.null(available) || length(available) == 0) {
        return(requested)
      }
      if (requested %in% available) {
        return(requested)
      }
      matched <- available[make.names(available) == make.names(requested)]
      if (length(matched) > 0) {
        return(matched[1])
      }
      requested
    }

    size_update_in_progress <- reactiveVal(FALSE)
    taxa_selected_memory <- reactiveVal(character(0))

    compute_taxa_comparison_auto_dims <- function(n_facet_cols, n_bars_per_facet, n_facet_rows) {
      width <- min(max(300, 220 + n_facet_cols * n_bars_per_facet * 20), 2200)
      height <- min(max(360, 220 + n_facet_rows * 200), 1600)
      list(width = width, height = height)
    }

    get_manual_facet_layout <- function(is_secondary, n_taxa) {
      if (isTRUE(is_secondary) || !isTRUE(input$manual_facet_layout)) {
        return(list(ncol = NULL, nrow = NULL))
      }
      manual_ncol <- input$facet_ncol
      manual_nrow <- input$facet_nrow
      manual_ncol <- if (!is.null(manual_ncol) && is.finite(manual_ncol) && manual_ncol > 0) as.integer(manual_ncol) else NULL
      manual_nrow <- if (!is.null(manual_nrow) && is.finite(manual_nrow) && manual_nrow > 0) as.integer(manual_nrow) else NULL
      if (!is.null(manual_ncol)) {
        manual_ncol <- min(manual_ncol, max(1L, as.integer(n_taxa)))
      }
      if (!is.null(manual_nrow)) {
        manual_nrow <- min(manual_nrow, max(1L, as.integer(n_taxa)))
      }
      # Prevent facet_wrap errors when manual nrow/ncol capacity is smaller than panel count.
      if (!is.null(manual_ncol) && !is.null(manual_nrow)) {
        capacity <- as.integer(manual_ncol * manual_nrow)
        if (capacity < as.integer(n_taxa)) {
          manual_nrow <- as.integer(ceiling(as.integer(n_taxa) / max(1L, manual_ncol)))
        }
      }
      list(ncol = manual_ncol, nrow = manual_nrow)
    }

    measure_facet_layout <- function(df, is_secondary, manual_ncol = NULL, manual_nrow = NULL) {
      if (nrow(df) == 0) {
        return(list(n_facet_cols = 1L, n_facet_rows = 1L))
      }
      if (isTRUE(is_secondary)) {
        return(list(
          n_facet_cols = max(1L, as.integer(dplyr::n_distinct(df$Taxa))),
          n_facet_rows = max(1L, as.integer(dplyr::n_distinct(df$PrimaryGroup)))
        ))
      }

      p_layout <- ggplot2::ggplot(df, ggplot2::aes(x = Group, y = AbundancePlot)) +
        ggplot2::geom_blank() +
        ggplot2::facet_wrap(~ Taxa, scales = "free_y", ncol = manual_ncol, nrow = manual_nrow)
      layout_df <- ggplot2::ggplot_build(p_layout)$layout$layout
      if (is.null(layout_df) || nrow(layout_df) == 0) {
        return(list(n_facet_cols = 1L, n_facet_rows = 1L))
      }
      list(
        n_facet_cols = max(1L, as.integer(max(layout_df$COL, na.rm = TRUE))),
        n_facet_rows = max(1L, as.integer(max(layout_df$ROW, na.rm = TRUE)))
      )
    }

    maybe_auto_adjust_taxa_plot_size <- function(width, height, force = FALSE) {
      if (size_update_in_progress()) {
        return()
      }
      current_width <- suppressWarnings(as.numeric(input$plot_width))
      current_height <- suppressWarnings(as.numeric(input$plot_height))
      if (!isTRUE(force) &&
          is.finite(current_width) && is.finite(current_height) &&
          as.integer(current_width) == as.integer(width) &&
          as.integer(current_height) == as.integer(height)) {
        return()
      }
      size_update_in_progress(TRUE)
      updateNumericInput(session, "plot_width", value = width)
      updateNumericInput(session, "plot_height", value = height)
      session$onFlushed(function() {
        size_update_in_progress(FALSE)
      }, once = TRUE)
    }

    apply_disambiguated_taxrank <- function(ps, tax_level) {
      if (is.null(ps) || tax_level == "ASV") {
        return(ps)
      }
      tt_obj <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
      if (is.null(tt_obj)) {
        return(ps)
      }
      tt <- as.data.frame(tt_obj, stringsAsFactors = FALSE)
      tax_cols <- colnames(tt)
      if (!tax_level %in% tax_cols) {
        return(ps)
      }

      # If target rank is "uncultured" / "Unassigned", prepend parent-rank label
      # so these entries are not merged across different higher taxa.
      target_raw <- as.character(tt[[tax_level]])
      target_norm <- tolower(trimws(target_raw))
      # Handle variants such as "g__uncultured", "uncultured bacterium", "k__Unassigned", etc.
      target_norm <- gsub("^[a-z]__", "", target_norm)
      idx_placeholder <- !is.na(target_norm) & grepl("(uncultured|unassigned)", target_norm)
      if (!any(idx_placeholder)) {
        return(ps)
      }

      tax_ranks <- phyloseq::rank_names(ps)
      rank_pos <- match(tax_level, tax_ranks)
      parent_candidates <- character(0)
      if (!is.na(rank_pos) && rank_pos > 1) {
        parent_candidates <- rev(tax_ranks[seq_len(rank_pos - 1)])
        parent_candidates <- parent_candidates[parent_candidates %in% tax_cols]
      }
      parent_val <- rep("UnclassifiedParent", nrow(tt))
      if (length(parent_candidates) > 0) {
        for (parent_rank in parent_candidates) {
          candidate_val <- as.character(tt[[parent_rank]])
          candidate_norm <- tolower(trimws(candidate_val))
          candidate_norm <- gsub("^[a-z]__", "", candidate_norm)
          candidate_is_placeholder <- is.na(candidate_norm) | !nzchar(candidate_norm) | grepl("(uncultured|unassigned)", candidate_norm)
          use_idx <- idx_placeholder & parent_val == "UnclassifiedParent" & !candidate_is_placeholder
          parent_val[use_idx] <- candidate_val[use_idx]
        }
      }

      tt[[tax_level]][idx_placeholder] <- parent_val[idx_placeholder]
      phyloseq::tax_table(ps) <- phyloseq::tax_table(as.matrix(tt))
      ps
    }

    is_active_tab <- reactive({
      if (is.null(active_tab)) {
        return(TRUE)
      }
      current_tab <- tryCatch(active_tab(), error = function(e) NULL)
      identical(current_tab, "Taxa Comparison")
    })
    
    output$group_selector <- renderUI({
      req(meta_cols())
      group_choices <- meta_cols()
      non_sampleid_choices <- setdiff(group_choices, "SampleID")
      selected_col <- if (length(non_sampleid_choices) > 0) non_sampleid_choices[1] else if (length(group_choices) > 0) group_choices[1] else NULL
      selectInput(
        session$ns("group_var"),
        "Primary Group:",
        choices = group_choices,
        selected = selected_col
      )
    })

    output$secondary_group_selector <- renderUI({
      req(meta_cols())
      primary_input <- input$group_var
      if (is.null(primary_input) || !nzchar(primary_input)) {
        base_choices <- setdiff(meta_cols(), "SampleID")
        primary_input <- if (length(base_choices) > 0) base_choices[1] else NULL
      }
      resolved_primary <- resolve_meta_colname(primary_input, meta_cols())
      group_choices <- setdiff(meta_cols(), c("SampleID", primary_input, resolved_primary))
      selectInput(
        session$ns("secondary_group_var"),
        "Secondary Group (Optional):",
        choices = c("(None)" = "none", group_choices),
        selected = "none"
      )
    })

    observe({
      req(meta_cols(), input$group_var)
      resolved_primary <- resolve_meta_colname(input$group_var, meta_cols())
      resolved_secondary <- resolve_meta_colname(input$secondary_group_var, meta_cols())
      longitudinal_choices <- setdiff(
        meta_cols(),
        c("SampleID", input$group_var, resolved_primary, input$secondary_group_var, resolved_secondary)
      )
      selected_subject <- isolate(input$subject_id_var)
      if (!is.null(selected_subject) && selected_subject %in% longitudinal_choices) {
        next_selected <- selected_subject
      } else if (length(longitudinal_choices) > 0) {
        next_selected <- longitudinal_choices[1]
      } else {
        next_selected <- character(0)
      }
      updateSelectInput(
        session,
        "subject_id_var",
        choices = longitudinal_choices,
        selected = next_selected
      )
    })

    observeEvent(input$show_trend_line, {
      if (isTRUE(input$show_trend_line) && isTRUE(input$show_p_val)) {
        updateCheckboxInput(session, "show_p_val", value = FALSE)
      }
      if (isTRUE(input$show_trend_line) && !identical(input$plot_type, "scatter")) {
        updateSelectInput(session, "plot_type", selected = "scatter")
      }
    }, ignoreInit = TRUE)

    taxa_long_data <- reactive({
      req(isTRUE(is_active_tab()))
      req(ps_obj(), input$group_var, input$tax_level)
      
      ps <- ps_obj()
      tax_level <- input$tax_level
      sample_cols <- colnames(as.data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE))
      primary_col <- resolve_meta_colname(input$group_var, sample_cols)
      secondary_col <- resolve_meta_colname(input$secondary_group_var, sample_cols)
      is_secondary <- !is.null(secondary_col) && secondary_col != "none"
      
      if (tax_level != "ASV") {
        validate(need(!is.null(phyloseq::tax_table(ps)), "Taxonomy table is missing."))
        tax_ranks <- phyloseq::rank_names(ps)
        validate(need(tax_level %in% tax_ranks, paste0("Level '", tax_level, "' not found.")))
        ps <- apply_disambiguated_taxrank(ps, tax_level)
        ps <- phyloseq::tax_glom(ps, taxrank = tax_level, NArm = FALSE)
      }
      
      abundance_mode <- input$abundance_mode
      if (is.null(abundance_mode) || !abundance_mode %in% c("clr", "relative", "log_tss")) {
        abundance_mode <- "clr"
      }

      if (abundance_mode %in% c("relative", "log_tss")) {
        ps <- phyloseq::transform_sample_counts(ps, function(x) {
          s <- sum(x)
          if (s == 0) x else x / s
        })
      } else {
        otu_mat <- as(phyloseq::otu_table(ps), "matrix")
        if (!phyloseq::taxa_are_rows(ps)) {
          otu_mat <- t(otu_mat)
        }
        
        otu_clr <- apply(otu_mat, 2, function(x) {
          x <- as.numeric(x) + 1
          log_x <- log(x)
          log_x - mean(log_x, na.rm = TRUE)
        })
        
        otu_clr <- as.matrix(otu_clr)
        rownames(otu_clr) <- rownames(otu_mat)
        colnames(otu_clr) <- colnames(otu_mat)
        
        otu_obj <- phyloseq::otu_table(otu_clr, taxa_are_rows = TRUE)
        sdata_obj <- phyloseq::sample_data(ps)
        tax_obj <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
        ps <- if (is.null(tax_obj)) {
          phyloseq::phyloseq(otu_obj, sdata_obj)
        } else {
          phyloseq::phyloseq(otu_obj, tax_obj, sdata_obj)
        }
      }
      
      df <- phyloseq::psmelt(ps)
      
      meta_df <- as.data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE)
      for (meta_col in colnames(meta_df)) {
        syntactic_col <- make.names(meta_col)
        if (!meta_col %in% names(df) && syntactic_col %in% names(df)) {
          names(df)[names(df) == syntactic_col] <- meta_col
        }
      }
      primary_col <- resolve_meta_colname(primary_col, names(df))
      secondary_col <- resolve_meta_colname(secondary_col, names(df))
      validate(need(primary_col %in% names(df), paste0("Primary grouping variable '", input$group_var, "' is not available.")))
      if (is_secondary) {
        validate(need(secondary_col %in% names(df), paste0("Secondary grouping variable '", input$secondary_group_var, "' is not available.")))
      }
      
      if (tax_level == "ASV") {
        df$Taxa <- as.character(df$OTU)
      } else {
        taxa_values <- as.character(df[[tax_level]])
        taxa_values[is.na(taxa_values) | taxa_values == ""] <- as.character(df$OTU[is.na(taxa_values) | taxa_values == ""])
        df$Taxa <- taxa_values
      }
      
      df$PrimaryGroup <- as.factor(df[[primary_col]])
      if (is_secondary) {
        df$SecondaryGroup <- as.factor(df[[secondary_col]])
        df <- df[!is.na(df$PrimaryGroup) & !is.na(df$SecondaryGroup), , drop = FALSE]
        df$Group <- df$SecondaryGroup
      } else {
        df <- df[!is.na(df$PrimaryGroup), , drop = FALSE]
        df$Group <- df$PrimaryGroup
      }
      
      if (abundance_mode == "relative") {
        df$AbundancePlot <- df$Abundance * 100
      } else if (abundance_mode == "log_tss") {
        df$AbundancePlot <- log10((df$Abundance * 100) + 1)
      } else {
        df$AbundancePlot <- df$Abundance
      }
      df <- df[is.finite(df$AbundancePlot), , drop = FALSE]
      
      df
    })

    taxa_stats <- reactive({
      req(isTRUE(is_active_tab()))
      req(taxa_long_data(), input$group_var)
      df <- taxa_long_data()
      subject_col <- resolve_meta_colname(input$subject_id_var, names(df))
      has_subject <- isTRUE(input$enable_longitudinal) &&
        !is.null(subject_col) &&
        nzchar(subject_col) &&
        subject_col %in% names(df)

      df %>%
        dplyr::group_by(Taxa) %>%
        dplyr::summarise(
          mean_abundance = mean(AbundancePlot, na.rm = TRUE),
          p_value = tryCatch({
            n_groups <- dplyr::n_distinct(Group)
            if (n_groups == 2) {
              if (!has_subject) {
                if (identical(input$stat_method, "t.test")) {
                  stats::t.test(AbundancePlot ~ Group)$p.value
                } else {
                  stats::wilcox.test(AbundancePlot ~ Group)$p.value
                }
              } else {
                sub_df <- dplyr::cur_data_all()
                valid_subject <- !is.na(sub_df[[subject_col]]) & nzchar(as.character(sub_df[[subject_col]]))
                sub_df <- sub_df[valid_subject, , drop = FALSE]
                if (nrow(sub_df) == 0) {
                  return(NA_real_)
                }
                g_levels <- unique(as.character(sub_df$Group))
                if (length(g_levels) != 2) {
                  return(NA_real_)
                }

                g1 <- sub_df[sub_df$Group == g_levels[1], c(subject_col, "AbundancePlot"), drop = FALSE]
                g2 <- sub_df[sub_df$Group == g_levels[2], c(subject_col, "AbundancePlot"), drop = FALSE]
                colnames(g1) <- c("Subject", "v1")
                colnames(g2) <- c("Subject", "v2")
                g1 <- g1 %>%
                  dplyr::group_by(Subject) %>%
                  dplyr::summarise(v1 = stats::median(v1, na.rm = TRUE), .groups = "drop")
                g2 <- g2 %>%
                  dplyr::group_by(Subject) %>%
                  dplyr::summarise(v2 = stats::median(v2, na.rm = TRUE), .groups = "drop")
                paired_df <- dplyr::inner_join(g1, g2, by = "Subject")
                paired_df <- paired_df[is.finite(paired_df$v1) & is.finite(paired_df$v2), , drop = FALSE]
                if (nrow(paired_df) <= 1) {
                  return(NA_real_)
                }

                if (identical(input$stat_method, "t.test")) {
                  stats::t.test(paired_df$v1, paired_df$v2, paired = TRUE)$p.value
                } else {
                  stats::wilcox.test(paired_df$v1, paired_df$v2, paired = TRUE)$p.value
                }
              }
            } else if (n_groups > 2) {
              stats::kruskal.test(AbundancePlot ~ Group)$p.value
            } else {
              NA_real_
            }
          }, error = function(e) NA_real_),
          .groups = "drop"
        )
    })
    
    tax_level_ps_raw <- reactive({
      req(isTRUE(is_active_tab()))
      req(ps_obj(), input$tax_level)
      ps <- ps_obj()
      tax_level <- input$tax_level
      
      if (tax_level != "ASV") {
        validate(
          need(!is.null(phyloseq::tax_table(ps)), "Taxonomy table is missing.")
        )
        tax_ranks <- phyloseq::rank_names(ps)
        validate(
          need(tax_level %in% tax_ranks, paste0("Taxonomic level '", tax_level, "' is not available in this dataset."))
        )
        ps <- apply_disambiguated_taxrank(ps, tax_level)
        ps <- phyloseq::tax_glom(ps, taxrank = tax_level, NArm = FALSE)
      }
      ps
    })
    
    tax_matrices <- reactive({
      req(isTRUE(is_active_tab()))
      req(tax_level_ps_raw(), input$tax_level)
      ps <- tax_level_ps_raw()
      tax_level <- input$tax_level
      
      otu_mat <- as(phyloseq::otu_table(ps), "matrix")
      if (!phyloseq::taxa_are_rows(ps)) {
        otu_mat <- t(otu_mat)
      }
      
      if (tax_level == "ASV") {
        taxa_labels <- rownames(otu_mat)
      } else {
        tax_tab <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE)
        taxa_labels <- as.character(tax_tab[[tax_level]])
        missing_idx <- is.na(taxa_labels) | taxa_labels == ""
        taxa_labels[missing_idx] <- rownames(otu_mat)[missing_idx]
      }
      taxa_labels <- make.unique(as.character(taxa_labels))
      rownames(otu_mat) <- taxa_labels
      
      raw_df <- as.data.frame(t(otu_mat), stringsAsFactors = FALSE, check.names = FALSE)
      raw_df <- tibble::rownames_to_column(raw_df, var = "SampleID")
      
      rel_mat <- apply(otu_mat, 2, function(x) {
        s <- sum(x, na.rm = TRUE)
        if (s == 0) rep(0, length(x)) else (x / s) * 100
      })
      rel_mat <- as.matrix(rel_mat)
      rownames(rel_mat) <- rownames(otu_mat)
      colnames(rel_mat) <- colnames(otu_mat)
      
      rel_df <- as.data.frame(t(rel_mat), stringsAsFactors = FALSE, check.names = FALSE)
      rel_df <- tibble::rownames_to_column(rel_df, var = "SampleID")
      
      list(raw = raw_df, rel = rel_df)
    })
    
    observeEvent(list(taxa_stats(), input$tax_level, input$show_p_lt_0_05_only), {
      req(isTRUE(is_active_tab()))
      taxa_summary <- taxa_stats()

      if (isTRUE(input$show_p_lt_0_05_only)) {
        taxa_summary <- taxa_summary %>%
          dplyr::filter(!is.na(p_value) & p_value < 0.05) %>%
          dplyr::arrange(p_value, dplyr::desc(mean_abundance))
      } else {
        taxa_summary <- taxa_summary %>%
          dplyr::arrange(dplyr::desc(mean_abundance))
      }
      
      taxa_choices <- as.character(taxa_summary$Taxa)
      taxa_choices <- taxa_choices[!is.na(taxa_choices) & nzchar(taxa_choices)]
      
      current_selected <- input$taxa_selected
      if (is.null(current_selected)) current_selected <- character(0)
      remembered_selected <- taxa_selected_memory()
      candidate_selected <- unique(c(as.character(current_selected), as.character(remembered_selected)))

      default_selected <- if (length(candidate_selected) > 0) {
        intersect(candidate_selected, taxa_choices)
      } else {
        head(taxa_choices, 1)
      }
      if (length(default_selected) == 0 && length(taxa_choices) > 0) {
        default_selected <- head(taxa_choices, 1)
      }
      if (length(default_selected) > 0) {
        taxa_selected_memory(default_selected)
      }
      
      updateSelectizeInput(
        session,
        "taxa_selected",
        choices = taxa_choices,
        selected = default_selected,
        server = TRUE
      )
    })

    observeEvent(input$taxa_selected, {
      selected_now <- input$taxa_selected
      if (!is.null(selected_now) && length(selected_now) > 0) {
        taxa_selected_memory(as.character(selected_now))
      }
    }, ignoreInit = FALSE)
    
    taxa_plot_data <- reactive({
      req(isTRUE(is_active_tab()))
      req(taxa_long_data(), input$group_var)
      df <- taxa_long_data()
      selected_taxa <- if (is.null(input$taxa_selected) || length(input$taxa_selected) == 0) {
        character(0)
      } else {
        input$taxa_selected
      }
      
      validate(need(length(selected_taxa) > 0, ""))
      
      df_sub <- df[df$Taxa %in% selected_taxa, , drop = FALSE]
      
      selected_order <- unique(as.character(selected_taxa))
      available_taxa <- unique(as.character(df_sub$Taxa))
      taxa_order <- selected_order[selected_order %in% available_taxa]
      if (length(taxa_order) == 0) {
        taxa_order <- available_taxa
      }
      df_sub$Taxa <- factor(df_sub$Taxa, levels = taxa_order)
      
      df_sub
    })

    observeEvent(taxa_plot_data(), {
      req(taxa_plot_data())
      df <- taxa_plot_data()
      if (nrow(df) == 0) {
        return()
      }
      n_bars_per_facet <- dplyr::n_distinct(df$Group)
      is_secondary <- !is.null(input$secondary_group_var) && input$secondary_group_var != "none"
      n_taxa <- dplyr::n_distinct(df$Taxa)
      manual_layout <- get_manual_facet_layout(is_secondary, n_taxa)
      measured_layout <- measure_facet_layout(df, is_secondary, manual_layout$ncol, manual_layout$nrow)
      n_facet_cols <- measured_layout$n_facet_cols
      n_facet_rows <- measured_layout$n_facet_rows
      dims <- compute_taxa_comparison_auto_dims(n_facet_cols, n_bars_per_facet, n_facet_rows)
      maybe_auto_adjust_taxa_plot_size(dims$width, dims$height)
    }, ignoreInit = FALSE)

    observeEvent(
      list(input$taxa_selected, input$group_var, input$secondary_group_var, input$manual_facet_layout, input$facet_ncol, input$facet_nrow),
      {
        req(taxa_plot_data())
        df <- taxa_plot_data()
        if (nrow(df) == 0) {
          return()
        }
        n_bars_per_facet <- dplyr::n_distinct(df$Group)
        is_secondary <- !is.null(input$secondary_group_var) && input$secondary_group_var != "none"
        n_taxa <- dplyr::n_distinct(df$Taxa)
        manual_layout <- get_manual_facet_layout(is_secondary, n_taxa)
        measured_layout <- measure_facet_layout(df, is_secondary, manual_layout$ncol, manual_layout$nrow)
        n_facet_cols <- measured_layout$n_facet_cols
        n_facet_rows <- measured_layout$n_facet_rows
        dims <- compute_taxa_comparison_auto_dims(n_facet_cols, n_bars_per_facet, n_facet_rows)
        maybe_auto_adjust_taxa_plot_size(dims$width, dims$height, force = TRUE)
      },
      ignoreInit = TRUE
    )
    
    taxa_plot_reactive <- reactive({
      req(isTRUE(is_active_tab()))
      req(taxa_plot_data())
      df <- taxa_plot_data()
      secondary_col <- input$secondary_group_var
      is_secondary <- !is.null(secondary_col) && secondary_col != "none"
      enable_longitudinal <- isTRUE(input$enable_longitudinal)
      subject_col <- resolve_meta_colname(input$subject_id_var, names(df))
      has_subject <- enable_longitudinal && !is.null(subject_col) && nzchar(subject_col) && subject_col %in% names(df)
      abundance_mode <- input$abundance_mode
      if (is.null(abundance_mode) || !abundance_mode %in% c("clr", "relative", "log_tss")) {
        abundance_mode <- "clr"
      }
      y_label <- switch(
        abundance_mode,
        "relative" = "Relative Abundance (%)",
        "log_tss" = "Log TSS (log10(% + 1))",
        "CLR Abundance"
      )
      base_size <- input$base_size
      if (is.null(base_size) || !is.finite(base_size)) {
        base_size <- 11
      }
      group_levels <- levels(factor(df$Group))
      num_groups <- length(group_levels)
      plot_type <- input$plot_type
      if (is.null(plot_type) || !plot_type %in% c("boxplot", "barplot", "scatter")) {
        plot_type <- "boxplot"
      }
      plot_title <- paste0(input$tax_level, "-Level Taxa Comparison")
      palette_key <- input$color_palette
      if (is.null(palette_key) || !palette_key %in% c("set2", "dark2", "paired", "gray")) {
        palette_key <- "set2"
      }
      use_pattern_requested <- isTRUE(input$use_ggpattern)
      use_pattern <- plot_type == "barplot" && use_pattern_requested && requireNamespace("ggpattern", quietly = TRUE)
      n_groups <- max(1, length(group_levels))
      fill_values <- switch(
        palette_key,
        "dark2" = grDevices::hcl.colors(n_groups, palette = "Dark 3"),
        "paired" = grDevices::hcl.colors(n_groups, palette = "Set 3"),
        "gray" = grDevices::gray.colors(n_groups, start = 0.2, end = 0.85),
        grDevices::hcl.colors(n_groups, palette = "Set 2")
      )
      names(fill_values) <- group_levels
      pattern_values <- rep(c("stripe", "stripe", "crosshatch", "none"), length.out = n_groups)
      names(pattern_values) <- group_levels
      pattern_angle_values <- rep(c(45, -45, 0, 0), length.out = n_groups)
      names(pattern_angle_values) <- group_levels
      y_plot_col <- "AbundancePlot"
      y_labels_fn <- scales::label_number(accuracy = 0.1, trim = TRUE)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = Group, y = .data[[y_plot_col]], fill = Group)) +
        ggplot2::scale_y_continuous(
          labels = y_labels_fn,
          expand = ggplot2::expansion(mult = c(0, 0.18))
        ) +
        ggplot2::theme_bw(base_size = base_size) +
        ggplot2::labs(title = plot_title,
                      x = if (is_secondary) secondary_col else input$group_var, y = y_label) +
        ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = base_size + 3, face = "bold"),
          axis.title.x = ggplot2::element_text(face = "bold", size = base_size + 1),
          axis.title.y = ggplot2::element_text(face = "bold", size = base_size + 1),
          axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = max(6, base_size - 1)),
          legend.position = "none",
          strip.background = ggplot2::element_rect(fill = "grey85", color = "grey20"),
          plot.margin = ggplot2::margin(t = 26, r = 8, b = 8, l = 8)
        ) +
        ggplot2::coord_cartesian(clip = "off")

      if (plot_type == "barplot") {
        facet_group_cols <- c("Taxa")
        if (is_secondary) {
          facet_group_cols <- c(facet_group_cols, "PrimaryGroup")
        }
        # facet_grid(..., scales = "free_y") shares y-scale by row when secondary grouping is used.
        # Use row-level minima so bar baselines align with the actual shared y-scale domain.
        scale_group_cols <- if (is_secondary) c("PrimaryGroup") else c("Taxa")
        summary_group_cols <- c(facet_group_cols, "Group")
        summary_df <- df %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(summary_group_cols))) %>%
          dplyr::summarise(
            mean_ab = mean(.data[[y_plot_col]], na.rm = TRUE),
            se_ab = stats::sd(.data[[y_plot_col]], na.rm = TRUE) / sqrt(dplyr::n()),
            .groups = "drop"
          )
        scale_baseline_df <- df %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(scale_group_cols))) %>%
          dplyr::summarise(
            scale_min = {
              vals <- .data[[y_plot_col]]
              vals <- vals[is.finite(vals)]
              if (length(vals) == 0) 0 else min(vals)
            },
            .groups = "drop"
          )
        summary_df <- dplyr::left_join(summary_df, scale_baseline_df, by = scale_group_cols)
        summary_df$mean_ab[!is.finite(summary_df$mean_ab)] <- NA_real_
        summary_df$se_ab[!is.finite(summary_df$se_ab)] <- 0
        summary_df$scale_min[!is.finite(summary_df$scale_min)] <- 0
        summary_df <- summary_df[is.finite(summary_df$mean_ab), , drop = FALSE]
        if (nrow(summary_df) > 0) {
          use_pattern_crossbar <- use_pattern &&
            "geom_crossbar_pattern" %in% getNamespaceExports("ggpattern")
          if (isTRUE(use_pattern_crossbar)) {
            p <- p +
              ggpattern::geom_crossbar_pattern(
                data = summary_df,
                mapping = ggplot2::aes(
                  x = Group,
                  y = mean_ab,
                  ymin = scale_min,
                  ymax = mean_ab,
                  fill = Group,
                  pattern = Group,
                  pattern_angle = Group
                ),
                inherit.aes = FALSE,
                alpha = 0.85,
                width = 0.7,
                color = "#333333",
                linewidth = 0.25,
                fatten = 0,
                pattern_fill = "white",
                pattern_colour = "#222222",
                pattern_density = 0.01,
                pattern_spacing = 0.03
              ) +
              ggpattern::scale_pattern_manual(values = pattern_values, drop = FALSE) +
              ggpattern::scale_pattern_angle_manual(values = pattern_angle_values, drop = FALSE) +
              ggplot2::guides(pattern = "none", pattern_angle = "none")
          } else {
            p <- p +
              ggplot2::geom_crossbar(
                data = summary_df,
                mapping = ggplot2::aes(
                  x = Group,
                  y = mean_ab,
                  ymin = scale_min,
                  ymax = mean_ab,
                  fill = Group
                ),
                inherit.aes = FALSE,
                alpha = 0.8,
                width = 0.7,
                color = "#333333",
                linewidth = 0.25,
                fatten = 0
              )
          }
          p <- p +
            ggplot2::geom_errorbar(
              data = summary_df,
              mapping = ggplot2::aes(
                x = Group,
                ymin = pmax(mean_ab - se_ab, scale_min),
                ymax = mean_ab + se_ab
              ),
              inherit.aes = FALSE,
              width = 0.2,
              linewidth = 0.4
            )
        }
        if (has_subject) {
          p <- p + ggplot2::geom_point(alpha = 0.6, size = 1.4)
        } else {
          p <- p + ggplot2::geom_point(
            position = ggplot2::position_jitter(width = 0.12, height = 0),
            alpha = 0.55,
            size = 1.4
          )
        }
      } else if (plot_type == "boxplot") {
        p <- p + ggplot2::geom_boxplot(alpha = 0.6, width = 0.65)
        if (has_subject) {
          p <- p + ggplot2::geom_point(alpha = 0.7, size = 1.5)
        } else {
          p <- p + ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.5)
        }
      } else {
        if (has_subject) {
          p <- p + ggplot2::geom_point(alpha = 0.7, size = 1.6)
        } else {
          p <- p + ggplot2::geom_jitter(width = 0.15, alpha = 0.65, size = 1.6)
        }
      }

      if (is_secondary) {
        p <- p + ggplot2::facet_grid(PrimaryGroup ~ Taxa, scales = "free_y")
      } else {
        n_taxa <- dplyr::n_distinct(df$Taxa)
        manual_layout <- get_manual_facet_layout(is_secondary = FALSE, n_taxa = n_taxa)
        if (!is.null(manual_layout$ncol) || !is.null(manual_layout$nrow)) {
          p <- p + ggplot2::facet_wrap(
            ~ Taxa,
            scales = "free_y",
            ncol = manual_layout$ncol,
            nrow = manual_layout$nrow
          )
        } else {
          p <- p + ggplot2::facet_wrap(~ Taxa, scales = "free_y")
        }
      }

      if (isTRUE(input$show_trend_line)) {
        trend_df <- df
        trend_df$GroupOrder <- as.numeric(factor(trend_df$Group, levels = group_levels))
        trend_df <- trend_df[is.finite(trend_df$GroupOrder) & is.finite(trend_df[[y_plot_col]]), , drop = FALSE]
        if (nrow(trend_df) >= 3) {
          trend_method <- input$trend_line_method
          if (is.null(trend_method) || !trend_method %in% c("spearman", "pearson", "lm")) {
            trend_method <- "spearman"
          }
          smooth_method <- if (identical(trend_method, "spearman")) "loess" else "lm"
          p <- p + ggplot2::geom_smooth(
            data = trend_df,
            mapping = ggplot2::aes(x = GroupOrder, y = .data[[y_plot_col]], group = 1),
            inherit.aes = FALSE,
            method = smooth_method,
            se = TRUE,
            color = "#d62728",
            linewidth = 0.8
          )

          facet_vars <- if (is_secondary) c("PrimaryGroup", "Taxa") else c("Taxa")
          facet_split <- interaction(trend_df[, facet_vars, drop = FALSE], drop = TRUE, lex.order = TRUE)
          trend_facet_data <- split(trend_df, facet_split)
          trend_p_rows <- vector("list", length(trend_facet_data))
          trend_idx <- 0L

          for (sub_data in trend_facet_data) {
            if (nrow(sub_data) < 3) next
            x_vals <- sub_data$GroupOrder
            y_vals <- sub_data[[y_plot_col]]
            ok <- is.finite(x_vals) & is.finite(y_vals)
            if (sum(ok) < 3) next

            p_val <- tryCatch({
              if (identical(trend_method, "spearman")) {
                suppressWarnings(stats::cor.test(x_vals[ok], y_vals[ok], method = "spearman", exact = FALSE)$p.value)
              } else if (identical(trend_method, "pearson")) {
                suppressWarnings(stats::cor.test(x_vals[ok], y_vals[ok], method = "pearson")$p.value)
              } else {
                lm_fit <- stats::lm(y_vals[ok] ~ x_vals[ok])
                lm_coef <- summary(lm_fit)$coefficients
                if (nrow(lm_coef) >= 2 && ncol(lm_coef) >= 4) as.numeric(lm_coef[2, 4]) else NA_real_
              }
            }, error = function(e) NA_real_)
            if (!is.finite(p_val)) next

            x_min <- suppressWarnings(min(x_vals[ok], na.rm = TRUE))
            x_max <- suppressWarnings(max(x_vals[ok], na.rm = TRUE))
            y_min <- suppressWarnings(min(y_vals[ok], na.rm = TRUE))
            y_max <- suppressWarnings(max(y_vals[ok], na.rm = TRUE))
            if (!is.finite(x_min) || !is.finite(x_max)) next
            if (!is.finite(y_min) || !is.finite(y_max)) next
            x_span <- max(1e-9, x_max - x_min)
            y_span <- max(1e-9, y_max - y_min)
            ann <- data.frame(
              x = x_min + 0.05 * x_span,
              y = y_max + 0.12 * y_span,
              p_label = paste0("p-value = ", formatC(p_val, format = "g", digits = 3)),
              stringsAsFactors = FALSE
            )
            for (fv in facet_vars) {
              ann[[fv]] <- sub_data[[fv]][1]
            }
            trend_idx <- trend_idx + 1L
            trend_p_rows[[trend_idx]] <- ann
          }

          trend_p_rows <- trend_p_rows[seq_len(trend_idx)]
          if (length(trend_p_rows) > 0) {
            trend_p_df <- do.call(rbind, trend_p_rows)
            p <- p + ggplot2::geom_text(
              data = trend_p_df,
              mapping = ggplot2::aes(x = x, y = y, label = p_label),
              inherit.aes = FALSE,
              hjust = 0,
              vjust = 0,
              size = max(3.4, base_size * 0.34),
              color = "#111111"
            )
          }
        }
      }

      if (has_subject) {
        df$LongitudinalSubject <- as.character(df[[subject_col]])
        group_levels <- levels(factor(df$Group))
        df$LongitudinalOrder <- as.numeric(factor(df$Group, levels = group_levels))

        line_df <- df %>%
          dplyr::filter(!is.na(LongitudinalSubject) & nzchar(LongitudinalSubject))
        if (is_secondary) {
          line_df <- line_df %>%
            dplyr::group_by(Taxa, PrimaryGroup, LongitudinalSubject) %>%
            dplyr::filter(dplyr::n_distinct(Group) >= 2) %>%
            dplyr::arrange(LongitudinalOrder, .by_group = TRUE) %>%
            dplyr::ungroup()
        } else {
          line_df <- line_df %>%
            dplyr::group_by(Taxa, LongitudinalSubject) %>%
            dplyr::filter(dplyr::n_distinct(Group) >= 2) %>%
            dplyr::arrange(LongitudinalOrder, .by_group = TRUE) %>%
            dplyr::ungroup()
        }

        if (nrow(line_df) > 0) {
          p <- p + ggplot2::geom_line(
            data = line_df,
            mapping = ggplot2::aes(
              x = Group,
              y = .data[[y_plot_col]],
              group = interaction(PrimaryGroup, Taxa, LongitudinalSubject)
            ),
            inherit.aes = FALSE,
            color = "#2F4F4F",
            linewidth = 0.35,
            alpha = 0.5
          )
        }
      }

      if (isTRUE(input$show_p_val) && num_groups >= 2) {
        all_comps <- utils::combn(group_levels, 2, simplify = FALSE)
        sig_only <- isTRUE(input$show_sig_only)
        sig_as_marks <- isTRUE(input$show_sig_as_marks)
        paired_mode <- has_subject
        p_adjust_method <- input$p_adjust_method
        if (is.null(p_adjust_method) || !p_adjust_method %in% c("none", "holm", "BH", "bonferroni")) {
          p_adjust_method <- "none"
        }
        pairwise_vjust <- if (sig_as_marks) 0.25 else -0.15

        compute_pairwise_p <- function(sub_data, comp_levels) {
          if (!all(comp_levels %in% unique(as.character(sub_data$Group)))) {
            return(NA_real_)
          }

          if (!paired_mode) {
            d1 <- sub_data$AbundancePlot[sub_data$Group == comp_levels[1]]
            d2 <- sub_data$AbundancePlot[sub_data$Group == comp_levels[2]]
            if (length(d1) <= 1 || length(d2) <= 1) {
              return(NA_real_)
            }
            return(
              tryCatch({
                if (input$stat_method == "wilcox.test") {
                  stats::wilcox.test(d1, d2)$p.value
                } else {
                  stats::t.test(d1, d2)$p.value
                }
              }, error = function(e) NA_real_)
            )
          }

          validate_subject <- !is.na(sub_data$LongitudinalSubject) & nzchar(sub_data$LongitudinalSubject)
          sub_data <- sub_data[validate_subject, , drop = FALSE]
          if (nrow(sub_data) == 0) {
            return(NA_real_)
          }

          g1 <- sub_data[sub_data$Group == comp_levels[1], c("LongitudinalSubject", "AbundancePlot"), drop = FALSE]
          g2 <- sub_data[sub_data$Group == comp_levels[2], c("LongitudinalSubject", "AbundancePlot"), drop = FALSE]
          if (nrow(g1) == 0 || nrow(g2) == 0) {
            return(NA_real_)
          }

          g1 <- g1 %>%
            dplyr::group_by(LongitudinalSubject) %>%
            dplyr::summarise(v1 = stats::median(AbundancePlot, na.rm = TRUE), .groups = "drop")
          g2 <- g2 %>%
            dplyr::group_by(LongitudinalSubject) %>%
            dplyr::summarise(v2 = stats::median(AbundancePlot, na.rm = TRUE), .groups = "drop")

          paired_df <- dplyr::inner_join(g1, g2, by = "LongitudinalSubject")
          paired_df <- paired_df[is.finite(paired_df$v1) & is.finite(paired_df$v2), , drop = FALSE]
          if (nrow(paired_df) <= 1) {
            return(NA_real_)
          }

          tryCatch({
            if (input$stat_method == "wilcox.test") {
              stats::wilcox.test(paired_df$v1, paired_df$v2, paired = TRUE)$p.value
            } else {
              stats::t.test(paired_df$v1, paired_df$v2, paired = TRUE)$p.value
            }
          }, error = function(e) NA_real_)
        }

        format_p_value <- function(pv) {
          if (!is.finite(pv)) return(NA_character_)
          if (pv < 0.001) return("< 0.001")
          sprintf("%.3f", pv)
        }
        signif_mark <- function(pv) {
          if (!is.finite(pv)) return("ns")
          if (pv < 0.001) return("***")
          if (pv < 0.01) return("**")
          if (pv < 0.05) return("*")
          "ns"
        }

        facet_vars <- if (is_secondary) c("PrimaryGroup", "Taxa") else c("Taxa")
        facet_split <- interaction(df[, facet_vars, drop = FALSE], drop = TRUE, lex.order = TRUE)
        facet_data <- split(df, facet_split)
        annotation_rows <- vector("list", length(facet_data))
        row_idx <- 0L

        for (sub_data in facet_data) {
          if (nrow(sub_data) == 0) next
          raw_p <- vapply(all_comps, function(this_comp) {
            compute_pairwise_p(sub_data, this_comp)
          }, numeric(1))

          adj_p <- raw_p
          if (p_adjust_method != "none") {
            ok <- is.finite(adj_p)
            if (any(ok)) {
              adj_p[ok] <- stats::p.adjust(adj_p[ok], method = p_adjust_method)
            }
          }

          valid_idx <- which(is.finite(adj_p))
          if (sig_only) valid_idx <- valid_idx[adj_p[valid_idx] < 0.05]
          if (length(valid_idx) == 0) next

          y_max <- suppressWarnings(max(sub_data[[y_plot_col]], na.rm = TRUE))
          if (!is.finite(y_max) || y_max <= 0) y_max <- 1
          base_mult <- if (sig_as_marks) 1.12 else 1.20
          step_mult <- if (sig_as_marks) 0.08 else 0.12

          ann_df <- data.frame(
            group1 = vapply(all_comps[valid_idx], `[[`, character(1), 1),
            group2 = vapply(all_comps[valid_idx], `[[`, character(1), 2),
            p.adj = adj_p[valid_idx],
            stringsAsFactors = FALSE
          )
          ann_df$p.signif <- vapply(ann_df$p.adj, signif_mark, character(1))
          ann_df$p.format <- vapply(ann_df$p.adj, format_p_value, character(1))
          ann_df$y.position <- y_max * (base_mult + step_mult * (seq_len(nrow(ann_df)) - 1))
          for (fv in facet_vars) {
            ann_df[[fv]] <- sub_data[[fv]][1]
          }

          row_idx <- row_idx + 1L
          annotation_rows[[row_idx]] <- ann_df
        }

        annotation_rows <- annotation_rows[seq_len(row_idx)]
        if (length(annotation_rows) > 0) {
          annotation_df <- do.call(rbind, annotation_rows)
          p <- p + ggpubr::stat_pvalue_manual(
            annotation_df,
            label = if (sig_as_marks) "p.signif" else "p.format",
            xmin = "group1",
            xmax = "group2",
            y.position = "y.position",
            tip.length = 0.01,
            size = max(2, base_size / 3.2),
            bracket.size = 0.3,
            vjust = pairwise_vjust,
            hide.ns = FALSE,
            inherit.aes = FALSE
          )
        }

        # Do not add the overall test label (e.g., "t-test, p=...") to avoid duplicate p-value text.
      }

      p
    })
    
    output$taxa_comparison_plot <- renderPlot({
      ps_now <- ps_obj()
      if (is.null(ps_now)) {
        graphics::plot.new()
        graphics::text(
          0.5, 0.5,
          "Applying selected samples and preparing taxa comparison.\nPlease wait...",
          cex = 0.85
        )
        return(invisible(NULL))
      }
      if (phyloseq::nsamples(ps_now) == 0) {
        graphics::plot.new()
        graphics::text(
          0.5, 0.5,
          "No samples are currently selected.\nPlease select at least one sample in Preprocessing.",
          cex = 0.85
        )
        return(invisible(NULL))
      }
      tryCatch(
        taxa_plot_reactive(),
        error = function(e) {
          graphics::plot.new()
          graphics::text(
            0.5, 0.5,
            paste0("Taxa comparison is not ready yet. Please wait...\n", conditionMessage(e)),
            cex = 0.85
          )
        }
      )
    }, height = function() { req(input$plot_height); input$plot_height },
       width = function() { req(input$plot_width); input$plot_width })

    output$comparison_status <- renderText({
      req(taxa_plot_data(), input$group_var, input$stat_method)
      df <- taxa_plot_data()
      secondary_col <- input$secondary_group_var
      is_secondary <- !is.null(secondary_col) && secondary_col != "none"
      subject_col <- resolve_meta_colname(input$subject_id_var, names(df))
      has_subject <- isTRUE(input$enable_longitudinal) &&
        !is.null(subject_col) &&
        nzchar(subject_col) &&
        subject_col %in% names(df)

      group_levels <- levels(factor(df$Group))
      n_groups <- length(group_levels)
      plot_type <- input$plot_type
      if (is.null(plot_type) || !plot_type %in% c("boxplot", "barplot", "scatter")) {
        plot_type <- "boxplot"
      }
      p_adjust_method <- input$p_adjust_method
      if (is.null(p_adjust_method) || !p_adjust_method %in% c("none", "holm", "BH", "bonferroni")) {
        p_adjust_method <- "none"
      }
      use_pattern_requested <- isTRUE(input$use_ggpattern)
      pattern_status <- if (plot_type == "barplot") {
        if (use_pattern_requested) {
          if (requireNamespace("ggpattern", quietly = TRUE)) {
            "On (stripe 45°, stripe -45°, crosshatch, solid; density = 0.01)"
          } else {
            "Requested, but ggpattern is not installed (fallback to standard fill)"
          }
        } else {
          "Off"
        }
      } else {
        "Not used"
      }
      method_label <- if (identical(input$stat_method, "t.test")) {
        "T-test"
      } else {
        "Wilcoxon test"
      }

      overall_label <- if (!has_subject) {
        if (n_groups > 2) {
          "Kruskal-Wallis (overall)"
        } else {
          paste0(method_label, " (overall)")
        }
      } else if (n_groups == 2) {
        paste0(method_label, " with paired = TRUE (overall)")
      } else {
        "Not shown in paired mode for >2 groups"
      }

      paste(
        c(
          paste0("Primary grouping variable: ", input$group_var),
          paste0("Secondary grouping variable: ", if (is_secondary) secondary_col else "(None)"),
          paste0("Plot type: ", if (plot_type == "barplot") "Bar plot" else if (plot_type == "scatter") "Scatter plot" else "Box plot"),
          paste0("ggpattern: ", pattern_status),
          paste0("Pairing condition axis: ", if (is_secondary) "Secondary group" else "Primary group"),
          paste0("Selected taxa count: ", dplyr::n_distinct(df$Taxa)),
          paste0("Group levels compared: ", paste(group_levels, collapse = ", ")),
          paste0("Within-subject pairing: ", if (isTRUE(input$enable_longitudinal)) "On" else "Off"),
          paste0("Subject ID variable: ", if (has_subject) input$subject_id_var else "Not applied"),
          paste0("Comparison mode: ", if (has_subject) "Paired / repeated" else "Independent groups"),
          paste0("Pairwise test: ", method_label, if (has_subject) " with paired = TRUE" else ""),
          paste0("Pairwise correction: ", p_adjust_method),
          paste0("Overall test: ", overall_label)
        ),
        collapse = "\n"
      )
    })

    output$comparison_status_box <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 400 else input$plot_width
      tags$div(
        class = "taxa-result-card",
        style = paste(sprintf("width: %spx;", box_width), "max-width: 100%; box-sizing: border-box;"),
        verbatimTextOutput(session$ns("comparison_status"))
      )
    })
    
    output$download_taxa_plot <- downloadHandler(
      filename = function() { paste0("taxa_plot_", Sys.Date(), ".png") },
      content = function(file) {
        ggplot2::ggsave(file, plot = taxa_plot_reactive(), device = "png",
                        width = input$plot_width/100, height = input$plot_height/100)
      }
    )
    
    output$download_taxa_data <- downloadHandler(
      filename = function() { paste0("taxa_data_", Sys.Date(), ".tsv") },
      content = function(file) {
        plot_df <- taxa_plot_data()
        secondary_col <- input$secondary_group_var
        is_secondary <- !is.null(secondary_col) && secondary_col != "none"
        abundance_col <- if (input$abundance_mode == "relative") {
          "Relative_Abundance_percent"
        } else if (input$abundance_mode == "log_tss") {
          "Log_TSS_log10_percent_plus1"
        } else {
          "CLR_Abundance"
        }

        if (is_secondary) {
          selected_cols <- c("Sample", "PrimaryGroup", "SecondaryGroup", "Taxa", "AbundancePlot")
          output_names <- c("SampleID", input$group_var, secondary_col, "Taxa", abundance_col)

          if (isTRUE(input$enable_longitudinal) && !is.null(input$subject_id_var) &&
              nzchar(input$subject_id_var)) {
            subject_col <- resolve_meta_colname(input$subject_id_var, names(plot_df))
            if (!is.null(subject_col) && nzchar(subject_col) && subject_col %in% names(plot_df)) {
              selected_cols <- c(selected_cols, subject_col)
              output_names <- c(output_names, input$subject_id_var)
            }
          }

          out_df <- plot_df[, selected_cols, drop = FALSE]
          colnames(out_df) <- output_names
        } else {
          out_df <- plot_df[, c("Sample", "Group", "Taxa", "AbundancePlot")]
          colnames(out_df) <- c("SampleID", input$group_var, "Taxa", abundance_col)
        }
        readr::write_tsv(out_df, file)
      }
    )
    
    output$download_raw_matrix <- downloadHandler(
      filename = function() {
        paste0("taxa_raw_matrix_", input$tax_level, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        req(tax_matrices())
        readr::write_tsv(tax_matrices()$raw, file)
      }
    )
    
    output$download_rel_matrix <- downloadHandler(
      filename = function() {
        paste0("taxa_relative_abundance_matrix_", input$tax_level, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        req(tax_matrices())
        readr::write_tsv(tax_matrices()$rel, file)
      }
    )

    output$taxa_legend_box <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 400 else input$plot_width
      tags$div(
        style = paste(
          "margin-top: 12px;",
          sprintf("width: %spx;", box_width),
          "max-width: 100%;",
          "padding: 12px 14px;",
          "border: 1px solid #e5e7eb;",
          "border-left: 4px solid #6b7280;",
          "border-radius: 8px;",
          "background: linear-gradient(180deg, #fcfcfd 0%, #f7f8fa 100%);",
          "box-shadow: 0 1px 2px rgba(0,0,0,0.04);",
          "box-sizing: border-box;"
        ),
        tags$div(
          style = "color: #1f2937; font-size: 12.5px; line-height: 1.55;",
          uiOutput(session$ns("taxa_figure_legend"))
        )
      )
    })

    output$taxa_status_separator <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 400 else input$plot_width
      tags$hr(
        style = paste(
          sprintf("width: %spx;", box_width),
          "max-width: 100%;",
          "margin: 14px 0 12px 0;",
          "border-top: 1px solid #d1d5db;"
        )
      )
    })

    output$taxa_figure_legend <- renderUI({
      req(input$tax_level, input$plot_type, input$abundance_mode)
      tax_level_label <- tolower(input$tax_level)
      figure_title <- paste0(input$tax_level, "-Level Taxa Comparison")
      plot_type_label <- if (identical(input$plot_type, "barplot")) {
        "bar plots summarize group means with standard error bars"
      } else if (identical(input$plot_type, "scatter")) {
        "scatter plots show individual sample points by group"
      } else {
        "box plots summarize group distributions with overlaid sample points"
      }
      abundance_label <- if (identical(input$abundance_mode, "relative")) {
        "relative abundance"
      } else if (identical(input$abundance_mode, "log_tss")) {
        "log-transformed relative abundance"
      } else {
        "centered log-ratio transformed abundance"
      }
      longitudinal_sentence <- if (isTRUE(input$enable_longitudinal) &&
        !is.null(input$subject_id_var) &&
        nzchar(input$subject_id_var)) {
        paste0(
          " Within-subject pairing is applied using ",
          input$subject_id_var,
          " as the subject identifier."
        )
      } else {
        ""
      }
      tags$div(
        tags$div(
          style = "font-weight: 600; margin-bottom: 4px;",
          figure_title
        ),
        tags$div(
          paste0(
            "This figure compares microbial taxa at ",
            tax_level_label,
            " level across groups; ",
            plot_type_label,
            ". The y-axis represents ",
            abundance_label,
            ". Pairwise significance is annotated above each comparison.",
            longitudinal_sentence
          )
        )
      )
    })
  })
}
