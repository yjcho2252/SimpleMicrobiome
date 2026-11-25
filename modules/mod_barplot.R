## UI
mod_barplot_ui <- function(id) {
  ns <- NS(id)
  tax_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
  
  palette_choices <- c(
    "Paired" = "Paired",
    "Custom (HCL 30+ Colors)" = "Custom_HCL",
    "Set1" = "Set1",
    "Set2" = "Set2",
    "Set3" = "Set3",
    "Dark2" = "Dark2",
    "Accent" = "Accent",
    "Pastel1" = "Pastel1",
    "Pastel2" = "Pastel2",
    additional_plot_palettes
  )
  
  tagList(
    tags$style(HTML("
      .well h4 { font-size: 16px; }
      .well h5 { font-size: 13px; }
      .well .control-label { font-size: 12px; }
      .well .checkbox label { font-size: 12px; }
      .well .form-control { font-size: 12px; }
      .well .btn { font-size: 11px; }
      .taxa-barplot-loading {
        display: none;
        margin: 12px 0;
        color: #475569;
      }
      .taxa-barplot-view:has(.plotly.html-widget-output:not(.shiny-output-error):not([data-taxa-ready])) > .taxa-barplot-loading,
      .taxa-barplot-view:has(.plotly.html-widget-output.recalculating:not(.shiny-output-error)) > .taxa-barplot-loading {
        display: block;
      }
    ")),
    tags$style(HTML(sprintf(
      "#%s .shiny-options-group label, #%s .shiny-options-group label { font-size: 12px; }",
      ns("top_n_basis"), ns("abundance_threshold_basis")
    ))),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        h4(icon("chart-bar"), "Taxa Bar Plot"),
        hr(),
        
        uiOutput(ns("local_group_selector")), 
        
        uiOutput(ns("secondary_group_selector")),

        uiOutput(ns("shared_sort_notice")),
        conditionalPanel(
          condition = paste0("input['", ns("plot_mode"), "'] == 'Sample'"),

          selectInput(ns("sort_method"), "Primary Sort Criterion (within Group):",
                      choices = c(
                        "Top Taxa Abundance" = "taxa_top1",
                        "Sample Name" = "sample_name",
                        "Metadata Variable" = "metadata_var"
                      ),
                      selected = "taxa_top1"),

          conditionalPanel(
            condition = paste0("input['", ns("sort_method"), "'] == 'metadata_var'"),
            hr(),
            uiOutput(ns("sort_by_metadata_ui")),
            selectInput(ns("secondary_sort"), "Secondary Sort (after Metadata):",
                        choices = c(
                          "None" = "none",
                          "Top Taxa Abundance" = "taxa_top1",
                          "Sample Name" = "sample_name_asc"
                        ),
                        selected = "taxa_top1")
          )
        ),
        hr(),
        h4(icon("sliders"), "Profile Settings"),
        selectInput(ns("plot_mode"), "Plot Display Mode:",
                    choices = c("Sample Level" = "Sample",
                                "Group Mean" = "Group_Mean"),
                    selected = "Sample"),        
        
        selectInput(ns("tax_level"), "Taxonomic Level:",
                    choices = tax_ranks,
                    selected = "Genus"),        
        selectInput(ns("taxa_display_mode"), "Taxa Display Mode:",
                    choices = c(
                      "Top N Taxa" = "top_n",
                      "Abundance Threshold" = "abundance_threshold"
                    ),
                    selected = "top_n"),
        conditionalPanel(
          condition = paste0("input['", ns("taxa_display_mode"), "'] == 'top_n'"),
          numericInput(ns("top_n_taxa"), "Taxa Numbers to Display:",
                       value = 15, min = 1, max = 50, step = 1),
          radioButtons(ns("top_n_basis"), "Select top taxa from",
                       choices = c("All selected samples" = "overall", "Each group" = "each_group"),
                       selected = "overall")
        ),
        conditionalPanel(
          condition = paste0("input['", ns("taxa_display_mode"), "'] == 'abundance_threshold'"),
          numericInput(ns("min_taxa_abundance_pct"), "Abundance Threshold (%):",
                       value = 1, min = 0, max = 100, step = 0.1),
          radioButtons(ns("abundance_threshold_basis"), "Threshold applied to",
                       choices = c("Mean across selected samples" = "mean",
                                   "Mean in at least one group" = "any_group"),
                       selected = "mean")
        ),
        selectInput(ns("name_display_mode"), "Taxa Rank to Display:",
                    choices = c("Full Hierarchy" = "full",
                                "Current Rank" = "current"),
                    selected = "current"),
        selectInput(ns("color_palette"), "Color Palette:",
                    choices = palette_choices,
                    selected = "Paired"),
        hr(),
        h4(icon("sliders"), "Display Options"),
        selectInput(ns("plot_size_mode"), "Plot Size", choices = c("Auto" = "auto", "Manual" = "manual"), selected = "auto"),
        conditionalPanel(
          "input.plot_size_mode === 'manual'", ns = ns,
          numericInput(ns("plot_width"), "Plot Width (px):",
                       value = 1000, min = 200, step = 50),
          numericInput(ns("plot_height"), "Plot Height (px):",
                       value = 400, min = 300, step = 50)
        ),
        numericInput(ns("base_size"), "Base Font Size:",
                     value = 12, min = 6, max = 30, step = 1),

        
        hr(),
        h5(icon("download"), "Download"),
        tags$div(
          style = "display: flex; gap: 4px; align-items: center; flex-wrap: nowrap;",
          downloadButton(
            ns("download_barplot"),
            "Plot (PNG)",
            style = "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
          ),
          downloadButton(
            ns("download_barplot_matrix"),
            "Matrix (TSV)",
            style = "font-size: 11px; padding: 3px 6px; width: calc(50% - 2px); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;"
          )
        )
      ),
      mainPanel(
        width = 9,
        h4("Taxa Bar Plot"),
        div(
          class = "taxa-barplot-view",
          style = "width: max-content;",
          uiOutput(ns("barplot_title")),
          div(class = "taxa-barplot-loading", role = "status",
              "Preparing the Taxa Bar Plot. Please wait..."),
          div(
            style = "display: flex; align-items: flex-start; gap: 12px;",
            plotly::plotlyOutput(ns("barplot_out"), width = "auto", height = "auto"),
            div(
              id = ns("barplot_taxa_legend"),
              style = "flex: 0 0 auto; align-self: center; min-width: 220px; margin: 4px 0; white-space: nowrap;"
            )
          )
        ),
        uiOutput(ns("barplot_legend_box"))
      )
    )
  )
}

## Server
mod_barplot_server <- function(id, ps_obj, meta_cols, plot_order = function() NULL) {
  moduleServer(id, function(input, output, session) {
    supported_taxonomic_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
    output$shared_sort_notice <- renderUI({
      samples <- plot_order()$samples
      active <- identical(input$plot_mode, "Sample") && length(samples) > 0L
      if (active) helpText("Sample order is set in Preprocessing.")
    })
    ps_barplot <- reactive({
      ps <- ps_obj()
      if (is.null(ps) || phyloseq::nsamples(ps) == 0L) {
        return(ps)
      }
      sample_depths <- phyloseq::sample_sums(ps)
      rarefaction_depth <- 2000L
      validate(need(all(sample_depths > 0), "Deselect samples with zero reads in Preprocessing."))
      counts <- as(phyloseq::otu_table(ps), "matrix")
      if (!phyloseq::taxa_are_rows(ps)) counts <- t(counts)
      sample_depths <- sample_depths[colnames(counts)]

      deep_samples <- which(sample_depths > rarefaction_depth)
      if (length(deep_samples) > 0L) {
        if (nrow(counts) == 1L) {
          counts[1L, deep_samples] <- rarefaction_depth
        } else {
          # vegan expects samples in rows and draws without replacement.
          rarefied <- withr::with_seed(
            42,
            vegan::rrarefy(t(counts[, deep_samples, drop = FALSE]), rarefaction_depth)
          )
          counts[, deep_samples] <- t(rarefied)
        }
      }
      phyloseq::otu_table(ps) <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)
      phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)
    })

    observeEvent(ps_obj(), {
      ps <- ps_obj()
      available_ranks <- tryCatch(phyloseq::rank_names(ps), error = function(e) character(0))
      rank_choices <- supported_taxonomic_ranks[supported_taxonomic_ranks %in% available_ranks]
      current_rank <- isolate(input$tax_level)
      selected_rank <- if (!is.null(current_rank) && current_rank %in% rank_choices) {
        current_rank
      } else if ("Genus" %in% rank_choices) {
        "Genus"
      } else if (length(rank_choices) > 0) {
        tail(rank_choices, 1)
      } else {
        NULL
      }
      updateSelectInput(session, "tax_level", choices = rank_choices, selected = selected_rank)
    }, ignoreInit = FALSE)

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

    auto_plot_width <- function(n_bars) {
      n_bars <- max(1, as.numeric(n_bars))
      per_bar_px <- min(26, max(6, 89.44 / sqrt(n_bars)))
      ceiling(min(max(300, 220 + n_bars * per_bar_px), 2200))
    }

    apply_disambiguated_taxrank <- function(ps, tax_level) {
      if (is.null(ps) || is.null(tax_level) || !nzchar(tax_level)) {
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

      target_raw <- as.character(tt[[tax_level]])
      target_norm <- tolower(trimws(target_raw))
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
    
    output$local_group_selector <- renderUI({
      req(meta_cols())
      group_choices <- meta_cols()
      non_sampleid_choices <- setdiff(group_choices, "SampleID")
      selected_col <- if (length(non_sampleid_choices) > 0) non_sampleid_choices[1] else group_choices[1]
      selectInput(session$ns("group_var"), "Primary Group (Facet):",
                  choices = c("None" = "__no_facet__", group_choices), selected = selected_col)
    })
    
    output$secondary_group_selector <- renderUI({
      req(meta_cols(), input$group_var)
      resolved_primary <- resolve_meta_colname(input$group_var, meta_cols())
      group_choices <- setdiff(meta_cols(), c(input$group_var, resolved_primary))
      current <- isolate(input$secondary_var)
      selected <- if (!is.null(current) && current %in% group_choices) current else "None"
      selectInput(session$ns("secondary_var"), "Secondary Group (Optional):",
                  choices = c("None" = "None", group_choices), 
                  selected = selected)
    })
    
    group_var <- reactive({
      req(input$group_var)
      if (identical(input$group_var, "__no_facet__")) return(NULL)
      resolve_meta_colname(input$group_var, meta_cols())
    })
    
    secondary_var <- reactive({
      req(input$secondary_var)
      if (input$secondary_var == "None") {
        return(NULL)
      } else {
        return(resolve_meta_colname(input$secondary_var, meta_cols()))
      }
    })

    sort_metadata_var_resolved <- reactive({
      req(input$sort_metadata_var)
      resolve_meta_colname(input$sort_metadata_var, meta_cols())
    })
    
    output$sort_by_metadata_ui <- renderUI({
      req(meta_cols(), input$group_var)
      selectInput(session$ns("sort_metadata_var"), "Metadata Variable for Primary Sort:",
                  choices = meta_cols(),
                  selected = group_var())
    })
    
    generate_taxa_colors <- function(n) {
      if (n <= 1) return(character(0))
      hues <- seq(15, 375, length.out = n + 1)[1:n]
      colors <- grDevices::hcl(h = hues, c = 100, l = 65)
      return(colors)
    }

    build_disambiguated_rank_labels <- function(df_taxa, current_rank, rank_order) {
      if (!current_rank %in% colnames(df_taxa)) {
        return(rep("Unassigned", nrow(df_taxa)))
      }
      rank_prefix_map <- c(
        Kingdom = "k__",
        Phylum = "p__",
        Class = "c__",
        Order = "o__",
        Family = "f__",
        Genus = "g__",
        Species = "s__",
        Strain = "t__"
      )
      get_rank_unassigned_label <- function(rank_name) {
        prefix <- rank_prefix_map[[rank_name]]
        if (is.null(prefix) || !nzchar(prefix)) {
          prefix <- paste0(substr(tolower(rank_name), 1, 1), "__")
        }
        paste0(prefix, "Unassigned")
      }

      is_unassigned_val <- function(x) {
        x <- as.character(x)
        x_trim <- trimws(x)
        is.na(x_trim) | !nzchar(x_trim) | grepl("unassigned|uncultured", x_trim, ignore.case = TRUE)
      }

      base_label <- as.character(df_taxa[[current_rank]])
      out <- base_label

      # For unassigned values at current rank, keep only the nearest parent rank label.
      parent_ranks <- rank_order[rank_order %in% colnames(df_taxa)]
      current_idx <- match(current_rank, parent_ranks)
      if (!is.na(current_idx) && current_idx > 1) {
        parent_ranks_above <- parent_ranks[seq_len(current_idx - 1)]
        is_unassigned_idx <- is_unassigned_val(base_label)
        if (any(is_unassigned_idx) && length(parent_ranks_above) > 0) {
          parent_vals <- rep("UnclassifiedParent", nrow(df_taxa))
          for (parent_rank in rev(parent_ranks_above)) {
            candidate_vals <- as.character(df_taxa[[parent_rank]])
            candidate_placeholder_idx <- is_unassigned_val(candidate_vals)
            use_idx <- is_unassigned_idx & parent_vals == "UnclassifiedParent" & !candidate_placeholder_idx
            parent_vals[use_idx] <- candidate_vals[use_idx]
          }
          out[is_unassigned_idx] <- parent_vals[is_unassigned_idx]
        }
      }

      # Remove "__Unassigned" placeholders for display
      out <- gsub("__Unassigned", "", out, fixed = TRUE)
      out[is.na(out) | !nzchar(out)] <- get_rank_unassigned_label(current_rank)

      # Handle duplicates with parent rank disambiguation
      parent_ranks <- rank_order[rank_order %in% colnames(df_taxa)]
      current_idx <- match(current_rank, parent_ranks)
      if (!is.na(current_idx) && current_idx > 1) {
        parent_ranks <- parent_ranks[seq_len(current_idx - 1)]
      } else {
        parent_ranks <- character(0)
      }

      # Default label scope: include parent ranks only up to Family (not broader than Order/Class/Phylum).
      if ("Family" %in% parent_ranks) {
        fam_idx <- match("Family", parent_ranks)
        parent_ranks <- parent_ranks[fam_idx:length(parent_ranks)]
      }

      dup_vals <- unique(out[duplicated(out) | duplicated(out, fromLast = TRUE)])
      if (length(dup_vals) > 0 && length(parent_ranks) > 0) {
        for (val in dup_vals) {
          idx <- which(out == val)
          sub_df <- df_taxa[idx, , drop = FALSE]
          suffix <- rep("", length(idx))

          for (r in rev(parent_ranks)) {
            rv <- as.character(sub_df[[r]])
            rv[is.na(rv) | !nzchar(rv)] <- get_rank_unassigned_label(r)
            if (length(unique(rv)) > 1) {
              suffix <- ifelse(nzchar(suffix), paste0(rv, ";", suffix), rv)
            }
          }

          if (any(nzchar(suffix))) {
            out[idx] <- paste0(val, " [", suffix, "]")
          }
        }
      }
      out
    }

    resolve_taxa_selection <- function(ps_glom, display_mode, top_n_taxa, min_taxa_abundance_pct,
                                       threshold_basis, top_n_basis, sample_data, group_columns) {
      taxa_ids <- names(sort(phyloseq::taxa_sums(ps_glom), decreasing = TRUE))
      if (is.null(taxa_ids) || length(taxa_ids) == 0) {
        return(character(0))
      }

      display_mode <- if (is.null(display_mode)) "top_n" else display_mode
      by_group <- if (identical(display_mode, "abundance_threshold")) {
        identical(threshold_basis, "any_group")
      } else {
        identical(top_n_basis, "each_group")
      }
      if (by_group) {
        group_means <- sample_data %>%
          dplyr::group_by(OTU, dplyr::across(dplyr::all_of(group_columns))) %>%
          dplyr::summarise(mean_abundance = mean(Abundance), .groups = "drop")
      }
      if (identical(display_mode, "abundance_threshold")) {
        threshold <- suppressWarnings(as.numeric(min_taxa_abundance_pct))
        if (length(threshold) != 1L || !is.finite(threshold) || threshold < 0) {
          threshold <- 0
        }
        threshold <- threshold / 100
        selected_ids <- if (by_group) {
          unique(group_means$OTU[group_means$mean_abundance >= threshold])
        } else {
          abundance <- phyloseq::taxa_sums(ps_glom) / phyloseq::nsamples(ps_glom)
          names(abundance)[abundance >= threshold]
        }
        if (length(selected_ids) == 0) {
          selected_ids <- taxa_ids[1]
        }
        return(selected_ids)
      }

      top_n <- suppressWarnings(as.integer(top_n_taxa))
      if (length(top_n) != 1L || !is.finite(top_n) || top_n < 1) {
        top_n <- 1L
      }
      if (by_group) {
        return(group_means %>%
          dplyr::filter(mean_abundance > 0) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(group_columns))) %>%
          dplyr::slice_max(mean_abundance, n = top_n, with_ties = FALSE) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(OTU) %>%
          dplyr::pull(OTU))
      }
      taxa_ids[seq_len(min(top_n, length(taxa_ids)))]
    }

    barplot_rank_data <- reactive({
      current_rank <- input$tax_level
      req(current_rank, ps_barplot())
      ps_rel <- phyloseq::transform_sample_counts(ps_barplot(), function(x) {
        total <- sum(x, na.rm = TRUE)
        if (!is.finite(total) || total <= 0) x else x / total
      })
      validate(
        need(
          current_rank %in% phyloseq::rank_names(ps_rel),
          paste0("Taxonomic rank '", current_rank, "' is not available in the loaded dataset.")
        )
      )
      ps_rel <- apply_disambiguated_taxrank(ps_rel, current_rank)
      ps_glom <- phyloseq::tax_glom(ps_rel, taxrank = current_rank)

      tax_df <- as.data.frame(phyloseq::tax_table(ps_glom))
      tax_ranks <- phyloseq::rank_names(ps_glom)
      rank_prefix_map <- c(
        Kingdom = "k__",
        Phylum = "p__",
        Class = "c__",
        Order = "o__",
        Family = "f__",
        Genus = "g__",
        Species = "s__",
        Strain = "t__"
      )
      tax_df$Full_Taxa_Name <- apply(tax_df, 1, function(taxa) {
        valid_ranks <- taxa[tax_ranks[tax_ranks %in% names(taxa)]]
        full_name_parts <- c()
        for (i in 1:length(valid_ranks)) {
          part <- valid_ranks[i]
          if (is.na(part) || part == "") {
            rk <- names(valid_ranks)[i]
            prefix <- rank_prefix_map[[rk]]
            if (is.null(prefix) || !nzchar(prefix)) {
              prefix <- paste0(substr(tolower(rk), 1, 1), "__")
            }
            part <- paste0(prefix, "Unassigned")
          }
          full_name_parts <- c(full_name_parts, part)
        }
        paste(full_name_parts, collapse = ";")
      })
      phyloseq::tax_table(ps_glom) <- phyloseq::tax_table(as.matrix(tax_df))

      df_sample <- phyloseq::psmelt(ps_glom)
      meta_df_base <- as.data.frame(phyloseq::sample_data(ps_rel), stringsAsFactors = FALSE)
      meta_cols <- colnames(meta_df_base)

      for (meta_col in meta_cols) {
        syntactic_col <- make.names(meta_col)
        if (!meta_col %in% names(df_sample) && syntactic_col %in% names(df_sample)) {
          names(df_sample)[names(df_sample) == syntactic_col] <- meta_col
        }
      }
      list(df_sample = df_sample, metadata = meta_df_base, ps_glom = ps_glom, tax_ranks = tax_ranks)
    })

    barplot_data <- reactive({
      req(ps_obj(), input$group_var, input$tax_level, input$name_display_mode,
          input$plot_mode, input$taxa_display_mode)

      current_rank <- input$tax_level
      primary_var <- group_var()
      secondary_var_val <- secondary_var()
      plot_mode <- input$plot_mode
      taxa_display_mode <- input$taxa_display_mode
      min_taxa_abundance_pct <- input$min_taxa_abundance_pct

      ranked <- barplot_rank_data()
      df_sample <- ranked$df_sample
      meta_df_base <- ranked$metadata
      ps_glom <- ranked$ps_glom
      tax_ranks <- ranked$tax_ranks
      primary_var <- resolve_meta_colname(primary_var, names(df_sample))
      secondary_var_val <- resolve_meta_colname(secondary_var_val, names(df_sample))
      has_primary <- !identical(input$group_var, "__no_facet__")
      if (!has_primary) {
        primary_var <- tail(make.unique(c(union(names(df_sample), names(meta_df_base)), ".barplot_group")), 1L)
        df_sample[[primary_var]] <- "All samples"
        meta_df_base[[primary_var]] <- "All samples"
      }

      if (is.null(primary_var) || !primary_var %in% names(df_sample)) {
        showNotification(
          paste0("Error: Primary grouping variable '", primary_var, "' was not found in melted data."),
          type = "error"
        )
        return(NULL)
      }
      order_config <- plot_order()
      df_sample[[primary_var]] <- plot_order_factor(df_sample[[primary_var]], primary_var, order_config)

      if (!is.null(secondary_var_val)) {
        if (!secondary_var_val %in% names(df_sample)) {
          showNotification(
            paste0("Error: Secondary grouping variable '", secondary_var_val, "' was not found in melted data."),
            type = "error"
          )
          return(NULL)
        }
        df_sample[[secondary_var_val]] <- plot_order_factor(df_sample[[secondary_var_val]], secondary_var_val, order_config)
      }

      if (input$name_display_mode == "full") {
        df_sample$Taxa_Name_Display <- df_sample$Full_Taxa_Name
      } else {
        df_sample$Taxa_Name_Display <- build_disambiguated_rank_labels(
          df_taxa = df_sample,
          current_rank = current_rank,
          rank_order = tax_ranks
        )
      }

      selected_otu_ids <- resolve_taxa_selection(
        ps_glom = ps_glom,
        display_mode = taxa_display_mode,
        top_n_taxa = input$top_n_taxa,
        min_taxa_abundance_pct = min_taxa_abundance_pct,
        threshold_basis = input$abundance_threshold_basis,
        top_n_basis = input$top_n_basis,
        sample_data = df_sample,
        group_columns = c(primary_var, secondary_var_val)
      )

      df_sample$Taxa_Group <- ifelse(df_sample$OTU %in% selected_otu_ids,
                                     as.character(df_sample$Taxa_Name_Display),
                                     "Others")

      taxa_mean_abundance <- df_sample %>%
        dplyr::group_by(Taxa_Group) %>%
        dplyr::summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(desc(Mean_Abundance))

      ordered_taxa_levels <- taxa_mean_abundance$Taxa_Group
      if ("Others" %in% ordered_taxa_levels) {
        ordered_taxa_levels <- c(setdiff(ordered_taxa_levels, "Others"), "Others")
      }
      df_sample$Taxa_Group <- factor(df_sample$Taxa_Group, levels = ordered_taxa_levels)

      sample_metadata <- meta_df_base
      if (!is.null(secondary_var_val)) {
        meta_df_base[[primary_var]] <- plot_order_factor(meta_df_base[[primary_var]], primary_var, order_config)
        meta_df_base[[secondary_var_val]] <- plot_order_factor(meta_df_base[[secondary_var_val]], secondary_var_val, order_config)

        meta_df_tibble <- tibble::as_tibble(meta_df_base)
        if (has_primary) {
          unique_groups <- meta_df_tibble %>%
            dplyr::arrange(.data[[primary_var]], .data[[secondary_var_val]]) %>%
            dplyr::mutate(Combined_Group = paste(.data[[primary_var]], .data[[secondary_var_val]], sep = "_")) %>%
            dplyr::select(Combined_Group) %>%
            unique() %>%
            dplyr::pull(Combined_Group)

          df_sample$Combined_Group <- paste(df_sample[[primary_var]], df_sample[[secondary_var_val]], sep = "_")
          df_sample$Combined_Group <- factor(df_sample$Combined_Group, levels = unique_groups)
        } else {
          df_sample$Combined_Group <- df_sample[[secondary_var_val]]
        }
      } else {
        df_sample$Combined_Group <- df_sample[[primary_var]]
      }

      if (plot_mode == "Group_Mean") {
        df_sample_grouped <- df_sample %>%
          dplyr::group_by(Sample, Combined_Group, .data[[primary_var]], Taxa_Group) %>%
          dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

        df_plot <- df_sample_grouped %>%
          dplyr::group_by(Combined_Group, .data[[primary_var]], Taxa_Group) %>%
          dplyr::summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")
      } else {
        df_plot <- df_sample
      }

      list(
        df_plot = df_plot, metadata = sample_metadata, primary_var = primary_var,
        secondary_var = secondary_var_val, taxa_levels = ordered_taxa_levels, has_primary = has_primary
      )
    })

    barplot_matrix_reactive <- reactive({
      data <- barplot_data()
      req(data)
      df_plot <- data$df_plot
      primary_var <- if (data$has_primary) data$primary_var
      secondary_var_val <- data$secondary_var
      if (!is.null(secondary_var_val) && secondary_var_val %in% names(df_plot)) {
        df_plot %>%
          dplyr::select(any_of(c("Sample", "Combined_Group", primary_var, secondary_var_val, "Taxa_Group", "Abundance")))
      } else {
        df_plot %>%
          dplyr::select(any_of(c("Sample", "Combined_Group", primary_var, "Taxa_Group", "Abundance")))
      }
    })

    manual_plot_size <- reactive(identical(input$plot_size_mode, "manual"))
    effective_plot_height <- reactive({
      if (manual_plot_size()) req(input$plot_height) else 400
    })

    effective_plot_width <- reactive({
      if (manual_plot_size()) return(req(input$plot_width))
      data <- barplot_data()
      req(data)
      df_plot <- data$df_plot
      x_var <- if (identical(input$plot_mode, "Group_Mean")) "Combined_Group" else "Sample"
      n_bars <- length(unique(as.character(df_plot[[x_var]])))
      width <- auto_plot_width(n_bars)
      available_width <- input$available_plot_width
      if (length(available_width) == 1L && is.finite(available_width)) {
        width <- min(width, max(150, floor(available_width)))
      }
      width
    })
    
    build_barplot <- function(interactive = FALSE) {
      data <- barplot_data()
      req(data, input$color_palette, input$base_size)
      current_rank <- input$tax_level
      primary_var <- data$primary_var
      secondary_var_val <- data$secondary_var
      ordered_taxa_levels <- data$taxa_levels
      df_plot <- data$df_plot
      order_config <- plot_order()
      plot_mode <- input$plot_mode
      base_size <- input$base_size
      subgroup_label_size <- base_size / 3.6
      subgroup_headroom <- 0
      fill_label <- if (input$name_display_mode == "full") {
        paste(current_rank, "Taxa (Full Hierarchy)")
      } else {
        current_rank
      }
      if (plot_mode == "Group_Mean") {
        x_var <- "Combined_Group"
        x_label <- paste(primary_var, if (!is.null(secondary_var_val)) secondary_var_val else "")
        
        df_plot[[x_var]] <- factor(df_plot[[x_var]])
        
      } else {
        x_var <- "Sample"
        x_label <- ""
        subgroup_rect_df <- NULL
        subgroup_boundary_df <- NULL
        subgroup_label_df <- NULL
        
        if (length(order_config$samples)) {
          ordered_sample_levels <- plot_order_levels(
            sort(unique(as.character(df_plot$Sample))), order_config$samples
          )
        } else {
          req(input$sort_method)
          sort_method <- input$sort_method
          ordered_sample_levels <- c()
          meta_df_base <- data$metadata
          if (!"Sample" %in% names(meta_df_base)) { meta_df_base$Sample <- rownames(meta_df_base) }

          required_vars <- c("Sample", primary_var)
          if (!is.null(secondary_var_val)) { required_vars <- c(required_vars, secondary_var_val) }

          if (sort_method == "metadata_var") {
            req(input$sort_metadata_var)
            sort_meta_var <- sort_metadata_var_resolved()
            if (sort_meta_var %in% names(meta_df_base)) {
              required_vars <- c(required_vars, sort_meta_var)
            } else {
              showNotification(paste0("Error: Metadata variable '", input$sort_metadata_var, "' not found in sample data."), type = "error")
              return(NULL)
            }
          }

          valid_vars <- required_vars[required_vars %in% names(meta_df_base)]
          meta_df_temp <- meta_df_base[, valid_vars, drop = FALSE]
          meta_df <- unique(meta_df_temp)
          meta_df <- tibble::as_tibble(meta_df)
          req(nrow(meta_df) > 0)
          for (variable in intersect(required_vars[-1L], names(order_config$groups))) {
            meta_df[[variable]] <- plot_order_factor(meta_df[[variable]], variable, order_config)
          }

          sort_args <- list(dplyr::sym(primary_var))
          if (!is.null(secondary_var_val)) {
            sort_args <- c(sort_args, dplyr::sym(secondary_var_val))
          }

          if (sort_method == "taxa_top1") {
            main_top_taxa_name <- as.character(ordered_taxa_levels[1])
            df_sort <- df_plot %>%
              dplyr::filter(Taxa_Group == main_top_taxa_name) %>%
              dplyr::group_by_at(c("Sample", primary_var, if(!is.null(secondary_var_val)) secondary_var_val)) %>%
              dplyr::summarise(Sort_Abundance = sum(Abundance, na.rm = TRUE), .groups = 'drop')

            sort_args_taxa <- c(sort_args, dplyr::quo(dplyr::desc(!!dplyr::sym("Sort_Abundance"))), dplyr::sym("Sample"))
            df_sort_ordered <- df_sort %>% dplyr::arrange(!!!sort_args_taxa)
            ordered_sample_levels <- unique(df_sort_ordered$Sample)

          } else if (sort_method == "sample_name") {
            sort_args_name <- c(sort_args, dplyr::sym("Sample"))
            ordered_sample_levels <- meta_df %>% dplyr::arrange(!!!sort_args_name) %>% dplyr::pull(Sample)

          } else if (sort_method == "metadata_var") {
            req(input$secondary_sort)
            sort_meta_var <- sort_metadata_var_resolved(); secondary_sort_crit <- input$secondary_sort; meta_df_sortable <- meta_df

            if (secondary_sort_crit == "taxa_top1") {
              main_top_taxa_name <- as.character(ordered_taxa_levels[1])
              taxa_abun_data <- df_plot %>% dplyr::filter(Taxa_Group == main_top_taxa_name) %>% dplyr::group_by(Sample) %>% dplyr::summarise(Sort_Abundance = sum(Abundance, na.rm = TRUE), .groups = 'drop')
              meta_df_sortable <- meta_df_sortable %>% dplyr::left_join(taxa_abun_data, by = "Sample")
            }

            sort_args_meta <- c(sort_args, dplyr::sym(sort_meta_var))
            if (secondary_sort_crit == "taxa_top1") {
              sort_args_meta <- c(sort_args_meta, dplyr::quo(dplyr::desc(!!dplyr::sym("Sort_Abundance"))))
            } else if (secondary_sort_crit == "sample_name_asc") {
              sort_args_meta <- c(sort_args_meta, dplyr::sym("Sample"))
            }

            meta_df_sorted <- meta_df_sortable %>% dplyr::arrange(!!!sort_args_meta)
            ordered_sample_levels <- meta_df_sorted %>% dplyr::pull(Sample) %>% unique()

          } else {
            sort_args_fallback <- c(sort_args, dplyr::sym("Sample"))
            ordered_sample_levels <- meta_df %>% dplyr::arrange(!!!sort_args_fallback) %>% dplyr::pull(Sample)
          }

        }
        df_plot$Sample <- factor(df_plot$Sample, levels = ordered_sample_levels)

        if (!is.null(secondary_var_val) && secondary_var_val %in% names(df_plot)) {
          sample_group_df <- df_plot %>%
            dplyr::distinct(.data[[primary_var]], Sample, .data[[secondary_var_val]]) %>%
            dplyr::mutate(
              Sample = as.character(Sample),
              sample_rank_global = as.numeric(factor(Sample, levels = levels(df_plot$Sample)))
            ) %>%
            dplyr::group_by(.data[[primary_var]]) %>%
            dplyr::arrange(sample_rank_global, .by_group = TRUE) %>%
            dplyr::mutate(
              x_pos = dplyr::row_number(),
              subgroup_block = dplyr::consecutive_id(.data[[secondary_var_val]])
            ) %>%
            dplyr::ungroup()

          subgroup_blocks <- sample_group_df %>%
            dplyr::mutate(subgroup_label = as.character(.data[[secondary_var_val]])) %>%
            dplyr::group_by(.data[[primary_var]], subgroup_label, subgroup_block) %>%
            dplyr::summarise(
              xmin = min(x_pos) - 0.5,
              xmax = max(x_pos) + 0.5,
              .groups = "drop"
            ) %>%
            dplyr::group_by(.data[[primary_var]]) %>%
            dplyr::arrange(xmin, .by_group = TRUE) %>%
            dplyr::ungroup()

          subgroup_rect_df <- subgroup_blocks %>%
            dplyr::group_by(.data[[primary_var]]) %>%
            dplyr::mutate(block_order = dplyr::row_number()) %>%
            dplyr::filter(block_order %% 2 == 1) %>%
            dplyr::ungroup()

          subgroup_boundary_df <- subgroup_blocks %>%
            dplyr::group_by(.data[[primary_var]]) %>%
            dplyr::arrange(xmin, .by_group = TRUE) %>%
            dplyr::mutate(xintercept = xmax) %>%
            dplyr::filter(dplyr::row_number() < dplyr::n()) %>%
            dplyr::ungroup()

          label_height_px <- subgroup_label_size * 96 / 25.4
          label_band_px <- 1.4 * label_height_px + 6
          panel_height_px <- max(100, effective_plot_height() - 120)
          subgroup_headroom <- label_band_px / (panel_height_px - label_band_px)
          subgroup_label_df <- subgroup_blocks %>%
            dplyr::mutate(
              xmid = (xmin + xmax) / 2,
              y = 1 + subgroup_headroom / 2
            )
        }
      }
      
      palette_choice <- input$color_palette
      num_top_taxa <- length(setdiff(ordered_taxa_levels, "Others"))
      top_taxa_names <- setdiff(ordered_taxa_levels, "Others")
      
      if (palette_choice %in% additional_plot_palettes) {
        taxa_base_colors <- group_plot_colors(num_top_taxa, palette_choice)
      } else if (palette_choice == "Custom_HCL") {
        taxa_base_colors <- generate_taxa_colors(num_top_taxa)
      } else {
        tryCatch({
          pal_info <- RColorBrewer::brewer.pal.info[palette_choice, ]
          max_colors <- pal_info$maxcolors
          n_colors_to_use <- min(num_top_taxa, max_colors)
          brewer_colors <- RColorBrewer::brewer.pal(n_colors_to_use, palette_choice)
          taxa_base_colors <- rep(brewer_colors, length.out = num_top_taxa)
        }, error = function(e) {
          showNotification(paste0("Warning: RColorBrewer function failed (", e$message, "). Using HCL fallback."), type = "warning", duration = 3)
          taxa_base_colors <- generate_taxa_colors(num_top_taxa)
        })
      }
      
      names(taxa_base_colors) <- top_taxa_names
      taxa_colors <- c(taxa_base_colors, "Others" = "gray80")
      missing_levels <- setdiff(ordered_taxa_levels, names(taxa_colors))
      if (length(missing_levels) > 0) {
        extra_colors <- rep("gray50", length(missing_levels))
        names(extra_colors) <- missing_levels
        taxa_colors <- c(taxa_colors, extra_colors)
      }
      
      fill_var <- "Taxa_Group"
      bar_mapping <- NULL
      if (interactive) {
        bar_columns <- unique(c(x_var, primary_var))
        df_plot <- df_plot %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(c(bar_columns, fill_var)))) %>%
          dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(bar_columns))) %>%
          dplyr::mutate(Displayed_Abundance = Abundance / sum(Abundance)) %>%
          dplyr::ungroup() %>%
          dplyr::filter(Abundance > 0)
        df_plot$Tooltip <- sprintf(
          "%s: %s<br>Taxon: %s<br>Relative abundance: %.2f%%",
          if (plot_mode == "Sample") "Sample" else "Group",
          htmltools::htmlEscape(as.character(df_plot[[x_var]])),
          htmltools::htmlEscape(as.character(df_plot[[fill_var]])),
          100 * df_plot$Displayed_Abundance
        )
        bar_mapping <- ggplot2::aes(text = Tooltip, group = Taxa_Group)
      }
      
      p <- ggplot2::ggplot(df_plot, ggplot2::aes_string(x = x_var, y = "Abundance", fill = fill_var))
      
      if (!is.null(secondary_var_val) && identical(plot_mode, "Sample") && !is.null(subgroup_rect_df) && nrow(subgroup_rect_df) > 0) {
        p <- p + ggplot2::geom_rect(
          data = subgroup_rect_df,
          ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
          inherit.aes = FALSE,
          fill = "grey70",
          alpha = 0.10
        )
      }
      
      legend_key_size <- grid::unit(6, "pt")
      draw_bar_key <- function(data, params, size) {
        grid::grobTree(
          ggplot2::draw_key_rect(data, params, size),
          vp = grid::viewport(width = legend_key_size, height = legend_key_size)
        )
      }
      p <- p +
        ggplot2::geom_bar(
          mapping = bar_mapping, stat = "identity", position = "fill",
          key_glyph = if (interactive) ggplot2::draw_key_rect else draw_bar_key
        ) +
        ggplot2::scale_fill_manual(values = taxa_colors) + 
        
        if (data$has_primary && (plot_mode == "Sample" || !is.null(secondary_var_val))) {
          ggplot2::facet_grid(
            cols = ggplot2::vars(.data[[primary_var]]),
            scales = "free_x",
            space = "free_x"
          )
        } else {
          NULL 
        }
      
      if (!is.null(secondary_var_val) && identical(plot_mode, "Sample") && !is.null(subgroup_boundary_df) && nrow(subgroup_boundary_df) > 0) {
        p <- p + ggplot2::geom_vline(
          data = subgroup_boundary_df,
          ggplot2::aes(xintercept = xintercept),
          inherit.aes = FALSE,
          linetype = "dotted",
          linewidth = 0.24,
          color = "grey30"
        )
      }

      if (!is.null(secondary_var_val) && identical(plot_mode, "Sample") && !is.null(subgroup_label_df) && nrow(subgroup_label_df) > 0) {
        p <- p + ggplot2::geom_text(
          data = subgroup_label_df,
          ggplot2::aes(x = xmid, y = y, label = subgroup_label),
          inherit.aes = FALSE,
          size = subgroup_label_size,
          fontface = "bold",
          color = "grey20",
          vjust = 0.5
        )
      }
      
      p <- p +
        ggplot2::theme_bw(base_size = base_size) +
        ggplot2::scale_y_continuous(
          expand = if (subgroup_headroom > 0) {
            ggplot2::expansion(add = c(0, subgroup_headroom / 2))
          } else {
            ggplot2::expansion(mult = c(0, 0.03))
          }
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = if (plot_mode == "Group_Mean") base_size else max(6, base_size - 2)
          ),
          axis.title.y = ggplot2::element_text(face = "bold", size = base_size + 1),
          panel.border = ggplot2::element_rect(colour = "#595959", fill = NA, linewidth = 0.5),
          strip.background = ggplot2::element_rect(colour = "#595959", fill = "grey85", linewidth = 0.5),
          strip.text = ggplot2::element_text(face = "bold"),
          panel.grid.major = ggplot2::element_line(color = "white", linewidth = 0.1),
          panel.grid.minor = ggplot2::element_line(color = "white", linewidth = 0.1),
          axis.title.x = ggplot2::element_blank(),
          legend.title = ggplot2::element_text(size = base_size),
          plot.margin = ggplot2::margin(6, 10, 4, 5.5)
        ) +
        ggplot2::labs(
          y = "Relative Abundance",
          fill = fill_label
        )
      if (interactive) {
        p <- p + ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.5))
      } else {
        p <- p + ggplot2::theme(
          legend.key.size = legend_key_size,
          legend.key.spacing.y = grid::unit(0, "pt"),
          legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 3))
        )
      }
      
      return(p)
    }

    barplot_reactive <- reactive(build_barplot())
    barplot_title <- reactive(paste0(req(input$tax_level), "-Level Microbial Community Composition"))

    output$barplot_title <- renderUI({
      req(barplot_data(), input$base_size)
      tags$h5(
        style = sprintf("font: bold %gpt sans-serif; margin: 4px 0 12px 6px; white-space: nowrap;", input$base_size + 3),
        barplot_title()
      )
    })

    output$barplot_out <- plotly::renderPlotly({
      ps_now <- ps_obj()
      validate(
        need(!is.null(ps_now), "Applying selected samples and preparing bar plot."),
        need(phyloseq::nsamples(ps_now) > 0, "Select at least one sample in Preprocessing.")
      )
      plot <- build_barplot(interactive = TRUE)
      plot_width <- effective_plot_width()
      widget <- plotly::ggplotly(
        plot, tooltip = "text", dynamicTicks = FALSE,
        width = plot_width, height = effective_plot_height()
      )
      layout <- widget$x$layout
      axis_names <- grep("^xaxis[0-9]*$", names(layout), value = TRUE)
      if (length(axis_names) > 1L) {
        # ggplotly gives free-width facets equal domains. Match domain width to the discrete axis span.
        axis_names <- axis_names[order(vapply(layout[axis_names], function(axis) axis$domain[1L], numeric(1)))]
        old_domains <- do.call(rbind, lapply(layout[axis_names], function(axis) axis$domain))
        spans <- vapply(layout[axis_names], function(axis) diff(axis$range), numeric(1))
        widths <- sum(old_domains[, 2L] - old_domains[, 1L]) * spans / sum(spans)
        gaps <- old_domains[-1L, 1L] - old_domains[-nrow(old_domains), 2L]
        starts <- old_domains[1L, 1L] + c(0, cumsum(head(widths, -1L) + gaps))
        new_domains <- cbind(starts, starts + widths)
        map_x <- stats::approxfun(c(t(old_domains)), c(t(new_domains)), ties = "ordered", rule = 2)
        for (i in seq_along(axis_names)) layout[[axis_names[i]]]$domain <- unname(new_domains[i, ])
        layout$shapes <- lapply(layout$shapes, function(shape) {
          if (identical(shape$xref, "paper")) {
            shape$x0 <- map_x(shape$x0)
            shape$x1 <- map_x(shape$x1)
          }
          shape
        })
        layout$annotations <- lapply(layout$annotations, function(annotation) {
          if (identical(annotation$xref, "paper") && !identical(annotation$annotationType, "axis")) {
            annotation$x <- map_x(annotation$x)
          }
          annotation
        })
      }
      if (inherits(plot$facet, "FacetGrid") || identical(input$group_var, "__no_facet__")) {
        # ggplotly keeps the strip background but drops its border width.
        strip_indices <- which(vapply(layout$shapes, function(shape) {
          identical(shape$type, "rect") && identical(shape$xref, "paper") &&
            identical(shape$yref, "paper") && identical(shape$ysizemode, "pixel")
        }, logical(1)))
        for (i in strip_indices) {
          strip <- layout$shapes[[i]]
          layout$shapes[[i]]$line <- list(color = "#595959", width = 1)
          label_indices <- which(vapply(layout$annotations, function(annotation) {
            identical(annotation$xref, "paper") && identical(annotation$yref, "paper") &&
              !identical(annotation$annotationType, "axis") &&
              identical(annotation$yanchor, "bottom") &&
              isTRUE(all.equal(annotation$y, strip$yanchor)) &&
              isTRUE(all.equal(annotation$x, mean(c(strip$x0, strip$x1))))
          }, logical(1)))
          for (j in label_indices) {
            layout$annotations[[j]]$yanchor <- "middle"
            layout$annotations[[j]]$yshift <- mean(c(strip$y0, strip$y1))
            if (!startsWith(layout$annotations[[j]]$text, "<b>")) {
              layout$annotations[[j]]$text <- paste0("<b>", layout$annotations[[j]]$text, "</b>")
            }
          }
        }
      }
      group_borders <- lapply(axis_names, function(axis_name) {
        axis <- layout[[axis_name]]
        y_axis_name <- sub("^y", "yaxis", axis$anchor)
        y_domain <- layout[[y_axis_name]]$domain
        list(
          type = "rect", xref = "paper", yref = "paper",
          x0 = axis$domain[1L], x1 = axis$domain[2L],
          y0 = y_domain[1L], y1 = y_domain[2L],
          line = list(color = "#595959", width = 1), layer = "above"
        )
      })
      layout$shapes <- c(layout$shapes, group_borders)
      layout$hovermode <- "closest"
      layout$showlegend <- FALSE
      layout$autosize <- FALSE
      if (input$plot_mode == "Sample" && !is.null(secondary_var())) {
        widget$x$data <- lapply(widget$x$data, function(trace) {
          if (identical(trace$mode, "text")) {
            plain <- !startsWith(trace$text, "<b>")
            trace$text[plain] <- paste0("<b>", trace$text[plain], "</b>")
          }
          trace
        })
      }
      if (!manual_plot_size()) {
        panel_width <- max(1, plot_width - layout$margin$l - layout$margin$r)
        widget$x$data <- lapply(widget$x$data, function(trace) {
          if (!identical(trace$type, "bar")) return(trace)
          axis_name <- sub("^x", "xaxis", trace$xaxis)
          axis <- layout[[axis_name]]
          pixels_per_unit <- panel_width * diff(axis$domain) / diff(axis$range)
          trace$width <- pmin(trace$width, 26 / pixels_per_unit)
          trace
        })
      }
      font_family <- "sans-serif"
      layout$font$family <- font_family
      layout$hoverlabel$font$family <- font_family
      for (axis_name in grep("^[xy]axis[0-9]*$", names(layout), value = TRUE)) {
        layout[[axis_name]]$tickfont$family <- font_family
      }
      layout$annotations <- lapply(layout$annotations, function(annotation) {
        annotation$font$family <- font_family
        annotation
      })
      widget$x$layout <- layout
      htmlwidgets::onRender(widget, "function(el, x, config) {
        el.setAttribute('data-taxa-ready', 'true');
        el.style.width = x.layout.width + 'px';
        el.style.flex = '0 0 auto';
        const legend = document.getElementById(config.id);
        if (!legend) return;
        legend.replaceChildren();
        legend.style.fontFamily = config.fontFamily;
        const title = document.createElement('strong');
        title.textContent = config.title;
        legend.appendChild(title);
        const items = document.createElement('div');
        items.style.cssText = 'display:flex;flex-direction:column;gap:1px;align-items:flex-start;';
        legend.appendChild(items);
        const buttons = new Map();
        let selected = null;
        function toggleTaxon(taxon) {
          selected = selected === taxon ? null : taxon;
          const indices = [];
          const opacity = [];
          el.data.forEach(function(item, index) {
            if (item.type !== 'bar') return;
            indices.push(index);
            opacity.push(selected === null || (item.legendgroup || item.name) === selected ? 1 : 0.2);
          });
          Plotly.restyle(el, {opacity: opacity}, indices);
          buttons.forEach(function(button, name) {
            button.setAttribute('aria-pressed', String(name === selected));
            button.style.opacity = selected === null || name === selected ? '1' : '0.45';
          });
        }
        el.data.forEach(function(trace) {
          if (trace.type !== 'bar') return;
          const taxon = trace.legendgroup || trace.name;
          if (buttons.has(taxon)) return;
          const button = document.createElement('button');
          button.type = 'button';
          button.setAttribute('aria-pressed', 'false');
          button.style.cssText = 'display:inline-flex;align-items:baseline;gap:3px;padding:0 2px;border:0;background:transparent;color:inherit;font:inherit;font-size:11px;line-height:1.15;text-align:left;white-space:nowrap;';
          const swatch = document.createElement('span');
          swatch.setAttribute('aria-hidden', 'true');
          swatch.style.cssText = 'display:inline-block;flex:0 0 8px;height:8px;';
          swatch.style.backgroundColor = Array.isArray(trace.marker.color) ? trace.marker.color[0] : trace.marker.color;
          const label = document.createElement('span');
          label.textContent = taxon;
          button.append(swatch, label);
          button.addEventListener('click', function() { toggleTaxon(taxon); });
          items.appendChild(button);
          buttons.set(taxon, button);
        });
        if (buttons.size > 1 && legend.getBoundingClientRect().height > x.layout.height) {
          items.style.cssText = 'display:grid;grid-auto-flow:column;gap:1px 10px;justify-items:start;';
          items.style.gridTemplateRows = 'repeat(' + Math.ceil(buttons.size / 2) + ', max-content)';
        }
        if (el.taxaWidthObserver) el.taxaWidthObserver.disconnect();
        if (!config.autoWidth) return;
        const canvas = el.parentElement;
        const panel = canvas.parentElement.parentElement;
        function reportWidth() {
          if (!panel.clientWidth) return;
          const panelStyle = getComputedStyle(panel);
          const gap = parseFloat(getComputedStyle(canvas).columnGap) || 0;
          const width = Math.floor(panel.clientWidth - parseFloat(panelStyle.paddingLeft) -
            parseFloat(panelStyle.paddingRight) - legend.getBoundingClientRect().width - gap);
          if (width === el.taxaAvailableWidth) return;
          el.taxaAvailableWidth = width;
          Shiny.setInputValue(config.widthInput, width);
        }
        el.taxaWidthObserver = new ResizeObserver(reportWidth);
        el.taxaWidthObserver.observe(panel);
        el.taxaWidthObserver.observe(legend);
        reportWidth();
      }", data = list(id = session$ns("barplot_taxa_legend"), title = plot$labels$fill,
                      fontFamily = font_family, autoWidth = !manual_plot_size(),
                      widthInput = session$ns("available_plot_width")))
    })
    
    output$download_barplot <- downloadHandler(
      filename = function() { paste0("taxa_barplot_", input$tax_level, "_", input$plot_mode, ".png") },
      content = function(file) {
        plot_width_in <- effective_plot_width() / 72
        plot_height_in <- effective_plot_height() / 72
        
        title_grob <- grid::textGrob(
          barplot_title(), x = grid::unit(6, "pt"), just = "left",
          gp = grid::gpar(fontfamily = "sans", fontface = "bold", fontsize = input$base_size + 3)
        )
        grDevices::pdf(NULL)
        sizes <- tryCatch({
          plot_grob <- ggplot2::ggplotGrob(barplot_reactive())
          legend_index <- match("guide-box-right", plot_grob$layout$name)
          legend_width <- if (!is.na(legend_index)) grid::grobWidth(plot_grob$grobs[[legend_index]]) else grid::unit(0, "in")
          chart_width <- plot_width_in + grid::convertWidth(legend_width, "in", valueOnly = TRUE)
          title_width <- grid::convertWidth(grid::grobWidth(title_grob) + grid::unit(16, "pt"), "in", valueOnly = TRUE)
          title_height <- grid::convertHeight(grid::grobHeight(title_grob) + grid::unit(18, "pt"), "in", valueOnly = TRUE)
          list(chart_width = chart_width, width = max(chart_width, title_width), title_height = title_height)
        }, finally = grDevices::dev.off())

        figure <- gtable::gtable(
          widths = grid::unit(c(sizes$chart_width, sizes$width - sizes$chart_width), "in"),
          heights = grid::unit(c(sizes$title_height, plot_height_in), "in")
        )
        figure <- gtable::gtable_add_grob(figure, title_grob, t = 1, l = 1, r = 2, clip = "off")
        figure <- gtable::gtable_add_grob(figure, plot_grob, t = 2, l = 1)
        ggplot2::ggsave(
          file, plot = figure, device = "png", width = sizes$width,
          height = plot_height_in + sizes$title_height, dpi = 300, bg = "white"
        )
      }
    )

    output$download_barplot_matrix <- downloadHandler(
      filename = function() { paste0("taxa_barplot_matrix_", input$tax_level, "_", input$plot_mode, ".tsv") },
      content = function(file) {
        req(barplot_matrix_reactive())
        readr::write_tsv(barplot_matrix_reactive(), file)
      }
    )

    output$barplot_legend_box <- renderUI({
      box_width <- effective_plot_width()
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
          uiOutput(session$ns("barplot_figure_legend"))
        )
      )
    })

    output$barplot_figure_legend <- renderUI({
      req(input$tax_level, input$plot_mode, input$name_display_mode, input$taxa_display_mode)
      secondary_group <- tryCatch(secondary_var(), error = function(e) NULL)
      secondary_label <- if (!is.null(secondary_group) && nzchar(secondary_group)) {
        secondary_group
      } else {
        "None"
      }
      tax_level_label <- tolower(input$tax_level)
      first_sentence <- if (identical(input$plot_mode, "Group_Mean")) {
        "Stacked bar plots show group-level mean relative abundance of microbial taxa"
      } else {
        "Stacked bar plots show the relative abundance of microbial taxa across samples"
      }
      taxa_name_label <- if (identical(input$name_display_mode, "full")) "full taxonomy hierarchy" else "current rank name"
      sort_sentence <- ""
      if (!identical(input$plot_mode, "Group_Mean")) {
        sort_method <- input$sort_method
        if (length(plot_order()$samples)) {
          sort_sentence <- "Sample order follows the Plot Order settings in Preprocessing. "
        } else if (identical(sort_method, "sample_name")) {
          sort_sentence <- "Within each group, sample order is alphabetical by sample name. "
        } else if (identical(sort_method, "metadata_var")) {
          secondary_sort <- input$secondary_sort
          if (identical(secondary_sort, "taxa_top1")) {
            sort_sentence <- "Within each group, sample order is determined first by the selected metadata variable and then by the relative abundance of the dominant taxon at the selected rank. "
          } else if (identical(secondary_sort, "sample_name_asc")) {
            sort_sentence <- "Within each group, sample order is determined first by the selected metadata variable and then alphabetically by sample name. "
          } else {
            sort_sentence <- "Within each group, sample order is determined by the selected metadata variable. "
          }
        } else {
          sort_sentence <- "Within each group, sample order is determined by the relative abundance of the dominant taxon at the selected rank. "
        }
        if (identical(input$group_var, "__no_facet__") && is.null(secondary_group)) {
          sort_sentence <- sub("^Within each group, sample", "Sample", sort_sentence)
        }
      }
      tags$div(
        tags$div(
          style = "font-weight: 600; margin-bottom: 4px;",
          "Taxonomic composition bar plot"
        ),
        tags$div(
          paste0(
            first_sentence,
            " at ",
            tax_level_label,
            " level. Each colored segment represents a taxon, and its height shows that taxon's share of the bar. ",
            if (identical(input$taxa_display_mode, "top_n")) {
              if (identical(input$top_n_basis, "each_group")) {
                paste0("The top ", input$top_n_taxa,
                       " taxa by mean abundance in each group share a common legend. Remaining taxa are grouped as Others")
              } else {
                paste0("Taxa outside the top ", input$top_n_taxa,
                       " most abundant taxa across selected samples are grouped as Others")
              }
            } else if (identical(input$abundance_threshold_basis, "any_group")) {
              paste0("Taxa with mean abundance below ", input$min_taxa_abundance_pct,
                     "% in every group are grouped as Others")
            } else {
              paste0("Taxa with mean abundance below ", input$min_taxa_abundance_pct,
                     "% across selected samples are grouped as Others")
            },
            if (!identical(secondary_label, "None")) paste0(" and sub-grouped by ", secondary_label) else "",
            ". ",
            sort_sentence,
            "The y-axis represents relative abundance, and taxa labels are shown as ",
            taxa_name_label,
            "; for unassigned entries, the closest annotated higher taxonomic rank is shown instead."            
          )
        )
      )
    })

  })
}

