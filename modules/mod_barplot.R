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
    "Pastel2" = "Pastel2"
  )
  
  tagList(
    tags$style(HTML("
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
        h4(icon("chart-bar"), "Taxa Bar plot"),
        hr(),
        
        uiOutput(ns("local_group_selector")), 
        
        uiOutput(ns("secondary_group_selector")),

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
        h4(icon("sliders"), "Plot Settings"),
        selectInput(ns("plot_mode"), "Plot Display Mode:",
                    choices = c("Sample Level" = "Sample",
                                "Group Mean" = "Group_Mean"),
                    selected = "Sample"),        
        
        selectInput(ns("tax_level"), "Taxonomic Level:",
                    choices = tax_ranks,
                    selected = "Genus"),        
        numericInput(ns("top_n_taxa"), "Taxa Numbers to Display:",
                     value = 15, min = 1, max = 50, step = 1),        
        selectInput(ns("name_display_mode"), "Taxa Rank to Display:",
                    choices = c("Full Hierarchy" = "full",
                                "Current Rank" = "current"),
                    selected = "current"),
        selectInput(ns("color_palette"), "Color Palette:",
                    choices = palette_choices,
                    selected = "Paired"),
        hr(),
        h4(icon("sliders"), "Plot Dimensions"),
        numericInput(ns("plot_width"), "Plot Width (px):",
                     value = 1000, min = 300, step = 50),
        numericInput(ns("plot_height"), "Plot Height (px):",
                     value = 600, min = 300, step = 50),
        numericInput(ns("base_size"), "Base Font Size:",
                     value = 11, min = 6, max = 30, step = 1),

        
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
        h4("Taxa Bar plot"),
        plotOutput(ns("barplot_out"), height = "auto"),
        uiOutput(ns("barplot_legend_box"))
      )
    )
  )
}

## Server
mod_barplot_server <- function(id, ps_obj, meta_cols) {
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

    compute_auto_plot_dims <- function(n_bars) {
      n_bars <- max(1, as.numeric(n_bars))
      per_bar_px <- min(26, max(6, 89.44 / sqrt(n_bars)))
      width <- ceiling(min(max(300, 220 + n_bars * per_bar_px), 2200))
      height <- 400
      list(width = width, height = height)
    }

    maybe_auto_adjust_plot_size <- function(width, height) {
      if (size_update_in_progress()) {
        return()
      }
      current_width <- suppressWarnings(as.numeric(input$plot_width))
      current_height <- suppressWarnings(as.numeric(input$plot_height))
      if (is.finite(current_width) && is.finite(current_height) &&
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
                  choices = group_choices, selected = selected_col)
    })
    
    output$secondary_group_selector <- renderUI({
      req(meta_cols(), input$group_var)
      resolved_primary <- resolve_meta_colname(input$group_var, meta_cols())
      group_choices <- setdiff(meta_cols(), c(input$group_var, resolved_primary))
      non_sampleid_choices <- setdiff(group_choices, "SampleID")
      selected_col <- if (length(non_sampleid_choices) > 0) non_sampleid_choices[1] else if (length(group_choices) > 0) group_choices[1] else NULL
      selectInput(session$ns("secondary_var"), "Secondary Group (Optional):",
                  choices = c("None" = "None", group_choices), 
                  selected = "None")
    })
    
    group_var <- reactive({
      req(input$group_var)
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
      req(meta_cols(), group_var())
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

    barplot_matrix_reactive <- reactive({
      req(ps_obj(), group_var(), input$tax_level, input$name_display_mode,
          input$plot_mode, input$top_n_taxa)

      current_rank <- input$tax_level
      primary_var <- group_var()
      secondary_var_val <- secondary_var()
      topN <- input$top_n_taxa
      plot_mode <- input$plot_mode

      ps_rel <- phyloseq::transform_sample_counts(ps_obj(), function(x) {
        total <- sum(x, na.rm = TRUE)
        if (!is.finite(total) || total <= 0) x else x / total
      })
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

      sorted_taxa <- names(sort(phyloseq::taxa_sums(ps_glom), decreasing = TRUE))
      top_otu_ids <- sorted_taxa[1:min(topN, length(sorted_taxa))]

      df_sample <- phyloseq::psmelt(ps_glom)
      meta_df_base <- as.data.frame(phyloseq::sample_data(ps_rel), stringsAsFactors = FALSE)
      meta_cols <- colnames(meta_df_base)

      for (meta_col in meta_cols) {
        syntactic_col <- make.names(meta_col)
        if (!meta_col %in% names(df_sample) && syntactic_col %in% names(df_sample)) {
          names(df_sample)[names(df_sample) == syntactic_col] <- meta_col
        }
      }
      primary_var <- resolve_meta_colname(primary_var, names(df_sample))
      secondary_var_val <- resolve_meta_colname(secondary_var_val, names(df_sample))

      if (!primary_var %in% names(df_sample)) {
        return(NULL)
      }
      df_sample[[primary_var]] <- factor(df_sample[[primary_var]])

      if (!is.null(secondary_var_val)) {
        if (!secondary_var_val %in% names(df_sample)) {
          return(NULL)
        }
        df_sample[[secondary_var_val]] <- factor(df_sample[[secondary_var_val]])
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

      df_sample$Taxa_Group <- ifelse(df_sample$OTU %in% top_otu_ids,
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

      if (!is.null(secondary_var_val)) {
        meta_df_base[[primary_var]] <- factor(meta_df_base[[primary_var]])
        meta_df_base[[secondary_var_val]] <- factor(meta_df_base[[secondary_var_val]])

        df_sample[[primary_var]] <- factor(df_sample[[primary_var]])
        df_sample[[secondary_var_val]] <- factor(df_sample[[secondary_var_val]])

        meta_df_tibble <- tibble::as_tibble(meta_df_base)
        unique_groups <- meta_df_tibble %>%
          dplyr::arrange(.data[[primary_var]], .data[[secondary_var_val]]) %>%
          dplyr::mutate(Combined_Group = paste(.data[[primary_var]], .data[[secondary_var_val]], sep = "_")) %>%
          dplyr::select(Combined_Group) %>%
          unique() %>%
          dplyr::pull(Combined_Group)

        df_sample$Combined_Group <- paste(df_sample[[primary_var]], df_sample[[secondary_var_val]], sep = "_")
        df_sample$Combined_Group <- factor(df_sample$Combined_Group, levels = unique_groups)
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

      if (!is.null(secondary_var_val) && secondary_var_val %in% names(df_plot)) {
        df_plot %>%
          dplyr::select(any_of(c("Sample", "Combined_Group", primary_var, secondary_var_val, "Taxa_Group", "Abundance")))
      } else {
        df_plot %>%
          dplyr::select(any_of(c("Sample", "Combined_Group", primary_var, "Taxa_Group", "Abundance")))
      }
    })

    observeEvent(list(barplot_matrix_reactive(), input$plot_mode), {
      req(barplot_matrix_reactive())
      df_plot <- barplot_matrix_reactive()
      x_var <- if (identical(input$plot_mode, "Group_Mean")) "Combined_Group" else "Sample"
      if (!x_var %in% names(df_plot) || nrow(df_plot) == 0) {
        return()
      }
      n_bars <- length(unique(as.character(df_plot[[x_var]])))
      dims <- compute_auto_plot_dims(n_bars)
      maybe_auto_adjust_plot_size(dims$width, dims$height)
    }, ignoreInit = FALSE)
    
    barplot_reactive <- reactive({
      req(ps_obj(), group_var(), input$tax_level, input$name_display_mode, 
          input$plot_mode, input$top_n_taxa, input$color_palette, input$base_size)
      
      current_rank <- input$tax_level
      primary_var <- group_var()
      secondary_var_val <- secondary_var()
      
      facet_var <- primary_var 
      topN <- input$top_n_taxa
      plot_mode <- input$plot_mode
      base_size <- input$base_size
      
      ps_rel <- phyloseq::transform_sample_counts(ps_obj(), function(x) {
        total <- sum(x, na.rm = TRUE)
        if (!is.finite(total) || total <= 0) x else x / total
      })
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
        return(paste(full_name_parts, collapse = ";"))
      })
      phyloseq::tax_table(ps_glom) <- phyloseq::tax_table(as.matrix(tax_df))
      
      sorted_taxa <- names(sort(phyloseq::taxa_sums(ps_glom), decreasing = TRUE))
      top_otu_ids <- sorted_taxa[1:min(topN, length(sorted_taxa))]
      
      df_sample <- phyloseq::psmelt(ps_glom)
      meta_df_base <- as.data.frame(phyloseq::sample_data(ps_rel), stringsAsFactors = FALSE)
      meta_cols <- colnames(meta_df_base)
      
      for (meta_col in meta_cols) {
        syntactic_col <- make.names(meta_col)
        if (!meta_col %in% names(df_sample) && syntactic_col %in% names(df_sample)) {
          names(df_sample)[names(df_sample) == syntactic_col] <- meta_col
        }
      }
      primary_var <- resolve_meta_colname(primary_var, names(df_sample))
      secondary_var_val <- resolve_meta_colname(secondary_var_val, names(df_sample))
      
      if (!primary_var %in% names(df_sample)) {
        showNotification(
          paste0("Error: Primary grouping variable '", primary_var, "' was not found in melted data."),
          type = "error"
        )
        return(NULL)
      }
      df_sample[[primary_var]] <- factor(df_sample[[primary_var]])
      
      if (!is.null(secondary_var_val)) {
        if (!secondary_var_val %in% names(df_sample)) {
          showNotification(
            paste0("Error: Secondary grouping variable '", secondary_var_val, "' was not found in melted data."),
            type = "error"
          )
          return(NULL)
        }
        df_sample[[secondary_var_val]] <- factor(df_sample[[secondary_var_val]])
      }
      
      if (input$name_display_mode == "full") {
        df_sample$Taxa_Name_Display <- df_sample$Full_Taxa_Name
        fill_label <- paste(current_rank, "Taxa (Full Hierarchy)")
      } else {
        df_sample$Taxa_Name_Display <- build_disambiguated_rank_labels(
          df_taxa = df_sample,
          current_rank = current_rank,
          rank_order = tax_ranks
        )
        fill_label <- current_rank
      }
      
      df_sample$Taxa_Group <- ifelse(df_sample$OTU %in% top_otu_ids,
                                     as.character(df_sample$Taxa_Name_Display),
                                     "Others")
      
      taxa_mean_abundance <- df_sample %>%
        dplyr::group_by(Taxa_Group) %>%
        dplyr::summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = 'drop') %>%
        dplyr::arrange(desc(Mean_Abundance))
      
      ordered_taxa_levels <- taxa_mean_abundance$Taxa_Group
      if ("Others" %in% ordered_taxa_levels) {
        ordered_taxa_levels <- c(setdiff(ordered_taxa_levels, "Others"), "Others")
      }
      df_sample$Taxa_Group <- factor(df_sample$Taxa_Group, levels = ordered_taxa_levels)
      
      if (!is.null(secondary_var_val)) {
        
        meta_df_base[[primary_var]] <- factor(meta_df_base[[primary_var]])
        meta_df_base[[secondary_var_val]] <- factor(meta_df_base[[secondary_var_val]])
        
        df_sample[[primary_var]] <- factor(df_sample[[primary_var]])
        df_sample[[secondary_var_val]] <- factor(df_sample[[secondary_var_val]])
        
        
        meta_df_tibble <- tibble::as_tibble(meta_df_base)
        
        unique_groups <- meta_df_tibble %>%
          dplyr::arrange(.data[[primary_var]], .data[[secondary_var_val]]) %>%
          dplyr::mutate(Combined_Group = paste(.data[[primary_var]], .data[[secondary_var_val]], sep = "_")) %>%
          dplyr::select(Combined_Group) %>% 
          unique() %>%
          dplyr::pull(Combined_Group)
        
        df_sample$Combined_Group <- paste(df_sample[[primary_var]], df_sample[[secondary_var_val]], sep = "_")
        df_sample$Combined_Group <- factor(df_sample$Combined_Group, levels = unique_groups)
        
        combined_title_suffix <- paste0(" by ", primary_var, " & ", secondary_var_val)    
      } else {
        df_sample$Combined_Group <- df_sample[[primary_var]]
        combined_title_suffix <- paste0(" by ", primary_var)
      }      
      
      if (plot_mode == "Group_Mean") {
        df_sample_grouped <- df_sample %>%
          dplyr::group_by(Sample, Combined_Group, .data[[primary_var]], Taxa_Group) %>%
          dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

        df_plot <- df_sample_grouped %>%
          dplyr::group_by(Combined_Group, .data[[primary_var]], Taxa_Group) %>% 
          dplyr::summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = 'drop')
        
        x_var <- "Combined_Group"
        x_label <- paste(primary_var, if (!is.null(secondary_var_val)) secondary_var_val else "")
        plot_title <- paste0(current_rank, "-Level Microbial Community Composition")
        
        df_plot[[x_var]] <- factor(df_plot[[x_var]])
        
      } else {
        df_plot <- df_sample
        x_var <- "Sample"
        x_label <- ""
        plot_title <- paste0(current_rank, "-Level Microbial Community Composition")
        subgroup_rect_df <- NULL
        subgroup_boundary_df <- NULL
        subgroup_label_df <- NULL
        
        req(input$sort_method)
        sort_method <- input$sort_method
        ordered_sample_levels <- c()
        meta_df_base <- as.data.frame(phyloseq::sample_data(ps_rel), stringsAsFactors = FALSE)
        if (!"Sample" %in% names(meta_df_base)) { meta_df_base$Sample <- rownames(meta_df_base) }
        
        sorting_vars <- c("Sample", primary_var)
        if (!is.null(secondary_var_val)) { sorting_vars <- c(sorting_vars, secondary_var_val) }
        
        
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
            dplyr::mutate(x_pos = dplyr::row_number()) %>%
            dplyr::ungroup()

          subgroup_blocks <- sample_group_df %>%
            dplyr::mutate(subgroup_label = as.character(.data[[secondary_var_val]])) %>%
            dplyr::group_by(.data[[primary_var]], subgroup_label) %>%
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

          subgroup_label_df <- subgroup_blocks %>%
            dplyr::mutate(
              xmid = (xmin + xmax) / 2,
              y = 1.005
            )
        }
      }
      
      palette_choice <- input$color_palette
      num_top_taxa <- length(setdiff(ordered_taxa_levels, "Others"))
      top_taxa_names <- setdiff(ordered_taxa_levels, "Others")
      
      if (palette_choice == "Custom_HCL") {
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
      
      p <- p +
        ggplot2::geom_bar(stat = "identity", position = "fill") +
        ggplot2::scale_fill_manual(values = taxa_colors) + 
        
        if (plot_mode == "Sample" || !is.null(secondary_var_val)) {
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
          linetype = "dashed",
          linewidth = 0.3,
          color = "grey30"
        )
      }

      if (!is.null(secondary_var_val) && identical(plot_mode, "Sample") && !is.null(subgroup_label_df) && nrow(subgroup_label_df) > 0) {
        p <- p + ggplot2::geom_text(
          data = subgroup_label_df,
          ggplot2::aes(x = xmid, y = y, label = subgroup_label),
          inherit.aes = FALSE,
          size = max(2.2, base_size / 3.6),
          fontface = "bold",
          color = "grey20",
          vjust = 0
        )
      }
      
      p <- p +
        ggplot2::theme_bw(base_size = base_size) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.03))) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = base_size + 3, face = "bold"),
          axis.text.x = ggplot2::element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = if (plot_mode == "Group_Mean") base_size else max(6, base_size - 2)
          ),
          axis.title.y = ggplot2::element_text(face = "bold", size = base_size + 1),
          panel.grid.major = ggplot2::element_line(color = "white", linewidth = 0.1),
          panel.grid.minor = ggplot2::element_line(color = "white", linewidth = 0.1),
          axis.title.x = ggplot2::element_blank(),
          legend.title = ggplot2::element_text(size = base_size),
          plot.margin = ggplot2::margin(6, 10, 4, 5.5)
        ) +
        ggplot2::labs(
          title = plot_title,
          y = "Relative Abundance",
          fill = fill_label
        )
      
      return(p)
    })
    
    output$barplot_out <- renderPlot({
      ps_now <- ps_obj()
      if (is.null(ps_now)) {
        graphics::plot.new()
        graphics::text(
          0.5, 0.5,
          "Applying selected samples and preparing bar plot.\nPlease wait...",
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
        barplot_reactive(),
        error = function(e) {
          graphics::plot.new()
          graphics::text(
            0.5, 0.5,
            paste0("Bar plot is not ready yet. Please wait...\n", conditionMessage(e)),
            cex = 0.85
          )
        }
      )
    },
    height = function() { input$plot_height },
    width = function() { input$plot_width }
    )
    
    output$download_barplot <- downloadHandler(
      filename = function() { paste0("taxa_barplot_", input$tax_level, "_", input$plot_mode, ".png") },
      content = function(file) {
        plot_width_in <- input$plot_width / 72
        plot_height_in <- input$plot_height / 72
        
        ggplot2::ggsave(file,
                        plot = barplot_reactive(),
                        device = "png",
                        width = plot_width_in,
                        height = plot_height_in,
                        dpi = 300
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
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 1000 else input$plot_width
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
      req(input$tax_level, input$top_n_taxa, input$plot_mode, input$name_display_mode)
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
        if (identical(sort_method, "sample_name")) {
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
            " level. Colors denote taxa, and taxa outside the top ",
            input$top_n_taxa,
            " most abundant taxa in the full displayed dataset are aggregated as others",
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

