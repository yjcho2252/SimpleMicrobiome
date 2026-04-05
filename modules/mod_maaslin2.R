## UI
mod_maaslin2_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4(icon("flask"), "MaAsLin2 Analysis"),
        hr(),

        selectInput(ns("group_var"), "1. Metadata grouping variable", choices = NULL),

        selectizeInput(
          ns("group_levels"),
          "2. Group levels to include",
          choices = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Select two or more levels",
            plugins = list("remove_button")
          )
        ),

        selectInput(ns("reference_level"), "3. Reference level", choices = NULL),
        verbatimTextOutput(ns("group_sample_counts")),
        hr(),

        selectInput(ns("tax_level"), "4. Taxonomic level",
                    choices = c("ASV", "Genus", "Species"), selected = "Genus"),

        selectInput(ns("volcano_y_axis"), "5. Statistical significance metric",
                    choices = c("FDR-adjusted p-value (q-value)" = "q_val",
                                "Raw p-value (p-value)" = "p_val"),
                    selected = "q_val"),

        numericInput(
          ns("prevalence_filter_pct"),
          "6. Prevalence filter cutoff (0-20%)",
          value = 5,
          min = 0,
          max = 20,
          step = 1
        ),

        tags$details(style = "margin-bottom: 12px;",
          tags$summary(strong("Advanced Options")),
          br(),
          selectizeInput(
            ns("fix_covariates"),
            "7. Additional covariates for fixed effects (optional)",
            choices = NULL,
            multiple = TRUE,
            options = list(
              placeholder = "Select metadata covariates to adjust for",
              plugins = list("remove_button")
            )
          ),
          verbatimTextOutput(ns("model_formula_preview"))
        ),

        numericInput(ns("plot_height"), "Plot height (px)", value = 600, min = 300, max = 2000, step = 50),
        numericInput(ns("plot_width"), "Plot width (px)", value = 900, min = 400, max = 2000, step = 50),
        actionButton(ns("run_maaslin2_btn"), "Run MaAsLin2", class = "btn-danger", style = "font-size: 12px;"),
        tags$script(HTML(
          "Shiny.addCustomMessageHandler('toggle-maaslin2-run-btn', function(msg) {
             var btn = document.getElementById(msg.id);
             if (!btn) return;
             btn.disabled = !!msg.disabled;
             if (msg.label) btn.textContent = msg.label;
           });"
        ))
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Table",
                   downloadButton(ns("download_maaslin2_table"), "Download Table (TSV)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
                   DTOutput(ns("maaslin2_table"))),
          tabPanel("Volcano Plot",
                   downloadButton(ns("download_volcano"), "Download Plot (PNG)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
                   plotOutput(ns("maaslin2_plot"), height = "auto")),
          tabPanel("Bar Plot",
                   downloadButton(ns("download_barplot"), "Download Plot (PNG)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
                   plotOutput(ns("maaslin2_barplot"), height = "auto"))
        )
      )
    )
  )
}

## Server
mod_maaslin2_server <- function(id, ps_obj) {
  moduleServer(id, function(input, output, session) {
    taxa_counts <- reactiveVal(list(before = NA_integer_, after = NA_integer_))

    observeEvent(ps_obj(), {
      req(ps_obj())
      meta_df <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
      meta_cols <- colnames(meta_df)
      candidate_groups <- setdiff(meta_cols, "SampleID")
      group_choices <- candidate_groups[vapply(candidate_groups, function(col) {
        level_values <- as.character(meta_df[[col]])
        level_values <- level_values[!is.na(level_values) & nzchar(level_values)]
        length(unique(level_values)) < 5
      }, logical(1))]
      selected_group <- if (length(group_choices) > 0) group_choices[1] else NULL

      updateSelectInput(session, "group_var",
                        choices = group_choices,
                        selected = selected_group)
    }, ignoreNULL = FALSE)

    observeEvent(list(ps_obj(), input$group_var), {
      req(ps_obj(), input$group_var)
      meta_df <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
      validate(
        need(input$group_var %in% colnames(meta_df), paste0("ERROR: Metadata variable '", input$group_var, "' not found."))
      )

      level_choices <- sort(unique(as.character(meta_df[[input$group_var]])))
      level_choices <- level_choices[!is.na(level_choices) & nzchar(level_choices)]

      selected_levels <- isolate(input$group_levels)
      if (is.null(selected_levels)) {
        selected_levels <- character(0)
      }
      selected_levels <- intersect(selected_levels, level_choices)
      if (length(selected_levels) == 0) {
        selected_levels <- head(level_choices, 5)
      }
      if (length(selected_levels) < 2 && length(level_choices) >= 2) {
        selected_levels <- head(level_choices, min(5, length(level_choices)))
      }

      updateSelectizeInput(session, "group_levels",
                           choices = level_choices,
                           selected = selected_levels,
                           server = TRUE)

      reference_level <- input$reference_level
      if (is.null(reference_level) || !reference_level %in% selected_levels) {
        reference_level <- if (length(selected_levels) > 0) selected_levels[1] else NULL
      }
      updateSelectInput(session, "reference_level", choices = selected_levels, selected = reference_level)

      covariate_choices <- setdiff(colnames(meta_df), c("SampleID", input$group_var))
      selected_covariates <- input$fix_covariates
      if (is.null(selected_covariates)) {
        selected_covariates <- character(0)
      }
      updateSelectizeInput(
        session,
        "fix_covariates",
        choices = covariate_choices,
        selected = intersect(selected_covariates, covariate_choices),
        server = TRUE
      )
    }, ignoreNULL = FALSE)

    group_selection_info <- reactive({
      req(ps_obj(), input$group_var)
      ps <- ps_obj()
      meta_df <- as.data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE)
      validate(
        need(input$group_var %in% colnames(meta_df), paste0("ERROR: Metadata variable '", input$group_var, "' not found in sample_data."))
      )

      selected_levels <- input$group_levels
      if (is.null(selected_levels)) {
        selected_levels <- character(0)
      }
      selected_levels <- selected_levels[nzchar(selected_levels)]

      group_values <- as.character(meta_df[[input$group_var]])
      sample_ids <- rownames(meta_df)
      selected_ids <- sample_ids[group_values %in% selected_levels]

      counts_by_level <- vapply(selected_levels, function(lvl) {
        sum(group_values == lvl, na.rm = TRUE)
      }, numeric(1))

      list(
        selected_levels = selected_levels,
        selected_ids = selected_ids,
        counts_by_level = counts_by_level
      )
    })

    output$group_sample_counts <- renderText({
      req(ps_obj(), input$group_var)
      info <- group_selection_info()
      counts <- taxa_counts()

      lines <- c(
        "Selected Sample Counts",
        paste0("Included levels: ", paste(info$selected_levels, collapse = ", ")),
        paste0("Total selected samples: ", length(info$selected_ids))
      )

      if (length(info$counts_by_level) > 0) {
        per_level <- paste0(names(info$counts_by_level), ": ", as.integer(info$counts_by_level), " sample(s)")
        lines <- c(lines, per_level)
      }

      if (!is.na(counts$before) && !is.na(counts$after)) {
        lines <- c(
          lines,
          paste0("Taxa before prevalence filter: ", counts$before),
          paste0("Taxa after prevalence filter: ", counts$after)
        )
      }

      paste(lines, collapse = "\n")
    })

    ps_filtered <- reactive({
      req(ps_obj(), input$group_var, input$group_levels, input$tax_level)

      info <- group_selection_info()
      selected_levels <- info$selected_levels

      validate(
        need(length(selected_levels) >= 2, "ERROR: Select at least two group levels.")
      )
      validate(
        need(!is.null(input$reference_level) && input$reference_level %in% selected_levels,
             "ERROR: Select a valid reference level from included levels.")
      )

      ps <- ps_obj()
      selected_ids <- info$selected_ids
      validate(
        need(length(selected_ids) > 0, "ERROR: No samples found for selected group levels.")
      )
      ps_sub <- phyloseq::prune_samples(selected_ids, ps)

      if (input$tax_level != "ASV") {
        tax_rank <- input$tax_level
        validate(
          need(!is.null(phyloseq::tax_table(ps_sub)), "ERROR: Taxonomic table missing for aggregation.")
        )
        tax_cols <- colnames(phyloseq::tax_table(ps_sub))
        validate(
          need(tax_rank %in% tax_cols, paste("ERROR: Taxonomic rank", tax_rank, "not found in taxonomy table."))
        )
        ps_sub <- phyloseq::tax_glom(ps_sub, taxrank = tax_rank, NArm = TRUE)

        validate(
          need(phyloseq::ntaxa(ps_sub) > 0, paste("ERROR: After aggregating to", tax_rank, ", no taxa remain."))
        )

        tax_table_sub <- phyloseq::tax_table(ps_sub)
        if (!is.null(tax_table_sub) && tax_rank %in% colnames(tax_table_sub)) {
          tax_df <- as.data.frame(tax_table_sub)
          rank_values <- as.character(tax_df[[tax_rank]])
          rank_values[is.na(rank_values) | rank_values == ""] <- phyloseq::taxa_names(ps_sub)[is.na(rank_values) | rank_values == ""]
          phyloseq::taxa_names(ps_sub) <- make.unique(rank_values)
        }
      }

      prevalence_cutoff <- input$prevalence_filter_pct
      if (is.null(prevalence_cutoff) || length(prevalence_cutoff) == 0 || is.na(prevalence_cutoff)) {
        prevalence_cutoff <- 5
      }
      prevalence_cutoff <- max(0, min(20, as.numeric(prevalence_cutoff)))

      taxa_before_filter <- phyloseq::ntaxa(ps_sub)

      otu_mat <- as.matrix(phyloseq::otu_table(ps_sub))
      if (!phyloseq::taxa_are_rows(phyloseq::otu_table(ps_sub))) {
        otu_mat <- t(otu_mat)
      }

      taxa_prevalence_pct <- rowMeans(otu_mat > 0) * 100
      keep_taxa <- taxa_prevalence_pct > prevalence_cutoff
      ps_sub <- phyloseq::prune_taxa(keep_taxa, ps_sub)

      taxa_after_filter <- phyloseq::ntaxa(ps_sub)
      taxa_counts(list(before = taxa_before_filter, after = taxa_after_filter))

      validate(
        need(
          phyloseq::ntaxa(ps_sub) > 0,
          paste0(
            "ERROR: No taxa remain after prevalence filtering (cutoff = ",
            prevalence_cutoff,
            "%). Lower the cutoff."
          )
        )
      )

      ps_sub
    })

    maaslin2_res <- eventReactive(input$run_maaslin2_btn, {
      req(ps_filtered(), input$group_var, input$reference_level)

      validate(
        need(requireNamespace("Maaslin2", quietly = TRUE),
             "MaAsLin2 package is not installed. Please install 'Maaslin2' first.")
      )

      current_group_var <- input$group_var

      session$sendCustomMessage("toggle-maaslin2-run-btn", list(
        id = session$ns("run_maaslin2_btn"),
        disabled = TRUE,
        label = "Running..."
      ))
      on.exit({
        session$sendCustomMessage("toggle-maaslin2-run-btn", list(
          id = session$ns("run_maaslin2_btn"),
          disabled = FALSE,
          label = "Run MaAsLin2"
        ))
      }, add = TRUE)

      result <- tryCatch({
        withProgress(message = "Running MaAsLin2...", value = 0, {
          ps_current <- ps_filtered()
          otu_mat <- as.matrix(phyloseq::otu_table(ps_current))
          if (!phyloseq::taxa_are_rows(phyloseq::otu_table(ps_current))) {
            otu_mat <- t(otu_mat)
          }
          input_data <- as.data.frame(t(otu_mat), stringsAsFactors = FALSE, check.names = FALSE)

          input_metadata <- as.data.frame(phyloseq::sample_data(ps_current), stringsAsFactors = FALSE)
          input_metadata$SampleID <- rownames(input_metadata)
          rownames(input_metadata) <- input_metadata$SampleID

          common_samples <- intersect(rownames(input_data), rownames(input_metadata))
          input_data <- input_data[common_samples, , drop = FALSE]
          input_metadata <- input_metadata[common_samples, , drop = FALSE]

          covariates <- input$fix_covariates
          if (is.null(covariates) || length(covariates) == 0) {
            covariates <- character(0)
          }
          covariates <- covariates[nzchar(covariates)]
          covariates <- unique(setdiff(covariates, current_group_var))
          fixed_effects <- c(current_group_var, covariates)

          validate(
            need(all(fixed_effects %in% colnames(input_metadata)),
                 "One or more fixed-effect metadata columns are missing.")
          )

          metadata_cols <- unique(c("SampleID", fixed_effects))
          input_metadata <- input_metadata[, metadata_cols, drop = FALSE]
          input_metadata <- as.data.frame(input_metadata, stringsAsFactors = FALSE, check.names = FALSE)
          colnames(input_metadata) <- make.unique(colnames(input_metadata))

          for (col in fixed_effects) {
            if (col %in% colnames(input_metadata)) {
              if (is.data.frame(input_metadata[[col]])) {
                input_metadata[[col]] <- input_metadata[[col]][, 1, drop = TRUE]
              }
              if (is.list(input_metadata[[col]])) {
                input_metadata[[col]] <- vapply(input_metadata[[col]], function(x) as.character(x)[1], character(1))
              }
              input_metadata[[col]] <- as.vector(input_metadata[[col]])
            }
          }

          make_reference <- function(var_name, ref_level) {
            paste0(var_name, ",", ref_level)
          }

          reference_terms <- character(0)
          reference_terms <- c(reference_terms, make_reference(current_group_var, input$reference_level))

          for (col in setdiff(fixed_effects, current_group_var)) {
            if (!col %in% colnames(input_metadata)) {
              next
            }
            col_vals <- input_metadata[[col]]
            if (is.numeric(col_vals)) {
              next
            }
            lvl <- unique(as.character(col_vals))
            lvl <- lvl[!is.na(lvl) & nzchar(lvl)]
            if (length(lvl) > 2) {
              ref_level <- sort(lvl)[1]
              reference_terms <- c(reference_terms, make_reference(col, ref_level))
            }
          }

          out_dir <- file.path(tempdir(), paste0("maaslin2_", session$token, "_", as.integer(Sys.time())))
          dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

          input_data_path <- file.path(out_dir, "maaslin2_input_data.tsv")
          input_metadata_path <- file.path(out_dir, "maaslin2_input_metadata.tsv")

          input_data_out <- tibble::rownames_to_column(input_data, var = "SampleID")
          readr::write_tsv(input_data_out, input_data_path)
          readr::write_tsv(input_metadata, input_metadata_path)

          fit <- Maaslin2::Maaslin2(
            input_data = input_data_path,
            input_metadata = input_metadata_path,
            output = out_dir,
            fixed_effects = fixed_effects,
            reference = reference_terms,
            random_effects = NULL,
            normalization = "TSS",
            transform = "LOG",
            analysis_method = "LM",
            standardize = FALSE,
            max_significance = 1,
            plot_heatmap = FALSE,
            plot_scatter = FALSE
          )

          res <- fit$results
          if (is.null(res) || nrow(res) == 0) {
            results_path <- file.path(out_dir, "all_results.tsv")
            if (file.exists(results_path)) {
              res <- readr::read_tsv(results_path, show_col_types = FALSE)
            }
          }

          validate(
            need(!is.null(res) && nrow(res) > 0, "MaAsLin2 returned no results.")
          )

          res <- as.data.frame(res, stringsAsFactors = FALSE)
          if ("metadata" %in% colnames(res)) {
            metadata_chr <- as.character(res$metadata)
            metadata_chr[is.na(metadata_chr)] <- ""
            metadata_keep <- metadata_chr == current_group_var
            if (!any(metadata_keep)) {
              metadata_keep <- endsWith(metadata_chr, current_group_var)
            }
            if (!any(metadata_keep)) {
              metadata_keep <- grepl(current_group_var, metadata_chr, fixed = TRUE)
            }
            res <- res[metadata_keep, , drop = FALSE]
          }
          if ("value" %in% colnames(res)) {
            selected_levels <- group_selection_info()$selected_levels
            contrast_levels <- setdiff(selected_levels, input$reference_level)
            value_chr <- as.character(res$value)
            value_chr[is.na(value_chr)] <- ""

            normalize_value_to_level <- function(v, levels, group_var_name) {
              if (!nzchar(v)) return(NA_character_)
              if (v %in% levels) return(v)

              prefix <- paste0(group_var_name)
              if (startsWith(v, prefix)) {
                candidate <- substring(v, nchar(prefix) + 1)
                if (candidate %in% levels) return(candidate)
              }

              suffix_hits <- levels[endsWith(v, levels)]
              if (length(suffix_hits) == 1) return(suffix_hits)
              if (length(suffix_hits) > 1) return(suffix_hits[which.max(nchar(suffix_hits))])

              contains_hits <- levels[grepl(v, levels, fixed = TRUE)]
              if (length(contains_hits) == 1) return(contains_hits)
              if (length(contains_hits) > 1) return(contains_hits[which.max(nchar(contains_hits))])

              NA_character_
            }

            value_norm <- vapply(
              value_chr,
              normalize_value_to_level,
              character(1),
              levels = contrast_levels,
              group_var_name = current_group_var
            )

            res$value <- value_norm
            res <- res[!is.na(res$value) & nzchar(res$value), , drop = FALSE]
          }
          validate(
            need(nrow(res) > 0, "No contrast rows were found for the selected grouping variable.")
          )

          res$feature_id <- if ("feature" %in% colnames(res)) as.character(res$feature) else rownames(res)

          tax_table_current <- phyloseq::tax_table(ps_current)
          if (!is.null(tax_table_current)) {
            tax_df <- as.data.frame(tax_table_current)
            tax_df$feature_id <- rownames(tax_df)
            res <- dplyr::left_join(res, tax_df, by = "feature_id")
          }

          preferred_order <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
          res$taxa_label <- NA_character_
          for (col in preferred_order) {
            if (col %in% colnames(res)) {
              idx <- is.na(res$taxa_label) & !is.na(res[[col]]) & res[[col]] != ""
              res$taxa_label[idx] <- as.character(res[[col]][idx])
            }
          }
          missing_idx <- is.na(res$taxa_label) | res$taxa_label == ""
          res$taxa_label[missing_idx] <- res$feature_id[missing_idx]
          row_key <- if ("value" %in% colnames(res)) {
            paste0(res$feature_id, "__", as.character(res$value))
          } else {
            as.character(res$feature_id)
          }
          rownames(res) <- make.unique(row_key)

          res
        })
      }, error = function(e) {
        showNotification(paste("MaAsLin2 Execution Failed:", e$message), type = "error", duration = NULL)
        NULL
      })

      validate(
        need(!is.null(result), "MaAsLin2 analysis failed. Check the error message above.")
      )

      result
    })

    output$model_formula_preview <- renderText({
      req(input$group_var)
      covariates <- input$fix_covariates
      if (is.null(covariates) || length(covariates) == 0) {
        covariates <- character(0)
      }
      covariates <- covariates[nzchar(covariates)]
      covariates <- unique(setdiff(covariates, input$group_var))
      paste0("Applied fixed_effects: ", paste(c(input$group_var, covariates), collapse = " + "))
    })

    round_numeric_columns <- function(df, digits = 2, exclude_cols = c("pval", "qval", "p_val", "q_val")) {
      is_num_col <- vapply(df, is.numeric, logical(1))
      round_target_cols <- setdiff(names(df)[is_num_col], exclude_cols)
      df[round_target_cols] <- lapply(df[round_target_cols], function(x) round(x, digits))
      df
    }

    maaslin2_table_data <- reactive({
      req(maaslin2_res())
      res <- maaslin2_res()

      if (!"feature" %in% colnames(res) && "feature_id" %in% colnames(res)) {
        res$feature <- res$feature_id
      }
      if ("coef" %in% colnames(res)) {
        coef_num <- suppressWarnings(as.numeric(res$coef))
        res$exp_coef <- exp(coef_num)
        res$lfc <- log2(res$exp_coef)
      } else {
        res$exp_coef <- NA_real_
        res$lfc <- NA_real_
      }

      core_cols <- c("feature", "metadata", "value", "coef", "exp_coef", "lfc", "stderr", "N", "pval", "qval")
      taxonomy_cols <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(res))
      id_cols <- intersect(c("taxa_label", "feature_id"), colnames(res))
      metric_cols <- intersect(core_cols, colnames(res))
      used_cols <- unique(c(taxonomy_cols, id_cols, metric_cols))
      remaining_cols <- setdiff(colnames(res), used_cols)
      res <- res[, c(used_cols, remaining_cols), drop = FALSE]
      round_numeric_columns(res, digits = 2)
    })

    maaslin2_processed <- reactive({
      req(maaslin2_res())
      res <- maaslin2_res()
      selected_levels <- group_selection_info()$selected_levels
      ref_level <- input$reference_level

      coef_values <- suppressWarnings(as.numeric(res$coef))
      lfc_values <- log2(exp(coef_values))
      se_values <- if ("stderr" %in% colnames(res)) suppressWarnings(as.numeric(res$stderr)) else rep(NA_real_, nrow(res))
      p_values <- if ("pval" %in% colnames(res)) suppressWarnings(as.numeric(res$pval)) else rep(NA_real_, nrow(res))
      q_values <- if ("qval" %in% colnames(res)) suppressWarnings(as.numeric(res$qval)) else p.adjust(p_values, method = "holm")

      taxonomy_cols <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(res))
      taxonomy_df <- if (length(taxonomy_cols) > 0) res[, taxonomy_cols, drop = FALSE] else NULL

      base_df <- data.frame(
        contrast = if ("value" %in% colnames(res)) as.character(res$value) else NA_character_,
        taxon = res$taxa_label,
        taxa_label = res$taxa_label,
        feature_id = res$feature_id,
        coef = coef_values,
        lfc = lfc_values,
        se = se_values,
        p_val = p_values,
        q_val = q_values,
        diff = ifelse(is.na(q_values), FALSE, q_values < 0.05),
        stringsAsFactors = FALSE
      )

      if (!is.null(taxonomy_df)) {
        base_df <- cbind(taxonomy_df, base_df)
      }
      facet_levels <- setdiff(selected_levels, ref_level)
      if (length(facet_levels) > 0) {
        base_df$contrast <- factor(base_df$contrast, levels = facet_levels)
      } else {
        base_df$contrast <- factor(base_df$contrast)
      }

      base_df$direction <- ifelse(base_df$lfc > 0,
                                  "Increase vs reference",
                                  "Decrease vs reference")
      base_df
    })

    selected_plot_data <- reactive({
      req(maaslin2_processed())
      maaslin2_processed()
    })

    output$maaslin2_table <- renderDT({
      req(maaslin2_table_data())
      datatable(maaslin2_table_data(), options = list(scrollX = TRUE))
    })

    output$download_maaslin2_table <- downloadHandler(
      filename = function() {
        group_tag <- gsub("[^A-Za-z0-9_]+", "_", input$group_var)
        ref_tag <- gsub("[^A-Za-z0-9_]+", "_", input$reference_level)
        paste0("maaslin2_table_", group_tag, "_ref_", ref_tag, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        req(maaslin2_table_data())
        readr::write_tsv(maaslin2_table_data(), file)
      }
    )

    volcano_plot_reactive <- reactive({
      req(selected_plot_data(), input$volcano_y_axis)
      res <- selected_plot_data()

      y_col_name <- input$volcano_y_axis
      y_axis_label <- if (y_col_name == "q_val") expression("-log"[10]*"(FDR-adjusted p-value)") else expression("-log"[10]*"(p-value)")
      y_data <- res[[y_col_name]]

      plot_df <- res[!is.na(res$lfc) & !is.na(y_data) & y_data > 0, , drop = FALSE]
      validate(
        need(nrow(plot_df) > 0, "No valid rows to plot. Check MaAsLin2 results.")
      )
      effect_cutoff <- 0.5

      plot_df$y_metric <- -log10(y_data[!is.na(res$lfc) & !is.na(y_data) & y_data > 0])
      is_significant <- !is.na(plot_df[[y_col_name]]) & is.finite(plot_df[[y_col_name]]) & (plot_df[[y_col_name]] < 0.05)
      plot_df$diffexpressed <- dplyr::case_when(
        is_significant & !is.na(plot_df$lfc) & is.finite(plot_df$lfc) & (plot_df$lfc >= effect_cutoff) ~ "Increased",
        is_significant & !is.na(plot_df$lfc) & is.finite(plot_df$lfc) & (plot_df$lfc <= -effect_cutoff) ~ "Decreased",
        TRUE ~ "Not significant"
      )
      plot_df$diffexpressed <- factor(plot_df$diffexpressed, levels = c("Decreased", "Not significant", "Increased"))
      plot_df$delabel <- ifelse(plot_df$diffexpressed == "Not significant", "", plot_df$taxa_label)
      y_vals <- plot_df$y_metric[is.finite(plot_df$y_metric) & !is.na(plot_df$y_metric)]
      y_upper <- if (length(y_vals) > 0) max(stats::quantile(y_vals, probs = 0.99, na.rm = TRUE), max(y_vals, na.rm = TRUE) * 0.1) + 1 else 10
      y_upper <- max(5, y_upper)

      x_vals <- plot_df$lfc[is.finite(plot_df$lfc) & !is.na(plot_df$lfc)]
      x_abs_max <- if (length(x_vals) > 0) as.numeric(stats::quantile(abs(x_vals), probs = 0.99, na.rm = TRUE)) + 0.5 else 2
      x_abs_max <- max(1, ceiling(x_abs_max))
      x_lim <- c(-x_abs_max, x_abs_max)
      x_step <- max(1, ceiling((2 * x_abs_max) / 10))
      x_breaks <- seq(from = -x_abs_max, to = x_abs_max, by = x_step)

      p <- ggplot(plot_df, aes(x = lfc, y = y_metric, col = diffexpressed, label = delabel)) +
        geom_vline(xintercept = c(-effect_cutoff, effect_cutoff), col = "gray", linetype = "dashed", linewidth = 0.4) +
        geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed", linewidth = 0.4) +
        geom_point(size = 2.8, alpha = 0.9) +
        scale_color_manual(
          values = c("Decreased" = "#00AFBB", "Not significant" = "grey", "Increased" = "#bb0c00"),
          labels = c("Decreased", "Not significant", "Increased")
        ) +
        coord_cartesian(ylim = c(0, y_upper), xlim = x_lim) +
        scale_x_continuous(breaks = x_breaks) +
        labs(title = paste("MaAsLin2 Volcano Plot (reference:", input$reference_level, ")"),
             x = expression("log"[2]*"FC"),
             y = y_axis_label,
             color = "Direction",
             subtitle = "Positive log2FC = increased vs reference") +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12)
        )

      if ("contrast" %in% colnames(plot_df)) {
        p <- p + facet_wrap(~ contrast, scales = "free_x", drop = FALSE)
      }

      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(max.overlaps = Inf, size = 3.5, show.legend = FALSE)
      } else {
        p <- p + geom_text(size = 3.5, vjust = -0.8, show.legend = FALSE, check_overlap = TRUE)
      }

      p
    })

    bar_plot_reactive <- reactive({
      req(selected_plot_data())
      res <- selected_plot_data()
      plot_df <- res[res$diff & !is.na(res$lfc), , drop = FALSE]
      if (nrow(plot_df) == 0) {
        plot_df <- res[!is.na(res$lfc), , drop = FALSE]
      }
      validate(
        need(nrow(plot_df) > 0, "No valid rows to plot for bar chart.")
      )

      top_n <- 10

      if ("contrast" %in% colnames(plot_df)) {
        res_top <- plot_df %>%
          dplyr::group_by(contrast) %>%
          dplyr::slice_max(order_by = abs(lfc), n = top_n, with_ties = FALSE) %>%
          dplyr::ungroup()
      } else {
        res_top <- plot_df[order(abs(plot_df$lfc), decreasing = TRUE), , drop = FALSE]
        res_top <- res_top[seq_len(min(top_n, nrow(res_top))), , drop = FALSE]
      }

      res_top$direction_flag <- res_top$lfc > 0

      p <- ggplot(res_top, aes(x = reorder(taxa_label, lfc), y = lfc, fill = direction_flag)) +
        geom_col(color = "black") +
        coord_flip() +
        scale_fill_manual(
          values = c("TRUE" = "red", "FALSE" = "blue"),
          labels = c(
            "TRUE" = "Increase vs reference",
            "FALSE" = "Decrease vs reference"
          )
        ) +
        labs(title = paste("Top", top_n, "Differential Taxa by log2FC (reference:", input$reference_level, ")"),
             x = input$tax_level,
             y = expression("log"[2]*"FC"),
             fill = "Direction") +
        theme_bw()

      if ("contrast" %in% colnames(res_top)) {
        p <- p + facet_wrap(~ contrast, scales = "free_y", drop = FALSE)
      }

      p
    })

    output$maaslin2_plot <- renderPlot(
      { volcano_plot_reactive() },
      height = function() { req(input$plot_height); input$plot_height },
      width = function() { req(input$plot_width); input$plot_width }
    )

    output$maaslin2_barplot <- renderPlot(
      { bar_plot_reactive() },
      height = function() { req(input$plot_height); input$plot_height },
      width = function() { req(input$plot_width); input$plot_width }
    )

    output$download_volcano <- downloadHandler(
      filename = function() {
        group_tag <- gsub("[^A-Za-z0-9_]+", "_", input$group_var)
        ref_tag <- gsub("[^A-Za-z0-9_]+", "_", input$reference_level)
        paste0("maaslin2_volcano_", group_tag, "_ref_", ref_tag, ".png")
      },
      content = function(file) {
        req(input$plot_height, input$plot_width)
        ggplot2::ggsave(file, plot = volcano_plot_reactive(), device = "png",
                        width = input$plot_width / 100, height = input$plot_height / 100)
      }
    )

    output$download_barplot <- downloadHandler(
      filename = function() {
        group_tag <- gsub("[^A-Za-z0-9_]+", "_", input$group_var)
        ref_tag <- gsub("[^A-Za-z0-9_]+", "_", input$reference_level)
        paste0("maaslin2_barplot_", group_tag, "_ref_", ref_tag, ".png")
      },
      content = function(file) {
        req(input$plot_height, input$plot_width)
        ggplot2::ggsave(file, plot = bar_plot_reactive(), device = "png",
                        width = input$plot_width / 100, height = input$plot_height / 100)
      }
    )
  })
}

