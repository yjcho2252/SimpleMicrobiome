## UI
mod_maaslin2_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML("
      .simple-result-card {
        width: 600px;
        max-width: 100%;
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
      .simple-result-card pre {
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
        h4(icon("flask"), "MaAsLin2"),
        hr(),

        selectInput(ns("group_var"), "1. Grouping variable", choices = NULL),

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

        selectInput(ns("tax_level"), "4. Taxonomic level",
                    choices = c("ASV", "Genus", "Species"), selected = "Genus"),

        selectInput(ns("volcano_y_axis"), "5. Statistical metric",
                    choices = c("FDR q-value" = "q_val",
                                "p-value" = "p_val"),
                    selected = "q_val"),

        numericInput(
          ns("prevalence_filter_pct"),
          "6. Prevalence filter (0-20%)",
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
          selectizeInput(
            ns("fix_interactions"),
            "8. Interaction terms for fixed effects (optional)",
            choices = NULL,
            multiple = TRUE,
            options = list(
              placeholder = "Select interaction terms among fixed effects",
              plugins = list("remove_button")
            )
          ),
          selectizeInput(
            ns("random_effects"),
            "9. Random effects (optional)",
            choices = NULL,
            multiple = TRUE,
            options = list(
              placeholder = "Select grouping variables for random effects",
              plugins = list("remove_button")
            )
          ),
          verbatimTextOutput(ns("model_formula_preview"))
        ),
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Plot Dimensions"),
        numericInput(ns("plot_width"), "Plot width (px)", value = 700, min = 400, max = 2000, step = 50),
        numericInput(ns("plot_height"), "Plot height (px)", value = 400, min = 300, max = 2000, step = 50),
        numericInput(ns("base_size"), "Base Font Size:", value = 11, min = 6, max = 30, step = 1),
        actionButton(ns("run_maaslin2_btn"), "Run MaAsLin2", class = "btn-danger", style = "font-size: 12px;"),
        tags$script(HTML(
          "Shiny.addCustomMessageHandler('toggle-maaslin2-run-btn', function(msg) {
             var btn = document.getElementById(msg.id);
             if (!btn) return;
             btn.disabled = !!msg.disabled;
             if (msg.label) btn.textContent = msg.label;
           });
           Shiny.addCustomMessageHandler('set-tab-container-width', function(msg) {
             var el = document.getElementById(msg.id);
             if (!el) return;
             el.style.width = msg.width;
             el.style.maxWidth = '100%';
           });"
        ))
      ),
      mainPanel(
        h4("MaAsLin2"),
        tags$div(
          id = ns("maaslin2_tab_container"),
          style = "max-width: 100%;",
          tabsetPanel(
            id = ns("maaslin2_active_tab"),
            tabPanel("Volcano Plot",
                     downloadButton(ns("download_volcano"), "Download Plot (PNG)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
                     tags$div(
                       style = "margin-top: 8px;",
                       plotOutput(ns("maaslin2_plot"), height = "auto")
                     )),
            tabPanel("Bar Plot",
                     downloadButton(ns("download_barplot"), "Download Plot (PNG)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
                     tags$div(
                       style = "margin-top: 8px;",
                       plotOutput(ns("maaslin2_barplot"), height = "auto")
                     )),
            tabPanel("Table",
                     downloadButton(ns("download_maaslin2_table"), "Download Table (TSV)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
                     DTOutput(ns("maaslin2_table")))
          )
        ),
        uiOutput(ns("maaslin2_legend_box")),
        uiOutput(ns("maaslin2_status_separator")),
        h5(icon("circle-info"), "MaAsLin2 Status"),
        uiOutput(ns("maaslin2_status_box"))
      )
    )
  )
}

## Server
mod_maaslin2_server <- function(id, ps_obj) {
  moduleServer(id, function(input, output, session) {
    draw_wait_message <- function(message_text) {
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      graphics::par(mar = c(0, 0, 0, 0))
      graphics::plot.new()
      graphics::text(0.5, 0.5, message_text)
    }

    observe({
      width_px <- suppressWarnings(as.integer(input$plot_width))
      if (!is.finite(width_px) || is.na(width_px) || width_px <= 0) width_px <- 900L
      session$sendCustomMessage(
        "set-tab-container-width",
        list(id = session$ns("maaslin2_tab_container"), width = paste0(width_px, "px"))
      )
    })

    taxa_counts <- reactiveVal(list(before = NA_integer_, after = NA_integer_))
    apply_disambiguated_taxrank <- function(ps, tax_level) {
      if (is.null(ps) || identical(tax_level, "ASV")) {
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
    group_var_resolved <- reactive({
      req(ps_obj(), input$group_var)
      meta_df <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
      resolve_meta_colname(input$group_var, colnames(meta_df))
    })
    
    build_interaction_choices <- function(group_var, fix_covariates) {
      if (is.null(group_var) || !nzchar(group_var)) {
        return(character(0))
      }
      covariates <- fix_covariates
      if (is.null(covariates) || length(covariates) == 0) {
        covariates <- character(0)
      }
      covariates <- covariates[nzchar(covariates)]
      covariates <- unique(setdiff(covariates, group_var))
      fixed_terms <- unique(c(group_var, covariates))
      if (length(fixed_terms) < 2) {
        return(character(0))
      }
      comb_mat <- utils::combn(fixed_terms, 2)
      apply(comb_mat, 2, function(x) paste(x[1], x[2], sep = ":"))
    }

    expand_interaction_columns <- function(metadata_df, interaction_terms) {
      sanitize_name <- function(x) {
        x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
        x <- gsub("^_+|_+$", "", x)
        if (!nzchar(x)) {
          return("x")
        }
        x
      }

      if (is.null(interaction_terms) || length(interaction_terms) == 0) {
        return(list(metadata = metadata_df, interaction_effects = character(0)))
      }

      out_metadata <- metadata_df
      created_effects <- character(0)
      interaction_terms <- interaction_terms[nzchar(interaction_terms)]

      for (term in interaction_terms) {
        parts <- strsplit(term, ":", fixed = TRUE)[[1]]
        if (length(parts) != 2) {
          next
        }
        vars <- trimws(parts)
        if (any(!vars %in% colnames(out_metadata))) {
          next
        }

        design_df <- data.frame(
          lhs = out_metadata[[vars[1]]],
          rhs = out_metadata[[vars[2]]],
          stringsAsFactors = FALSE
        )

        interaction_mat <- tryCatch(
          stats::model.matrix(~ 0 + lhs:rhs, data = design_df),
          error = function(e) NULL
        )
        if (is.null(interaction_mat) || ncol(interaction_mat) == 0) {
          next
        }

        keep_col <- vapply(seq_len(ncol(interaction_mat)), function(i) {
          vals <- as.numeric(interaction_mat[, i])
          any(!is.na(vals)) && stats::sd(vals, na.rm = TRUE) > 0
        }, logical(1))
        interaction_mat <- interaction_mat[, keep_col, drop = FALSE]
        if (ncol(interaction_mat) == 0) {
          next
        }

        for (j in seq_len(ncol(interaction_mat))) {
          mm_name <- colnames(interaction_mat)[j]
          mm_label <- gsub("^lhs", paste0(vars[1], "_"), mm_name)
          mm_label <- gsub(":rhs", paste0("_X_", vars[2], "_"), mm_label, fixed = TRUE)
          mm_label <- gsub("^rhs", paste0(vars[2], "_"), mm_label)
          mm_label <- sanitize_name(mm_label)
          base_name <- paste0(
            "int_",
            sanitize_name(vars[1]),
            "_X_",
            sanitize_name(vars[2]),
            "__",
            mm_label
          )
          new_name <- base_name
          suffix <- 1L
          while (new_name %in% colnames(out_metadata)) {
            suffix <- suffix + 1L
            new_name <- paste0(base_name, "_", suffix)
          }
          out_metadata[[new_name]] <- as.numeric(interaction_mat[, j])
          created_effects <- c(created_effects, new_name)
        }
      }

      list(metadata = out_metadata, interaction_effects = created_effects)
    }

    observeEvent(ps_obj(), {
      req(ps_obj())
      meta_df <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
      meta_cols <- colnames(meta_df)
      candidate_groups <- setdiff(meta_cols, "SampleID")
      group_choices <- candidate_groups[vapply(candidate_groups, function(col) {
        level_values <- as.character(meta_df[[col]])
        level_values <- level_values[!is.na(level_values) & nzchar(level_values)]
        length(unique(level_values)) <= 10
      }, logical(1))]
      selected_group <- if (length(group_choices) > 0) group_choices[1] else NULL

      updateSelectInput(session, "group_var",
                        choices = group_choices,
                        selected = selected_group)
    }, ignoreNULL = FALSE)

    observeEvent(list(ps_obj(), input$group_var), {
      req(ps_obj(), input$group_var)
      meta_df <- as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)
      group_var <- group_var_resolved()
      validate(
        need(group_var %in% colnames(meta_df), paste0("ERROR: Metadata variable '", input$group_var, "' not found."))
      )

      level_choices <- sort(unique(as.character(meta_df[[group_var]])))
      level_choices <- level_choices[!is.na(level_choices) & nzchar(level_choices)]

      selected_levels <- isolate(input$group_levels)
      if (is.null(selected_levels)) {
        selected_levels <- character(0)
      }
      selected_levels <- intersect(selected_levels, level_choices)
      if (length(selected_levels) == 0) {
        selected_levels <- head(level_choices, 2)
      }
      if (length(selected_levels) < 2 && length(level_choices) >= 2) {
        selected_levels <- head(level_choices, 2)
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

      covariate_choices <- setdiff(colnames(meta_df), c("SampleID", group_var))
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

      interaction_choices <- build_interaction_choices(group_var, selected_covariates)
      selected_interactions <- input$fix_interactions
      if (is.null(selected_interactions)) {
        selected_interactions <- character(0)
      }
      updateSelectizeInput(
        session,
        "fix_interactions",
        choices = interaction_choices,
        selected = intersect(selected_interactions, interaction_choices),
        server = TRUE
      )

      random_choices <- setdiff(colnames(meta_df), "SampleID")
      selected_random <- input$random_effects
      if (is.null(selected_random)) {
        selected_random <- character(0)
      }
      updateSelectizeInput(
        session,
        "random_effects",
        choices = random_choices,
        selected = intersect(selected_random, random_choices),
        server = TRUE
      )
    }, ignoreNULL = FALSE)
    
    observeEvent(list(input$group_var, input$fix_covariates), {
      interaction_choices <- build_interaction_choices(group_var_resolved(), input$fix_covariates)
      selected_interactions <- input$fix_interactions
      if (is.null(selected_interactions)) {
        selected_interactions <- character(0)
      }
      updateSelectizeInput(
        session,
        "fix_interactions",
        choices = interaction_choices,
        selected = intersect(selected_interactions, interaction_choices),
        server = TRUE
      )
    }, ignoreNULL = FALSE)

    group_selection_info <- reactive({
      req(ps_obj(), input$group_var)
      ps <- ps_obj()
      meta_df <- as.data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE)
      group_var <- group_var_resolved()
      validate(
        need(group_var %in% colnames(meta_df), paste0("ERROR: Metadata variable '", input$group_var, "' not found in sample_data."))
      )

      selected_levels <- input$group_levels
      if (is.null(selected_levels)) {
        selected_levels <- character(0)
      }
      selected_levels <- selected_levels[nzchar(selected_levels)]

      group_values <- as.character(meta_df[[group_var]])
      sample_ids <- rownames(meta_df)
      selected_ids <- sample_ids[group_values %in% selected_levels]

      counts_by_level <- vapply(selected_levels, function(lvl) {
        sum(group_values == lvl, na.rm = TRUE)
      }, numeric(1))

      list(
        group_var = group_var,
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
        ps_sub <- apply_disambiguated_taxrank(ps_sub, tax_rank)
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

    maaslin2_running <- reactiveVal(FALSE)

    maaslin2_res <- eventReactive(input$run_maaslin2_btn, {
      req(ps_filtered(), input$group_var, input$reference_level)
      maaslin2_running(TRUE)

      validate(
        need(requireNamespace("Maaslin2", quietly = TRUE),
             "MaAsLin2 package is not installed. Please install 'Maaslin2' first.")
      )

      current_group_var <- group_var_resolved()

      session$sendCustomMessage("toggle-maaslin2-run-btn", list(
        id = session$ns("run_maaslin2_btn"),
        disabled = TRUE,
        label = "Running..."
      ))
      on.exit({
        maaslin2_running(FALSE)
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
          interaction_terms <- input$fix_interactions
          if (is.null(interaction_terms) || length(interaction_terms) == 0) {
            interaction_terms <- character(0)
          }
          interaction_terms <- interaction_terms[nzchar(interaction_terms)]
          valid_interactions <- build_interaction_choices(current_group_var, covariates)
          interaction_terms <- intersect(interaction_terms, valid_interactions)
          random_effects <- input$random_effects
          if (is.null(random_effects) || length(random_effects) == 0) {
            random_effects <- character(0)
          }
          random_effects <- random_effects[nzchar(random_effects)]
          random_effects <- unique(setdiff(random_effects, "SampleID"))
          main_effects <- c(current_group_var, covariates)

          validate(
            need(all(main_effects %in% colnames(input_metadata)),
                 "One or more fixed-effect metadata columns are missing.")
          )
          validate(
            need(all(random_effects %in% colnames(input_metadata)),
                 "One or more random-effect metadata columns are missing.")
          )

          metadata_cols <- unique(c("SampleID", main_effects, random_effects))
          input_metadata <- input_metadata[, metadata_cols, drop = FALSE]
          input_metadata <- as.data.frame(input_metadata, stringsAsFactors = FALSE, check.names = FALSE)
          colnames(input_metadata) <- make.unique(colnames(input_metadata))

          for (col in unique(c(main_effects, random_effects))) {
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

          expanded_interactions <- expand_interaction_columns(input_metadata, interaction_terms)
          input_metadata <- expanded_interactions$metadata
          interaction_effects <- expanded_interactions$interaction_effects
          fixed_effects <- unique(c(main_effects, interaction_effects))
          random_effects_arg <- if (length(random_effects) > 0) random_effects else NULL

          if (length(interaction_terms) > 0 && length(interaction_effects) == 0) {
            showNotification(
              "Selected interaction terms could not be encoded and were skipped.",
              type = "warning",
              duration = 8
            )
          }

          make_reference <- function(var_name, ref_level) {
            paste0(var_name, ",", ref_level)
          }

          reference_terms <- character(0)
          reference_terms <- c(reference_terms, make_reference(current_group_var, input$reference_level))

          for (col in setdiff(main_effects, current_group_var)) {
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
            random_effects = random_effects_arg,
            normalization = "TSS",
            transform = "LOG",
            analysis_method = "LM",
            standardize = FALSE,
            cores = 4,
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
      group_var <- group_var_resolved()
      covariates <- input$fix_covariates
      if (is.null(covariates) || length(covariates) == 0) {
        covariates <- character(0)
      }
      covariates <- covariates[nzchar(covariates)]
      covariates <- unique(setdiff(covariates, group_var))
      interaction_terms <- input$fix_interactions
      if (is.null(interaction_terms) || length(interaction_terms) == 0) {
        interaction_terms <- character(0)
      }
      interaction_terms <- interaction_terms[nzchar(interaction_terms)]
      valid_interactions <- build_interaction_choices(group_var, covariates)
      interaction_terms <- intersect(interaction_terms, valid_interactions)
      random_effects <- input$random_effects
      if (is.null(random_effects) || length(random_effects) == 0) {
        random_effects <- character(0)
      }
      random_effects <- random_effects[nzchar(random_effects)]
      random_effects <- unique(setdiff(random_effects, "SampleID"))
      fixed_effects <- unique(c(group_var, covariates))
      preview <- paste0("Applied fixed_effects: ", paste(fixed_effects, collapse = " + "))
      if (length(interaction_terms) > 0) {
        preview <- paste0(
          preview,
          "\nSelected interactions (expanded into derived metadata columns at run time): ",
          paste(interaction_terms, collapse = ", ")
        )
      }
      preview <- paste0(
        preview,
        "\nApplied random_effects: ",
        if (length(random_effects) > 0) paste(random_effects, collapse = " + ") else "NULL"
      )
      preview
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
        coef_num_raw <- suppressWarnings(as.numeric(res$coef))
        coef_num <- coef_num_raw
        res$coef <- coef_num
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

      coef_values_raw <- suppressWarnings(as.numeric(res$coef))
      coef_values <- coef_values_raw
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

      contrast_label <- as.character(base_df$contrast)
      contrast_label[is.na(contrast_label) | !nzchar(contrast_label)] <- "selected comparison levels"
      base_df$direction <- ifelse(
        base_df$lfc > 0,
        paste0("Increase in ", contrast_label),
        paste0("Decrease in ", contrast_label)
      )
      base_df
    })

    selected_plot_data <- reactive({
      req(maaslin2_processed())
      maaslin2_processed()
    })

    output$maaslin2_table <- renderDT({
      req(maaslin2_table_data())
      res <- maaslin2_table_data()
      datatable(
        res,
        options = list(scrollX = TRUE),
        class = "compact stripe hover cell-border",
        style = "bootstrap"
      ) %>%
        formatStyle(
          columns = colnames(res),
          `font-size` = "11px",
          `padding` = "4px"
        )
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

    resolve_current_taxa_labels <- function(df, tax_level) {
      is_unassigned_taxa <- function(x) {
        x <- as.character(x)
        x_trim <- trimws(x)
        is.na(x_trim) | !nzchar(x_trim) | grepl("unassigned|uncultured", x_trim, ignore.case = TRUE)
      }

      is_raw_placeholder_taxa <- function(x) {
        x <- as.character(x)
        x_trim <- trimws(x)
        x_norm <- tolower(gsub("^[a-z]__", "", x_trim))
        is.na(x_trim) | !nzchar(x_trim) | grepl("unassigned|uncultured", x_norm, ignore.case = TRUE)
      }

      # Extract current rank only from potential full hierarchy
      extract_current_rank <- function(x, rank_name) {
        x <- as.character(x)
        
        # If contains ";", extract the part corresponding to current rank
        if (grepl(";", x, fixed = TRUE)) {
          parts <- strsplit(x, ";", fixed = TRUE)[[1]]
          parts <- trimws(parts)
          # Return the last non-empty part (assumed to be current rank)
          parts <- parts[nzchar(parts)]
          if (length(parts) > 0) return(parts[length(parts)])
        }
        
        # Keep original rank prefix (e.g., f__, g__) if present.
        x
      }

      rank_order <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      out <- if (identical(tax_level, "ASV")) {
        if ("feature_id" %in% colnames(df)) as.character(df$feature_id) else rep(NA_character_, nrow(df))
      } else {
        label_vec <- rep(NA_character_, nrow(df))
        if (tax_level %in% rank_order) {
          selected_idx <- match(tax_level, rank_order)
          raw_rank_vals <- as.character(df[[tax_level]])
          current_rank_vals <- sapply(raw_rank_vals, function(x) extract_current_rank(x, tax_level), USE.NAMES = FALSE)
          
          # Initialize with extracted current rank value
          label_vec <- current_rank_vals

          # For unassigned current taxa, use the nearest parent rank label.
          is_unassigned_idx <- is_unassigned_taxa(label_vec)
          if (any(is_unassigned_idx)) {
            # Get parent ranks in order
            parent_ranks <- rev(rank_order[seq_len(selected_idx - 1)])
            parent_ranks <- parent_ranks[parent_ranks %in% colnames(df)]
            if (length(parent_ranks) > 0) {
              parent_vals <- rep("UnclassifiedParent", nrow(df))
              for (parent_rank in parent_ranks) {
                candidate_vals <- as.character(df[[parent_rank]])
                candidate_vals <- sapply(candidate_vals, function(x) extract_current_rank(x, parent_rank), USE.NAMES = FALSE)
                candidate_placeholder_idx <- is_unassigned_taxa(candidate_vals)
                use_idx <- is_unassigned_idx & parent_vals == "UnclassifiedParent" & !candidate_placeholder_idx
                parent_vals[use_idx] <- candidate_vals[use_idx]
              }
              label_vec[is_unassigned_idx] <- parent_vals[is_unassigned_idx]
            }
          }
        }
        label_vec
      }

      fallback_taxa <- if ("taxa_label" %in% colnames(df)) as.character(df$taxa_label) else rep("", nrow(df))
      fallback_feature <- if ("feature_id" %in% colnames(df)) as.character(df$feature_id) else rep("", nrow(df))
      out[is.na(out) | !nzchar(out) | is_raw_placeholder_taxa(out)] <- fallback_taxa[is.na(out) | !nzchar(out) | is_raw_placeholder_taxa(out)]
      out[is.na(out) | !nzchar(out) | is_raw_placeholder_taxa(out)] <- fallback_feature[is.na(out) | !nzchar(out) | is_raw_placeholder_taxa(out)]
      out
    }

    volcano_plot_reactive <- reactive({
      req(selected_plot_data(), input$volcano_y_axis, input$base_size)
      res <- selected_plot_data()
      selected_levels <- group_selection_info()$selected_levels
      target_levels <- setdiff(selected_levels, input$reference_level)
      target_label <- if (length(target_levels) == 1) target_levels[1] else "selected comparison levels"

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
      label_by_rank <- resolve_current_taxa_labels(plot_df, input$tax_level)
      plot_df$delabel <- ifelse(plot_df$diffexpressed == "Not significant", "", label_by_rank)
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
          breaks = c("Decreased", "Not significant", "Increased"),
          labels = c(
            "Decreased" = "Decreased",
            "Not significant" = "Not significant",
            "Increased" = "Increased"
          ),
          drop = FALSE
        ) +
        coord_cartesian(ylim = c(0, y_upper), xlim = x_lim) +
        scale_x_continuous(breaks = x_breaks) +
        labs(title = paste0(input$tax_level, "-Level MaAsLin2 Volcano Plot (reference: ", input$reference_level, ")"),
             x = expression("log"[2]*"FC"),
             y = y_axis_label,
             color = "Direction",
             subtitle = paste0("Positive log2FC = increased in ", target_label, " (vs ", input$reference_level, ")")) +
        theme_minimal(base_size = input$base_size) +
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
      req(selected_plot_data(), input$volcano_y_axis, input$base_size)
      res <- selected_plot_data()
      selected_levels <- group_selection_info()$selected_levels
      target_levels <- setdiff(selected_levels, input$reference_level)
      target_label <- if (length(target_levels) == 1) target_levels[1] else "selected comparison levels"
      sig_col <- input$volcano_y_axis
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

      sig_values <- suppressWarnings(as.numeric(res_top[[sig_col]]))
      res_top$stat_value_label <- ifelse(
        is.na(sig_values),
        "NA",
        formatC(sig_values, format = "e", digits = 2)
      )
      res_top$taxa_label_current <- make.unique(resolve_current_taxa_labels(res_top, input$tax_level))

      res_top$direction <- ifelse(res_top$lfc > 0, "Increased", "Decreased")
      res_top$direction <- factor(res_top$direction, levels = c("Decreased", "Increased"))
      y_abs_max <- max(abs(res_top$lfc), na.rm = TRUE)
      if (!is.finite(y_abs_max) || y_abs_max <= 0) {
        y_abs_max <- 1
      }
      label_offset <- y_abs_max * 0.08
      label_right <- max(0, max(res_top$lfc, na.rm = TRUE)) + label_offset
      res_top$label_y <- label_right
      res_top$label_hjust <- 0
      y_limits <- c(min(res_top$lfc, na.rm = TRUE), label_right + (label_offset * 4))

      p <- ggplot(res_top, aes(x = reorder(taxa_label_current, lfc), y = lfc, fill = direction)) +
        geom_col(color = "black") +
        geom_text(
          aes(y = label_y, label = stat_value_label, hjust = label_hjust),
          size = 3
        ) +
        coord_flip(clip = "off") +
        scale_y_continuous(limits = y_limits) +
        scale_fill_manual(
          values = c("Decreased" = "blue", "Increased" = "red"),
          breaks = c("Decreased", "Increased"),
          labels = c(
            "Decreased" = paste0("Decrease in ", target_label),
            "Increased" = paste0("Increase in ", target_label)
          ),
          drop = FALSE
        ) +
        labs(title = paste0("Top ", top_n, " ", input$tax_level, "-Level Differential Taxa by log2FC (reference: ", input$reference_level, ")"),
             x = input$tax_level,
             y = expression("log"[2]*"FC"),
             fill = "Direction") +
        theme_bw(base_size = input$base_size) +
        theme(plot.margin = ggplot2::margin(6, 90, 6, 90))

      if ("contrast" %in% colnames(res_top)) {
        p <- p + facet_wrap(~ contrast, scales = "free_y", drop = FALSE)
      }

      p
    })

    output$maaslin2_plot <- renderPlot(
      {
        if (is.null(input$run_maaslin2_btn) || input$run_maaslin2_btn < 1) {
          draw_wait_message("Click 'Run MaAsLin2' to start analysis.")
          return(invisible(NULL))
        }
        if (isTRUE(maaslin2_running())) {
          draw_wait_message("MaAsLin2 is running. Please wait...")
          return(invisible(NULL))
        }
        volcano_plot_reactive()
      },
      height = function() {
        if (is.null(input$run_maaslin2_btn) || input$run_maaslin2_btn < 1) {
          200
        } else {
          h <- suppressWarnings(as.numeric(input$plot_height))
          if (!is.finite(h) || is.na(h) || h <= 0) 400 else h
        }
      },
      width = function() {
        w <- suppressWarnings(as.numeric(input$plot_width))
        if (!is.finite(w) || is.na(w) || w <= 0) 700 else w
      }
    )

    output$maaslin2_barplot <- renderPlot(
      {
        if (is.null(input$run_maaslin2_btn) || input$run_maaslin2_btn < 1) {
          draw_wait_message("Click 'Run MaAsLin2' to start analysis.")
          return(invisible(NULL))
        }
        if (isTRUE(maaslin2_running())) {
          draw_wait_message("MaAsLin2 is running. Please wait...")
          return(invisible(NULL))
        }
        bar_plot_reactive()
      },
      height = function() {
        if (is.null(input$run_maaslin2_btn) || input$run_maaslin2_btn < 1) {
          200
        } else {
          h <- suppressWarnings(as.numeric(input$plot_height))
          if (!is.finite(h) || is.na(h) || h <= 0) 400 else h
        }
      },
      width = function() {
        w <- suppressWarnings(as.numeric(input$plot_width))
        if (!is.finite(w) || is.na(w) || w <= 0) 700 else w
      }
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

    output$maaslin2_legend_box <- renderUI({
      req(input$plot_width)
      active_tab <- if (is.null(input$maaslin2_active_tab) || !nzchar(input$maaslin2_active_tab)) "Volcano Plot" else input$maaslin2_active_tab
      if (identical(active_tab, "Table")) {
        return(NULL)
      }
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 900 else input$plot_width
      tags$div(
        style = paste(
          "margin-top: 12px;",
          "clear: both;",
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
          uiOutput(session$ns("maaslin2_figure_legend"))
        )
      )
    })

    output$maaslin2_status_box <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 900 else input$plot_width
      tags$div(
        class = "simple-result-card",
        style = paste0("width: ", box_width, "px; max-width: 100%;"),
        verbatimTextOutput(session$ns("group_sample_counts"))
      )
    })

    output$maaslin2_status_separator <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 900 else input$plot_width
      tags$hr(
        style = paste(
          sprintf("width: %spx;", box_width),
          "max-width: 100%;",
          "margin: 14px 0 12px 0;",
          "border-top: 1px solid #d1d5db;"
        )
      )
    })

    output$maaslin2_figure_legend <- renderUI({
      req(input$tax_level, input$reference_level)
      active_tab <- if (is.null(input$maaslin2_active_tab) || !nzchar(input$maaslin2_active_tab)) "Volcano Plot" else input$maaslin2_active_tab
      y_metric_label <- if (identical(input$volcano_y_axis, "q_val")) {
        "FDR-adjusted p-value (q-value)"
      } else {
        "raw p-value"
      }
      legend_title <- if (identical(active_tab, "Bar Plot")) {
        "MaAsLin2 differential taxa bar plot"
      } else if (identical(active_tab, "Table")) {
        "MaAsLin2 results table"
      } else {
        "MaAsLin2 volcano plot"
      }
      legend_body <- if (identical(active_tab, "Bar Plot")) {
        paste0(
          "Bars show the top 10 taxa ranked by absolute log2 fold-change (log2FC) from MaAsLin2 at ",
          tolower(input$tax_level),
          " level against reference group ",
          input$reference_level,
          ". Red indicates taxa increased in the comparison group and blue indicates taxa decreased."
        )
      } else if (identical(active_tab, "Table")) {
        paste0(
          "This table reports MaAsLin2 coefficient estimates, transformed effect sizes, and significance values for selected contrasts at ",
          tolower(input$tax_level),
          " level. Reference group is ",
          input$reference_level,
          "."
        )
      } else {
        paste0(
          "Points represent taxa modeled by MaAsLin2 at ",
          tolower(input$tax_level),
          " level. The x-axis is log2FC relative to reference group ",
          input$reference_level,
          ", and the y-axis is -log10(",
          y_metric_label,
          "). Dashed lines indicate |log2FC| = 0.5 and p/q = 0.05 thresholds."
        )
      }
      tags$div(
        tags$div(
          style = "font-weight: 600; margin-bottom: 4px;",
          legend_title
        ),
        tags$div(legend_body)
      )
    })
  })
}
