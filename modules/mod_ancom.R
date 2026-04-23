library(shiny)
library(ggplot2)
library(dplyr)
library(phyloseq)

## UI
mod_ancom_ui <- function(id) {
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
        h4(icon("vial-circle-check"), "ANCOM-BC2"),
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
                    choices = c("q-value (FDR)" = "q_val",
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
            "7. Additional covariates for fix_formula (optional)",
            choices = NULL,
            multiple = TRUE,
            options = list(
              placeholder = "Select metadata covariates to adjust for",
              plugins = list("remove_button")
            )
          ),
          selectizeInput(
            ns("fix_interactions"),
            "8. Interaction terms for fix_formula (optional)",
            choices = NULL,
            multiple = TRUE,
            options = list(
              placeholder = "Select interaction terms among fixed effects",
              plugins = list("remove_button")
            )
          ),
          selectizeInput(
            ns("rand_covariates"),
            "9. Random-effect grouping variables (optional)",
            choices = NULL,
            multiple = TRUE,
            options = list(
              placeholder = "Select grouping variables for random intercepts",
              plugins = list("remove_button")
            )
          ),
          verbatimTextOutput(ns("formula_preview"))
        ),
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Plot Dimensions"),
        numericInput(ns("plot_width"), "Plot width (px)", value = 700, min = 400, max = 2000, step = 50),
        numericInput(ns("plot_height"), "Plot height (px)", value = 400, min = 300, max = 2000, step = 50),
        numericInput(ns("base_size"), "Base Font Size:", value = 11, min = 6, max = 30, step = 1),
        actionButton(ns("run_ancom_btn"), "Run ANCOM-BC2", class = "btn-danger", style = "font-size: 12px;"),
        tags$script(HTML(
          "Shiny.addCustomMessageHandler('toggle-ancom-run-btn', function(msg) {
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
        h4("ANCOM-BC2"),
        tags$div(
          id = ns("ancom_tab_container"),
          style = "max-width: 100%;",
          tabsetPanel(
            id = ns("ancom_active_tab"),
            tabPanel("Volcano plot",
                     downloadButton(ns("download_volcano"), "Download plot (PNG)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
                     tags$div(
                       style = "margin-top: 8px;",
                       plotOutput(ns("ancom_plot"), height = "auto")
                     )),
            tabPanel("Bar plot",
                     downloadButton(ns("download_barplot"), "Download plot (PNG)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
                     tags$div(
                       style = "margin-top: 8px;",
                       plotOutput(ns("ancom_barplot"), height = "auto")
                     )),
            tabPanel("Table",
                     downloadButton(ns("download_ancom_table"), "Download Table (TSV)", style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"),
                     DTOutput(ns("ancom_table")))
          )
        ),
        uiOutput(ns("ancom_legend_box")),
        uiOutput(ns("ancom_status_separator")),
        h5(icon("circle-info"), "ANCOM-BC2 Status"),
        uiOutput(ns("ancom_status_box"))
      )
    )
  )
}

## Server
mod_ancom_server <- function(id, ps_obj) {
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
        list(id = session$ns("ancom_tab_container"), width = paste0(width_px, "px"))
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
    
    build_interaction_choices <- function(group_var, all_metadata_cols) {
      if (is.null(group_var) || !nzchar(group_var)) {
        return(character(0))
      }
      all_vars <- setdiff(all_metadata_cols, c("SampleID", group_var))
      
      if (length(all_vars) == 0) {
        return(character(0))
      }
      
      # Single variables (for group_var * variable interactions)
      single_interactions <- all_vars
      
      # Multi-variable interactions
      if (length(all_vars) >= 2) {
        comb_mat <- utils::combn(all_vars, 2)
        multi_interactions <- apply(comb_mat, 2, function(x) paste(x[1], x[2], sep = ":"))
        return(c(single_interactions, multi_interactions))
      } else {
        return(single_interactions)
      }
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
      if (length(selected_levels) > 5) {
        selected_levels <- selected_levels[1:5]
        showNotification("You can select up to 5 group levels. Keeping the first 5 selected levels.", type = "warning")
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

      selected_rand_covariates <- input$rand_covariates
      if (is.null(selected_rand_covariates)) {
        selected_rand_covariates <- character(0)
      }
      updateSelectizeInput(
        session,
        "rand_covariates",
        choices = covariate_choices,
        selected = intersect(selected_rand_covariates, covariate_choices),
        server = TRUE
      )

      interaction_choices <- build_interaction_choices(group_var, colnames(meta_df))
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
    
    observeEvent(list(input$group_var, input$fix_covariates), {
      interaction_choices <- build_interaction_choices(group_var_resolved(), colnames(as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE)))
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

      ref_level <- input$reference_level
      non_ref_levels <- setdiff(selected_levels, ref_level)
      ordered_levels <- c(ref_level, non_ref_levels)
      group_var <- info$group_var
      group_vals <- as.character(phyloseq::sample_data(ps_sub)[[group_var]])
      phyloseq::sample_data(ps_sub)[[group_var]] <- factor(group_vals, levels = ordered_levels)

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

    ancom_running <- reactiveVal(FALSE)

    ancom_res <- eventReactive(input$run_ancom_btn, {
      req(ps_filtered(), input$group_var, input$reference_level)
      current_group_var <- group_var_resolved()
      ancom_running(TRUE)

      session$sendCustomMessage("toggle-ancom-run-btn", list(
        id = session$ns("run_ancom_btn"),
        disabled = TRUE,
        label = "Running..."
      ))
      on.exit({
        ancom_running(FALSE)
        session$sendCustomMessage("toggle-ancom-run-btn", list(
          id = session$ns("run_ancom_btn"),
          disabled = FALSE,
          label = "Run ANCOM-BC2"
        ))
      }, add = TRUE)

      result <- tryCatch({
        withProgress(message = "Running ANCOM-BC2...", value = 0, {
          ps_current <- ps_filtered()
          otu_for_var <- as.matrix(phyloseq::otu_table(ps_current))
          if (!phyloseq::taxa_are_rows(phyloseq::otu_table(ps_current))) {
            otu_for_var <- t(otu_for_var)
          }
          taxa_var <- apply(otu_for_var, 1, stats::var, na.rm = TRUE)
          keep_taxa <- rownames(otu_for_var)[is.finite(taxa_var) & taxa_var > 0]
          removed_taxa <- setdiff(phyloseq::taxa_names(ps_current), keep_taxa)
          if (length(removed_taxa) > 0) {
            ps_current <- phyloseq::prune_taxa(keep_taxa, ps_current)
            preview_removed <- paste(utils::head(removed_taxa, 5), collapse = ", ")
            more_suffix <- if (length(removed_taxa) > 5) paste0(" ... (+", length(removed_taxa) - 5, " more)") else ""
            showNotification(
              paste0(
                "Removed ",
                length(removed_taxa),
                " zero-variance taxa before ANCOM-BC2: ",
                preview_removed,
                more_suffix
              ),
              type = "warning",
              duration = 8
            )
          }
          validate(
            need(phyloseq::ntaxa(ps_current) > 0, "No taxa remain after removing zero-variance taxa.")
          )
          tax_level_arg <- if (input$tax_level == "ASV") NULL else input$tax_level
          selected_levels <- group_selection_info()$selected_levels
          use_multi_group_tests <- length(selected_levels) >= 3

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
          all_metadata_cols <- colnames(as.data.frame(phyloseq::sample_data(ps_current), stringsAsFactors = FALSE))
          valid_interactions <- build_interaction_choices(current_group_var, all_metadata_cols)
          interaction_terms <- intersect(interaction_terms, valid_interactions)
          
          # Process interaction terms: single variables and multi-variable interactions
          interaction_terms_formatted <- sapply(interaction_terms, function(term) {
            if (grepl(":", term)) {
              # Multi-variable interaction: Age:BMI → Age*BMI
              gsub(":", "*", term)
            } else {
              # Single variable: Age → group_var*Age
              paste0(current_group_var, "*", term)
            }
          }, USE.NAMES = FALSE)
          
          # Check if any interaction term already contains current_group_var
          has_group_interaction <- any(grepl(paste0("^", current_group_var, "\\*"), interaction_terms_formatted))
          
          # Build formula without duplicating current_group_var
          if (has_group_interaction) {
            fix_formula_str <- paste(c(covariates, interaction_terms_formatted), collapse = " + ")
          } else {
            fix_formula_str <- paste(c(current_group_var, covariates, interaction_terms_formatted), collapse = " + ")
          }

          rand_covariates <- input$rand_covariates
          if (is.null(rand_covariates) || length(rand_covariates) == 0) {
            rand_covariates <- character(0)
          }
          rand_covariates <- rand_covariates[nzchar(rand_covariates)]
          rand_covariates <- unique(setdiff(rand_covariates, current_group_var))
          rand_covariates <- setdiff(rand_covariates, covariates)
          rand_formula_str <- if (length(rand_covariates) > 0) {
            paste0("(1|", rand_covariates, ")", collapse = " + ")
          } else {
            NULL
          }

          out <- ANCOMBC::ancombc2(
            data = ps_current,
            tax_level = tax_level_arg,
            fix_formula = fix_formula_str,
            rand_formula = rand_formula_str,
            group = current_group_var,
            struc_zero = TRUE,
            neg_lb = TRUE,
            global = use_multi_group_tests,
            pairwise = use_multi_group_tests,
            trend = use_multi_group_tests,
            alpha = 0.05,
            n_cl = 4
          )

          res <- out$res
          res$feature_id <- res$taxon

          tax_table_current <- phyloseq::tax_table(ps_current)
          if (!is.null(tax_table_current)) {
            tax_df <- as.data.frame(tax_table_current)
            tax_df$feature_id <- rownames(tax_df)
            res <- dplyr::left_join(res, tax_df, by = "feature_id")
          }

          taxonomy_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          rank_default_unassigned <- c(
            Kingdom = "k__Unassigned",
            Phylum = "p__Unassigned",
            Class = "c__Unassigned",
            Order = "o__Unassigned",
            Family = "f__Unassigned",
            Genus = "g__Unassigned",
            Species = "s__Unassigned"
          )
          for (rk in taxonomy_ranks) {
            if (!rk %in% colnames(res)) {
              res[[rk]] <- NA_character_
            }
          }

          extract_rank_from_feature <- function(feature_vec, rank_name) {
            prefix_map <- list(
              Kingdom = c("k__", "d__", "kingdom__"),
              Phylum = c("p__", "phylum__"),
              Class = c("c__", "class__"),
              Order = c("o__", "order__"),
              Family = c("f__", "family__"),
              Genus = c("g__", "genus__"),
              Species = c("s__", "species__")
            )
            prefixes <- prefix_map[[rank_name]]
            out <- rep(NA_character_, length(feature_vec))
            if (is.null(prefixes) || length(prefixes) == 0) return(out)

            for (i in seq_along(feature_vec)) {
              fid <- as.character(feature_vec[i])
              if (is.na(fid) || !nzchar(trimws(fid))) next
              fid_norm <- trimws(fid)
              if (!grepl(";", fid_norm, fixed = TRUE) && grepl("__", fid_norm, fixed = TRUE)) {
                fid_norm <- gsub("(?<!^)([A-Za-z][A-Za-z0-9]*__)", "; \\1", fid_norm, perl = TRUE)
              }
              parts <- trimws(unlist(strsplit(fid_norm, ";", fixed = TRUE)))
              if (length(parts) == 0) next
              lower_parts <- tolower(parts)
              hit_mask <- Reduce(`|`, lapply(prefixes, function(px) {
                startsWith(lower_parts, tolower(px))
              }))
              hit_idx <- which(hit_mask)
              if (length(hit_idx) > 0) {
                token <- parts[hit_idx[1]]
                token <- gsub("_+$", "", token)
                out[i] <- token
              }
            }
            out
          }

          for (rk in taxonomy_ranks) {
            vals <- as.character(res[[rk]])
            vals <- gsub("_+$", "", vals)
            missing_idx <- is.na(vals) | !nzchar(trimws(vals))
            if (any(missing_idx) && "feature_id" %in% colnames(res)) {
              parsed_vals <- extract_rank_from_feature(res$feature_id, rk)
              vals[missing_idx] <- parsed_vals[missing_idx]
            }
            still_missing <- is.na(vals) | !nzchar(trimws(vals))
            vals[still_missing] <- rank_default_unassigned[[rk]]
            res[[rk]] <- vals
          }

          current_labels <- resolve_current_taxa_labels(res, input$tax_level)
          res$taxa_label <- current_labels
          res$taxon <- current_labels

          rownames(res) <- res$feature_id
          res
        })
      }, error = function(e) {
        showNotification(paste("ANCOM-BC2 Execution Failed:", e$message), type = "error", duration = NULL)
        NULL
      })

      validate(
        need(!is.null(result), "ANCOM-BC2 analysis failed. Check the error message above.")
      )

      result
    })

    output$formula_preview <- renderText({
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
      all_metadata_cols_preview <- colnames(as.data.frame(phyloseq::sample_data(ps_obj()), stringsAsFactors = FALSE))
      valid_interactions <- build_interaction_choices(group_var, all_metadata_cols_preview)
      interaction_terms <- intersect(interaction_terms, valid_interactions)
      
      # Process interaction terms for display: single variables and multi-variable interactions
      interaction_terms_formatted <- sapply(interaction_terms, function(term) {
        if (grepl(":", term)) {
          # Multi-variable interaction: Age:BMI → Age*BMI
          gsub(":", "*", term)
        } else {
          # Single variable: Age → group_var*Age
          paste0(group_var, "*", term)
        }
      }, USE.NAMES = FALSE)
      
      # Check if any interaction term already contains group_var
      has_group_interaction <- any(grepl(paste0("^", group_var, "\\*"), interaction_terms_formatted))
      
      # Build formula without duplicating group_var
      if (has_group_interaction) {
        fix_formula_str <- paste(c(covariates, interaction_terms_formatted), collapse = " + ")
      } else {
        fix_formula_str <- paste(c(group_var, covariates, interaction_terms_formatted), collapse = " + ")
      }
      
      rand_covariates <- input$rand_covariates
      if (is.null(rand_covariates) || length(rand_covariates) == 0) {
        rand_covariates <- character(0)
      }
      rand_covariates <- rand_covariates[nzchar(rand_covariates)]
      rand_covariates <- unique(setdiff(rand_covariates, group_var))
      
      rand_formula_str <- if (length(rand_covariates) > 0) {
        paste0("(1|", rand_covariates, ")", collapse = " + ")
      } else {
        "None"
      }

      paste(
        c(
          paste0("Applied fix_formula: ", fix_formula_str),
          paste0("Applied rand_formula: ", rand_formula_str)
        ),
        collapse = "\n"
      )
    })

    round_numeric_columns <- function(df, digits = 2, exclude_cols = c("p_val", "q_val")) {
      is_num_col <- vapply(df, is.numeric, logical(1))
      round_target_cols <- setdiff(names(df)[is_num_col], exclude_cols)
      df[round_target_cols] <- lapply(df[round_target_cols], function(x) round(x, digits))
      df
    }

    ancom_processed <- reactive({
      req(ancom_res(), input$group_var, input$reference_level)
      res <- ancom_res()
      current_group_var <- group_var_resolved()

      selected_levels <- group_selection_info()$selected_levels
      ref_level <- input$reference_level
      contrast_levels <- setdiff(selected_levels, ref_level)

      validate(
        need(length(contrast_levels) > 0, "ERROR: At least one non-reference group level is required.")
      )

      to_numeric <- function(values) suppressWarnings(as.numeric(values))
      pick_col <- function(df, names_to_try) {
        hit <- names_to_try[names_to_try %in% colnames(df)]
        if (length(hit) == 0) return(NA_character_)
        hit[1]
      }

      normalize_token <- function(x) {
        x <- as.character(x)
        x <- gsub("[^A-Za-z0-9]", "", x)
        tolower(x)
      }

      group_norm <- normalize_token(current_group_var)

      build_suffix_candidates <- function(level_name) {
        level_raw <- as.character(level_name)
        level_clean <- gsub("[^A-Za-z0-9]", "", level_raw)
        level_norm <- normalize_token(level_raw)
        unique(c(
          paste0(current_group_var, level_raw),
          paste0(current_group_var, "_", level_raw),
          paste0(current_group_var, level_clean),
          paste0(current_group_var, "_", level_clean),
          paste0(group_norm, level_norm),
          paste0(group_norm, "_", level_norm),
          level_raw,
          level_clean,
          level_norm,
          paste0(level_raw, "V"),
          paste0(level_clean, "V"),
          paste0(level_norm, "V")
        ))
      }

      all_cols_norm <- setNames(normalize_token(colnames(res)), colnames(res))

      get_col_for_level <- function(prefixes, suffix_candidates, fallbacks = character()) {
        direct_candidates <- unique(c(
          unlist(lapply(prefixes, function(pref) paste0(pref, suffix_candidates))),
          fallbacks
        ))
        direct_hit <- pick_col(res, direct_candidates)
        if (!is.na(direct_hit)) {
          return(direct_hit)
        }

        for (pref in prefixes) {
          pref_norm <- normalize_token(pref)
          suffix_norm <- unique(normalize_token(suffix_candidates))
          match_names <- names(all_cols_norm)[startsWith(all_cols_norm, pref_norm)]
          if (length(match_names) == 0) {
            next
          }
          tail_norm <- substring(all_cols_norm[match_names], nchar(pref_norm) + 1)
          hit_idx <- which(tail_norm %in% suffix_norm)
          if (length(hit_idx) > 0) {
            return(match_names[hit_idx[1]])
          }
        }

        NA_character_
      }

      taxonomy_cols <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(res))
      taxonomy_df <- if (length(taxonomy_cols) > 0) res[, taxonomy_cols, drop = FALSE] else NULL

      contrast_results <- lapply(contrast_levels, function(level_name) {
        suffix_candidates <- build_suffix_candidates(level_name)

        lfc_col <- get_col_for_level(c("lfc_"), suffix_candidates, fallbacks = c("lfc"))
        se_col <- get_col_for_level(c("se_"), suffix_candidates)
        w_col <- get_col_for_level(c("W_"), suffix_candidates, fallbacks = c("W"))
        p_col <- get_col_for_level(c("p_"), suffix_candidates, fallbacks = c("p_val", "pvalue"))
        q_col <- get_col_for_level(c("q_", "q_val_"), suffix_candidates, fallbacks = c("q_val", "padj"))
        diff_col <- get_col_for_level(c("diff_"), suffix_candidates, fallbacks = c("diff", "diff_abn"))

        if (is.na(lfc_col) || (is.na(p_col) && is.na(q_col))) {
          return(NULL)
        }

        lfc_values_raw <- to_numeric(res[[lfc_col]])
        lfc_values <- lfc_values_raw
        se_values <- if (!is.na(se_col)) to_numeric(res[[se_col]]) else rep(NA_real_, nrow(res))
        w_values <- if (!is.na(w_col)) to_numeric(res[[w_col]]) else rep(NA_real_, nrow(res))
        p_values <- if (!is.na(p_col)) to_numeric(res[[p_col]]) else rep(NA_real_, nrow(res))
        q_values <- if (!is.na(q_col)) to_numeric(res[[q_col]]) else p.adjust(p_values, method = "holm")

        diff_values <- if (!is.na(diff_col) && diff_col %in% colnames(res)) res[[diff_col]] else res[["diff_abn"]]
        if (is.null(diff_values)) {
          diff_values <- rep(FALSE, nrow(res))
        }
        diff_values <- as.logical(diff_values)
        diff_values[is.na(diff_values)] <- FALSE

        base_df <- data.frame(
          contrast = as.character(level_name),
          taxon = if ("taxon" %in% colnames(res)) res$taxon else res$taxa_label,
          taxa_label = res$taxa_label,
          feature_id = res$feature_id,
          lfc = lfc_values,
          se = se_values,
          W = w_values,
          p_val = p_values,
          q_val = q_values,
          diff = diff_values,
          stringsAsFactors = FALSE
        )

        if (!is.null(taxonomy_df)) {
          base_df <- cbind(taxonomy_df, base_df)
        }

        base_df
      })

      contrast_results <- Filter(Negate(is.null), contrast_results)
      validate(
        need(length(contrast_results) > 0,
             "No contrast-specific ANCOM-BC2 columns were found for selected levels and reference.")
      )

      base_df <- dplyr::bind_rows(contrast_results)
      base_df$contrast <- factor(base_df$contrast, levels = contrast_levels)
      contrast_label <- as.character(base_df$contrast)
      contrast_label[is.na(contrast_label) | !nzchar(contrast_label)] <- "selected comparison levels"
      base_df$direction <- ifelse(
        base_df$lfc > 0,
        paste0("Increase in ", contrast_label),
        paste0("Decrease in ", contrast_label)
      )
      base_df
    })

    output$ancom_table <- renderDT({
      req(ancom_processed())
      res <- ancom_processed()
      leading_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "contrast", "taxon", "taxa_label", "feature_id")
      existing_leading <- intersect(leading_cols, colnames(res))
      remaining_cols <- setdiff(colnames(res), existing_leading)
      res <- res[, c(existing_leading, remaining_cols), drop = FALSE]
      res <- round_numeric_columns(res, digits = 2)
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

    output$download_ancom_table <- downloadHandler(
      filename = function() {
        group_tag <- gsub("[^A-Za-z0-9_]+", "_", input$group_var)
        ref_tag <- gsub("[^A-Za-z0-9_]+", "_", input$reference_level)
        paste0("ancombc2_table_", group_tag, "_ref_", ref_tag, "_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        req(ancom_processed())
        rounded_res <- round_numeric_columns(ancom_processed(), digits = 2)
        readr::write_tsv(rounded_res, file)
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

      # Final fallback: use feature_id only (avoid falling back to full hierarchy labels).
      fallback_feature <- if ("feature_id" %in% colnames(df)) as.character(df$feature_id) else rep("", nrow(df))
      use_idx <- is.na(out) | !nzchar(out) | is_raw_placeholder_taxa(out)
      out[use_idx] <- fallback_feature[use_idx]
      out
    }

    volcano_plot_reactive <- reactive({
      req(ancom_processed(), input$volcano_y_axis, input$base_size)
      res <- ancom_processed()
      selected_levels <- group_selection_info()$selected_levels
      target_levels <- setdiff(selected_levels, input$reference_level)
      target_label <- if (length(target_levels) == 1) target_levels[1] else "selected comparison levels"

      y_col_name <- input$volcano_y_axis
      y_axis_label <- if (y_col_name == "q_val") expression("-log"[10]*"(FDR-adjusted p-value)") else expression("-log"[10]*"(p-value)")
      y_data <- res[[y_col_name]]

      plot_df <- res[!is.na(res$lfc) & !is.na(y_data) & y_data > 0, , drop = FALSE]

      validate(
        need(nrow(plot_df) > 0, "No valid rows to plot. Check ANCOM results.")
      )

      effect_cutoff <- 0.5
      plot_df$y_metric <- -log10(y_data[!is.na(res$lfc) & !is.na(y_data) & y_data > 0])
      is_significant <- !is.na(plot_df[[y_col_name]]) & is.finite(plot_df[[y_col_name]]) & (plot_df[[y_col_name]] < 0.05)
      plot_df$diffexpressed <- dplyr::case_when(
        is_significant & !is.na(plot_df$lfc) & is.finite(plot_df$lfc) & (plot_df$lfc >= effect_cutoff) ~ "Increased",
        is_significant & !is.na(plot_df$lfc) & is.finite(plot_df$lfc) & (plot_df$lfc <= -effect_cutoff) ~ "Decreased",
        TRUE ~ "Not significant"
      )
      plot_df$diffexpressed <- factor(
        plot_df$diffexpressed,
        levels = c("Decreased", "Not significant", "Increased")
      )
      label_by_rank <- resolve_current_taxa_labels(plot_df, input$tax_level)
      plot_df$delabel <- ifelse(plot_df$diffexpressed == "Not significant", "", label_by_rank)

      y_vals <- plot_df$y_metric[is.finite(plot_df$y_metric) & !is.na(plot_df$y_metric)]
      y_upper <- if (length(y_vals) > 0) {
        max(stats::quantile(y_vals, probs = 0.99, na.rm = TRUE), max(y_vals, na.rm = TRUE) * 0.1) + 1
      } else {
        10
      }
      y_upper <- max(5, y_upper)

      x_vals <- plot_df$lfc[is.finite(plot_df$lfc) & !is.na(plot_df$lfc)]
      x_abs_max <- if (length(x_vals) > 0) {
        as.numeric(stats::quantile(abs(x_vals), probs = 0.99, na.rm = TRUE)) + 0.5
      } else {
        2
      }
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
        labs(title = paste0(input$tax_level, "-Level ANCOM-BC2 Volcano plot (reference: ", input$reference_level, ")"),
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
      req(ancom_processed(), input$volcano_y_axis, input$base_size)
      res <- ancom_processed()
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

      res_top$direction_flag <- res_top$lfc > 0
      y_abs_max <- max(abs(res_top$lfc), na.rm = TRUE)
      if (!is.finite(y_abs_max) || y_abs_max <= 0) {
        y_abs_max <- 1
      }
      label_offset <- y_abs_max * 0.08
      label_right <- max(0, max(res_top$lfc, na.rm = TRUE)) + label_offset
      res_top$label_y <- label_right
      res_top$label_hjust <- 0
      y_limits <- c(min(res_top$lfc, na.rm = TRUE), label_right + (label_offset * 4))

      p <- ggplot(res_top, aes(x = reorder(taxa_label_current, lfc), y = lfc, fill = direction_flag)) +
        geom_col(color = "black") +
        geom_text(
          aes(y = label_y, label = stat_value_label, hjust = label_hjust),
          size = 3
        ) +
        coord_flip(clip = "off") +
        scale_y_continuous(limits = y_limits) +
        scale_fill_manual(
          values = c("TRUE" = "red", "FALSE" = "blue"),
          labels = c(
            "TRUE" = paste0("Increase in ", target_label),
            "FALSE" = paste0("Decrease in ", target_label)
          )
        ) +
        labs(title = paste0("Top ", top_n, " ", input$tax_level, "-Level Differential Taxa by LFC (reference: ", input$reference_level, ")"),
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

    output$ancom_plot <- renderPlot(
      {
        if (is.null(input$run_ancom_btn) || input$run_ancom_btn < 1) {
          draw_wait_message("Click 'Run ANCOM-BC2' to start analysis.")
          return(invisible(NULL))
        }
        if (isTRUE(ancom_running())) {
          draw_wait_message("ANCOM-BC2 is running. Please wait...")
          return(invisible(NULL))
        }
        volcano_plot_reactive()
      },
      height = function() {
        if (is.null(input$run_ancom_btn) || input$run_ancom_btn < 1) {
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
    output$ancom_barplot <- renderPlot(
      {
        if (is.null(input$run_ancom_btn) || input$run_ancom_btn < 1) {
          draw_wait_message("Click 'Run ANCOM-BC2' to start analysis.")
          return(invisible(NULL))
        }
        if (isTRUE(ancom_running())) {
          draw_wait_message("ANCOM-BC2 is running. Please wait...")
          return(invisible(NULL))
        }
        bar_plot_reactive()
      },
      height = function() {
        if (is.null(input$run_ancom_btn) || input$run_ancom_btn < 1) {
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
        paste0("ancombc2_volcano_", group_tag, "_ref_", ref_tag, ".png")
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
        paste0("ancombc2_barplot_", group_tag, "_ref_", ref_tag, ".png")
      },
      content = function(file) {
        req(input$plot_height, input$plot_width)
        ggplot2::ggsave(file, plot = bar_plot_reactive(), device = "png",
                        width = input$plot_width / 100, height = input$plot_height / 100)
      }
    )

    output$ancom_legend_box <- renderUI({
      req(input$plot_width)
      active_tab <- if (is.null(input$ancom_active_tab) || !nzchar(input$ancom_active_tab)) "Volcano plot" else input$ancom_active_tab
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
          uiOutput(session$ns("ancom_figure_legend"))
        )
      )
    })

    output$ancom_status_box <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 900 else input$plot_width
      tags$div(
        class = "simple-result-card",
        style = paste0("width: ", box_width, "px; max-width: 100%;"),
        verbatimTextOutput(session$ns("group_sample_counts"))
      )
    })

    output$ancom_status_separator <- renderUI({
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

    output$ancom_figure_legend <- renderUI({
      req(input$tax_level, input$reference_level)
      active_tab <- if (is.null(input$ancom_active_tab) || !nzchar(input$ancom_active_tab)) "Volcano plot" else input$ancom_active_tab
      y_metric_label <- if (identical(input$volcano_y_axis, "q_val")) {
        "FDR-adjusted p-value (q-value)"
      } else {
        "raw p-value"
      }
      legend_title <- if (identical(active_tab, "Bar plot")) {
        "ANCOM-BC2 differential taxa bar plot"
      } else if (identical(active_tab, "Table")) {
        "ANCOM-BC2 results table"
      } else {
        "ANCOM-BC2 volcano plot"
      }
      legend_body <- if (identical(active_tab, "Bar plot")) {
        paste0(
          "Bars show the top 10 taxa ranked by absolute log2 fold-change (log2FC) from ANCOM-BC2 at ",
          tolower(input$tax_level),
          " level against reference group ",
          input$reference_level,
          ". Red indicates taxa increased in the comparison group and blue indicates taxa decreased."
        )
      } else if (identical(active_tab, "Table")) {
        paste0(
          "This table reports ANCOM-BC2 coefficient estimates, test statistics, and significance values for selected contrasts at ",
          tolower(input$tax_level),
          " level. Reference group is ",
          input$reference_level,
          "."
        )
      } else {
        paste0(
          "Points represent taxa tested by ANCOM-BC2 at ",
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
