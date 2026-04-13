## UI
mod_beta_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML("
      .beta-result-card {
        width: 600px;
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
      .beta-result-card pre {
        margin: 0;
        padding: 0;
        border: 0;
        background: transparent;
        font-size: 12px;
        line-height: 1.45;
        white-space: pre-wrap;
      }
      .beta-result-card p {
        margin-bottom: 0.4rem;
      }
    ")),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        h4(icon("project-diagram"), "Beta Diversity"),
        hr(),
        
        uiOutput(ns("primary_group_selector")),       
        
        uiOutput(ns("secondary_group_selector")),
        selectInput(
          ns("norm_method"),
          "Normalization method:",
          choices = c("TSS" = "tss", "CLR" = "clr"),
          selected = "tss"
        ),
        selectInput(
          ns("beta_tax_level"),
          "Taxonomic level:",
          choices = c("ASV" = "ASV", "Species" = "Species", "Genus" = "Genus"),
          selected = "ASV"
        ),
        hr(),
        
        h4(icon("sliders"), "Plot Settings"),
        numericInput(ns("plot_width_px"), "Plot Width (pixels):", value = 600, min = 300, max = 1500, step = 50),
        numericInput(ns("plot_height_px"), "Plot Height (pixels):", value = 450, min = 300, max = 1500, step = 50),
        numericInput(ns("dot_size"), "Dot Size (point size):", value = 4, min = 0.5, max = 10, step = 0.5),
        checkboxInput(ns("show_dot_outline"), "Show Dot Outline", value = TRUE),
        checkboxInput(ns("show_ellipses"), "Show Group Ellipses", value = FALSE),
        checkboxInput(ns("show_sample_names"), "Show Sample Names", value = FALSE),
        uiOutput(ns("primary_color_controls")),
        tags$details(style = "margin-top: 8px; margin-bottom: 8px;",
          tags$summary(strong("Axis Limits")),
          br(),
          helpText("Leave blank to use automatic limits."),
          numericInput(ns("x_min"), "X-axis minimum:", value = NA),
          numericInput(ns("x_max"), "X-axis maximum:", value = NA),
          numericInput(ns("y_min"), "Y-axis minimum:", value = NA),
          numericInput(ns("y_max"), "Y-axis maximum:", value = NA)
        ),
        hr(),

        tags$details(style = "margin-top: 8px; margin-bottom: 8px;",
          tags$summary(strong("Clustering")),
          br(),
          selectInput(
            ns("cluster_mode"),
            "Clustering mode:",
            choices = c("Manual k" = "manual", "Auto k (Silhouette)" = "auto"),
            selected = "manual"
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'manual'", ns("cluster_mode")),
            numericInput(ns("cluster_k"), "Manual k:", value = 3, min = 2, max = 20, step = 1)
          ),
          conditionalPanel(
            condition = sprintf("input['%s'] == 'auto'", ns("cluster_mode")),
            numericInput(ns("cluster_k_min"), "Auto k min:", value = 2, min = 2, max = 20, step = 1),
            numericInput(ns("cluster_k_max"), "Auto k max:", value = 10, min = 3, max = 30, step = 1)
          ),
          checkboxInput(ns("show_cluster_labels"), "Show cluster labels on plot", value = TRUE),
          tags$div(
            style = "display: flex; gap: 8px; align-items: center; flex-wrap: wrap;",
            actionButton(ns("run_clustering"), "Run Clustering", class = "btn-secondary btn-sm", style = "font-size: 12px;"),
            downloadButton(ns("download_clustering_table"), "Download Cluster Assignment", class = "btn-secondary btn-sm", style = "font-size: 12px;")
          ),
          tags$script(HTML(
            "Shiny.addCustomMessageHandler('toggle-cluster-run-btn', function(msg) {
               var btn = document.getElementById(msg.id);
               if (!btn) return;
               btn.disabled = !!msg.disabled;
               if (msg.label) btn.textContent = msg.label;
            });"
          ))
        ),
        tags$details(style = "margin-top: 8px; margin-bottom: 8px;",
          tags$summary(strong("EnvFit")),
          br(),
          uiOutput(ns("envfit_var_selector")),
          uiOutput(ns("envfit_taxa_selector")),
          checkboxInput(ns("envfit_only_sig"), "Show only significant vectors", value = TRUE),
          numericInput(ns("envfit_p_cutoff"), "p-value cutoff:", value = 0.05, min = 0.001, max = 1, step = 0.01),
          checkboxInput(ns("show_envfit_vectors"), "Overlay EnvFit vectors on plots", value = TRUE),
          actionButton(ns("run_envfit"), "Run EnvFit", class = "btn-secondary btn-sm", style = "font-size: 12px;"),
          tags$script(HTML(
            "Shiny.addCustomMessageHandler('toggle-envfit-run-btn', function(msg) {
               var btn = document.getElementById(msg.id);
               if (!btn) return;
               btn.disabled = !!msg.disabled;
               if (msg.label) btn.textContent = msg.label;
             });"
          ))
        )
      ),
      
      mainPanel(
        width = 9,
        tabsetPanel(
          id = ns("beta_tabs"),
          tabPanel("PCoA",
                   h4(icon("chart-line"), "PCoA Ordination", style = "margin-top: 12px;"),
                   uiOutput(ns("pcoa_plot_ui")),
                   tags$div(style = "margin-top: 12px;",
                            downloadButton(
                              ns("download_pcoa"),
                              "Download PCoA Plot (PNG)",
                              style = "width: 240px; height: 40px; font-size: 12px; display: flex; align-items: center; justify-content: center;"
                            ))
          ),
          tabPanel("NMDS",
                   h4(icon("chart-line"), "NMDS Ordination", style = "margin-top: 12px;"),
                   uiOutput(ns("nmds_plot_ui")),
                   tags$div(
                     style = "margin-top: 12px; display: flex; gap: 10px; align-items: center;",
                     downloadButton(
                       ns("download_nmds"),
                       "Download NMDS Plot (PNG)",
                       style = "width: 240px; height: 40px; font-size: 12px; display: flex; align-items: center; justify-content: center;"
                     ),
                     actionButton(
                       ns("redraw_nmds"),
                       "Redraw Plot",
                       style = "width: 240px; height: 40px; font-size: 12px;"
                     )
                   )
          )
        ),
        
        hr(),
        h4(icon("flask"), "Statistical Test: PERMANOVA Results"),
        div(
          class = "beta-result-card",
          htmlOutput(ns("permanova_results_out"))
        ),
        h4(icon("object-group"), "Clustering Results"),
        div(
          class = "beta-result-card",
          verbatimTextOutput(ns("clustering_results_out"))
        ),
        h4(icon("arrows-to-dot"), "EnvFit Results"),
        div(
          class = "beta-result-card",
          verbatimTextOutput(ns("envfit_results_out"))
        )
      )
    )
  )
}

## Server
mod_beta_server <- function(id, ps_obj, meta_cols) {
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
    
    output$primary_group_selector <- renderUI({
      req(meta_cols())
      group_choices <- setdiff(meta_cols(), "SampleID")
      selected_col <- if (length(group_choices) > 0) group_choices[1] else NULL
      selectInput(session$ns("primary_group_var"), "Select Primary Variable (Color):",
                  choices = group_choices, selected = selected_col)
    })
    
    output$secondary_group_selector <- renderUI({
      req(meta_cols(), input$primary_group_var)
      resolved_primary <- resolve_meta_colname(input$primary_group_var, meta_cols())
      group_choices <- c("None", setdiff(meta_cols(), c("SampleID", input$primary_group_var, resolved_primary)))
      selected_col <- if (length(group_choices) > 1) group_choices[2] else group_choices[1]
      selectInput(session$ns("secondary_group_var"), "Select Secondary Variable (Shape):",
                  choices = group_choices, selected = "None")
    })
    
    primary_group_var <- reactive({
      req(input$primary_group_var, ps_obj())
      metadata <- as(phyloseq::sample_data(ps_obj()), "data.frame")
      resolve_meta_colname(input$primary_group_var, colnames(metadata))
    })
    secondary_group_var <- reactive({
      req(ps_obj())
      if (is.null(input$secondary_group_var) || identical(input$secondary_group_var, "None")) {
        return("None")
      }
      metadata <- as(phyloseq::sample_data(ps_obj()), "data.frame")
      resolve_meta_colname(input$secondary_group_var, colnames(metadata))
    })

    coerce_group_vars_to_factor <- function(ps_data, primary_var, secondary_var = "None") {
      if (is.null(ps_data)) return(ps_data)
      sdata <- phyloseq::sample_data(ps_data)
      sdata_df <- as(sdata, "data.frame")

      if (!is.null(primary_var) && primary_var %in% colnames(sdata_df)) {
        sdata_df[[primary_var]] <- as.factor(as.character(sdata_df[[primary_var]]))
      }
      if (!is.null(secondary_var) && !identical(secondary_var, "None") && secondary_var %in% colnames(sdata_df)) {
        sdata_df[[secondary_var]] <- as.factor(as.character(sdata_df[[secondary_var]]))
      }

      phyloseq::sample_data(ps_data) <- phyloseq::sample_data(sdata_df)
      ps_data
    }
    
    default_group_colors <- reactive({
      req(ps_obj(), primary_group_var())
      metadata <- as(phyloseq::sample_data(ps_obj()), "data.frame")
      validate(
        need(primary_group_var() %in% colnames(metadata), "Selected primary variable was not found in metadata.")
      )
      level_values <- unique(as.character(metadata[[primary_group_var()]]))
      level_values <- level_values[!is.na(level_values) & nzchar(level_values)]
      if (length(level_values) == 0) {
        return(list(levels = character(0), colors = character(0)))
      }
      default_colors <- grDevices::hcl.colors(length(level_values), palette = "Set 2")
      names(default_colors) <- level_values
      list(levels = level_values, colors = default_colors)
    })
    
    output$primary_color_controls <- renderUI({
      group_info <- default_group_colors()
      req(length(group_info$levels) > 0)
      ns <- session$ns
      
      controls <- lapply(seq_along(group_info$levels), function(i) {
        lvl <- group_info$levels[i]
        input_id <- paste0("primary_color_", i)
        input_id_ns <- ns(input_id)
        shiny_set_value_js <- sprintf(
          "Shiny.setInputValue('%s', this.value, {priority: 'event'})",
          input_id_ns
        )
        
        tagList(
          tags$label(`for` = input_id_ns, paste0("Color - ", lvl)),
          tags$input(
            id = input_id_ns,
            type = "color",
            value = group_info$colors[[lvl]],
            style = "width: 100%; height: 40px; padding: 2px;",
            oninput = shiny_set_value_js,
            onchange = shiny_set_value_js
          )
        )
      })
      
      tags$div(        
        tags$details(
          style = "margin-top: 6px; margin-bottom: 6px;",
          tags$summary(strong("Primary Group Colors")),
          br(),
          do.call(tagList, controls)
        )        
      )
    })
    
    primary_color_map <- reactive({
      group_info <- default_group_colors()
      req(length(group_info$levels) > 0)
      color_values <- vapply(seq_along(group_info$levels), function(i) {
        input_val <- input[[paste0("primary_color_", i)]]
        if (is.null(input_val) || !grepl("^#([A-Fa-f0-9]{6})$", input_val)) {
          group_info$colors[i]
        } else {
          input_val
        }
      }, FUN.VALUE = character(1))
      names(color_values) <- group_info$levels
      color_values
    })
    
    plot_shape_var <- reactive({
      if (secondary_group_var() == "None") { NULL } else { secondary_group_var() }
    })
    
    dist_method <- reactive({
      req(input$norm_method)
      if (identical(input$norm_method, "clr")) {
        "euclidean"
      } else {
        "bray"
      }
    })

    output$envfit_var_selector <- renderUI({
      req(data_prepared())
      meta_df <- data_prepared()$metadata
      candidate_vars <- setdiff(colnames(meta_df), "SampleID")
      selectizeInput(
        session$ns("envfit_vars"),
        "EnvFit variables (numeric + categorical):",
        choices = candidate_vars,
        selected = head(candidate_vars, min(1, length(candidate_vars))),
        multiple = TRUE,
        options = list(
          placeholder = "Select one or more metadata variables",
          plugins = list("remove_button")
        )
      )
    })

    output$envfit_taxa_selector <- renderUI({
      req(ps_normalized_local())
      otu_mat <- as(phyloseq::otu_table(ps_normalized_local()), "matrix")
      if (!phyloseq::taxa_are_rows(ps_normalized_local())) {
        otu_mat <- t(otu_mat)
      }
      taxa_choices <- rownames(otu_mat)
      taxa_choices <- taxa_choices[order(rowSums(otu_mat, na.rm = TRUE), decreasing = TRUE)]

      selectizeInput(
        session$ns("envfit_taxa_vars"),
        "Additional taxa abundance variables (optional):",
        choices = taxa_choices,
        selected = character(0),
        multiple = TRUE,
        options = list(
          placeholder = "Select taxa to include as EnvFit variables",
          plugins = list("remove_button")
        )
      )
    })

    observeEvent(data_prepared(), {
      meta_df <- data_prepared()$metadata
      candidate_vars <- setdiff(colnames(meta_df), "SampleID")
      current_sel <- input$envfit_vars
      if (is.null(current_sel)) current_sel <- character(0)
      current_sel <- intersect(current_sel, candidate_vars)
      if (length(current_sel) == 0 && length(candidate_vars) > 0) {
        current_sel <- candidate_vars[1]
      }
      updateSelectizeInput(
        session,
        "envfit_vars",
        choices = candidate_vars,
        selected = current_sel,
        server = TRUE
      )
    }, ignoreInit = FALSE)

    distance_label <- reactive({
      if (identical(dist_method(), "euclidean")) {
        "Aitchison distance (CLR-transformed data)"
      } else {
        "Bray-Curtis"
      }
    })

    ps_normalized_local <- reactive({
      req(ps_obj(), input$norm_method, input$beta_tax_level)
      ps_data <- ps_obj()
      validate(
        need(phyloseq::nsamples(ps_data) > 0, "No samples available for beta diversity."),
        need(phyloseq::ntaxa(ps_data) > 0, "No taxa available for beta diversity.")
      )

      tax_level <- input$beta_tax_level
      if (!identical(tax_level, "ASV")) {
        validate(
          need(!is.null(phyloseq::tax_table(ps_data)), "Taxonomy table is required for taxonomic aggregation.")
        )
        tax_cols <- colnames(phyloseq::tax_table(ps_data))
        validate(
          need(tax_level %in% tax_cols, paste0("Taxonomic rank '", tax_level, "' is not available in taxonomy table."))
        )
        ps_data <- phyloseq::tax_glom(ps_data, taxrank = tax_level, NArm = FALSE)

        tax_df <- as.data.frame(phyloseq::tax_table(ps_data), stringsAsFactors = FALSE)
        rank_values <- as.character(tax_df[[tax_level]])
        missing_idx <- is.na(rank_values) | rank_values == ""
        rank_values[missing_idx] <- phyloseq::taxa_names(ps_data)[missing_idx]
        phyloseq::taxa_names(ps_data) <- make.unique(rank_values)
      }

      if (identical(input$norm_method, "clr")) {
        phyloseq::transform_sample_counts(ps_data, function(x) {
          x_num <- as.numeric(x)
          names(x_num) <- names(x)
          log_x <- log(x_num + 1)
          clr_x <- log_x - mean(log_x)
          names(clr_x) <- names(x)
          clr_x
        })
      } else {
        phyloseq::transform_sample_counts(ps_data, function(x) {
          total <- sum(x)
          if (total <= 0) x else x / total
        })
      }
    })

    data_prepared <- reactive({
      req(ps_normalized_local(), dist_method(), primary_group_var())
      ps_data <- ps_normalized_local()
      ps_data <- coerce_group_vars_to_factor(ps_data, primary_group_var(), secondary_group_var())
      dist_mat <- phyloseq::distance(ps_data, method = dist_method())
      metadata <- as(phyloseq::sample_data(ps_data), "data.frame")
      list(ps = ps_data, dist_mat = dist_mat, metadata = metadata, dist_method = dist_method())
    })

    cluster_result_val <- reactiveVal(list(
      result = NULL,
      note = "Clustering has not been run yet."
    ))

    observeEvent(input$run_clustering, {
      if (!requireNamespace("cluster", quietly = TRUE)) {
        cluster_result_val(list(
          result = NULL,
          note = "The 'cluster' package is required for silhouette-based clustering."
        ))
        return()
      }

      if (is.null(data_prepared())) {
        cluster_result_val(list(
          result = NULL,
          note = "Data is not ready for clustering yet."
        ))
        return()
      }

      session$sendCustomMessage("toggle-cluster-run-btn", list(
        id = session$ns("run_clustering"),
        disabled = TRUE,
        label = "Running..."
      ))
      on.exit({
        session$sendCustomMessage("toggle-cluster-run-btn", list(
          id = session$ns("run_clustering"),
          disabled = FALSE,
          label = "Run Clustering"
        ))
      }, add = TRUE)

      run_result <- tryCatch({
        dmat <- as.matrix(data_prepared()$dist_mat)
        n <- nrow(dmat)

        validate(
          need(n >= 3, "At least 3 samples are required for clustering."),
          need(!is.null(rownames(dmat)), "Sample IDs are missing in distance matrix.")
        )

        run_pam <- function(k_val) {
          pam_fit <- cluster::pam(dmat, k = k_val, diss = TRUE)
          sil <- cluster::silhouette(pam_fit$clustering, dist = as.dist(dmat))
          mean_sil <- mean(sil[, "sil_width"], na.rm = TRUE)
          list(fit = pam_fit, mean_sil = mean_sil)
        }

        if (identical(input$cluster_mode, "auto")) {
          k_min <- max(2, as.integer(input$cluster_k_min))
          k_max <- max(k_min, as.integer(input$cluster_k_max))
          k_max <- min(k_max, n - 1)
          k_min <- min(k_min, k_max)

          candidate_k <- seq(k_min, k_max)
          validate(
            need(length(candidate_k) > 0, "No valid k range for auto clustering. Adjust k min/max.")
          )

          candidate_res <- lapply(candidate_k, run_pam)
          sil_values <- vapply(candidate_res, function(x) x$mean_sil, numeric(1))
          best_idx <- which.max(sil_values)
          best_k <- candidate_k[best_idx]
          best <- candidate_res[[best_idx]]

          list(
            mode = "auto",
            k = best_k,
            silhouette = best$mean_sil,
            cluster = as.integer(best$fit$clustering),
            sample_ids = names(best$fit$clustering),
            k_min = k_min,
            k_max = k_max
          )
        } else {
          k_val <- max(2, min(as.integer(input$cluster_k), n - 1))
          one <- run_pam(k_val)
          list(
            mode = "manual",
            k = k_val,
            silhouette = one$mean_sil,
            cluster = as.integer(one$fit$clustering),
            sample_ids = names(one$fit$clustering)
          )
        }
      }, error = function(e) {
        cluster_result_val(list(
          result = NULL,
          note = paste0("Clustering failed: ", conditionMessage(e))
        ))
        NULL
      })

      if (!is.null(run_result)) {
        cluster_result_val(list(
          result = run_result,
          note = "Clustering completed successfully."
        ))
      }
    }, ignoreInit = TRUE)

    cluster_lookup <- reactive({
      cluster_obj <- cluster_result_val()
      res <- cluster_obj$result
      req(res)
      data.frame(
        SampleID = as.character(res$sample_ids),
        Cluster = as.factor(res$cluster),
        stringsAsFactors = FALSE
      )
    })

    envfit_result_val <- reactiveVal(list(
      pcoa = data.frame(),
      nmds = data.frame(),
      note = "EnvFit has not been run yet."
    ))

    observeEvent(input$run_envfit, {
      if (!requireNamespace("vegan", quietly = TRUE)) {
        envfit_result_val(list(
          pcoa = data.frame(),
          nmds = data.frame(),
          note = "The 'vegan' package is required for EnvFit."
        ))
        return()
      }

      if (is.null(data_prepared())) {
        envfit_result_val(list(
          pcoa = data.frame(),
          nmds = data.frame(),
          note = "Data is not ready for EnvFit yet."
        ))
        return()
      }

      vars_selected <- input$envfit_vars
      if (is.null(vars_selected)) vars_selected <- character(0)
      vars_selected <- vars_selected[nzchar(vars_selected)]
      taxa_selected <- input$envfit_taxa_vars
      if (is.null(taxa_selected)) taxa_selected <- character(0)
      taxa_selected <- taxa_selected[nzchar(taxa_selected)]

      if (length(vars_selected) == 0 && length(taxa_selected) == 0) {
        meta_df_default <- data_prepared()$metadata
        fallback_vars <- setdiff(colnames(meta_df_default), "SampleID")
        fallback_vars <- fallback_vars[vapply(fallback_vars, function(vn) {
          vv <- meta_df_default[[vn]]
          if (is.numeric(vv)) {
            is.finite(stats::sd(vv, na.rm = TRUE)) && stats::sd(vv, na.rm = TRUE) > 0
          } else {
            ux <- unique(as.character(vv))
            ux <- ux[!is.na(ux) & nzchar(ux)]
            length(ux) > 1
          }
        }, logical(1))]
        vars_selected <- head(fallback_vars, min(1, length(fallback_vars)))
        if (length(vars_selected) > 0) {
          updateSelectizeInput(session, "envfit_vars", selected = vars_selected, server = TRUE)
        } else {
          envfit_result_val(list(
            pcoa = data.frame(),
            nmds = data.frame(),
            note = "No valid metadata variables selected for EnvFit."
          ))
          return()
        }
      }

      session$sendCustomMessage("toggle-envfit-run-btn", list(
        id = session$ns("run_envfit"),
        disabled = TRUE,
        label = "Running..."
      ))
      on.exit({
        session$sendCustomMessage("toggle-envfit-run-btn", list(
          id = session$ns("run_envfit"),
          disabled = FALSE,
          label = "Run EnvFit"
        ))
      }, add = TRUE)

      ps_for_envfit <- tryCatch(ps_plot_obj(), error = function(e) data_prepared()$ps)

      build_envfit_df <- function(ord_obj, axis_x, axis_y) {
        if (is.null(ord_obj)) {
          return(data.frame())
        }
        ord_df <- phyloseq::plot_ordination(ps_for_envfit, ord_obj, justDF = TRUE)
        if (!all(c(axis_x, axis_y, "SampleID") %in% colnames(ord_df))) {
          return(data.frame())
        }

        ord_df$SampleID <- as.character(ord_df$SampleID)
        ord_df <- ord_df[!is.na(ord_df[[axis_x]]) & !is.na(ord_df[[axis_y]]), , drop = FALSE]

        meta_df <- data_prepared()$metadata
        meta_df$SampleID <- rownames(meta_df)
        common_samples <- intersect(ord_df$SampleID, meta_df$SampleID)
        if (length(common_samples) < 3) {
          return(data.frame())
        }

        ord_sub <- ord_df[match(common_samples, ord_df$SampleID), , drop = FALSE]
        env_sub <- if (length(vars_selected) > 0) {
          meta_df[match(common_samples, meta_df$SampleID), vars_selected, drop = FALSE]
        } else {
          data.frame(row.names = common_samples)
        }
        env_sub <- as.data.frame(env_sub, stringsAsFactors = FALSE)

        if (length(taxa_selected) > 0) {
          otu_mat <- as(phyloseq::otu_table(ps_normalized_local()), "matrix")
          if (!phyloseq::taxa_are_rows(ps_normalized_local())) {
            otu_mat <- t(otu_mat)
          }
          taxa_selected_safe <- intersect(taxa_selected, rownames(otu_mat))
          if (length(taxa_selected_safe) > 0) {
            taxa_df <- as.data.frame(t(otu_mat[taxa_selected_safe, , drop = FALSE]), stringsAsFactors = FALSE)
            taxa_df$SampleID <- rownames(taxa_df)
            taxa_df <- taxa_df[match(common_samples, taxa_df$SampleID), , drop = FALSE]
            taxa_df$SampleID <- NULL
            colnames(taxa_df) <- paste0("Taxa::", colnames(taxa_df))
            env_sub <- cbind(env_sub, taxa_df)
          }
        }

        primary_var <- primary_group_var()
        secondary_var <- secondary_group_var()
        if (!is.null(primary_var) && primary_var %in% colnames(env_sub)) {
          env_sub[[primary_var]] <- as.factor(as.character(env_sub[[primary_var]]))
        }
        if (!is.null(secondary_var) && !identical(secondary_var, "None") && secondary_var %in% colnames(env_sub)) {
          env_sub[[secondary_var]] <- as.factor(as.character(env_sub[[secondary_var]]))
        }

        for (cn in colnames(env_sub)) {
          if (!is.numeric(env_sub[[cn]])) {
            env_sub[[cn]] <- as.factor(as.character(env_sub[[cn]]))
          }
        }

        keep_valid <- vapply(env_sub, function(x) {
          if (is.numeric(x)) {
            is.finite(stats::sd(x, na.rm = TRUE)) && stats::sd(x, na.rm = TRUE) > 0
          } else {
            ux <- unique(as.character(x))
            ux <- ux[!is.na(ux) & nzchar(ux)]
            length(ux) > 1
          }
        }, logical(1))
        env_sub <- env_sub[, keep_valid, drop = FALSE]
        if (ncol(env_sub) == 0) {
          return(data.frame())
        }

        complete_idx <- stats::complete.cases(env_sub)
        env_sub <- env_sub[complete_idx, , drop = FALSE]
        ord_sub <- ord_sub[complete_idx, , drop = FALSE]
        if (nrow(env_sub) < 3) {
          return(data.frame())
        }

        ord_mat <- as.matrix(ord_sub[, c(axis_x, axis_y), drop = FALSE])
        rownames(ord_mat) <- ord_sub$SampleID

        fit <- vegan::envfit(ord_mat, env_sub, permutations = 999, na.rm = TRUE)
        vec_num <- tryCatch(vegan::scores(fit, display = "vectors"), error = function(e) NULL)
        vec_fac <- tryCatch(vegan::scores(fit, display = "factors"), error = function(e) NULL)

        vec_num <- if (is.null(vec_num)) data.frame() else as.data.frame(vec_num)
        vec_fac <- if (is.null(vec_fac)) data.frame() else as.data.frame(vec_fac)

        out_parts <- list()

        if (nrow(vec_num) > 0) {
          vec_num$variable <- rownames(vec_num)
          vec_num$r2 <- as.numeric(fit$vectors$r[vec_num$variable])
          vec_num$p_value <- as.numeric(fit$vectors$pvals[vec_num$variable])
          names(vec_num)[1:2] <- c("v1", "v2")
          vec_num$type <- "numeric"
          out_parts[[length(out_parts) + 1]] <- vec_num
        }

        if (nrow(vec_fac) > 0) {
          vec_fac$level <- rownames(vec_fac)
          vec_fac$variable <- sub("^([^:]+):.*$", "\\1", vec_fac$level)
          missing_var <- !(vec_fac$variable %in% names(fit$factors$pvals))
          if (any(missing_var)) {
            vec_fac$variable[missing_var] <- sub("^(.*?)\\..*$", "\\1", vec_fac$level[missing_var])
          }
          vec_fac$r2 <- as.numeric(fit$factors$r[vec_fac$variable])
          vec_fac$p_value <- as.numeric(fit$factors$pvals[vec_fac$variable])
          names(vec_fac)[1:2] <- c("v1", "v2")
          vec_fac$type <- "factor"
          out_parts[[length(out_parts) + 1]] <- vec_fac
        }

        if (length(out_parts) == 0) {
          return(data.frame())
        }

        vec <- dplyr::bind_rows(out_parts)

        if (isTRUE(input$envfit_only_sig)) {
          cutoff <- as.numeric(input$envfit_p_cutoff)
          if (is.na(cutoff)) cutoff <- 0.05
          vec <- vec[!is.na(vec$p_value) & vec$p_value <= cutoff, , drop = FALSE]
        }

        if (nrow(vec) == 0) {
          return(vec)
        }

        vec$is_taxa <- startsWith(as.character(vec$variable), "Taxa::")
        scale_ref <- vec[!vec$is_taxa, , drop = FALSE]
        if (nrow(scale_ref) == 0) {
          scale_ref <- vec
        }

        xr <- range(ord_sub[[axis_x]], na.rm = TRUE)
        yr <- range(ord_sub[[axis_y]], na.rm = TRUE)
        max_v1 <- max(abs(scale_ref$v1), na.rm = TRUE)
        max_v2 <- max(abs(scale_ref$v2), na.rm = TRUE)
        max_v1 <- ifelse(is.finite(max_v1) && max_v1 > 0, max_v1, 1)
        max_v2 <- ifelse(is.finite(max_v2) && max_v2 > 0, max_v2, 1)
        scale_factor <- min(diff(xr) / (2 * max_v1), diff(yr) / (2 * max_v2))
        if (!is.finite(scale_factor) || scale_factor <= 0) scale_factor <- 1
        scale_factor <- scale_factor * 0.8

        vec$x <- 0
        vec$y <- 0
        vec$xend <- vec$v1 * scale_factor
        vec$yend <- vec$v2 * scale_factor
        vec$label <- ifelse(
          vec$type == "factor" & "level" %in% colnames(vec),
          paste0(vec$level, " (p=", format(round(vec$p_value, 3), nsmall = 3), ")"),
          paste0(vec$variable, " (p=", format(round(vec$p_value, 3), nsmall = 3), ")")
        )
        vec
      }

      pcoa_ord <- tryCatch(pcoa_ordination_reactive(), error = function(e) NULL)
      if (is.null(pcoa_ord)) {
        pcoa_ord <- tryCatch(
          phyloseq::ordinate(ps_for_envfit, method = "PCoA", distance = dist_method()),
          error = function(e) NULL
        )
      }

      nmds_ord <- tryCatch(nmds_ordination_reactive(), error = function(e) NULL)
      if (is.null(nmds_ord)) {
        nmds_ord <- tryCatch(
          phyloseq::ordinate(ps_for_envfit, method = "NMDS", distance = dist_method()),
          error = function(e) NULL
        )
      }

      out <- tryCatch({
        pcoa_vec <- build_envfit_df(pcoa_ord, "Axis.1", "Axis.2")
        nmds_vec <- build_envfit_df(nmds_ord, "NMDS1", "NMDS2")
        list(pcoa = pcoa_vec, nmds = nmds_vec, note = NULL)
      }, error = function(e) {
        list(
          pcoa = data.frame(),
          nmds = data.frame(),
          note = paste0("EnvFit execution issue: ", conditionMessage(e))
        )
      })

      envfit_result_val(out)
    }, ignoreInit = TRUE)

    output$download_clustering_table <- downloadHandler(
      filename = function() {
        paste0("beta_clustering_assignments_", Sys.Date(), ".tsv")
      },
      content = function(file) {
        cluster_obj <- cluster_result_val()
        req(cluster_obj$result)
        res <- cluster_obj$result
        assignment_df <- data.frame(
          SampleID = as.character(res$sample_ids),
          Cluster = as.integer(res$cluster),
          stringsAsFactors = FALSE
        )
        assignment_df <- assignment_df[order(assignment_df$Cluster, assignment_df$SampleID), , drop = FALSE]
        assignment_df$Mode <- res$mode
        assignment_df$Selected_k <- res$k
        assignment_df$Average_Silhouette <- round(res$silhouette, 6)
        readr::write_tsv(assignment_df, file)
      }
    )
    
    ps_plot_obj <- reactive({
      req(data_prepared())
      
      ps_data <- data_prepared()$ps      
      shape_var <- plot_shape_var()
      
      if (!is.null(shape_var)) {
        sdata <- phyloseq::sample_data(ps_data)
        
        if (shape_var %in% names(sdata)) {
          current_col <- sdata[[shape_var]]
          
          if (!is.factor(current_col) && (is.numeric(current_col) || is.character(current_col))) {
            sdata[[shape_var]] <- as.factor(current_col)
            phyloseq::sample_data(ps_data) <- sdata
          }
        }
      }
      return(ps_data)
    })
    
    output$pcoa_plot_ui <- renderUI({ 
      req(input$plot_width_px, input$plot_height_px)
      plotOutput(session$ns("pcoa_plot_out"),
                 height = paste0(input$plot_height_px, "px"),
                 width = paste0(input$plot_width_px, "px"))
    })
    
    output$nmds_plot_ui <- renderUI({ 
      req(input$plot_width_px, input$plot_height_px)
      plotOutput(session$ns("nmds_plot_out"),
                 height = paste0(input$plot_height_px, "px"),
                 width = paste0(input$plot_width_px, "px"))
    })
    
    axis_limits <- reactive({
      build_limits <- function(min_val, max_val) {
        if (all(is.na(c(min_val, max_val)))) { NULL } else { c(min_val, max_val) }
      }
      list(x = build_limits(input$x_min, input$x_max),
           y = build_limits(input$y_min, input$y_max))
    })
    
    
    pcoa_ordination_reactive <- reactive({
      req(ps_plot_obj())
      phyloseq::ordinate(ps_plot_obj(), method = "PCoA", distance = dist_method())
    })
    
    pcoa_plot_reactive <- reactive({
      req(ps_plot_obj(), primary_group_var(), ord <- pcoa_ordination_reactive(), input$dot_size)
      ps_data_for_plot <- ps_plot_obj()
      
      if (!is.null(ord$values) && "Eigenvalues" %in% colnames(ord$values)) {
        eigenvalues <- ord$values$Eigenvalues
        total_var <- sum(eigenvalues)
        pc1_percent <- round(eigenvalues[1] / total_var * 100, 2)
        pc2_percent <- round(eigenvalues[2] / total_var * 100, 2)
        x_axis_label <- paste0("PC1 (", pc1_percent, "%)")
        y_axis_label <- paste0("PC2 (", pc2_percent, "%)")
      } else {
        x_axis_label <- "PC1"
        y_axis_label <- "PC2"
      }
      
      if (is.null(plot_shape_var())) {
        p <- phyloseq::plot_ordination(ps_data_for_plot, ord, color = primary_group_var())
      } else {
        p <- phyloseq::plot_ordination(ps_data_for_plot, ord, color = primary_group_var(), shape = plot_shape_var())
      }
      if (isTRUE(input$show_dot_outline)) {
        outline_size <- input$dot_size + 0.4
        p <- p + ggplot2::geom_point(
          color = "black",
          size = outline_size,
          show.legend = FALSE
        )
      }
      p <- p + ggplot2::geom_point(size = input$dot_size)
      p <- p + ggplot2::scale_color_manual(values = primary_color_map())
      
      if (input$show_ellipses) {
        p <- p + ggplot2::stat_ellipse(ggplot2::aes_string(group = primary_group_var()))
      }
      
      p <- p + ggplot2::theme_bw() +
        ggplot2::labs(title = paste("PCoA -", distance_label()),
                      color = primary_group_var(),
                      shape = if (is.null(plot_shape_var())) NULL else plot_shape_var(),
                      x = x_axis_label,
                      y = y_axis_label)
      
      plot_data <- phyloseq::plot_ordination(ps_data_for_plot, ord, justDF = TRUE)
      
      if (input$show_sample_names) {
        p <- p + ggplot2::geom_text(data = plot_data, mapping = ggplot2::aes(x = Axis.1, y = Axis.2, label = SampleID), color = "black", size = 3, hjust = 0.2, vjust = 0.5)
      }

      if (isTRUE(input$show_cluster_labels) && !is.null(cluster_result_val()$result)) {
        clu_df <- tryCatch(cluster_lookup(), error = function(e) NULL)
        if (!is.null(clu_df)) {
        plot_data$SampleID <- as.character(plot_data$SampleID)
        plot_data <- dplyr::left_join(plot_data, clu_df, by = "SampleID")
        if ("Cluster" %in% colnames(plot_data)) {
          plot_data$ClusterLabel <- ifelse(is.na(plot_data$Cluster), "", paste0("C", plot_data$Cluster))
          p <- p + ggplot2::geom_text(
            data = plot_data,
            mapping = ggplot2::aes(x = Axis.1, y = Axis.2, label = ClusterLabel),
            color = "black", size = 2.7, vjust = -0.8, show.legend = FALSE
          )
        }
        }
      }

      if (isTRUE(input$show_envfit_vectors) && isTRUE(input$run_envfit > 0)) {
        vec_obj <- envfit_result_val()
        vec_df <- if (is.null(vec_obj)) NULL else vec_obj$pcoa
        if (!is.null(vec_df) && nrow(vec_df) > 0) {
          x_rng <- range(plot_data$Axis.1, na.rm = TRUE)
          y_rng <- range(plot_data$Axis.2, na.rm = TRUE)
          x_pad <- ifelse(is.finite(diff(x_rng)) && diff(x_rng) > 0, diff(x_rng) * 0.02, 0)
          y_pad <- ifelse(is.finite(diff(y_rng)) && diff(y_rng) > 0, diff(y_rng) * 0.02, 0)
          x_min <- x_rng[1] + x_pad
          x_max <- x_rng[2] - x_pad
          y_min <- y_rng[1] + y_pad
          y_max <- y_rng[2] - y_pad

          taxa_flag <- if ("is_taxa" %in% colnames(vec_df)) as.logical(vec_df$is_taxa) else rep(FALSE, nrow(vec_df))
          taxa_flag[is.na(taxa_flag)] <- FALSE
          vec_arrow <- vec_df[!taxa_flag, , drop = FALSE]
          vec_point <- vec_df[taxa_flag, , drop = FALSE]

          if (nrow(vec_arrow) > 0) {
          vec_arrow$xend_plot <- pmin(pmax(vec_arrow$xend, x_min), x_max)
          vec_arrow$yend_plot <- pmin(pmax(vec_arrow$yend, y_min), y_max)
          p <- p +
            ggplot2::geom_segment(
              data = vec_arrow,
              ggplot2::aes(x = x, y = y, xend = xend_plot, yend = yend_plot),
              inherit.aes = FALSE,
              arrow = grid::arrow(length = grid::unit(0.18, "cm")),
              color = "#2E7D32",
              linewidth = 0.6
            )
          if (requireNamespace("ggrepel", quietly = TRUE)) {
            p <- p + ggrepel::geom_text_repel(
              data = vec_arrow,
              ggplot2::aes(x = xend_plot, y = yend_plot, label = label),
              inherit.aes = FALSE,
              color = "#1B5E20",
              size = 3,
              min.segment.length = 0,
              box.padding = 0.2,
              point.padding = 0.15,
              max.overlaps = Inf
            )
          } else {
            p <- p + ggplot2::geom_text(
              data = vec_arrow,
              ggplot2::aes(x = xend_plot, y = yend_plot, label = label),
              inherit.aes = FALSE,
              color = "#1B5E20",
              size = 3,
              hjust = -0.05,
              check_overlap = TRUE
            )
          }
          }

          if (nrow(vec_point) > 0) {
            vec_point$x_plot <- pmin(pmax(vec_point$xend, x_min), x_max)
            vec_point$y_plot <- pmin(pmax(vec_point$yend, y_min), y_max)

            p <- p +
              ggplot2::geom_point(
                data = vec_point,
                ggplot2::aes(x = x_plot, y = y_plot),
                inherit.aes = FALSE,
                shape = 21,
                fill = "#7B1FA2",
                color = "white",
                stroke = 0.4,
                size = 2.8
              )
            if (requireNamespace("ggrepel", quietly = TRUE)) {
              p <- p + ggrepel::geom_text_repel(
                data = vec_point,
                ggplot2::aes(x = x_plot, y = y_plot, label = label),
                inherit.aes = FALSE,
                color = "#4A148C",
                size = 2.8,
                min.segment.length = 0,
                box.padding = 0.2,
                point.padding = 0.15,
                max.overlaps = Inf
              )
            } else {
              p <- p + ggplot2::geom_text(
                data = vec_point,
                ggplot2::aes(x = x_plot, y = y_plot, label = label),
                inherit.aes = FALSE,
                color = "#4A148C",
                size = 2.8,
                hjust = -0.05,
                vjust = -0.2,
                check_overlap = TRUE
              )
            }
          }
        }
      }
      
      
      limits <- axis_limits()
      if (!is.null(limits$x) || !is.null(limits$y)) {
        p <- p + ggplot2::coord_cartesian(xlim = limits$x, ylim = limits$y)
      }
      
      return(p)
    })
    
    output$pcoa_plot_out <- renderPlot({ pcoa_plot_reactive() })
    
    output$download_pcoa <- downloadHandler(
      filename = function() { paste0("PCoA_", gsub("-", "", distance_label()), "_", Sys.Date(), ".png") },
      content = function(file) {
        dpi_val <- 300
        width_in <- input$plot_width_px / dpi_val
        height_in <- input$plot_height_px / dpi_val
        ggplot2::ggsave(file, plot = pcoa_plot_reactive(), device = "png",
                        width = width_in, height = height_in, units = "in", dpi = dpi_val)
      }
    )
    
    nmds_ordination_reactive <- eventReactive(input$redraw_nmds, {
      req(ps_plot_obj())
      isolate(phyloseq::ordinate(ps_plot_obj(), method = "NMDS", distance = dist_method()))
    }, ignoreNULL = FALSE)
    
    nmds_plot_reactive <- reactive({
      req(ps_plot_obj(), primary_group_var(), ord <- nmds_ordination_reactive(), input$dot_size)
      ps_data_for_plot <- ps_plot_obj()
      
      if (is.null(plot_shape_var())) {
        p <- phyloseq::plot_ordination(ps_data_for_plot, ord, color = primary_group_var())
      } else {
        p <- phyloseq::plot_ordination(ps_data_for_plot, ord, color = primary_group_var(), shape = plot_shape_var())
      }
      if (isTRUE(input$show_dot_outline)) {
        outline_size <- input$dot_size + 0.4
        p <- p + ggplot2::geom_point(
          color = "black",
          size = outline_size,
          show.legend = FALSE
        )
      }
      p <- p + ggplot2::geom_point(size = input$dot_size)
      p <- p + ggplot2::scale_color_manual(values = primary_color_map())
      
      if (input$show_ellipses) {
        p <- p + ggplot2::stat_ellipse(ggplot2::aes_string(group = primary_group_var()))
      }
      
      p <- p + ggplot2::theme_bw() +
        ggplot2::labs(title = paste0("NMDS - ", distance_label(), " (Stress: ", round(ord$stress, 3), ")"),
                      color = primary_group_var(),
                      shape = if (is.null(plot_shape_var())) NULL else plot_shape_var())
      
      plot_data <- phyloseq::plot_ordination(ps_data_for_plot, ord, justDF = TRUE)
      
      if (input$show_sample_names) {
        p <- p + ggplot2::geom_text(data = plot_data, mapping = ggplot2::aes(x = NMDS1, y = NMDS2, label = SampleID), color = "black", size = 3, hjust = -0.1, vjust = 0.5)
      }

      if (isTRUE(input$show_cluster_labels) && !is.null(cluster_result_val()$result)) {
        clu_df <- tryCatch(cluster_lookup(), error = function(e) NULL)
        if (!is.null(clu_df)) {
        plot_data$SampleID <- as.character(plot_data$SampleID)
        plot_data <- dplyr::left_join(plot_data, clu_df, by = "SampleID")
        if ("Cluster" %in% colnames(plot_data)) {
          plot_data$ClusterLabel <- ifelse(is.na(plot_data$Cluster), "", paste0("C", plot_data$Cluster))
          p <- p + ggplot2::geom_text(
            data = plot_data,
            mapping = ggplot2::aes(x = NMDS1, y = NMDS2, label = ClusterLabel),
            color = "black", size = 2.7, vjust = -0.8, show.legend = FALSE
          )
        }
        }
      }

      if (isTRUE(input$show_envfit_vectors) && isTRUE(input$run_envfit > 0)) {
        vec_obj <- envfit_result_val()
        vec_df <- if (is.null(vec_obj)) NULL else vec_obj$nmds
        if (!is.null(vec_df) && nrow(vec_df) > 0) {
          x_rng <- range(plot_data$NMDS1, na.rm = TRUE)
          y_rng <- range(plot_data$NMDS2, na.rm = TRUE)
          x_pad <- ifelse(is.finite(diff(x_rng)) && diff(x_rng) > 0, diff(x_rng) * 0.02, 0)
          y_pad <- ifelse(is.finite(diff(y_rng)) && diff(y_rng) > 0, diff(y_rng) * 0.02, 0)
          x_min <- x_rng[1] + x_pad
          x_max <- x_rng[2] - x_pad
          y_min <- y_rng[1] + y_pad
          y_max <- y_rng[2] - y_pad

          taxa_flag <- if ("is_taxa" %in% colnames(vec_df)) as.logical(vec_df$is_taxa) else rep(FALSE, nrow(vec_df))
          taxa_flag[is.na(taxa_flag)] <- FALSE
          vec_arrow <- vec_df[!taxa_flag, , drop = FALSE]
          vec_point <- vec_df[taxa_flag, , drop = FALSE]

          if (nrow(vec_arrow) > 0) {
          vec_arrow$xend_plot <- pmin(pmax(vec_arrow$xend, x_min), x_max)
          vec_arrow$yend_plot <- pmin(pmax(vec_arrow$yend, y_min), y_max)
          p <- p +
            ggplot2::geom_segment(
              data = vec_arrow,
              ggplot2::aes(x = x, y = y, xend = xend_plot, yend = yend_plot),
              inherit.aes = FALSE,
              arrow = grid::arrow(length = grid::unit(0.18, "cm")),
              color = "#2E7D32",
              linewidth = 0.6
            )
          if (requireNamespace("ggrepel", quietly = TRUE)) {
            p <- p + ggrepel::geom_text_repel(
              data = vec_arrow,
              ggplot2::aes(x = xend_plot, y = yend_plot, label = label),
              inherit.aes = FALSE,
              color = "#1B5E20",
              size = 3,
              min.segment.length = 0,
              box.padding = 0.2,
              point.padding = 0.15,
              max.overlaps = Inf
            )
          } else {
            p <- p + ggplot2::geom_text(
              data = vec_arrow,
              ggplot2::aes(x = xend_plot, y = yend_plot, label = label),
              inherit.aes = FALSE,
              color = "#1B5E20",
              size = 3,
              hjust = -0.05,
              check_overlap = TRUE
            )
          }
          }

          if (nrow(vec_point) > 0) {
            vec_point$x_plot <- pmin(pmax(vec_point$xend, x_min), x_max)
            vec_point$y_plot <- pmin(pmax(vec_point$yend, y_min), y_max)

            p <- p +
              ggplot2::geom_point(
                data = vec_point,
                ggplot2::aes(x = x_plot, y = y_plot),
                inherit.aes = FALSE,
                shape = 21,
                fill = "#7B1FA2",
                color = "white",
                stroke = 0.4,
                size = 2.8
              )
            if (requireNamespace("ggrepel", quietly = TRUE)) {
              p <- p + ggrepel::geom_text_repel(
                data = vec_point,
                ggplot2::aes(x = x_plot, y = y_plot, label = label),
                inherit.aes = FALSE,
                color = "#4A148C",
                size = 2.8,
                min.segment.length = 0,
                box.padding = 0.2,
                point.padding = 0.15,
                max.overlaps = Inf
              )
            } else {
              p <- p + ggplot2::geom_text(
                data = vec_point,
                ggplot2::aes(x = x_plot, y = y_plot, label = label),
                inherit.aes = FALSE,
                color = "#4A148C",
                size = 2.8,
                hjust = -0.05,
                vjust = -0.2,
                check_overlap = TRUE
              )
            }
          }
        }
      }
      
      
      limits <- axis_limits()
      if (!is.null(limits$x) || !is.null(limits$y)) {
        p <- p + ggplot2::coord_cartesian(xlim = limits$x, ylim = limits$y)
      }
      
      return(p)
    })
    
    output$nmds_plot_out <- renderPlot({ nmds_plot_reactive() })
    
    output$download_nmds <- downloadHandler(
      filename = function() { paste0("NMDS_", gsub("-", "", distance_label()), "_", Sys.Date(), ".png") },
      content = function(file) {
        dpi_val <- 300
        width_in <- input$plot_width_px / dpi_val
        height_in <- input$plot_height_px / dpi_val
        ggplot2::ggsave(file, plot = nmds_plot_reactive(), device = "png",
                        width = width_in, height = height_in, units = "in", dpi = dpi_val)
      }
    )
    
    permanova_reactive <- reactive({
      req(data_prepared(), primary_group_var())
      data <- data_prepared()
      
      if (!requireNamespace("vegan", quietly = TRUE)) {
        stop("The 'vegan' package is required for PERMANOVA. Please install it.")
      }
      
      permanova_formula <- as.formula(paste("data$dist_mat ~", primary_group_var()))
      
      permanova_result <- vegan::adonis2(permanova_formula,
                                         data = data$metadata,
                                         permutations = 999)
      
      R2 <- permanova_result$R2[[1]]
      Pval <- permanova_result$`Pr(>F)`[[1]]
      
      list(R2 = R2, Pval = Pval)
    })
    
    output$permanova_results_out <- renderUI({
      req(permanova_result <- permanova_reactive(), primary_group_var())
      HTML(paste0(
        "<strong>Grouping Variable:</strong> ", primary_group_var(), "<br>",
        "<strong>Distance Metric:</strong> ", distance_label(), "<br>",
        "---<br>",
        "<strong>R-squared:</strong> ", format(permanova_result$R2, digits = 3, nsmall = 3), "<br>",
        "<strong>P-value:</strong> ", format(permanova_result$Pval, digits = 3)
      ))
    })

    output$clustering_results_out <- renderText({
      cluster_obj <- cluster_result_val()
      res <- cluster_obj$result
      if (is.null(res)) {
        return(cluster_obj$note)
      }
      cluster_counts <- table(res$cluster)
      cluster_lines <- paste0("Cluster ", names(cluster_counts), ": ", as.integer(cluster_counts), " sample(s)")

      mode_line <- if (identical(res$mode, "auto")) {
        paste0("Mode: Auto k (Silhouette), searched k=", res$k_min, " to ", res$k_max)
      } else {
        "Mode: Manual k"
      }

      paste(
        c(
          mode_line,
          paste0("Selected k: ", res$k),
          paste0("Average silhouette width: ", format(round(res$silhouette, 4), nsmall = 4)),
          "Cluster sizes:",
          cluster_lines
        ),
        collapse = "\n"
      )
    })

    output$envfit_results_out <- renderText({
      if (!isTRUE(input$run_envfit > 0)) {
        return("EnvFit has not been run yet.\nOpen Advanced Options and click 'Run EnvFit'.")
      }
      env_obj <- envfit_result_val()
      if (is.null(env_obj)) {
        return("EnvFit is running or did not return results yet. Please click 'Run EnvFit' again.")
      }
      pcoa_n <- if (is.null(env_obj$pcoa)) 0 else nrow(env_obj$pcoa)
      nmds_n <- if (is.null(env_obj$nmds)) 0 else nrow(env_obj$nmds)
      vars_selected <- input$envfit_vars
      if (is.null(vars_selected)) vars_selected <- character(0)
      taxa_selected <- input$envfit_taxa_vars
      if (is.null(taxa_selected)) taxa_selected <- character(0)
      env_note <- env_obj$note

      pcoa_taxa <- if (!is.null(env_obj$pcoa) && nrow(env_obj$pcoa) > 0) {
        env_obj$pcoa[startsWith(as.character(env_obj$pcoa$variable), "Taxa::"), c("variable", "r2", "p_value"), drop = FALSE]
      } else data.frame()
      nmds_taxa <- if (!is.null(env_obj$nmds) && nrow(env_obj$nmds) > 0) {
        env_obj$nmds[startsWith(as.character(env_obj$nmds$variable), "Taxa::"), c("variable", "r2", "p_value"), drop = FALSE]
      } else data.frame()

      fmt_taxa_lines <- function(df, title) {
        if (is.null(df) || nrow(df) == 0) {
          return(c(paste0(title, ": None")))
        }
        df <- df[order(df$r2, decreasing = TRUE), , drop = FALSE]
        top_df <- head(df, 10)
        lines <- paste0(
          "- ",
          gsub("^Taxa::", "", as.character(top_df$variable)),
          " | r2=",
          format(round(as.numeric(top_df$r2), 4), nsmall = 4),
          " | p=",
          format(round(as.numeric(top_df$p_value), 4), nsmall = 4)
        )
        c(paste0(title, " (top ", nrow(top_df), " by r2):"), lines)
      }

      taxa_effect_lines <- c(
        fmt_taxa_lines(pcoa_taxa, "Taxa effect size in PCoA"),
        fmt_taxa_lines(nmds_taxa, "Taxa effect size in NMDS")
      )

      paste(
        c(
          paste0("Selected metadata variables: ", if (length(vars_selected) == 0) "None" else paste(vars_selected, collapse = ", ")),
          paste0("Selected taxa variables: ", if (length(taxa_selected) == 0) "None" else paste(taxa_selected, collapse = ", ")),
          if (!is.null(env_note) && nzchar(env_note)) paste0("Note: ", env_note) else NULL,
          paste0("Only significant vectors: ", ifelse(isTRUE(input$envfit_only_sig), "Yes", "No")),
          paste0("p-value cutoff: ", input$envfit_p_cutoff),
          paste0("PCoA vectors shown: ", pcoa_n),
          paste0("NMDS vectors shown: ", nmds_n),
          "",
          taxa_effect_lines
        ),
        collapse = "\n"
      )
    })
    
  })
}
