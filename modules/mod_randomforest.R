## UI
mod_randomforest_ui <- function(id) {
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
        h4(icon("tree"), "Random Forest"),
        hr(),
        selectInput(ns("group_var"), "1. Primary grouping variable", choices = NULL),
        selectInput(ns("primary_level"), "2. Primary level to include", choices = NULL),
        selectInput(ns("outcome_var"), "3. Secondary grouping variable", choices = NULL),
        selectInput(
          ns("outcome_type"),
          "4. Outcome type",
          choices = c("Auto", "Classification", "Regression"),
          selected = "Classification"
        ),
        selectizeInput(
          ns("outcome_levels"),
          "5. Levels to include (optional)",
          choices = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Select levels to include",
            plugins = list("remove_button")
          )
        ),
        selectInput(ns("shap_target_class"), "6. SHAP target class", choices = NULL),
        selectInput(
          ns("tax_level"),
          "7. Taxonomic level",
          choices = c("ASV", "Genus", "Species"),
          selected = "Genus"
        ),
        selectInput(
          ns("transform_method"),
          "8. Feature transform",
          choices = c("TSS", "CLR", "Presence/Absence", "log"),
          selected = "TSS"
        ),
        numericInput(ns("prevalence_filter_pct"), "9. Feature prevalence cutoff (0-20%)", value = 5, min = 0, max = 20, step = 1),
        sliderInput(ns("train_ratio"), "10. Train ratio", min = 0.6, max = 0.9, value = 0.8, step = 0.05),
        tags$details(
          style = "margin-bottom: 10px;",
          tags$summary("Advanced Options"),
          br(),
          numericInput(ns("top_n_features"), "11. Top N features by mean abundance", value = 100, min = 10, max = 5000, step = 10),
          numericInput(ns("ntree"), "12. Number of trees (ntree)", value = 500, min = 100, max = 5000, step = 100),
          numericInput(ns("mtry"), "13. mtry (0 = auto)", value = 0, min = 0, max = 10000, step = 1),
          selectInput(
            ns("validation_mode"),
            "14. Validation mode",
            choices = c("Holdout split", "K-fold CV"),
            selected = "Holdout split"
          ),
          numericInput(ns("cv_folds"), "15. K for K-fold CV", value = 5, min = 3, max = 10, step = 1),
          numericInput(ns("seed"), "16. Random seed", value = 1234, min = 1, max = 999999, step = 1)
        ),
        hr(),
        h4(icon("up-right-and-down-left-from-center"), "Plot Dimensions"),
        numericInput(ns("plot_width"), "Plot width (px)", value = 1160, min = 600, max = 3200, step = 50),
        numericInput(ns("plot_height"), "Plot height (px)", value = 500, min = 300, max = 2400, step = 50),
        numericInput(ns("base_size"), "Base Font Size:", value = 11, min = 6, max = 30, step = 1),
        tags$div(
          style = "display: flex; align-items: center; gap: 8px; flex-wrap: wrap;",
          actionButton(ns("run_rf"), "Run Random Forest", class = "btn-danger", style = "font-size: 12px;"),
          tags$span("May take a long time.", style = "font-size: 11px; color: #b94a48;")
        ),
        tags$script(HTML(
          "Shiny.addCustomMessageHandler('toggle-rf-run-btn', function(msg) {
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
        width = 10,
        h4("Random Forest"),
        tags$div(
          id = ns("rf_tab_container"),
          style = "max-width: 100%;",
          tabsetPanel(
            id = ns("rf_active_tab"),
            tabPanel(
              "RF Table",
              downloadButton(
                ns("download_rf_table"),
                "Download Table (TSV)",
                style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"
              ),
              DTOutput(ns("importance_table"))
            ),
            tabPanel(
              "RF Bar plot",
              downloadButton(
                ns("download_rf_barplot"),
                "Download Plot (PNG)",
                style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"
              ),
              div(
                style = "display: flex; justify-content: center; margin-top: 8px;",
                uiOutput(ns("importance_plot_ui"))
              )
            ),
            tabPanel(
              "ROC / AUC",
              downloadButton(
                ns("download_rf_rocplot"),
                "Download ROC (PNG)",
                style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"
              ),
              div(
                style = "display: flex; justify-content: center; margin-top: 8px;",
                plotOutput(ns("roc_plot"), width = "1160px", height = "500px")
              ),
              verbatimTextOutput(ns("roc_auc_text"))
            ),
            tabPanel(
              "SHAP Table",
              downloadButton(
                ns("download_rf_shap_summary"),
                "Download SHAP Summary (TSV)",
                style = "width: 220px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"
              ),
              DTOutput(ns("shap_summary_table"))
            ),
            tabPanel(
              "SHAP Bar plot",
              downloadButton(
                ns("download_rf_shap_plot"),
                "Download SHAP Plot (PNG)",
                style = "width: 210px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"
              ),
              div(
                style = "display: flex; justify-content: center; margin-top: 8px;",
                uiOutput(ns("shap_summary_plot_ui"))
              )
            )
          )
        ),
        uiOutput(ns("rf_legend_box")),
        uiOutput(ns("rf_results_separator")),
        h4(icon("square-poll-vertical"), "Result"),
        uiOutput(ns("rf_metrics_box"))
      )
    )
  )
}

## Server
mod_randomforest_server <- function(id, ps_obj_filtered_raw) {
  moduleServer(id, function(input, output, session) {
    observe({
      width_px <- suppressWarnings(as.integer(input$plot_width))
      if (!is.finite(width_px) || is.na(width_px) || width_px <= 0) width_px <- 1160L
      session$sendCustomMessage(
        "set-tab-container-width",
        list(id = session$ns("rf_tab_container"), width = paste0(width_px, "px"))
      )
    })

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

    model_result <- reactiveVal(NULL)
    status_text <- reactiveVal("Waiting for model run.")
    rf_running <- reactiveVal(FALSE)
    filtered_outcome_vars <- reactiveVal(character(0))
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
    outcome_var_resolved <- reactive({
      req(ps_obj_filtered_raw(), input$outcome_var)
      if (identical(input$outcome_var, "None")) {
        return(group_var_resolved())
      }
      md <- data.frame(phyloseq::sample_data(ps_obj_filtered_raw()))
      resolve_meta_colname(input$outcome_var, colnames(md))
    })
    group_var_resolved <- reactive({
      req(ps_obj_filtered_raw(), input$group_var)
      md <- data.frame(phyloseq::sample_data(ps_obj_filtered_raw()))
      resolve_meta_colname(input$group_var, colnames(md))
    })
    log_rf <- function(msg) {
      cat(
        sprintf(
          "[%s] [RF][%s] %s\n",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          session$ns("run_rf"),
          msg
        ),
        file = stderr()
      )
    }

    is_identifier_like <- function(vec) {
      v <- as.character(vec)
      v <- v[!is.na(v) & nzchar(v)]
      if (length(v) == 0) return(TRUE)
      uniq_ratio <- length(unique(v)) / length(v)
      uniq_ratio > 0.9
    }

    get_valid_outcome_vars <- function(md, outcome_mode) {
      if (is.null(md) || ncol(md) == 0) return(character(0))
      vars <- colnames(md)
      vars <- vars[!tolower(vars) %in% c("sampleid", "sample_id")]
      vars <- vars[!grepl("(^|_)id$", tolower(vars))]

      keep <- vapply(vars, function(vn) {
        x <- md[[vn]]
        x_non_na <- x[!is.na(x)]
        if (length(x_non_na) < 3) return(FALSE)
        if (length(unique(x_non_na)) < 2) return(FALSE)
        if (is_identifier_like(x_non_na)) return(FALSE)

        if (identical(outcome_mode, "Classification")) {
          lv <- table(as.character(x_non_na))
          return(length(lv) >= 2 && all(lv >= 3))
        }
        if (identical(outcome_mode, "Auto")) {
          if (is.factor(x) || is.character(x)) {
            lv <- table(as.character(x_non_na))
            return(length(lv) >= 2 && all(lv >= 3))
          }
        }
        TRUE
      }, logical(1))

      vars[keep]
    }

    observe({
      req(ps_obj_filtered_raw(), input$outcome_type)
      ps <- ps_obj_filtered_raw()
      md <- data.frame(phyloseq::sample_data(ps))
      vars <- get_valid_outcome_vars(md, input$outcome_type)
      filtered_outcome_vars(vars)
      group_choices <- setdiff(colnames(md), c("SampleID", "sample_id"))
      current_group <- input$group_var
      if (is.null(current_group) || !current_group %in% group_choices) {
        current_group <- if (length(group_choices) > 0) group_choices[1] else NULL
      }
      updateSelectInput(session, "group_var", choices = group_choices, selected = current_group)

      primary_level_choices <- if (!is.null(current_group) && current_group %in% colnames(md)) {
        v <- sort(unique(as.character(md[[current_group]])))
        v[!is.na(v) & nzchar(v)]
      } else {
        character(0)
      }
      current_primary_level <- input$primary_level
      if (is.null(current_primary_level) || !current_primary_level %in% primary_level_choices) {
        current_primary_level <- if (length(primary_level_choices) > 0) primary_level_choices[1] else NULL
      }
      updateSelectInput(session, "primary_level", choices = primary_level_choices, selected = current_primary_level)

      current <- input$outcome_var
      choices_outcome <- c("None", vars)
      if (is.null(current) || !current %in% choices_outcome) {
        current <- "None"
      }
      updateSelectInput(session, "outcome_var", choices = choices_outcome, selected = current)
    })

    observeEvent(list(ps_obj_filtered_raw(), input$outcome_var, input$outcome_type, input$group_var), {
      req(ps_obj_filtered_raw(), input$outcome_var, input$outcome_type, input$group_var)
      ps <- ps_obj_filtered_raw()
      md <- data.frame(phyloseq::sample_data(ps))
      if (!identical(input$outcome_var, "None")) {
        validate(
          need(input$outcome_var %in% filtered_outcome_vars(), "Selected outcome variable is filtered out by automatic validity checks.")
        )
      }
      outcome_var <- outcome_var_resolved()
      validate(need(outcome_var %in% colnames(md), "Selected outcome is not available."))

      v <- md[[outcome_var]]
      mode_selected <- input$outcome_type
      is_classification_mode <- if (identical(mode_selected, "Auto")) {
        is.factor(v) || is.character(v)
      } else {
        identical(mode_selected, "Classification")
      }
      if (is_classification_mode) {
        choices <- sort(unique(as.character(v)))
        choices <- choices[!is.na(choices) & nzchar(choices)]
      } else {
        choices <- character(0)
      }
      selected <- isolate(input$outcome_levels)
      if (is.null(selected)) selected <- character(0)
      selected <- selected[selected %in% choices]
      updateSelectizeInput(session, "outcome_levels", choices = choices, selected = selected, server = TRUE)
      shap_selected <- isolate(input$shap_target_class)
      if (is.null(shap_selected) || !shap_selected %in% choices) {
        shap_selected <- if (length(choices) >= 2) choices[2] else if (length(choices) == 1) choices[1] else character(0)
      }
      updateSelectInput(session, "shap_target_class", choices = choices, selected = shap_selected)
    }, ignoreInit = FALSE)

    observeEvent(input$run_rf, {
      req(ps_obj_filtered_raw(), input$outcome_var)

      rf_running(TRUE)
      status_text("Running Random Forest...")
      model_result(NULL)
      session$sendCustomMessage("toggle-rf-run-btn", list(
        id = session$ns("run_rf"),
        disabled = TRUE,
        label = "Running..."
      ))
      on.exit({
        rf_running(FALSE)
        session$sendCustomMessage("toggle-rf-run-btn", list(
          id = session$ns("run_rf"),
          disabled = FALSE,
          label = "Run Random Forest"
        ))
      }, add = TRUE)
      log_rf("run started")
      tryCatch({
        withCallingHandlers({
          ps <- ps_obj_filtered_raw()
          md <- data.frame(phyloseq::sample_data(ps))
          outcome_var <- outcome_var_resolved()
          validate(need(outcome_var %in% colnames(md), "Outcome variable is not available."))
          primary_var <- group_var_resolved()
          validate(need(primary_var %in% colnames(md), "Primary grouping variable is not available."))

          y_raw <- md[[outcome_var]]
          taxonomy_info <- NULL
          feature_map <- NULL

          if (!identical(input$tax_level, "ASV")) {
            validate(
              need(!is.null(phyloseq::tax_table(ps)), "Taxonomy table is required for taxonomic aggregation.")
            )
            tax_cols <- colnames(phyloseq::tax_table(ps))
            validate(
              need(input$tax_level %in% tax_cols, paste0("Selected taxonomic rank '", input$tax_level, "' is not available."))
            )
            ps <- apply_disambiguated_taxrank(ps, input$tax_level)
            ps <- phyloseq::tax_glom(ps, taxrank = input$tax_level, NArm = FALSE)
          }

          otu <- as(phyloseq::otu_table(ps), "matrix")
          if (phyloseq::taxa_are_rows(ps)) otu <- t(otu)
          x <- as.data.frame(otu, check.names = FALSE)
          feature_ids <- colnames(x)
          feature_model_names <- make.names(feature_ids, unique = TRUE)
          colnames(x) <- feature_model_names

          if (!is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
            tax_df <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE)
            tax_df$FeatureID <- rownames(tax_df)
            taxonomy_info <- tax_df
          }
          feature_map <- data.frame(
            FeatureModel = feature_model_names,
            FeatureID = feature_ids,
            stringsAsFactors = FALSE
          )

          common_samples <- intersect(rownames(md), rownames(x))
          validate(need(length(common_samples) >= 10, "At least 10 samples are required."))

          md <- md[common_samples, , drop = FALSE]
          x <- x[common_samples, , drop = FALSE]
          if (!identical(input$outcome_var, "None")) {
            keep_primary <- as.character(md[[primary_var]]) == as.character(input$primary_level)
            keep_primary[is.na(keep_primary)] <- FALSE
            md <- md[keep_primary, , drop = FALSE]
            x <- x[keep_primary, , drop = FALSE]
          }
          y_raw <- md[[outcome_var]]

          keep_non_na <- !is.na(y_raw)
          y_raw <- y_raw[keep_non_na]
          x <- x[keep_non_na, , drop = FALSE]

          if ((is.factor(y_raw) || is.character(y_raw)) && length(input$outcome_levels) > 0) {
            level_chr <- as.character(y_raw)
            keep_levels <- level_chr %in% input$outcome_levels
            x <- x[keep_levels, , drop = FALSE]
            y_raw <- y_raw[keep_levels]
          }

          validate(need(nrow(x) >= 10, "Not enough samples after filtering."))

          x_counts <- x
          prevalence_cutoff <- as.numeric(input$prevalence_filter_pct) / 100
          prevalence <- colMeans(x_counts > 0, na.rm = TRUE)
          keep_features <- prevalence >= prevalence_cutoff
          x_counts <- x_counts[, keep_features, drop = FALSE]
          validate(need(ncol(x_counts) >= 2, "Not enough features after prevalence filtering."))

          top_n <- min(as.integer(input$top_n_features), ncol(x_counts))
          mean_abund <- colMeans(x_counts, na.rm = TRUE)
          top_features <- names(sort(mean_abund, decreasing = TRUE))[seq_len(top_n)]
          x_counts <- x_counts[, top_features, drop = FALSE]

          transform_method <- input$transform_method
          x <- switch(
            transform_method,
            "TSS" = {
              rs <- rowSums(x_counts, na.rm = TRUE)
              rs[rs <= 0 | is.na(rs)] <- 1
              sweep(x_counts, 1, rs, "/")
            },
            "CLR" = {
              pseudo <- 1
              x_pos <- x_counts + pseudo
              log_x <- log(x_pos)
              gm <- exp(rowMeans(log_x, na.rm = TRUE))
              gm[gm <= 0 | is.na(gm)] <- 1
              log_x - log(gm)
            },
            "Presence/Absence" = {
              (x_counts > 0) * 1
            },
            "log" = {
              log1p(x_counts)
            },
            x_counts
          )
          x <- as.data.frame(x, check.names = FALSE)
          log_rf(sprintf("matrix prepared: n=%d, p=%d, size=%s", nrow(x), ncol(x), format(utils::object.size(x), units = "MB")))

          set.seed(as.integer(input$seed))

          mode_selected <- input$outcome_type
          is_classification <- if (identical(mode_selected, "Auto")) {
            is.factor(y_raw) || is.character(y_raw)
          } else {
            identical(mode_selected, "Classification")
          }
          compute_binary_roc <- function(y_true_bin, y_score) {
            keep <- !is.na(y_true_bin) & !is.na(y_score)
            y_true_bin <- y_true_bin[keep]
            y_score <- y_score[keep]
            if (length(y_true_bin) < 2) return(NULL)
            n_pos <- sum(y_true_bin == 1)
            n_neg <- sum(y_true_bin == 0)
            if (n_pos == 0 || n_neg == 0) return(NULL)

            ord <- order(y_score, decreasing = TRUE)
            y_sorted <- y_true_bin[ord]
            score_sorted <- y_score[ord]
            tp <- cumsum(y_sorted == 1)
            fp <- cumsum(y_sorted == 0)
            tpr <- tp / n_pos
            fpr <- fp / n_neg
            change_idx <- which(c(diff(score_sorted) != 0, TRUE))
            tpr_u <- c(0, tpr[change_idx], 1)
            fpr_u <- c(0, fpr[change_idx], 1)
            auc <- sum(diff(fpr_u) * (head(tpr_u, -1) + tail(tpr_u, -1)) / 2)
            list(curve = data.frame(FPR = fpr_u, TPR = tpr_u), auc = as.numeric(auc))
          }

          if (is_classification) {
            y <- factor(as.character(y_raw))
            validate(need(length(levels(y)) >= 2, "Classification requires at least two classes."))
          } else {
            y <- as.numeric(y_raw)
            validate(need(!all(is.na(y)), "Outcome has only missing values."))
          }

          validation_mode <- if (!is.null(input$validation_mode)) input$validation_mode else "Holdout split"
          cv_folds <- as.integer(input$cv_folds)
          if (is.na(cv_folds) || cv_folds < 3) cv_folds <- 5L
          cv_folds <- min(cv_folds, 10L)
          cv_folds <- min(cv_folds, nrow(x))
          if (cv_folds < 3) validation_mode <- "Holdout split"

          get_mtry <- function(p) {
            mtry_val <- as.integer(input$mtry)
            if (is.na(mtry_val) || mtry_val <= 0) mtry_val <- floor(sqrt(p))
            max(1L, min(mtry_val, p))
          }

          roc_payload <- NULL
          if (identical(validation_mode, "K-fold CV")) {
            create_folds_classification <- function(y_fac, k) {
              folds <- vector("list", k)
              by_class <- split(seq_along(y_fac), y_fac)
              for (idx in by_class) {
                idx <- sample(idx, length(idx))
                for (j in seq_along(idx)) {
                  fold_id <- ((j - 1) %% k) + 1
                  folds[[fold_id]] <- c(folds[[fold_id]], idx[j])
                }
              }
              lapply(folds, sort)
            }
            create_folds_regression <- function(n, k) {
              idx <- sample(seq_len(n), n)
              split(idx, cut(seq_along(idx), breaks = k, labels = FALSE))
            }

            folds <- if (is_classification) create_folds_classification(y, cv_folds) else create_folds_regression(nrow(x), cv_folds)
            if (is_classification) {
              lv <- levels(y)
              oof_pred <- rep(NA_character_, length(y))
              oof_prob <- matrix(NA_real_, nrow = length(y), ncol = length(lv), dimnames = list(NULL, lv))
              acc_vec <- numeric(length(folds))
              macro_f1_vec <- numeric(length(folds))
              for (i in seq_along(folds)) {
                valid_idx <- as.integer(folds[[i]])
                train_idx <- setdiff(seq_len(nrow(x)), valid_idx)
                x_tr <- x[train_idx, , drop = FALSE]
                y_tr <- y[train_idx]
                x_va <- x[valid_idx, , drop = FALSE]
                y_va <- y[valid_idx]
                mtry_fold <- get_mtry(ncol(x_tr))
                rf_fold <- randomForest::randomForest(
                  x = x_tr, y = y_tr, ntree = as.integer(input$ntree), mtry = mtry_fold, importance = FALSE
                )
                pred_fold <- stats::predict(rf_fold, newdata = x_va)
                pred_prob_fold <- as.data.frame(stats::predict(rf_fold, newdata = x_va, type = "prob"), check.names = FALSE)
                oof_pred[valid_idx] <- as.character(pred_fold)
                for (cl in colnames(pred_prob_fold)) {
                  if (cl %in% colnames(oof_prob)) oof_prob[valid_idx, cl] <- as.numeric(pred_prob_fold[[cl]])
                }
                acc_vec[i] <- mean(pred_fold == y_va)
                f1_each <- vapply(lv, function(cl) {
                  tp <- sum(pred_fold == cl & y_va == cl)
                  fp <- sum(pred_fold == cl & y_va != cl)
                  fn <- sum(pred_fold != cl & y_va == cl)
                  precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
                  recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
                  if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
                }, numeric(1))
                macro_f1_vec[i] <- mean(f1_each)
              }

              oof_pred_fac <- factor(oof_pred, levels = levels(y))
              cm <- table(Predicted = oof_pred_fac, Actual = y)
              roc_list <- lapply(levels(y), function(cl) {
                y_bin <- as.integer(y == cl)
                score <- as.numeric(oof_prob[, cl])
                roc_res <- compute_binary_roc(y_bin, score)
                if (is.null(roc_res)) return(NULL)
                roc_res$class <- cl
                roc_res
              })
              names(roc_list) <- levels(y)
              roc_list <- roc_list[!vapply(roc_list, is.null, logical(1))]
              if (length(roc_list) > 0) {
                auc_values <- vapply(roc_list, function(xx) xx$auc, numeric(1))
                roc_payload <- list(curves = roc_list, auc_by_class = auc_values, macro_auc = mean(auc_values, na.rm = TRUE))
              }

              metrics_text <- paste0(
                "Task: Classification\n",
                "Outcome type mode: ", mode_selected, "\n",
                "Taxonomic level: ", input$tax_level, "\n",
                "Transform: ", transform_method, "\n",
                "Validation: K-fold CV (K=", cv_folds, ")\n",
                "Samples: ", nrow(x), "\n",
                "Features used: ", ncol(x), "\n",
                "Accuracy (mean±sd): ", sprintf("%.4f ± %.4f", mean(acc_vec), stats::sd(acc_vec)), "\n",
                "Macro-F1 (mean±sd): ", sprintf("%.4f ± %.4f", mean(macro_f1_vec), stats::sd(macro_f1_vec)), "\n\n",
                "OOF Confusion Matrix:\n",
                paste(capture.output(print(cm)), collapse = "\n")
              )
            } else {
              rmse_vec <- numeric(length(folds))
              mae_vec <- numeric(length(folds))
              r2_vec <- numeric(length(folds))
              for (i in seq_along(folds)) {
                valid_idx <- as.integer(folds[[i]])
                train_idx <- setdiff(seq_len(nrow(x)), valid_idx)
                x_tr <- x[train_idx, , drop = FALSE]
                y_tr <- y[train_idx]
                x_va <- x[valid_idx, , drop = FALSE]
                y_va <- y[valid_idx]
                mtry_fold <- get_mtry(ncol(x_tr))
                rf_fold <- randomForest::randomForest(
                  x = x_tr, y = y_tr, ntree = as.integer(input$ntree), mtry = mtry_fold, importance = FALSE
                )
                pred_fold <- as.numeric(stats::predict(rf_fold, newdata = x_va))
                rmse_vec[i] <- sqrt(mean((pred_fold - y_va)^2))
                mae_vec[i] <- mean(abs(pred_fold - y_va))
                sst <- sum((y_va - mean(y_va))^2)
                sse <- sum((y_va - pred_fold)^2)
                r2_vec[i] <- if (sst == 0) NA_real_ else 1 - sse / sst
              }
              metrics_text <- paste0(
                "Task: Regression\n",
                "Outcome type mode: ", mode_selected, "\n",
                "Taxonomic level: ", input$tax_level, "\n",
                "Transform: ", transform_method, "\n",
                "Validation: K-fold CV (K=", cv_folds, ")\n",
                "Samples: ", nrow(x), "\n",
                "Features used: ", ncol(x), "\n",
                "RMSE (mean±sd): ", sprintf("%.4f ± %.4f", mean(rmse_vec, na.rm = TRUE), stats::sd(rmse_vec, na.rm = TRUE)), "\n",
                "MAE (mean±sd): ", sprintf("%.4f ± %.4f", mean(mae_vec, na.rm = TRUE), stats::sd(mae_vec, na.rm = TRUE)), "\n",
                "R-squared (mean±sd): ", sprintf("%.4f ± %.4f", mean(r2_vec, na.rm = TRUE), stats::sd(r2_vec, na.rm = TRUE))
              )
            }

            mtry_val <- get_mtry(ncol(x))
            rf_options_text <- paste0(
              "RF options:",
              "\n- package: randomForest",
              "\n- validation mode: K-fold CV",
              "\n- K folds: ", as.integer(cv_folds),
              "\n- seed: ", as.integer(input$seed),
              "\n- prevalence cutoff (%): ", as.numeric(input$prevalence_filter_pct),
              "\n- top_n_features before transform: ", as.integer(input$top_n_features),
              "\n- ntree: ", as.integer(input$ntree),
              "\n- mtry used: ", as.integer(mtry_val),
              "\n- selected outcome levels: ", if (length(input$outcome_levels) > 0) paste(input$outcome_levels, collapse = ", ") else "All"
            )
            log_rf(sprintf("fit final model after CV: ntree=%d, mtry=%d", as.integer(input$ntree), mtry_val))
            rf_fit <- randomForest::randomForest(
              x = x, y = y, ntree = as.integer(input$ntree), mtry = mtry_val, importance = TRUE
            )
            log_rf("final fit completed")

            x_train <- x
            y_train <- y
            explain_n <- min(nrow(x), max(20L, floor(nrow(x) * 0.2)))
            explain_idx <- sort(sample(seq_len(nrow(x)), explain_n))
            x_test <- x[explain_idx, , drop = FALSE]
            y_test <- y[explain_idx]
          } else {
            if (is_classification) {
              split_by_class <- split(seq_len(length(y)), y)
              train_idx <- unlist(lapply(split_by_class, function(idx) {
                n_tr <- max(1, floor(length(idx) * input$train_ratio))
                sample(idx, size = n_tr)
              }), use.names = FALSE)
              train_idx <- sort(unique(train_idx))
              test_idx <- setdiff(seq_len(nrow(x)), train_idx)
              validate(need(length(test_idx) > 0, "Test set is empty. Lower train ratio or add samples."))
            } else {
              train_n <- floor(nrow(x) * input$train_ratio)
              train_idx <- sample(seq_len(nrow(x)), size = train_n)
              test_idx <- setdiff(seq_len(nrow(x)), train_idx)
              validate(need(length(test_idx) > 0, "Test set is empty. Lower train ratio or add samples."))
            }

            x_train <- x[train_idx, , drop = FALSE]
            y_train <- y[train_idx]
            x_test <- x[test_idx, , drop = FALSE]
            y_test <- y[test_idx]

            mtry_val <- get_mtry(ncol(x_train))
            rf_options_text <- paste0(
              "RF options:",
              "\n- package: randomForest",
              "\n- split strategy: ", if (is_classification) "stratified by class" else "random split",
              "\n- validation mode: Holdout split",
              "\n- train ratio: ", sprintf("%.2f", as.numeric(input$train_ratio)),
              "\n- seed: ", as.integer(input$seed),
              "\n- prevalence cutoff (%): ", as.numeric(input$prevalence_filter_pct),
              "\n- top_n_features before transform: ", as.integer(input$top_n_features),
              "\n- ntree: ", as.integer(input$ntree),
              "\n- mtry used: ", as.integer(mtry_val),
              "\n- selected outcome levels: ", if (length(input$outcome_levels) > 0) paste(input$outcome_levels, collapse = ", ") else "All"
            )

            log_rf(sprintf("fit start: ntree=%d, mtry=%d", as.integer(input$ntree), mtry_val))
            rf_fit <- randomForest::randomForest(
              x = x_train,
              y = y_train,
              ntree = as.integer(input$ntree),
              mtry = mtry_val,
              importance = TRUE
            )
            log_rf("fit completed")

            pred <- stats::predict(rf_fit, newdata = x_test)
            if (is_classification) {
              cm <- table(Predicted = pred, Actual = y_test)
              accuracy <- mean(pred == y_test)
              pred_prob <- stats::predict(rf_fit, newdata = x_test, type = "prob")
              lv <- levels(y)
              f1_each <- vapply(lv, function(cl) {
                tp <- sum(pred == cl & y_test == cl)
                fp <- sum(pred == cl & y_test != cl)
                fn <- sum(pred != cl & y_test == cl)
                precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
                recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
                if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
              }, numeric(1))
              macro_f1 <- mean(f1_each)
              roc_list <- lapply(lv, function(cl) {
                y_bin <- as.integer(y_test == cl)
                score <- as.numeric(pred_prob[, cl])
                roc_res <- compute_binary_roc(y_bin, score)
                if (is.null(roc_res)) return(NULL)
                roc_res$class <- cl
                roc_res
              })
              names(roc_list) <- lv
              roc_list <- roc_list[!vapply(roc_list, is.null, logical(1))]
              if (length(roc_list) > 0) {
                auc_values <- vapply(roc_list, function(x1) x1$auc, numeric(1))
                macro_auc <- mean(auc_values, na.rm = TRUE)
                roc_payload <- list(
                  curves = roc_list,
                  auc_by_class = auc_values,
                  macro_auc = macro_auc
                )
              }
              metrics_text <- paste0(
                "Task: Classification\n",
                "Outcome type mode: ", mode_selected, "\n",
                "Taxonomic level: ", input$tax_level, "\n",
                "Transform: ", transform_method, "\n",
                "Validation: Holdout split\n",
                "Samples: ", nrow(x), " (Train: ", length(train_idx), ", Test: ", length(test_idx), ")\n",
                "Features used: ", ncol(x), "\n",
                "Accuracy: ", round(accuracy, 4), "\n",
                "Macro-F1: ", round(macro_f1, 4), "\n\n",
                "Confusion Matrix:\n",
                paste(capture.output(print(cm)), collapse = "\n")
              )
            } else {
              rmse <- sqrt(mean((pred - y_test)^2))
              mae <- mean(abs(pred - y_test))
              sst <- sum((y_test - mean(y_test))^2)
              sse <- sum((y_test - pred)^2)
              r2 <- if (sst == 0) NA_real_ else 1 - sse / sst
              metrics_text <- paste0(
                "Task: Regression\n",
                "Outcome type mode: ", mode_selected, "\n",
                "Taxonomic level: ", input$tax_level, "\n",
                "Transform: ", transform_method, "\n",
                "Validation: Holdout split\n",
                "Samples: ", nrow(x), " (Train: ", length(train_idx), ", Test: ", length(test_idx), ")\n",
                "Features used: ", ncol(x), "\n",
                "RMSE: ", round(rmse, 4), "\n",
                "MAE: ", round(mae, 4), "\n",
                "R-squared: ", round(r2, 4)
              )
            }
          }

          perm_raw <- randomForest::importance(rf_fit, type = 1, scale = TRUE)
          perm_df <- as.data.frame(perm_raw)
          perm_df$Feature <- rownames(perm_df)

          if ("MeanDecreaseAccuracy" %in% colnames(perm_df)) {
            perm_col <- "MeanDecreaseAccuracy"
          } else if ("%IncMSE" %in% colnames(perm_df)) {
            perm_col <- "%IncMSE"
          } else {
            perm_col <- colnames(perm_df)[1]
          }

          imp_tbl <- perm_df %>%
            dplyr::select(Feature, PermutationImportance = all_of(perm_col)) %>%
            dplyr::arrange(dplyr::desc(PermutationImportance))

          imp_tbl <- imp_tbl %>%
            dplyr::left_join(feature_map, by = c("Feature" = "FeatureModel"))

          if (!is.null(taxonomy_info)) {
            taxonomy_cols <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(taxonomy_info))
            if (length(taxonomy_cols) > 0) {
              imp_tbl <- imp_tbl %>%
                dplyr::left_join(
                  taxonomy_info %>% dplyr::select(FeatureID, dplyr::all_of(taxonomy_cols)),
                  by = "FeatureID"
                )
            }
          }

          if (is_classification) {
            y_group <- factor(as.character(y))
            if (nlevels(y_group) == 2) {
              g1 <- levels(y_group)[1]
              g2 <- levels(y_group)[2]
              g1_idx <- which(y_group == g1)
              g2_idx <- which(y_group == g2)
              g1_mean <- colMeans(x[g1_idx, , drop = FALSE], na.rm = TRUE)
              g2_mean <- colMeans(x[g2_idx, , drop = FALSE], na.rm = TRUE)
              diff_vec <- g2_mean - g1_mean
              mean_diff_df <- data.frame(
                Feature = names(diff_vec),
                Group1Mean = as.numeric(g1_mean[names(diff_vec)]),
                Group2Mean = as.numeric(g2_mean[names(diff_vec)]),
                Group1Median = as.numeric(apply(x[g1_idx, , drop = FALSE], 2, stats::median, na.rm = TRUE)[names(diff_vec)]),
                Group2Median = as.numeric(apply(x[g2_idx, , drop = FALSE], 2, stats::median, na.rm = TRUE)[names(diff_vec)]),
                MeanDifference = as.numeric(diff_vec),
                stringsAsFactors = FALSE
              )
              fc_vec <- (g2_mean + 1e-09) / (g1_mean + 1e-09)
              pval_vec <- vapply(names(diff_vec), function(feat) {
                suppressWarnings(tryCatch(
                  stats::wilcox.test(x[g1_idx, feat], x[g2_idx, feat])$p.value,
                  error = function(e) NA_real_
                ))
              }, numeric(1))
              direction_vec <- ifelse(diff_vec > 0, paste0("Higher in ", g2), paste0("Higher in ", g1))
              direction_vec[abs(diff_vec) < .Machine$double.eps] <- "No difference"
              mean_diff_df$FoldChange <- as.numeric(fc_vec[names(diff_vec)])
              mean_diff_df$Log2FoldChange <- log2(mean_diff_df$FoldChange)
              mean_diff_df$UnivariatePValue <- as.numeric(pval_vec[names(diff_vec)])
              mean_diff_df$Direction <- as.character(direction_vec[names(diff_vec)])
              colnames(mean_diff_df)[2] <- paste0("Mean_", g1)
              colnames(mean_diff_df)[3] <- paste0("Mean_", g2)
              colnames(mean_diff_df)[4] <- paste0("Median_", g1)
              colnames(mean_diff_df)[5] <- paste0("Median_", g2)
              colnames(mean_diff_df)[6] <- paste0("Diff_", g2, "_minus_", g1)
              imp_tbl <- imp_tbl %>% dplyr::left_join(mean_diff_df, by = "Feature")
            } else {
              group_means <- sapply(levels(y_group), function(gr) {
                colMeans(x[y_group == gr, , drop = FALSE], na.rm = TRUE)
              }, simplify = "matrix")
              if (is.null(dim(group_means))) {
                group_means <- matrix(group_means, ncol = 1)
                rownames(group_means) <- colnames(x)
              }
              max_mean <- apply(group_means, 1, max, na.rm = TRUE)
              min_mean <- apply(group_means, 1, min, na.rm = TRUE)
              diff_vec <- max_mean - min_mean
              max_group_idx <- apply(group_means, 1, function(v) {
                if (all(is.na(v))) return(NA_integer_)
                idx <- which(v == max(v, na.rm = TRUE))
                if (length(idx) == 1) idx else NA_integer_
              })
              direction_vec <- rep("No dominant group", length(max_group_idx))
              valid_idx <- which(!is.na(max_group_idx))
              if (length(valid_idx) > 0) {
                direction_vec[valid_idx] <- paste0("Higher in ", colnames(group_means)[max_group_idx[valid_idx]])
              }
              mean_diff_df <- data.frame(
                Feature = names(diff_vec),
                MaxGroupMean = as.numeric(max_mean[names(diff_vec)]),
                MinGroupMean = as.numeric(min_mean[names(diff_vec)]),
                MeanDifference_MaxMinusMin = as.numeric(diff_vec),
                Direction = as.character(direction_vec[names(diff_vec)]),
                UnivariatePValue = as.numeric(vapply(names(diff_vec), function(feat) {
                  suppressWarnings(tryCatch(
                    stats::kruskal.test(x[, feat] ~ y_group)$p.value,
                    error = function(e) NA_real_
                  ))
                }, numeric(1))),
                stringsAsFactors = FALSE
              )
              group_mean_df <- as.data.frame(group_means, stringsAsFactors = FALSE)
              group_mean_df$Feature <- rownames(group_mean_df)
              colnames(group_mean_df) <- c(paste0("Mean_", colnames(group_means)), "Feature")
              mean_diff_df <- mean_diff_df %>% dplyr::left_join(group_mean_df, by = "Feature")
              imp_tbl <- imp_tbl %>% dplyr::left_join(mean_diff_df, by = "Feature")
            }
          } else {
            pred_all <- as.numeric(stats::predict(rf_fit, newdata = x))
            corr_vec <- vapply(seq_len(ncol(x)), function(j) {
              suppressWarnings(stats::cor(x[[j]], pred_all, use = "complete.obs", method = "spearman"))
            }, numeric(1))
            corr_df <- data.frame(
              Feature = colnames(x),
              FeaturePredictionCorrelation = as.numeric(corr_vec),
              stringsAsFactors = FALSE
            )
            imp_tbl <- imp_tbl %>% dplyr::left_join(corr_df, by = "Feature")
          }

          leading_cols <- c("Feature", "FeatureID")
          taxonomy_cols_present <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(imp_tbl))
          metric_cols <- setdiff(colnames(imp_tbl), c(leading_cols, taxonomy_cols_present))
          imp_tbl <- imp_tbl[, c(leading_cols, taxonomy_cols_present, metric_cols), drop = FALSE]
          if (length(taxonomy_cols_present) > 0) {
            rank_prefix_map <- c(
              Kingdom = "k__",
              Phylum = "p__",
              Class = "c__",
              Order = "o__",
              Family = "f__",
              Genus = "g__",
              Species = "s__"
            )
            get_rank_unassigned_label <- function(rank_name) {
              prefix <- rank_prefix_map[[rank_name]]
              if (is.null(prefix) || !nzchar(prefix)) {
                prefix <- paste0(substr(tolower(rank_name), 1, 1), "__")
              }
              paste0(prefix, "Unassigned")
            }
            label_rank_order <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
            if (!identical(input$tax_level, "ASV") && input$tax_level %in% label_rank_order) {
              label_rank_order <- label_rank_order[seq_len(match(input$tax_level, label_rank_order))]
            }
            label_ranks_present <- intersect(label_rank_order, taxonomy_cols_present)
            tax_mat <- as.data.frame(imp_tbl[, label_ranks_present, drop = FALSE], stringsAsFactors = FALSE)
            for (cc in colnames(tax_mat)) tax_mat[[cc]] <- as.character(tax_mat[[cc]])
            taxa_label <- apply(tax_mat, 1, function(taxa_row) {
              parts <- vapply(seq_along(taxa_row), function(i) {
                v <- taxa_row[i]
                rk <- names(taxa_row)[i]
                if (is.na(v) || !nzchar(v)) {
                  get_rank_unassigned_label(rk)
                } else {
                  v
                }
              }, character(1))
              paste(parts, collapse = ";")
            })
            taxa_label[is.na(taxa_label) | !nzchar(taxa_label)] <- imp_tbl$FeatureID[is.na(taxa_label) | !nzchar(taxa_label)]
            imp_tbl$taxa_label <- taxa_label
          } else {
            imp_tbl$taxa_label <- imp_tbl$FeatureID
          }

          shap_values_long <- NULL
          shap_summary_tbl <- NULL
          shap_class_summary <- NULL
          shap_local_default <- NULL
          shap_target_class_used <- NULL
          shap_baseline <- NA_real_
          shap_status <- "SHAP is available for classification only."
          shap_approach_used <- NA_character_
          shap_iterative_used <- TRUE
          shap_max_n_coalitions_used <- 256L
          shap_fallback_used <- FALSE
          shap_model_scope <- "Full RF model features"
          shap_boot_n <- 100L
          shap_perm_n <- 200L

          if (is_classification) {
            if (!requireNamespace("shapr", quietly = TRUE)) {
              shap_status <- "Package 'shapr' is not installed. Install it to see SHAP results."
            } else {
              shap_target_class_used <- input$shap_target_class
              if (is.null(shap_target_class_used) || !shap_target_class_used %in% levels(y)) {
                shap_target_class_used <- if (nlevels(y) >= 2) levels(y)[2] else levels(y)[1]
              }

              build_predict_target_prob <- function(feature_frame) {
                function(model, newdata) {
                  nd <- as.data.frame(newdata, check.names = FALSE)
                  missing_cols <- setdiff(colnames(feature_frame), colnames(nd))
                  if (length(missing_cols) > 0) {
                    fill_vals <- colMeans(feature_frame, na.rm = TRUE)
                    for (cc in missing_cols) nd[[cc]] <- as.numeric(fill_vals[[cc]])
                  }
                  nd <- nd[, colnames(feature_frame), drop = FALSE]
                  prob <- stats::predict(model, newdata = nd, type = "prob")
                  prob <- as.data.frame(prob, check.names = FALSE)
                  as.numeric(prob[[shap_target_class_used]])
                }
              }
              build_model_specs <- function(feature_frame) {
                function(model) {
                  feature_labels <- colnames(feature_frame)
                  feature_classes <- vapply(feature_frame, function(col) class(col)[1], character(1))
                  factor_levels <- lapply(feature_frame, function(col) if (is.factor(col)) levels(col) else character(0))
                  list(labels = feature_labels, classes = feature_classes[feature_labels], factor_levels = factor_levels[feature_labels])
                }
              }
              run_shapr_explain <- function(model_obj, x_train_df, x_test_df) {
                predict_target_prob <- build_predict_target_prob(x_train_df)
                get_model_specs_rf <- build_model_specs(x_train_df)
                shap_baseline <<- mean(predict_target_prob(model_obj, x_train_df), na.rm = TRUE)
                shapr_ns <- asNamespace("shapr")
                explain_fn <- get("explain", envir = shapr_ns)
                call_candidates <- list(
                  list(x_explain = x_test_df, model = model_obj, x_train = x_train_df, predict_model = predict_target_prob),
                  list(x_explain = x_test_df, model = model_obj, x = x_train_df, predict_model = predict_target_prob)
                )
                approaches <- c("empirical", "independence")
                explanation <- NULL
                explain_errors <- character(0)
                for (args_i in call_candidates) {
                  for (ap in approaches) {
                    args_now <- args_i
                    args_now$approach <- ap
                    args_now$iterative <- TRUE
                    args_now$max_n_coalitions <- 256
                    args_now$get_model_specs <- get_model_specs_rf
                    args_now$phi0 <- shap_baseline
                    args_now$prediction_zero <- shap_baseline
                    call_res <- tryCatch(do.call(explain_fn, args_now), error = function(e) e)
                    if (!inherits(call_res, "error")) {
                      explanation <- call_res
                      shap_approach_used <<- ap
                      break
                    }
                    msg <- gsub("\u001b\\[[0-9;]*m", "", conditionMessage(call_res))
                    explain_errors <- c(explain_errors, paste0("[", ap, "] ", msg))
                  }
                  if (!is.null(explanation)) break
                }
                if (is.null(explanation)) {
                  stop(paste("No compatible shapr::explain() signature found.", paste(unique(explain_errors), collapse = " | ")))
                }
                explanation
              }
              parse_shap_output <- function(explanation, x_test_df, class_vec) {
                shap_matrix <- NULL
                if (!is.null(explanation$shapley_values)) {
                  shap_matrix <- as.data.frame(explanation$shapley_values, check.names = FALSE)
                } else if (!is.null(explanation$dt) && is.data.frame(explanation$dt)) {
                  dt_df <- as.data.frame(explanation$dt, check.names = FALSE)
                  maybe_cols <- intersect(colnames(x_test_df), colnames(dt_df))
                  if (length(maybe_cols) == ncol(x_test_df)) {
                    shap_matrix <- dt_df[, maybe_cols, drop = FALSE]
                  }
                }
                if (is.null(shap_matrix) && is.data.frame(explanation)) {
                  exp_df <- as.data.frame(explanation, check.names = FALSE)
                  maybe_cols <- intersect(colnames(x_test_df), colnames(exp_df))
                  if (length(maybe_cols) == ncol(x_test_df)) {
                    shap_matrix <- exp_df[, maybe_cols, drop = FALSE]
                  }
                }
                if (is.null(shap_matrix) && is.list(explanation)) {
                  for (nm in names(explanation)) {
                    obj <- explanation[[nm]]
                    if (is.null(obj) || !is.data.frame(obj)) next
                    obj_df <- as.data.frame(obj, check.names = FALSE)
                    maybe_cols <- intersect(colnames(x_test_df), colnames(obj_df))
                    if (length(maybe_cols) == ncol(x_test_df)) {
                      shap_matrix <- obj_df[, maybe_cols, drop = FALSE]
                      break
                    }
                  }
                }
                if (is.null(shap_matrix) && is.list(explanation) && !is.null(explanation$dt) && is.data.frame(explanation$dt)) {
                  dt <- as.data.frame(explanation$dt, check.names = FALSE)
                  long_candidates <- c("feature", "Feature", "variable", "name")
                  value_candidates <- c("phi", "Phi", "shapley_value", "value", "contribution")
                  feat_col <- long_candidates[long_candidates %in% colnames(dt)]
                  val_col <- value_candidates[value_candidates %in% colnames(dt)]
                  if (length(feat_col) > 0 && length(val_col) > 0) {
                    feat_col <- feat_col[1]
                    val_col <- val_col[1]
                    id_candidates <- c("id", "index", "row_id", "sample_id", "SampleID")
                    id_col <- id_candidates[id_candidates %in% colnames(dt)]
                    if (length(id_col) > 0) {
                      id_col <- id_col[1]
                      dt_sub <- dt[, c(id_col, feat_col, val_col), drop = FALSE]
                      colnames(dt_sub) <- c("SampleID", "Feature", "SHAP")
                      dt_sub$SampleID <- as.character(dt_sub$SampleID)
                      wide <- tryCatch(
                        tidyr::pivot_wider(dt_sub, names_from = Feature, values_from = SHAP),
                        error = function(e) NULL
                      )
                      if (!is.null(wide)) {
                        wide <- as.data.frame(wide, check.names = FALSE)
                        rn <- as.character(wide$SampleID)
                        wide$SampleID <- NULL
                        maybe_cols <- intersect(colnames(x_test_df), colnames(wide))
                        if (length(maybe_cols) > 0) {
                          shap_matrix <- wide[, maybe_cols, drop = FALSE]
                          rownames(shap_matrix) <- rn
                        }
                      }
                    }
                  }
                }
                if (is.null(shap_matrix)) {
                  top_names <- if (is.list(explanation)) paste(names(explanation), collapse = ", ") else class(explanation)[1]
                  stop(paste0("Could not parse SHAP values from shapr output. Available fields: ", top_names))
                }
                shap_matrix$SampleID <- rownames(x_test_df)
                if (is.null(shap_matrix$SampleID)) shap_matrix$SampleID <- paste0("Sample_", seq_len(nrow(shap_matrix)))
                shap_values_long <- shap_matrix %>%
                  tidyr::pivot_longer(cols = -SampleID, names_to = "Feature", values_to = "SHAP") %>%
                  dplyr::left_join(data.frame(SampleID = rownames(x_test_df), ActualClass = as.character(class_vec), stringsAsFactors = FALSE), by = "SampleID")
                list(shap_values_long = shap_values_long)
              }
              shap_try <- tryCatch({
                x_train_df <- as.data.frame(x_train, check.names = FALSE)
                x_test_df <- as.data.frame(x_test, check.names = FALSE)
                explanation <- run_shapr_explain(rf_fit, x_train_df, x_test_df)
                parsed <- parse_shap_output(explanation, x_test_df, y_test)
                shap_values_long <- parsed$shap_values_long
                shap_summary_tbl <- build_shap_summary_metrics(
                  shap_values_long = shap_values_long,
                  feature_map = feature_map,
                  imp_tbl = imp_tbl,
                  seed_val = as.integer(input$seed),
                  n_boot = shap_boot_n,
                  n_perm = shap_perm_n
                )
                shap_class_summary <- shap_values_long %>%
                  dplyr::group_by(ActualClass, Feature) %>%
                  dplyr::summarise(MeanAbsSHAP = mean(abs(SHAP), na.rm = TRUE), MeanSHAP = mean(SHAP, na.rm = TRUE), .groups = "drop") %>%
                  dplyr::arrange(ActualClass, dplyr::desc(MeanAbsSHAP))
                first_sample <- shap_values_long$SampleID[1]
                shap_local_default <- shap_values_long %>%
                  dplyr::filter(SampleID == first_sample) %>%
                  dplyr::left_join(feature_map, by = c("Feature" = "FeatureModel")) %>%
                  dplyr::arrange(dplyr::desc(abs(SHAP)))
                NULL
              }, error = function(e) conditionMessage(e))
              if (!is.null(shap_try) && grepl("solve\\(\\): solution not found", shap_try)) {
                shap_try <- tryCatch({
                  shap_fallback_used <- TRUE
                  top_n <- min(30L, ncol(x_train))
                  shap_model_scope <- paste0("Fallback RF model (top ", top_n, " features)")
                  top_feats <- imp_tbl %>% dplyr::slice_max(order_by = PermutationImportance, n = top_n, with_ties = FALSE) %>% dplyr::pull(Feature)
                  x_train_small <- as.data.frame(x_train[, top_feats, drop = FALSE], check.names = FALSE)
                  x_test_small <- as.data.frame(x_test[, top_feats, drop = FALSE], check.names = FALSE)
                  rf_shap <- randomForest::randomForest(
                    x = x_train_small, y = y_train,
                    ntree = max(200L, floor(as.integer(input$ntree) / 2)),
                    mtry = max(1L, floor(sqrt(ncol(x_train_small)))),
                    importance = FALSE
                  )
                  explanation <- run_shapr_explain(rf_shap, x_train_small, x_test_small)
                  parsed <- parse_shap_output(explanation, x_test_small, y_test)
                  shap_values_long <- parsed$shap_values_long
                  shap_summary_tbl <- build_shap_summary_metrics(
                    shap_values_long = shap_values_long,
                    feature_map = feature_map,
                    imp_tbl = imp_tbl,
                    seed_val = as.integer(input$seed),
                    n_boot = shap_boot_n,
                    n_perm = shap_perm_n
                  )
                  shap_class_summary <- shap_values_long %>%
                    dplyr::group_by(ActualClass, Feature) %>%
                    dplyr::summarise(MeanAbsSHAP = mean(abs(SHAP), na.rm = TRUE), MeanSHAP = mean(SHAP, na.rm = TRUE), .groups = "drop") %>%
                    dplyr::arrange(ActualClass, dplyr::desc(MeanAbsSHAP))
                  first_sample <- shap_values_long$SampleID[1]
                  shap_local_default <- shap_values_long %>%
                    dplyr::filter(SampleID == first_sample) %>%
                    dplyr::left_join(feature_map, by = c("Feature" = "FeatureModel")) %>%
                    dplyr::arrange(dplyr::desc(abs(SHAP)))
                  shap_status <<- paste0("SHAP computed with fallback model (top ", top_n, " features).")
                  NULL
                }, error = function(e) conditionMessage(e))
              }

              if (is.null(shap_try)) {
                shap_status <- paste0("SHAP computed on test samples: ", length(unique(shap_values_long$SampleID)))
              } else {
                shap_status <- paste0("SHAP calculation failed: ", shap_try)
              }
            }
          }

          shap_section_text <- paste0(
            "SHAP status: ", shap_status,
            if (!is.null(shap_target_class_used) && nzchar(shap_target_class_used)) {
              paste0("\nSHAP target class: ", shap_target_class_used)
            } else "",
            if (is.finite(shap_baseline)) {
              paste0("\nSHAP baseline probability: ", sprintf("%.6f", shap_baseline))
            } else "",
            if (is_classification) {
              paste0(
                "\nSHAP options:",
                "\n- package: shapr",
                "\n- approach used: ", ifelse(is.na(shap_approach_used), "NA", shap_approach_used),
                "\n- iterative: ", ifelse(isTRUE(shap_iterative_used), "TRUE", "FALSE"),
                "\n- max_n_coalitions: ", as.integer(shap_max_n_coalitions_used),
                "\n- model scope: ", shap_model_scope,
                "\n- fallback used: ", ifelse(isTRUE(shap_fallback_used), "TRUE", "FALSE"),
                "\n- stability bootstrap reps: ", as.integer(shap_boot_n),
                "\n- empirical p-value permutations: ", as.integer(shap_perm_n)
              )
            } else ""
          )

          metrics_text <- paste0(
            metrics_text,
            "\n\n--------------------\n[RF]\n",
            rf_options_text,
            "\n\n--------------------\n[SHAP]\n",
            shap_section_text
          )

          model_result(list(
            metrics = metrics_text,
            importance = imp_tbl,
            score_col = "PermutationImportance",
            roc = roc_payload,
            shap_values = shap_values_long,
            shap_summary = shap_summary_tbl,
            shap_class_summary = shap_class_summary,
            shap_local_default = shap_local_default,
            shap_target_class = shap_target_class_used,
            shap_baseline = shap_baseline
          ))

          status_text("Model completed.")
          log_rf("run completed")
        }, warning = function(w) {
          log_rf(paste("warning:", conditionMessage(w)))
          invokeRestart("muffleWarning")
        })
      }, error = function(e) {
        log_rf(paste("error:", conditionMessage(e)))
        status_text(paste("Error:", conditionMessage(e)))
      })
    })

    output$rf_metrics <- renderText({
      res <- model_result()
      if (is.null(res)) return("Run the model to see results.")
      res$metrics
    })

    output$rf_metrics_box <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 1160 else input$plot_width
      tags$div(
        class = "simple-result-card",
        style = paste0("width: ", box_width, "px; max-width: 100%;"),
        verbatimTextOutput(session$ns("rf_metrics"))
      )
    })

    output$importance_table <- renderDT({
      res <- model_result()
      req(res)
      datatable(
        res$importance,
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE),
        class = "compact stripe hover cell-border",
        style = "bootstrap"
      ) %>%
        formatStyle(
          columns = colnames(res$importance),
          `font-size` = "11px",
          `padding` = "4px"
        )
    })

    output$download_rf_table <- downloadHandler(
      filename = function() {
        paste0("random_forest_table_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv")
      },
      content = function(file) {
        res <- model_result()
        req(res)
        utils::write.table(
          res$importance,
          file = file,
          sep = "\t",
          row.names = FALSE,
          quote = FALSE,
          na = ""
        )
      }
    )

    output$importance_plot_ui <- renderUI({
      w <- as.integer(input$plot_width)
      h <- as.integer(input$plot_height)
      w <- ifelse(is.na(w) || w < 100, 1160L, w)
      h <- ifelse(is.na(h) || h < 100, 500L, h)
      plotOutput(
        session$ns("importance_plot"),
        width = paste0(w, "px"),
        height = paste0(h, "px")
      )
    })

    output$shap_summary_plot_ui <- renderUI({
      w <- as.integer(input$plot_width)
      h <- as.integer(input$plot_height)
      w <- ifelse(is.na(w) || w < 100, 1160L, w)
      h <- ifelse(is.na(h) || h < 100, 500L, h)
      plotOutput(
        session$ns("shap_summary_plot"),
        width = paste0(w, "px"),
        height = paste0(h, "px")
      )
    })

    derive_shap_summary <- function(res) {
      if (!is.null(res$shap_summary) && nrow(res$shap_summary) > 0) {
        return(res$shap_summary)
      }
      if (is.null(res$shap_values) || nrow(res$shap_values) == 0) {
        return(NULL)
      }
      res$shap_values %>%
        dplyr::group_by(Feature) %>%
        dplyr::summarise(
          MeanAbsSHAP = mean(abs(SHAP), na.rm = TRUE),
          MeanSHAP = mean(SHAP, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::left_join(
          res$importance %>%
            dplyr::select(
              Feature,
              dplyr::any_of(c("FeatureID", "taxa_label", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
            ),
          by = "Feature"
        ) %>%
        dplyr::arrange(dplyr::desc(MeanAbsSHAP))
    }
    build_shap_summary_metrics <- function(shap_values_long, feature_map, imp_tbl, seed_val = 1234L, n_boot = 50L, n_perm = 200L) {
      if (is.null(shap_values_long) || nrow(shap_values_long) == 0) return(NULL)
      set.seed(as.integer(seed_val))
      by_feature <- split(shap_values_long$SHAP, shap_values_long$Feature)
      feats <- names(by_feature)
      metric_list <- lapply(feats, function(ft) {
        v <- as.numeric(by_feature[[ft]])
        v <- v[is.finite(v)]
        if (length(v) == 0) {
          return(data.frame(
            Feature = ft, N = 0L, MeanAbsSHAP = NA_real_, MeanSHAP = NA_real_,
            ImpactScore = NA_real_, DirectionScore = NA_real_,
            StabilityCV = NA_real_, StabilityScore = NA_real_,
            BootMeanAbs_LCL = NA_real_, BootMeanAbs_UCL = NA_real_,
            EmpiricalP_MeanSHAP = NA_real_, WilcoxonP_MeanSHAP = NA_real_,
            stringsAsFactors = FALSE
          ))
        }
        mean_abs <- mean(abs(v), na.rm = TRUE)
        mean_raw <- mean(v, na.rm = TRUE)
        direction_score <- if (is.na(mean_abs) || mean_abs <= 0) NA_real_ else mean_raw / mean_abs
        boot_means <- replicate(as.integer(n_boot), {
          idx <- sample.int(length(v), size = length(v), replace = TRUE)
          mean(abs(v[idx]), na.rm = TRUE)
        })
        boot_means <- as.numeric(boot_means)
        st_cv <- if (is.na(mean(boot_means)) || mean(boot_means) == 0) NA_real_ else stats::sd(boot_means) / mean(boot_means)
        st_score <- if (is.na(st_cv)) NA_real_ else 1 / (1 + st_cv)
        boot_lcl <- stats::quantile(boot_means, probs = 0.025, na.rm = TRUE, names = FALSE)
        boot_ucl <- stats::quantile(boot_means, probs = 0.975, na.rm = TRUE, names = FALSE)
        obs <- abs(mean_raw)
        perm_vals <- replicate(as.integer(n_perm), {
          sgn <- sample(c(-1, 1), size = length(v), replace = TRUE)
          abs(mean(v * sgn, na.rm = TRUE))
        })
        perm_vals <- as.numeric(perm_vals)
        emp_p <- (1 + sum(perm_vals >= obs, na.rm = TRUE)) / (as.integer(n_perm) + 1)
        wil_p <- suppressWarnings(tryCatch(stats::wilcox.test(v, mu = 0)$p.value, error = function(e) NA_real_))
        data.frame(
          Feature = ft,
          N = length(v),
          MeanAbsSHAP = mean_abs,
          MeanSHAP = mean_raw,
          ImpactScore = mean_abs,
          DirectionScore = direction_score,
          StabilityCV = st_cv,
          StabilityScore = st_score,
          BootMeanAbs_LCL = as.numeric(boot_lcl),
          BootMeanAbs_UCL = as.numeric(boot_ucl),
          EmpiricalP_MeanSHAP = as.numeric(emp_p),
          WilcoxonP_MeanSHAP = as.numeric(wil_p),
          stringsAsFactors = FALSE
        )
      })
      summary_df <- dplyr::bind_rows(metric_list)
      summary_df$EmpiricalFDR_MeanSHAP <- stats::p.adjust(summary_df$EmpiricalP_MeanSHAP, method = "BH")
      summary_df$WilcoxonFDR_MeanSHAP <- stats::p.adjust(summary_df$WilcoxonP_MeanSHAP, method = "BH")
      summary_df %>%
        dplyr::left_join(feature_map, by = c("Feature" = "FeatureModel")) %>%
        dplyr::left_join(
          imp_tbl %>% dplyr::select(
            Feature,
            dplyr::any_of(c("taxa_label", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
          ),
          by = "Feature"
        ) %>%
        dplyr::arrange(dplyr::desc(ImpactScore))
    }

    output$shap_summary_table <- renderDT({
      res <- model_result()
      req(res)
      ss <- derive_shap_summary(res)
      validate(need(!is.null(ss) && nrow(ss) > 0, "SHAP summary is not available."))
      display <- ss %>%
        dplyr::mutate(
          TargetClass = res$shap_target_class,
          BaselineProbability = res$shap_baseline
        )
      front_cols <- c(
        "Feature", "taxa_label",
        "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
        "TargetClass"
      )
      # Hide FeatureID from SHAP table display while keeping it in underlying data objects.
      display <- display %>%
        dplyr::select(-dplyr::any_of("FeatureID"))
      ordered_cols <- c(front_cols[front_cols %in% colnames(display)], setdiff(colnames(display), front_cols))
      display <- display[, ordered_cols, drop = FALSE]

      datatable(
        display,
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE),
        class = "compact stripe hover cell-border",
        style = "bootstrap"
      ) %>%
        formatStyle(columns = colnames(display), `font-size` = "11px", `padding` = "4px")
    })

    compute_plot_geom <- function(plot_w, plot_h, n_rows) {
      w <- as.integer(plot_w)
      h <- as.integer(plot_h)
      if (is.na(w) || w < 100) w <- 1160L
      if (is.na(h) || h < 100) h <- 500L
      nr <- max(1L, as.integer(n_rows))
      cell_mm <- max(3, min(14, (h * 0.72) / (nr * 3.78)))
      bar_width_pt <- max(90, min(360, w * 0.16 * 0.75))
      list(cell_mm = cell_mm, bar_width_pt = bar_width_pt)
    }

    build_shap_plot <- function(res, top_n = 20, taxa_fontsize = 8, cell_mm = 4.4, bar_width_pt = 120, base_size = 11) {
      ss <- derive_shap_summary(res)
      validate(need(!is.null(ss) && nrow(ss) > 0, "SHAP summary is not available."))
      ss <- as.data.frame(ss, check.names = FALSE)
      if (!"taxa_label" %in% colnames(ss)) {
        ss <- ss %>%
          dplyr::left_join(
            res$importance %>% dplyr::select(Feature, taxa_label),
            by = "Feature"
          )
      }
      if (!"taxa_label" %in% colnames(ss)) ss$taxa_label <- NA_character_

      top_shap <- ss %>%
        dplyr::mutate(
          Feature = as.character(Feature),
          taxa_label = as.character(taxa_label),
          taxa_label_plot = dplyr::if_else(is.na(taxa_label) | !nzchar(taxa_label), Feature, taxa_label),
          MeanAbsSHAP = as.numeric(MeanAbsSHAP)
        ) %>%
        dplyr::arrange(dplyr::desc(MeanAbsSHAP)) %>%
        utils::head(top_n)

      dup_label <- duplicated(top_shap$taxa_label_plot) | duplicated(top_shap$taxa_label_plot, fromLast = TRUE)
      if (any(dup_label)) {
        top_shap$taxa_label_plot[dup_label] <- paste0(
          top_shap$taxa_label_plot[dup_label],
          " [",
          top_shap$FeatureID[dup_label],
          "]"
        )
      }

      max_imp <- suppressWarnings(max(top_shap$MeanAbsSHAP, na.rm = TRUE))
      if (!is.finite(max_imp) || max_imp <= 0) max_imp <- 1
      bar_x_max <- max(0.1, ceiling(max_imp * 10) / 10)
      feat_levels <- rev(top_shap$taxa_label_plot[order(top_shap$MeanAbsSHAP)])
      top_shap <- top_shap[match(feat_levels, top_shap$taxa_label_plot), , drop = FALSE]

      if (!is.null(res$shap_class_summary) && nrow(res$shap_class_summary) > 0) {
        class_df <- res$shap_class_summary %>%
          dplyr::filter(Feature %in% top_shap$Feature) %>%
          dplyr::group_by(Feature, ActualClass) %>%
          dplyr::summarise(MeanAbsSHAP = mean(MeanAbsSHAP, na.rm = TRUE), .groups = "drop")
        class_wide <- class_df %>%
          tidyr::pivot_wider(names_from = ActualClass, values_from = MeanAbsSHAP, values_fill = 0)
        class_wide <- as.data.frame(class_wide, check.names = FALSE)
        rownames(class_wide) <- class_wide$Feature
        class_wide$Feature <- NULL
        class_wide <- class_wide[top_shap$Feature, , drop = FALSE]
        if (ncol(class_wide) > 0) {
          class_names <- colnames(class_wide)
          fill_key <- sapply(seq_len(nrow(class_wide)), function(i) {
            row_vals <- as.numeric(class_wide[i, ])
            if (all(!is.finite(row_vals))) return(rep("NA", length(class_names)))
            winner_idx <- which(row_vals == max(row_vals, na.rm = TRUE))
            if (length(winner_idx) == 0) return(rep("Not dominant", length(class_names)))
            out <- rep("Not dominant", length(class_names))
            out[winner_idx[1]] <- class_names[winner_idx[1]]
            out
          })
          hm_mat <- t(fill_key)
          rownames(hm_mat) <- top_shap$taxa_label_plot
          colnames(hm_mat) <- class_names
          dom_classes <- sort(unique(as.vector(hm_mat)))
          dom_classes <- dom_classes[!dom_classes %in% c("Not dominant", "NA")]
          hm_col <- setNames(
            grDevices::hcl.colors(max(1, length(dom_classes)), palette = "Dark 3")[seq_along(dom_classes)],
            dom_classes
          )
          hm_col <- c(hm_col, "Not dominant" = "grey92", "NA" = "grey70")
        } else {
          hm_mat <- matrix("Not available", nrow = nrow(top_shap), ncol = 1)
          rownames(hm_mat) <- top_shap$taxa_label_plot
          colnames(hm_mat) <- "Class"
          hm_col <- c("Not available" = "grey92")
        }
      } else {
        hm_mat <- matrix("Not available", nrow = nrow(top_shap), ncol = 1)
        rownames(hm_mat) <- top_shap$taxa_label_plot
        colnames(hm_mat) <- "Class"
        hm_col <- c("Not available" = "grey92")
      }

      row_label_width <- ComplexHeatmap::max_text_width(
        rownames(hm_mat),
        gp = grid::gpar(fontsize = taxa_fontsize)
      ) + grid::unit(3, "mm")

      hm <- ComplexHeatmap::Heatmap(
        hm_mat,
        name = "Class",
        col = hm_col,
        rect_gp = grid::gpar(type = "none"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          side_mm <- min(
            grid::convertWidth(w, "mm", valueOnly = TRUE),
            grid::convertHeight(h, "mm", valueOnly = TRUE)
          ) * 0.84
          side <- grid::unit(side_mm, "mm")
          grid::grid.rect(
            x = x, y = y, width = side, height = side,
            gp = grid::gpar(fill = fill, col = "#1a1a1a", lwd = 0.5)
          )
        },
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_heatmap_legend = FALSE,
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = taxa_fontsize),
        row_names_max_width = row_label_width,
        column_names_gp = grid::gpar(fontsize = 9),
        column_names_rot = 90,
        width = grid::unit(cell_mm * ncol(hm_mat), "mm"),
        height = grid::unit(cell_mm * nrow(hm_mat), "mm")
      )

      bar_anno <- ComplexHeatmap::rowAnnotation(
        `Mean(|SHAP|)` = ComplexHeatmap::anno_barplot(
          top_shap$MeanAbsSHAP,
          gp = grid::gpar(fill = "#1f78b4", col = NA),
          border = FALSE,
          axis_param = list(
            at = c(0, bar_x_max / 2, bar_x_max),
            labels = sprintf("%.2f", c(0, bar_x_max / 2, bar_x_max)),
            side = "bottom",
            gp = grid::gpar(fontsize = 8)
          ),
          ylim = c(0, bar_x_max),
          bar_width = 0.7,
          width = grid::unit(bar_width_pt, "pt")
        )
      )
      hm + bar_anno
    }

    output$shap_summary_plot <- renderPlot(
      {
        if (is.null(input$run_rf) || input$run_rf < 1) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Click 'Run Random Forest' to start analysis.", cex = 0.85)
          return(invisible(NULL))
        }
        if (isTRUE(rf_running())) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Random Forest is running. Please wait...", cex = 0.85)
          return(invisible(NULL))
        }
        res <- model_result()
        if (is.null(res)) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Random Forest results are not ready yet. Please wait...", cex = 0.85)
          return(invisible(NULL))
        }
        ss <- derive_shap_summary(res)
        validate(need(!is.null(ss) && nrow(ss) > 0, "SHAP plot is not available."))
        geom <- compute_plot_geom(input$plot_width, input$plot_height, n_rows = 20L)
        p <- build_shap_plot(
          res,
          top_n = 20,
          taxa_fontsize = 8,
          cell_mm = geom$cell_mm,
          bar_width_pt = geom$bar_width_pt,
          base_size = input$base_size
        )
        ComplexHeatmap::draw(
          p,
          merge_legend = TRUE,
          padding = grid::unit(c(6, 2, 6, 130), "mm")
        )
      },
      width = function() {
        req(input$plot_width)
        as.integer(input$plot_width)
      },
      height = function() {
        req(input$plot_height)
        as.integer(input$plot_height)
      }
    )

    output$download_rf_shap_summary <- downloadHandler(
      filename = function() {
        paste0("random_forest_shap_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv")
      },
      content = function(file) {
        res <- model_result()
        req(res)
        ss <- derive_shap_summary(res)
        validate(need(!is.null(ss) && nrow(ss) > 0, "SHAP summary is not available."))
        ss_out <- ss %>% dplyr::select(-dplyr::any_of("FeatureID"))
        utils::write.table(ss_out, file = file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
      }
    )

    output$download_rf_shap_plot <- downloadHandler(
      filename = function() {
        paste0("random_forest_shap_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
      },
      content = function(file) {
        res <- model_result()
        req(res)
        ss <- derive_shap_summary(res)
        validate(need(!is.null(ss) && nrow(ss) > 0, "SHAP plot is not available."))
        w <- as.integer(input$plot_width)
        h <- as.integer(input$plot_height)
        w <- ifelse(is.na(w) || w < 100, 1160L, w)
        h <- ifelse(is.na(h) || h < 100, 500L, h)
        grDevices::png(filename = file, width = w, height = h, res = 72)
        geom <- compute_plot_geom(input$plot_width, input$plot_height, n_rows = 25L)
        p <- build_shap_plot(
          res,
          top_n = 25,
          taxa_fontsize = 7,
          cell_mm = geom$cell_mm,
          bar_width_pt = geom$bar_width_pt,
          base_size = input$base_size
        )
        ComplexHeatmap::draw(
          p,
          merge_legend = TRUE,
          padding = grid::unit(c(6, 2, 6, 130), "mm")
        )
        grDevices::dev.off()
      }
    )

    build_importance_plot <- function(res, taxa_fontsize = 8, cell_mm = 4.4, bar_width_pt = 120, base_size = 11) {
      top_imp <- head(res$importance, 20)
      top_imp$taxa_label_plot <- as.character(top_imp$taxa_label)
      dup_label <- duplicated(top_imp$taxa_label_plot) | duplicated(top_imp$taxa_label_plot, fromLast = TRUE)
      if (any(dup_label)) {
        top_imp$taxa_label_plot[dup_label] <- paste0(
          top_imp$taxa_label_plot[dup_label],
          " [",
          top_imp$FeatureID[dup_label],
          "]"
        )
      }
      max_imp <- suppressWarnings(max(top_imp$PermutationImportance, na.rm = TRUE))
      if (!is.finite(max_imp) || max_imp <= 0) max_imp <- 1
      bar_x_max <- max(1, ceiling(max_imp * 10) / 10)
      feat_levels <- rev(top_imp$taxa_label_plot[order(top_imp$PermutationImportance)])
      top_imp <- top_imp[match(feat_levels, top_imp$taxa_label_plot), , drop = FALSE]

      mean_cols <- grep("^Mean_", colnames(top_imp), value = TRUE)
      if (length(mean_cols) == 0) {
        hm_mat <- matrix("Not available", nrow = nrow(top_imp), ncol = 1)
        rownames(hm_mat) <- top_imp$taxa_label_plot
        colnames(hm_mat) <- "Group"
        hm_col <- c("Not available" = "grey92")
      } else {
        mean_mat <- as.matrix(top_imp[, mean_cols, drop = FALSE])
        feature_winner <- apply(mean_mat, 1, function(v) {
          if (all(is.na(v))) return(NA_character_)
          idx <- which(v == max(v, na.rm = TRUE))
          if (length(idx) >= 1) mean_cols[idx[1]] else NA_character_
        })
        fill_key <- sapply(seq_len(nrow(mean_mat)), function(i) {
          row_vals <- mean_mat[i, ]
          w <- unname(feature_winner[i])
          sapply(seq_along(row_vals), function(j) {
            col_name <- mean_cols[j]
            if (is.na(row_vals[j])) return("NA")
            if (!is.na(w) && col_name == w) return(sub("^Mean_", "", col_name))
            "Not dominant"
          })
        })
        hm_mat <- t(fill_key)
        rownames(hm_mat) <- top_imp$taxa_label_plot
        colnames(hm_mat) <- sub("^Mean_", "", mean_cols)

        dominant_groups <- unique(as.vector(hm_mat))
        dominant_groups <- dominant_groups[!dominant_groups %in% c("Not dominant", "NA")]
        dominant_groups <- sort(dominant_groups)
        hm_col <- setNames(
          grDevices::hcl.colors(max(1, length(dominant_groups)), palette = "Dark 3")[seq_along(dominant_groups)],
          dominant_groups
        )
        hm_col <- c(hm_col, "Not dominant" = "grey92", "NA" = "grey70")
      }

      row_label_width <- ComplexHeatmap::max_text_width(
        rownames(hm_mat),
        gp = grid::gpar(fontsize = taxa_fontsize)
      ) + grid::unit(3, "mm")

      hm <- ComplexHeatmap::Heatmap(
        hm_mat,
        name = "Group",
        col = hm_col,
        rect_gp = grid::gpar(type = "none"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          side_mm <- min(
            grid::convertWidth(w, "mm", valueOnly = TRUE),
            grid::convertHeight(h, "mm", valueOnly = TRUE)
          ) * 0.84
          side <- grid::unit(side_mm, "mm")
          grid::grid.rect(
            x = x,
            y = y,
            width = side,
            height = side,
            gp = grid::gpar(fill = fill, col = "#1a1a1a", lwd = 0.5)
          )
        },
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_heatmap_legend = FALSE,
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = taxa_fontsize),
        row_names_max_width = row_label_width,
        column_names_gp = grid::gpar(fontsize = 9),
        column_names_rot = 90,
        width = grid::unit(cell_mm * ncol(hm_mat), "mm"),
        height = grid::unit(cell_mm * nrow(hm_mat), "mm")
      )

      bar_anno <- ComplexHeatmap::rowAnnotation(
        `Permutation importance` = ComplexHeatmap::anno_barplot(
          top_imp$PermutationImportance,
          gp = grid::gpar(fill = "#2c7fb8", col = NA),
          border = FALSE,
          axis_param = list(
            at = c(0, bar_x_max / 2, bar_x_max),
            labels = sprintf("%.2f", c(0, bar_x_max / 2, bar_x_max)),
            side = "bottom",
            gp = grid::gpar(fontsize = 8)
          ),
          ylim = c(0, bar_x_max),
          bar_width = 0.7,
          width = grid::unit(bar_width_pt, "pt")
        )
      )

      hm + bar_anno
    }

    output$importance_plot <- renderPlot(
      {
        if (is.null(input$run_rf) || input$run_rf < 1) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Click 'Run Random Forest' to start analysis.", cex = 0.85)
          return(invisible(NULL))
        }
        if (isTRUE(rf_running())) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Random Forest is running. Please wait...", cex = 0.85)
          return(invisible(NULL))
        }
        res <- model_result()
        if (is.null(res)) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Random Forest results are not ready yet. Please wait...", cex = 0.85)
          return(invisible(NULL))
        }
        geom <- compute_plot_geom(input$plot_width, input$plot_height, n_rows = 20L)
        p <- build_importance_plot(
          res,
          taxa_fontsize = 8,
          cell_mm = geom$cell_mm,
          bar_width_pt = geom$bar_width_pt,
          base_size = input$base_size
        )
        ComplexHeatmap::draw(
          p,
          merge_legend = TRUE,
          padding = grid::unit(c(6, 2, 6, 130), "mm")
        )
      },
      width = function() {
        req(input$plot_width)
        as.integer(input$plot_width)
      },
      height = function() {
        req(input$plot_height)
        as.integer(input$plot_height)
      }
    )

    output$download_rf_barplot <- downloadHandler(
      filename = function() {
        paste0("random_forest_barplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
      },
      content = function(file) {
        res <- model_result()
        req(res)
        w <- as.integer(input$plot_width)
        h <- as.integer(input$plot_height)
        w <- ifelse(is.na(w) || w < 100, 1160L, w)
        h <- ifelse(is.na(h) || h < 100, 500L, h)
        grDevices::png(filename = file, width = w, height = h, res = 72)
        geom <- compute_plot_geom(input$plot_width, input$plot_height, n_rows = 20L)
        p <- build_importance_plot(
          res,
          taxa_fontsize = 7,
          cell_mm = geom$cell_mm,
          bar_width_pt = geom$bar_width_pt,
          base_size = input$base_size
        )
        ComplexHeatmap::draw(
          p,
          merge_legend = TRUE,
          padding = grid::unit(c(6, 2, 6, 130), "mm")
        )
        grDevices::dev.off()
      }
    )

    output$roc_auc_text <- renderText({
      res <- model_result()
      req(res)
      if (is.null(res$roc) || length(res$roc$auc_by_class) == 0) {
        return("ROC/AUC is available only for classification results with valid test-set probabilities.")
      }
      auc_lines <- paste0(
        names(res$roc$auc_by_class),
        ": ",
        sprintf("%.4f", as.numeric(res$roc$auc_by_class))
      )
      paste0(
        "AUC by class (one-vs-rest):\n",
        paste(auc_lines, collapse = "\n"),
        "\n\nMacro AUC: ",
        sprintf("%.4f", as.numeric(res$roc$macro_auc))
      )
    })

    build_roc_plot <- function(res, base_size = 13) {
      curve_names <- names(res$roc$curves)
      curve_df <- dplyr::bind_rows(lapply(curve_names, function(cl) {
        df <- res$roc$curves[[cl]]$curve
        df$Class <- cl
        df
      }))
      auc_df <- data.frame(
        Class = curve_names,
        AUC = as.numeric(res$roc$auc_by_class[curve_names]),
        stringsAsFactors = FALSE
      )
      legend_labels <- setNames(
        paste0(auc_df$Class, " (AUC=", sprintf("%.3f", auc_df$AUC), ")"),
        auc_df$Class
      )

      ggplot2::ggplot(curve_df, ggplot2::aes(x = FPR, y = TPR, color = Class)) +
        ggplot2::geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "dashed",
          linewidth = 0.6,
          color = "grey60"
        ) +
        ggplot2::geom_path(linewidth = 1.2, alpha = 0.95) +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
        ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
        ggplot2::scale_color_manual(
          values = setNames(grDevices::hcl.colors(length(curve_names), "Dark 3"), curve_names),
          labels = legend_labels,
          breaks = curve_names
        ) +
        ggplot2::labs(
          title = "ROC Curves (One-vs-Rest)",
          subtitle = paste0("Macro AUC = ", sprintf("%.3f", as.numeric(res$roc$macro_auc))),
          x = "False Positive Rate",
          y = "True Positive Rate",
          color = "Class"
        ) +
        ggplot2::theme_minimal(base_size = base_size) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_line(color = "#e6e6e6", linewidth = 0.4),
          legend.position = "right",
          legend.title = ggplot2::element_text(face = "bold"),
          legend.text = ggplot2::element_text(size = base_size - 2)
        )
    }

    output$roc_plot <- renderPlot(
      {
        if (is.null(input$run_rf) || input$run_rf < 1) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Click 'Run Random Forest' to start analysis.", cex = 0.85)
          return(invisible(NULL))
        }
        if (isTRUE(rf_running())) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Random Forest is running. Please wait...", cex = 0.85)
          return(invisible(NULL))
        }
        res <- model_result()
        if (is.null(res)) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "Random Forest results are not ready yet. Please wait...", cex = 0.85)
          return(invisible(NULL))
        }
        if (is.null(res$roc) || length(res$roc$curves) == 0) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "ROC curve is available for classification results only.", cex = 0.85)
          return(invisible(NULL))
        }
        p <- build_roc_plot(res, base_size = input$base_size)
        print(p)
      },
      width = 1160,
      height = 500
    )

    output$download_rf_rocplot <- downloadHandler(
      filename = function() {
        paste0("random_forest_roc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
      },
      content = function(file) {
        res <- model_result()
        req(res)
        grDevices::png(filename = file, width = 3600, height = 1800, res = 300)
        if (is.null(res$roc) || length(res$roc$curves) == 0) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "ROC curve is available for classification results only.", cex = 0.85)
          grDevices::dev.off()
          return(invisible(NULL))
        }
        p <- build_roc_plot(res, base_size = input$base_size)
        print(p)
        grDevices::dev.off()
      }
    )

    output$rf_legend_box <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 1160 else input$plot_width
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
          uiOutput(session$ns("rf_figure_legend"))
        )
      )
    })

    output$rf_results_separator <- renderUI({
      req(input$plot_width)
      box_width <- if (is.null(input$plot_width) || !is.finite(input$plot_width)) 1160 else input$plot_width
      tags$hr(
        style = paste(
          sprintf("width: %spx;", box_width),
          "max-width: 100%;",
          "margin: 14px 0 12px 0;",
          "border-top: 1px solid #d1d5db;"
        )
      )
    })

    output$rf_figure_legend <- renderUI({
      req(input$outcome_var, input$tax_level, input$transform_method)
      active_tab <- if (is.null(input$rf_active_tab) || !nzchar(input$rf_active_tab)) "RF Bar plot" else input$rf_active_tab
      model_mode <- if (identical(input$outcome_type, "Auto")) "Auto-detected" else input$outcome_type
      validation_label <- if (is.null(input$validation_mode) || !nzchar(input$validation_mode)) "Holdout split" else input$validation_mode

      legend_title <- if (identical(active_tab, "ROC / AUC")) {
        "Random Forest ROC/AUC plot"
      } else if (identical(active_tab, "SHAP Bar plot")) {
        "Random Forest SHAP importance plot"
      } else if (identical(active_tab, "SHAP Table")) {
        "Random Forest SHAP summary table"
      } else if (identical(active_tab, "RF Table")) {
        "Random Forest permutation importance table"
      } else {
        "Random Forest permutation importance plot"
      }

      legend_body <- if (identical(active_tab, "ROC / AUC")) {
        "ROC curves summarize one-vs-rest classification performance for each class. AUC values in the legend and macro AUC indicate overall discriminative performance."
      } else if (identical(active_tab, "SHAP Bar plot")) {
        "The SHAP plot shows the top taxa ranked by mean absolute SHAP value. Larger bars indicate stronger contribution of each taxon to model predictions."
      } else if (identical(active_tab, "SHAP Table")) {
        "This table summarizes SHAP statistics per feature, including class-specific contributions when available."
      } else if (identical(active_tab, "RF Table")) {
        "This table reports feature-level permutation importance and associated abundance summaries used to interpret model drivers."
      } else {
        "The RF bar plot ranks taxa by permutation importance. Higher values indicate larger performance drop when that feature is permuted, implying stronger predictive contribution."
      }

      tags$div(
        tags$div(
          style = "font-weight: 600; margin-bottom: 4px;",
          legend_title
        ),
        tags$div(
          paste0(
            legend_body,
            " Model setup: outcome = ",
            input$outcome_var,
            ", outcome type = ",
            model_mode,
            ", taxonomic level = ",
            input$tax_level,
            ", transform = ",
            input$transform_method,
            ", validation = ",
            validation_label,
            "."
          )
        )
      )
    })
  })
}
