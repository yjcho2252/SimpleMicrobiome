## UI
mod_randomforest_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        width = 2,
        h4(icon("tree"), "Random Forest"),
        hr(),
        selectInput(ns("outcome_var"), "1. Outcome variable", choices = NULL),
        selectInput(
          ns("outcome_type"),
          "2. Outcome type",
          choices = c("Auto", "Classification", "Regression"),
          selected = "Auto"
        ),
        selectInput(
          ns("tax_level"),
          "3. Taxonomic level",
          choices = c("ASV", "Genus", "Species"),
          selected = "Genus"
        ),
        selectizeInput(
          ns("outcome_levels"),
          "4. Levels to include (optional)",
          choices = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Select levels to include",
            plugins = list("remove_button")
          )
        ),
        selectInput(
          ns("transform_method"),
          "5. Feature transform",
          choices = c("TSS", "CLR", "log"),
          selected = "TSS"
        ),
        numericInput(ns("prevalence_filter_pct"), "6. Feature prevalence cutoff (0-20%)", value = 5, min = 0, max = 20, step = 1),
        numericInput(ns("top_n_features"), "7. Top N features by mean abundance", value = 100, min = 10, max = 5000, step = 10),
        sliderInput(ns("train_ratio"), "8. Train ratio", min = 0.6, max = 0.9, value = 0.8, step = 0.05),
        tags$details(
          style = "margin-bottom: 10px;",
          tags$summary("Advanced Options"),
          br(),
          numericInput(ns("ntree"), "9. Number of trees (ntree)", value = 500, min = 100, max = 5000, step = 100),
          numericInput(ns("mtry"), "10. mtry (0 = auto)", value = 0, min = 0, max = 10000, step = 1),
          numericInput(ns("seed"), "11. Random seed", value = 1234, min = 1, max = 999999, step = 1)
        ),
        actionButton(ns("run_rf"), "Run Random Forest", class = "btn-danger", style = "font-size: 12px;"),        
        tags$script(HTML(
          "Shiny.addCustomMessageHandler('toggle-rf-run-btn', function(msg) {
             var btn = document.getElementById(msg.id);
             if (!btn) return;
             btn.disabled = !!msg.disabled;
             if (msg.label) btn.textContent = msg.label;
           });"
        ))
      ),
      mainPanel(
        width = 10,
        tabsetPanel(
          tabPanel(
            "Table",
            downloadButton(
              ns("download_rf_table"),
              "Download Table (TSV)",
              style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"
            ),
            DTOutput(ns("importance_table"))
          ),
          tabPanel(
            "Bar Plot",
            downloadButton(
              ns("download_rf_barplot"),
              "Download Plot (PNG)",
              style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"
            ),
            div(
              style = "display: flex; justify-content: center;",
              plotOutput(ns("importance_plot"), width = "1160px", height = "500px")
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
              style = "display: flex; justify-content: center;",
              plotOutput(ns("roc_plot"), width = "1160px", height = "500px")
            ),
            verbatimTextOutput(ns("roc_auc_text"))
          )
        ),
        hr(),
        h4(icon("square-poll-vertical"), "Result"),
        verbatimTextOutput(ns("rf_metrics"))
      )
    )
  )
}

## Server
mod_randomforest_server <- function(id, ps_obj_filtered_raw) {
  moduleServer(id, function(input, output, session) {
    model_result <- reactiveVal(NULL)
    status_text <- reactiveVal("Waiting for model run.")
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

    observe({
      req(ps_obj_filtered_raw())
      ps <- ps_obj_filtered_raw()
      md <- data.frame(phyloseq::sample_data(ps))
      vars <- colnames(md)
      current <- input$outcome_var
      if (is.null(current) || !current %in% vars) {
        current <- if (length(vars) > 0) vars[1] else NULL
      }
      updateSelectInput(session, "outcome_var", choices = vars, selected = current)
    })

    observeEvent(list(ps_obj_filtered_raw(), input$outcome_var, input$outcome_type), {
      req(ps_obj_filtered_raw(), input$outcome_var, input$outcome_type)
      ps <- ps_obj_filtered_raw()
      md <- data.frame(phyloseq::sample_data(ps))
      validate(need(input$outcome_var %in% colnames(md), "Selected outcome is not available."))

      v <- md[[input$outcome_var]]
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
    }, ignoreInit = FALSE)

    observeEvent(input$run_rf, {
      req(ps_obj_filtered_raw(), input$outcome_var)

      status_text("Running Random Forest...")
      model_result(NULL)
      session$sendCustomMessage("toggle-rf-run-btn", list(
        id = session$ns("run_rf"),
        disabled = TRUE,
        label = "Running..."
      ))
      on.exit({
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
          validate(need(input$outcome_var %in% colnames(md), "Outcome variable is not available."))

          y_raw <- md[[input$outcome_var]]
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
          y_raw <- md[[input$outcome_var]]

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
          if (is_classification) {
            y <- factor(as.character(y_raw))
            validate(need(length(levels(y)) >= 2, "Classification requires at least two classes."))

            split_by_class <- split(seq_len(length(y)), y)
            train_idx <- unlist(lapply(split_by_class, function(idx) {
              n_tr <- max(1, floor(length(idx) * input$train_ratio))
              sample(idx, size = n_tr)
            }), use.names = FALSE)
            train_idx <- sort(unique(train_idx))
            test_idx <- setdiff(seq_len(nrow(x)), train_idx)
            validate(need(length(test_idx) > 0, "Test set is empty. Lower train ratio or add samples."))
          } else {
            y <- as.numeric(y_raw)
            validate(need(!all(is.na(y)), "Outcome has only missing values."))
            train_n <- floor(nrow(x) * input$train_ratio)
            train_idx <- sample(seq_len(nrow(x)), size = train_n)
            test_idx <- setdiff(seq_len(nrow(x)), train_idx)
            validate(need(length(test_idx) > 0, "Test set is empty. Lower train ratio or add samples."))
          }

          x_train <- x[train_idx, , drop = FALSE]
          y_train <- y[train_idx]
          x_test <- x[test_idx, , drop = FALSE]
          y_test <- y[test_idx]

          mtry_val <- as.integer(input$mtry)
          if (is.na(mtry_val) || mtry_val <= 0) {
            mtry_val <- floor(sqrt(ncol(x_train)))
          }
          mtry_val <- max(1, min(mtry_val, ncol(x_train)))

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

      roc_payload <- NULL
      if (is_classification) {
        cm <- table(Predicted = pred, Actual = y_test)
        accuracy <- mean(pred == y_test)
        pred_prob <- stats::predict(rf_fit, newdata = x_test, type = "prob")

        compute_binary_roc <- function(y_true_bin, y_score) {
          keep <- !is.na(y_true_bin) & !is.na(y_score)
          y_true_bin <- y_true_bin[keep]
          y_score <- y_score[keep]
          if (length(y_true_bin) < 2) {
            return(NULL)
          }
          n_pos <- sum(y_true_bin == 1)
          n_neg <- sum(y_true_bin == 0)
          if (n_pos == 0 || n_neg == 0) {
            return(NULL)
          }

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

          list(
            curve = data.frame(FPR = fpr_u, TPR = tpr_u),
            auc = as.numeric(auc)
          )
        }

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
          auc_values <- vapply(roc_list, function(x) x$auc, numeric(1))
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
          "Samples: ", nrow(x), " (Train: ", length(train_idx), ", Test: ", length(test_idx), ")\n",
          "Features used: ", ncol(x), "\n",
          "RMSE: ", round(rmse, 4), "\n",
          "MAE: ", round(mae, 4), "\n",
          "R-squared: ", round(r2, 4)
        )
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
              paste0(rk, "__Unassigned")
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

          model_result(list(
            metrics = metrics_text,
            importance = imp_tbl,
            score_col = "PermutationImportance",
            roc = roc_payload
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

    build_importance_plot <- function(res, taxa_fontsize = 8) {
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
        width = grid::unit(4.4 * ncol(hm_mat), "mm"),
        height = grid::unit(4.4 * nrow(hm_mat), "mm")
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
          width = grid::unit(120, "pt")
        )
      )

      hm + bar_anno
    }

    output$importance_plot <- renderPlot(
      {
        res <- model_result()
        req(res)
        p <- build_importance_plot(res, taxa_fontsize = 8)
        ComplexHeatmap::draw(
          p,
          merge_legend = TRUE,
          padding = grid::unit(c(6, 2, 6, 130), "mm")
        )
      },
      width = 1160,
      height = 500
    )

    output$download_rf_barplot <- downloadHandler(
      filename = function() {
        paste0("random_forest_barplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
      },
      content = function(file) {
        res <- model_result()
        req(res)
        grDevices::png(filename = file, width = 3600, height = 1800, res = 300)
        p <- build_importance_plot(res, taxa_fontsize = 7)
        ComplexHeatmap::draw(
          p,
          merge_legend = TRUE,
          padding = grid::unit(c(6, 2, 6, 70), "mm")
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
        res <- model_result()
        req(res)
        if (is.null(res$roc) || length(res$roc$curves) == 0) {
          graphics::plot.new()
          graphics::text(0.5, 0.5, "ROC curve is available for classification results only.")
          return(invisible(NULL))
        }
        p <- build_roc_plot(res, base_size = 13)
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
          graphics::text(0.5, 0.5, "ROC curve is available for classification results only.")
          grDevices::dev.off()
          return(invisible(NULL))
        }
        p <- build_roc_plot(res, base_size = 18)
        print(p)
        grDevices::dev.off()
      }
    )
  })
}
