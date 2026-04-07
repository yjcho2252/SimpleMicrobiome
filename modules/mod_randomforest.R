library(shiny)
library(phyloseq)
library(dplyr)
library(DT)
library(ggplot2)

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
        numericInput(ns("prevalence_filter_pct"), "6. Feature prevalence cutoff (%)", value = 5, min = 0, max = 100, step = 1),
        numericInput(ns("top_n_features"), "7. Top N features by mean abundance", value = 100, min = 10, max = 5000, step = 10),
        sliderInput(ns("train_ratio"), "8. Train ratio", min = 0.6, max = 0.9, value = 0.8, step = 0.05),
        numericInput(ns("ntree"), "9. Number of trees (ntree)", value = 500, min = 100, max = 5000, step = 100),
        numericInput(ns("mtry"), "10. mtry (0 = auto)", value = 0, min = 0, max = 10000, step = 1),
        numericInput(ns("seed"), "11. Random seed", value = 1234, min = 1, max = 999999, step = 1),
        actionButton(ns("run_rf"), "Run Random Forest", class = "btn-danger", style = "font-size: 12px;")
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
            plotOutput(ns("importance_plot"), height = "500px")
          )
        ),
        hr(),
        h4(icon("square-poll-vertical"), "Result"),
        verbatimTextOutput(ns("rf_status")),
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

    observeEvent(list(ps_obj_filtered_raw(), input$outcome_var), {
      req(ps_obj_filtered_raw(), input$outcome_var)
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
      selected <- selected[selected %in% choices]
      updateSelectizeInput(session, "outcome_levels", choices = choices, selected = selected, server = TRUE)
    }, ignoreInit = FALSE)

    observeEvent(input$run_rf, {
      req(ps_obj_filtered_raw(), input$outcome_var)

      status_text("Running Random Forest...")
      model_result(NULL)

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

      rf_fit <- randomForest::randomForest(
        x = x_train,
        y = y_train,
        ntree = as.integer(input$ntree),
        mtry = mtry_val,
        importance = TRUE
      )

      pred <- stats::predict(rf_fit, newdata = x_test)

      if (is_classification) {
        cm <- table(Predicted = pred, Actual = y_test)
        accuracy <- mean(pred == y_test)

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
          mean_diff_df <- data.frame(
            Feature = names(diff_vec),
            MaxGroupMean = as.numeric(max_mean[names(diff_vec)]),
            MinGroupMean = as.numeric(min_mean[names(diff_vec)]),
            MeanDifference_MaxMinusMin = as.numeric(diff_vec),
            UnivariatePValue = as.numeric(vapply(names(diff_vec), function(feat) {
              suppressWarnings(tryCatch(
                stats::kruskal.test(x[, feat] ~ y_group)$p.value,
                error = function(e) NA_real_
              ))
            }, numeric(1))),
            stringsAsFactors = FALSE
          )
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

      model_result(list(
        metrics = metrics_text,
        importance = imp_tbl,
        score_col = "PermutationImportance"
      ))

      status_text("Model completed.")
    })

    output$rf_status <- renderText({
      status_text()
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

    output$importance_plot <- renderPlot({
      res <- model_result()
      req(res)
      top_imp <- head(res$importance, 20)
      ggplot(top_imp, aes(x = reorder(Feature, PermutationImportance), y = PermutationImportance)) +
        geom_col(fill = "#2c7fb8") +
        coord_flip() +
        labs(x = "Feature", y = "Permutation importance", title = "Top 20 Feature Importance") +
        theme_minimal(base_size = 12)
    })
  })
}
