mod_maaslin3_ui <- function(id) {
  ns <- NS(id)
  plot_scroll_layout(
    sidebarPanel(
      width = 2,
      h4(icon("flask"), "MaAsLin3"),
      hr(),
      selectInput(ns("primary"), "1. Primary Explanatory Variable", choices = NULL),
      uiOutput(ns("primary_type_ui")),
      checkboxInput(ns("use_subgroup"), "Select Subgroup of Primary Explanatory Variable", value = FALSE),
      uiOutput(ns("subgroup_ui")),
      uiOutput(ns("comparison_ui")),
      uiOutput(ns("reference_ui")),
      selectInput(ns("rank"), "4. Taxonomic Level", choices = NULL),
      numericInput(ns("q"), "5. Significance Threshold (q)", 0.1, min = 0.001, max = 1, step = 0.01),
      numericInput(ns("prevalence"), "6. Prevalence Filter (%)", 5, min = 0, max = 100, step = 1),
      tags$details(
        style = "margin-bottom: 12px;",
        tags$summary(strong("Advanced Options")),
        br(),
        selectizeInput(ns("covariates"), "7. Additional Covariates (Optional)", choices = NULL,
                       multiple = TRUE, options = list(plugins = list("remove_button"))),
        selectizeInput(ns("interactions"), "8. Interaction Terms (Optional)", choices = NULL,
                       multiple = TRUE, options = list(plugins = list("remove_button"))),
        selectizeInput(ns("subject"), "9. Random Intercepts (Optional)", choices = NULL,
                       multiple = TRUE, options = list(plugins = list("remove_button"))),
        checkboxInput(ns("depth"), "10. Adjust for Input Read Depth", TRUE),
        selectizeInput(ns("group"), "11. Group Tests (Optional)", choices = NULL,
                       multiple = TRUE, options = list(plugins = list("remove_button"))),
        selectInput(ns("ordered"), "12. Ordered Contrasts (Optional)", choices = c("None" = "")),
        uiOutput(ns("ordered_levels_ui")),
        selectInput(ns("strata"), "13. Matched Strata (Optional)", choices = c("None" = "")),
        helpText("Matched strata cannot be combined with random intercepts."),
        verbatimTextOutput(ns("model"))
      ),
      hr(),
      h4(icon("up-right-and-down-left-from-center"), "Display Options"),
      numericInput(ns("summary_width"), "Plot Width (px)", 700, min = 400, max = 4000, step = 50),
      numericInput(ns("summary_height"), "Plot Height (px)", 1000, min = 400, max = 4000, step = 50),
      numericInput(ns("summary_font"), "Base Font Size (pt)", 10, min = 6, max = 30, step = 1),
      actionButton(ns("run"), "Run MaAsLin3", class = "btn-danger maaslin3-run-btn",
                   style = "font-size: 12px;"),
      tags$script(HTML("
        function setMaaslin3RunButton(button, running) {
          button.disabled = running;
          button.textContent = running ? 'Running...' : 'Run MaAsLin3';
          if (running) {
            const spinner = document.createElement('span');
            spinner.className = 'spinner-border spinner-border-sm me-1';
            spinner.setAttribute('aria-hidden', 'true');
            button.prepend(spinner);
          }
        }
        Shiny.addCustomMessageHandler('toggle-maaslin3-run-btn', function(message) {
          const button = document.getElementById(message.id);
          if (button) setMaaslin3RunButton(button, message.running);
        });
        $(document).on('click.maaslin3', '.maaslin3-run-btn', function() {
          setMaaslin3RunButton(this, true);
        });
      "))
    ),
    mainPanel(
      width = 10,
      h4("MaAsLin3"),
      analysis_run_view(ns, "run",
      uiOutput(ns("results_width_style")),
      div(
        id = ns("results_content"),
        style = "width: var(--results-width, 700px); max-width: 100%;",
      tabsetPanel(
          id = ns("result_tab"),
          tabPanel("Summary Plot", uiOutput(ns("summary"))),
          tabPanel("Volcano Plot",
            checkboxInput(ns("volcano_primary_only"), "Show only the primary variable", TRUE),
            downloadButton(
              ns("volcano_download"),
              "Download volcano PNG",
              style = "width: 200px; height: 34px; font-size: 11px; display: flex; align-items: center; justify-content: center; margin: 10px 0 12px 0;"
            ),
            div(style = "overflow-x: auto;", plotOutput(ns("volcano"), height = "auto"))
          ),
          tabPanel(
            "Results Table",
            div(
              style = "margin-top: 12px; margin-bottom: 18px;",
              p("Abundance effects are relative fold changes among detected samples. Prevalence effects are detection odds ratios. Individual q-values apply to each model. Joint q-values combine evidence across both models."),
              fluidRow(
                column(6,
                  selectInput(ns("view"), "Model", c("Abundance", "Prevalence", "All results"), width = "100%"),
                  checkboxInput(ns("primary_only"), "Show only the primary variable", TRUE)
                ),
                column(6,
                  selectInput(ns("test"), "Significance test", c("Individual", "Joint"), width = "100%"),
                  checkboxInput(ns("significant"), "Significant results only", FALSE)
                )
              ),
              downloadButton(ns("tsv"), "Download full TSV"),
              p("Select a result row to inspect its abundance and prevalence."),
              DT::DTOutput(ns("results")),
              uiOutput(ns("association_detail"))
            )
          )
        ),
        div(
          style = "clear: both; margin-top: 16px;",
          hr(),
          h5(icon("circle-info"), "MaAsLin3 Status"),
          div(class = "app-result-card app-status-card", verbatimTextOutput(ns("status")))
        )
      )
      )
    )
  )
}

mod_maaslin3_server <- function(id, ps_obj, plot_order = function() NULL) {
  moduleServer(id, function(input, output, session) {
    output$run_prompt <- analysis_run_prompt("Run MaAsLin3")
    result <- reactiveVal(NULL)
    status <- reactiveVal("")
    running <- reactiveVal(FALSE)
    next_run_at <- -Inf

    draw_wait_message <- function(message_text) {
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      graphics::par(mar = c(0, 0, 0, 0))
      graphics::plot.new()
      graphics::text(0.5, 0.5, message_text, cex = 0.85)
    }

    metadata <- reactive({
      req(ps_obj())
      as(phyloseq::sample_data(ps_obj()), "data.frame")
    })

    observeEvent(ps_obj(), {
      taxonomy <- phyloseq::tax_table(ps_obj(), errorIfNULL = FALSE)
      ranks <- if (!is.null(taxonomy)) colnames(taxonomy)
      meta <- metadata()
      variables <- setdiff(names(meta), "SampleID")
      rank_choices <- c("Feature level" = "__feature_level__", stats::setNames(ranks, ranks))
      current_rank <- isolate(input$rank)
      selected_rank <- if (isTruthy(current_rank) && current_rank %in% unname(rank_choices)) {
        current_rank
      } else if ("Genus" %in% ranks) {
        "Genus"
      } else {
        "__feature_level__"
      }
      selected_primary <- if (length(variables) > 0L) variables[1L] else ""
      updateSelectInput(session, "rank",
                        choices = rank_choices,
                        selected = selected_rank)
      updateSelectInput(session, "primary", choices = variables, selected = selected_primary)
      updateSelectizeInput(session, "covariates", choices = variables, server = TRUE)
      updateSelectizeInput(session, "subject", choices = variables, server = TRUE)
      updateSelectInput(session, "strata", choices = c("None" = "", variables))
    })

    numeric_primary <- reactive({
      isTruthy(input$primary) && !is.null(ps_obj()) &&
        is.numeric(metadata()[[input$primary]])
    })

    output$primary_type_ui <- renderUI({
      req(numeric_primary())
      selectInput(session$ns("primary_type"), "Treat as",
                  choices = c("Continuous" = "continuous", "Categorical" = "categorical"))
    })

    categorical_primary <- reactive({
      !numeric_primary() || identical(input$primary_type, "categorical")
    })

    output$subgroup_ui <- renderUI({
      req(input$primary, input$primary %in% names(metadata()))
      if (!isTRUE(input$use_subgroup)) return(NULL)
      if (!categorical_primary()) {
        return(helpText("Treat the primary explanatory variable as categorical to select a subgroup."))
      }

      primary_values <- as.character(metadata()[[input$primary]])
      primary_levels <- levels(factor(primary_values[!is.na(primary_values) & nzchar(primary_values)]))
      categorical_variables <- setdiff(
        names(metadata())[!vapply(metadata(), is.numeric, logical(1))],
        c("SampleID", input$primary)
      )
      tagList(
        selectInput(session$ns("primary_level"), "Primary Explanatory Level to Include",
                    choices = primary_levels,
                    selected = if (isTruthy(isolate(input$primary_level)) &&
                                  isolate(input$primary_level) %in% primary_levels) {
                      isolate(input$primary_level)
                    } else {
                      head(primary_levels, 1L)
                    }),
        selectInput(session$ns("secondary"), "Secondary Grouping Variable",
                    choices = categorical_variables,
                    selected = if (isTruthy(isolate(input$secondary)) &&
                                  isolate(input$secondary) %in% categorical_variables) {
                      isolate(input$secondary)
                    } else {
                      head(categorical_variables, 1L)
                    })
      )
    })

    subgroup_active <- reactive({
      isTRUE(input$use_subgroup) && categorical_primary()
    })

    analysis_primary <- reactive({
      if (subgroup_active() && isTruthy(input$secondary)) input$secondary else input$primary
    })

    analysis_numeric <- reactive({
      variable <- analysis_primary()
      isTruthy(variable) && variable %in% names(metadata()) && is.numeric(metadata()[[variable]])
    })

    analysis_categorical <- reactive({
      if (subgroup_active()) TRUE else categorical_primary()
    })

    observe({
      req(ps_obj())
      excluded <- if (subgroup_active()) input$primary else character()
      choices <- setdiff(setdiff(names(metadata()), "SampleID"), excluded)
      updateSelectizeInput(session, "covariates", choices = choices,
                           selected = intersect(isolate(input$covariates), choices), server = TRUE)
      updateSelectizeInput(session, "subject", choices = choices,
                           selected = intersect(isolate(input$subject), choices), server = TRUE)
      current_strata <- isolate(input$strata)
      updateSelectInput(session, "strata", choices = c("None" = "", choices),
                        selected = if (isTruthy(current_strata) && current_strata %in% choices) current_strata else "")
    })

    special_effect_choices <- reactive({
      meta <- metadata()
      fixed <- intersect(names(meta), unique(c(analysis_primary(), input$covariates)))
      categorical <- !vapply(meta[fixed], is.numeric, logical(1))
      fixed[categorical | (fixed %in% analysis_primary() & analysis_categorical())]
    })

    observe({
      req(ps_obj())
      choices <- special_effect_choices()
      group <- intersect(isolate(input$group), choices)
      ordered <- intersect(isolate(input$ordered), choices)
      updateSelectizeInput(session, "group", choices = choices, selected = group, server = TRUE)
      updateSelectInput(session, "ordered", choices = c("None" = "", choices),
                        selected = if (length(ordered)) ordered else "")
    })

    primary_levels <- reactive({
      primary <- analysis_primary()
      req(primary, primary %in% names(metadata()))
      values_raw <- metadata()[[primary]]
      if (analysis_categorical()) {
        values <- as.character(values_raw)
        levels(factor(values[!is.na(values) & nzchar(values)]))
      } else {
        character()
      }
    })

    output$comparison_ui <- renderUI({
      if (is.null(ps_obj()) || !isTruthy(input$primary)) {
        return(selectizeInput(
          session$ns("comparison_groups"),
          "2. Explanatory Variable Levels to Compare",
          choices = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Select two or more levels to compare",
            plugins = list("remove_button")
          )
        ))
      }
      available <- primary_levels()
      if (length(available)) {
        tagList(
          selectizeInput(session$ns("comparison_groups"), "2. Explanatory Variable Levels to Compare",
                         choices = available, selected = available, multiple = TRUE,
                         options = list(plugins = list("remove_button"))),
          if (length(available) > 10L) {
            p(class = "text-warning-emphasis",
              sprintf(
                "This variable has %d observed levels. Review sample counts and keep the groups relevant to your analysis.",
                length(available)
              ))
          }
        )
      }
    })

    output$ordered_levels_ui <- renderUI({
      req(input$ordered, input$ordered %in% special_effect_choices())
      available <- levels(factor(metadata()[[input$ordered]]))
      if (identical(input$ordered, analysis_primary())) {
        available <- intersect(available, input$comparison_groups)
      }
      selectizeInput(session$ns("ordered_levels"), "Select Levels from Lowest to Highest",
                     choices = available, selected = intersect(isolate(input$ordered_levels), available),
                     multiple = TRUE,
                     options = list(plugins = list("remove_button")))
    })

    interaction_pairs <- reactive({
      variables <- names(metadata())
      fixed <- intersect(variables, unique(c(analysis_primary(), input$covariates)))
      if (length(fixed) < 2L) return(list())
      pairs <- utils::combn(fixed, 2, simplify = FALSE)
      keys <- vapply(pairs, function(pair) paste(match(pair, variables), collapse = ":"), character(1))
      setNames(pairs, keys)
    })

    observe({
      pairs <- if (is.null(ps_obj())) list() else interaction_pairs()
      labels <- vapply(pairs, function(pair) paste(pair, collapse = " × "), character(1))
      updateSelectizeInput(session, "interactions",
                           choices = setNames(as.character(names(pairs)), labels),
                           selected = intersect(isolate(input$interactions), names(pairs)),
                           server = TRUE)
    })

    output$reference_ui <- renderUI({
      if (is.null(ps_obj()) || !isTruthy(input$primary)) {
        return(selectInput(session$ns("reference"), "3. Reference Level", choices = NULL))
      }
      primary <- analysis_primary()
      req(primary, primary %in% names(metadata()))
      if (analysis_categorical() && !primary %in% c(input$group, input$ordered)) {
        choices <- intersect(primary_levels(), input$comparison_groups)
        previous <- isolate(input$reference)
        selectInput(session$ns("reference"), "3. Reference Level", choices = choices,
                    selected = if (isTruthy(previous) && previous %in% choices) previous else head(choices, 1L))
      } else {
        selectInput(session$ns("reference"), "3. Reference Level", choices = NULL)
      }
    })

    settings <- reactive(list(
      ps = ps_obj(), rank = input$rank, primary = analysis_primary(),
      use_subgroup = subgroup_active(),
      secondary = if (subgroup_active()) input$secondary else "",
      filter_variable = if (subgroup_active()) input$primary else "",
      filter_level = if (subgroup_active()) input$primary_level else "",
      numeric_as_category = analysis_numeric() && analysis_categorical(),
      comparison_groups = if (isTruthy(analysis_primary()) && !is.null(ps_obj())) {
        intersect(primary_levels(), input$comparison_groups)
      },
      reference = input$reference, covariates = input$covariates,
      interactions = if (!is.null(ps_obj())) {
        pairs <- interaction_pairs()
        unname(pairs[intersect(input$interactions, names(pairs))])
      },
      subject = input$subject, depth = input$depth, q = input$q,
      group = input$group, ordered = input$ordered, ordered_levels = input$ordered_levels,
      strata = input$strata,
      prevalence = if (is.null(input$prevalence)) 5 else input$prevalence
    ))

    prepare_model <- function(config) {
      ps <- config$ps
      if (is.null(ps) || !isTruthy(config$primary)) stop("Load data and select a primary variable.")
      if (isTRUE(config$use_subgroup) && !isTruthy(config$secondary)) {
        stop("Select a secondary grouping variable for the subgroup analysis.")
      }
      meta <- as(phyloseq::sample_data(ps), "data.frame")
      if (config$numeric_as_category) {
        meta[[config$primary]] <- as.character(meta[[config$primary]])
      }
      fixed <- unique(c(config$primary, config$covariates))
      subject <- unique(config$subject[nzchar(config$subject)])
      group <- unique(config$group)
      ordered <- if (isTruthy(config$ordered)) config$ordered else character()
      strata <- if (isTruthy(config$strata)) config$strata else character()
      special <- c(group, ordered)
      variables <- unique(c(fixed, subject, strata))
      if (!all(variables %in% names(meta))) stop("Select metadata variables from the current dataset.")
      if (any(subject %in% fixed)) stop("Choose a random intercept variable that is not a fixed effect.")
      if (length(strata) && length(subject)) stop("Matched strata and random intercepts cannot be used together.")
      if (any(strata %in% fixed)) stop("Choose a matched strata variable that is not a fixed effect.")
      if (anyDuplicated(special) || !all(special %in% fixed) ||
          any(vapply(meta[special], is.numeric, logical(1)))) {
        stop("Choose distinct categorical fixed effects for group tests and ordered contrasts.")
      }
      if (any(unlist(config$interactions, use.names = FALSE) %in% special)) {
        stop("Remove interaction terms involving variables used for group tests or ordered contrasts.")
      }

      ids <- phyloseq::sample_names(ps)
      depth <- phyloseq::sample_sums(ps)[ids]
      in_subgroup <- rep(TRUE, length(ids))
      if (isTruthy(config$filter_variable)) {
        if (!config$filter_variable %in% names(meta) || !isTruthy(config$filter_level)) {
          stop("Select a valid primary explanatory variable subgroup.")
        }
        in_subgroup <- as.character(meta[ids, config$filter_variable]) == config$filter_level
        in_subgroup[is.na(in_subgroup)] <- FALSE
      }
      meta <- meta[ids, variables, drop = FALSE]
      in_comparison <- rep(TRUE, nrow(meta))
      if (!is.numeric(meta[[config$primary]])) {
        if (length(config$comparison_groups) < 2L) stop("Select at least two explanatory variable levels to compare.")
        in_comparison <- meta[[config$primary]] %in% config$comparison_groups
      }
      keep <- in_subgroup & in_comparison & stats::complete.cases(meta) & is.finite(depth) & depth > 0
      if (sum(keep) < 3L) stop("At least three samples with complete metadata and positive read counts are required.")
      meta <- meta[keep, , drop = FALSE]
      for (name in variables) {
        if (!is.numeric(meta[[name]]) || name %in% c(subject, strata)) {
          meta[[name]] <- factor(meta[[name]])
        }
        if (length(unique(meta[[name]])) < 2L) {
          stop(sprintf("'%s' has fewer than two values after sample filtering.", name))
        }
      }
      for (name in subject) {
        if (nlevels(meta[[name]]) >= nrow(meta)) {
          stop(sprintf("Random intercept '%s' requires repeated observations.", name))
        }
      }
      if (length(strata) && !anyDuplicated(meta[[strata]])) {
        stop("Matched strata require samples that share a stratum.")
      }
      if (length(strata)) {
        within_strata <- split(meta[[config$primary]], meta[[strata]])
        if (!any(vapply(within_strata, function(x) length(unique(x)) > 1L, logical(1)))) {
          stop("The primary variable must vary within at least one matched stratum.")
        }
      }
      if (length(ordered)) {
        present <- levels(meta[[ordered]])
        level_order <- config$ordered_levels
        if (anyDuplicated(level_order) || !all(present %in% level_order)) {
          stop("Select every observed level of the ordered variable from lowest to highest.")
        }
        meta[[ordered]] <- factor(meta[[ordered]], levels = intersect(level_order, present))
      }
      primary <- meta[[config$primary]]
      use_reference <- is.factor(primary) && !config$primary %in% special
      if (use_reference) {
        if (!isTruthy(config$reference) || !config$reference %in% levels(primary)) {
          stop("Select a reference level represented in the analysis samples.")
        }
        meta[[config$primary]] <- stats::relevel(primary, ref = config$reference)
      }

      safe_names <- setNames(make.names(variables, unique = TRUE), variables)
      names(meta) <- unname(safe_names)
      fixed_terms <- unname(safe_names[setdiff(fixed, special)])
      adjusted_depth <- isTRUE(config$depth) && length(unique(depth[keep])) > 1L
      if (adjusted_depth) {
        depth_name <- tail(make.unique(c(names(meta), "Input_read_depth")), 1)
        meta[[depth_name]] <- unname(depth[keep])
        fixed_terms <- c(fixed_terms, depth_name)
      }
      interaction_terms <- vapply(config$interactions, function(pair) {
        paste(safe_names[pair], collapse = ":")
      }, character(1))
      terms <- c(fixed_terms, interaction_terms,
                 if (length(group)) sprintf("group(%s)", safe_names[group]),
                 if (length(ordered)) sprintf("ordered(%s)", safe_names[ordered]),
                 if (length(strata)) sprintf("strata(%s)", safe_names[strata]),
                 if (length(subject)) sprintf("(1 | %s)", safe_names[subject]))
      list(metadata = meta, formula = paste("~", paste(terms, collapse = " + ")),
           primary = unname(safe_names[config$primary]), names = safe_names,
           special = unname(safe_names[special]),
           group = unname(safe_names[group]), ordered = unname(safe_names[ordered]),
           fixed_effects = setdiff(fixed, special), random_effects = subject, strata = strata,
           reference = if (use_reference) config$reference else "",
           filter_variable = config$filter_variable, filter_level = config$filter_level,
           excluded = sum(!keep), adjusted_depth = adjusted_depth)
    }

    model <- reactive(prepare_model(settings()))
    output$model <- renderText({
      req(ps_obj(), input$primary)
      tryCatch({
        spec <- model()
        config <- settings()
        format_terms <- function(terms) {
          if (length(terms)) paste(terms, collapse = " + ") else "NULL"
        }
        lines <- c(
          if (isTruthy(spec$filter_variable)) {
            paste0("Primary explanatory variable filter: ", spec$filter_variable, " = ", spec$filter_level)
          },
          paste(if (isTruthy(spec$filter_variable)) "Secondary grouping variable:" else "Primary explanatory variable:",
                config$primary),
          paste("Applied fixed_effects:", format_terms(spec$fixed_effects)),
          if (length(config$comparison_groups)) {
            paste("Explanatory variable levels to compare:", paste(config$comparison_groups, collapse = ", "))
          },
          if (length(config$interactions)) {
            interactions <- vapply(config$interactions, paste, character(1), collapse = " × ")
            paste("Selected interactions:", paste(interactions, collapse = ", "))
          },
          paste("Applied random_effects:", format_terms(spec$random_effects)),
          if (length(config$group)) paste("Group tests:", format_terms(config$group)),
          if (length(spec$ordered)) {
            sprintf("Ordered contrasts: %s (%s)", config$ordered,
                    paste(levels(spec$metadata[[spec$ordered]]), collapse = " → "))
          },
          if (length(spec$strata)) paste("Matched strata:", spec$strata),
          paste("Read depth adjustment:", if (spec$adjusted_depth) "on" else "off"),
          sprintf("Samples: %d (%d excluded)", nrow(spec$metadata), spec$excluded)
        )
        paste(lines, collapse = "\n")
      }, error = conditionMessage)
    })

    observeEvent(input$run, {
      if (isTRUE(running())) {
        showNotification("MaAsLin3 is already running.", type = "message", duration = 3)
        return(invisible(NULL))
      }
      update_run_button <- function(active) {
        session$sendCustomMessage(
          "toggle-maaslin3-run-btn", list(id = session$ns("run"), running = active)
        )
      }
      on.exit({
        running(FALSE)
        update_run_button(FALSE)
      }, add = TRUE)
      now <- as.numeric(Sys.time())
      if (now < next_run_at) {
        showNotification(sprintf("Please wait %d second(s) before running again.",
                                 ceiling(next_run_at - now)), type = "message", duration = 3)
        return(invisible(NULL))
      }
      tryCatch({
        config <- settings()
        if (!isTruthy(config$rank)) stop("Select a taxonomic level.")
        spec <- model()
        if (length(config$q) != 1L || !is.finite(config$q) || config$q <= 0 || config$q > 1) {
          stop("Enter a q-value threshold greater than 0 and at most 1.")
        }
        if (length(config$prevalence) != 1L || !is.finite(config$prevalence) ||
            config$prevalence < 0 || config$prevalence > 100) {
          stop("Enter a prevalence filter between 0% and 100%.")
        }
        if (!requireNamespace("maaslin3", quietly = TRUE)) {
          stop("The maaslin3 package is unavailable. Restore the project dependencies and restart the app.")
        }
        next_run_at <<- as.numeric(Sys.time()) + 10
        running(TRUE)
        update_run_button(TRUE)
        result(NULL)
        ps <- phyloseq::prune_samples(rownames(spec$metadata), config$ps)
        if (!identical(config$rank, "__feature_level__")) {
          ps <- phyloseq::tax_glom(ps, taxrank = config$rank, NArm = FALSE)
        }
        counts <- as(phyloseq::otu_table(ps), "matrix")
        if (phyloseq::taxa_are_rows(ps)) counts <- t(counts)
        counts <- counts[rownames(spec$metadata), colSums(counts) > 0, drop = FALSE]
        if (!ncol(counts) || any(rowSums(counts) <= 0)) stop("No usable count table at this rank.")
        if (!identical(config$rank, "__feature_level__")) {
          labels <- as(phyloseq::tax_table(ps), "matrix")[colnames(counts), config$rank]
          missing <- is.na(labels) | !nzchar(labels)
          labels[missing] <- colnames(counts)[missing]
          colnames(counts) <- make.unique(labels)
        }

        out <- tempfile("maaslin3-")
        on.exit(unlink(out, recursive = TRUE), add = TRUE)
        withProgress(message = "Running MaAsLin3", value = 0.5, {
          maaslin3::maaslin3(
            input_data = as.data.frame(counts, check.names = FALSE),
            input_metadata = spec$metadata, output = out, formula = spec$formula,
            normalization = "TSS", transform = "LOG", standardize = TRUE,
            median_comparison_abundance = TRUE, median_comparison_prevalence = FALSE,
            min_prevalence = config$prevalence / 100,
            max_significance = config$q, plot_summary_plot = TRUE,
            plot_associations = FALSE, save_plots_rds = TRUE, cores = 4
          )
        })
        tsv <- file.path(out, "all_results.tsv")
        summary_rds <- file.path(out, "figures", "summary_plot_gg.RDS")
        read_bytes <- function(path) readBin(path, "raw", n = file.size(path))
        read_features <- function(filename) {
          utils::read.delim(file.path(out, "features", filename), row.names = 1,
                           check.names = FALSE, comment.char = "")
        }
        result(list(
          settings = config, model = spec,
          table = utils::read.delim(tsv, check.names = FALSE, comment.char = ""),
          filtered = read_features("filtered_data.tsv"),
          transformed = read_features("data_transformed.tsv"),
          tsv = read_bytes(tsv),
          summary_plot = if (file.exists(summary_rds)) readRDS(summary_rds)
        ))
        status("Analysis complete.")
      }, error = function(e) {
        message <- paste("MaAsLin3:", conditionMessage(e))
        status(message)
        showNotification(message, type = "error")
      })
    })

    current <- reactive({
      run <- result()
      req(run, identical(run$settings, settings()))
      run
    })

    output$status <- renderText({
      req(ps_obj())
      run <- result()
      message <- status()
      if (startsWith(message, "MaAsLin3:")) return(message)
      if (!is.null(run) && !identical(run$settings, settings())) {
        return("Data or model settings changed. Run MaAsLin3 again to see updated results.")
      }
      if (is.null(run)) {
        return(if (nzchar(message)) message else "Ready to run MaAsLin3.")
      }
      rank <- run$settings$rank
      if (identical(rank, "__feature_level__")) rank <- "Feature level"
      lines <- c(
        message,
        if (isTruthy(run$model$filter_variable)) {
          paste0("Primary explanatory variable filter: ", run$model$filter_variable,
                 " = ", run$model$filter_level)
        },
        paste(if (isTruthy(run$model$filter_variable)) "Secondary grouping variable:" else "Primary explanatory variable:",
              run$settings$primary),
        if (length(run$settings$comparison_groups)) {
          paste("Explanatory variable levels to compare:",
                paste(run$settings$comparison_groups, collapse = ", "))
        },
        paste("Samples used:", nrow(run$model$metadata)),
        paste("Taxonomic level:", rank),
        sprintf("Prevalence filter: %g%%", run$settings$prevalence),
        "",
        switch(
          if (is.null(input$result_tab)) "Summary Plot" else input$result_tab,
          "Volcano Plot" = sprintf(
            "Points show model effects relative to the tested null value. The dashed line marks the individual q-value cutoff (%g). Red and blue mark positive and negative effects at or below that cutoff.",
            run$settings$q
          ),
          "Results Table" = "Abundance effects describe fold changes among samples where the taxon was detected. Prevalence effects describe odds ratios for detection. Joint q-values combine evidence from both models.",
          "The summary plot shows associations from the abundance and prevalence models."
        )
      )
      paste(lines, collapse = "\n")
    })

    styled_summary <- reactive({
      plot <- current()$summary_plot
      req(plot)
      width <- input$summary_width
      height <- input$summary_height
      font <- input$summary_font
      valid <- function(value, lower, upper) {
        length(value) == 1L && is.finite(value) && value >= lower && value <= upper
      }
      validate(
        need(valid(width, 400, 4000) && valid(height, 400, 4000) &&
               width == floor(width) && height == floor(height),
             "Enter whole-number dimensions between 400 and 4000 px."),
        need(valid(font, 6, 30), "Enter a font size between 6 and 30 pt.")
      )
      font_theme <- ggplot2::theme(
        text = ggplot2::element_text(size = font),
        axis.title = ggplot2::element_text(size = font + 2),
        axis.text = ggplot2::element_text(size = font, colour = "black"),
        axis.text.x = ggplot2::element_text(size = font, colour = "black"),
        legend.title = ggplot2::element_text(size = font + 2),
        legend.text = ggplot2::element_text(size = font),
        strip.text = ggplot2::element_text(size = font)
      )
      style_panel <- function(panel) {
        fill_scale <- panel$scales$get_scales("fill")
        if (!is.null(fill_scale) && identical(fill_scale$name, "Beta coefficient")) {
          panel$scales <- panel$scales$clone()
          fill_scale <- panel$scales$get_scales("fill")
          fill_scale$labels <- function(breaks) {
            # MaAsLin3 assigns boundary values to the interval starting at that boundary.
            labels <- sub("^\\((?!-Inf,)(?=[^,]+,)", "[", breaks, perl = TRUE)
            sub("\\]$", ")", labels)
          }
        }
        # The heatmap hides its y-axis labels when a coefficient panel is present.
        if (!inherits(panel$theme$axis.text.y, "element_blank")) {
          panel <- panel + ggplot2::theme(axis.text.y = ggplot2::element_text(size = font, colour = "black"))
        }
        panel + font_theme
      }
      if (inherits(plot, "patchwork")) {
        for (i in seq_len(length(plot))) {
          plot[[i]] <- style_panel(plot[[i]])
        }
      } else {
        plot <- style_panel(plot)
      }
      list(plot = plot, width = width, height = height)
    })

    display_dpi <- 96
    export_dpi <- 300
    render_summary_png <- function(dpi) {
      summary <- styled_summary()
      file <- tempfile(fileext = ".png")
      on.exit(unlink(file), add = TRUE)
      ggplot2::ggsave(file, plot = summary$plot, device = ragg::agg_png, units = "in",
                      width = summary$width / display_dpi, height = summary$height / display_dpi,
                      dpi = dpi, bg = "white")
      readBin(file, "raw", n = file.size(file))
    }
    summary_png <- reactive({
      render_summary_png(display_dpi)
    })
    output$summary <- renderUI({
      if (is.null(result())) {
        return(plotOutput(session$ns("summary_wait"), height = "200px"))
      }
      run <- current()
      if (is.null(run$summary_plot)) {
        return(p(sprintf("No summary plot was produced at q ≤ %g. Check the results table and fit status.",
                         run$settings$q)))
      }
      tagList(
        downloadButton(session$ns("summary_download"), "Download summary PNG"),
        downloadButton(session$ns("summary_pdf"), "Download summary PDF"),
        div(style = "overflow-x: auto;",
            imageOutput(session$ns("summary_image"), height = "auto"))
      )
    })
    output$summary_wait <- renderPlot({
      if (isTRUE(running())) {
        draw_wait_message("MaAsLin3 is running. Please wait...")
      } else {
        draw_wait_message("Click 'Run MaAsLin3' to start analysis.")
      }
    }, height = 200)
    output$summary_image <- renderImage({
      png <- summary_png()
      file <- tempfile(fileext = ".png")
      writeBin(png, file)
      list(src = file, contentType = "image/png", alt = "MaAsLin3 summary plot")
    }, deleteFile = TRUE)
    output$summary_download <- downloadHandler("maaslin3-summary.png", function(file) {
      writeBin(render_summary_png(export_dpi), file)
    }, contentType = "image/png")
    output$summary_pdf <- downloadHandler("maaslin3-summary.pdf", function(file) {
      summary <- styled_summary()
      ggplot2::ggsave(file, plot = summary$plot, device = "pdf", units = "in",
                      width = summary$width / display_dpi, height = summary$height / display_dpi,
                      bg = "white")
    }, contentType = "application/pdf")
    output$tsv <- downloadHandler("maaslin3-all-results.tsv", function(file) {
      writeBin(current()$tsv, file)
    }, contentType = "text/tab-separated-values")

    volcano_plot <- reactive({
      run <- current()
      tab <- run$table
      if (isTRUE(input$volcano_primary_only)) {
        tab <- tab[tab$metadata %in% run$model$primary, , drop = FALSE]
      }
      models <- c(abundance = "Abundance", prevalence = "Prevalence")
      tab$model_label <- factor(unname(models[tolower(tab$model)]), levels = unname(models))
      tab$effect <- tab$coef - tab$null_hypothesis
      q <- tab$qval_individual
      valid <- !is.na(tab$model_label) & !tab$metadata %in% run$model$group &
        is.finite(tab$effect) & is.finite(q) & q >= 0 & q <= 1
      tab <- tab[valid, , drop = FALSE]
      validate(need(nrow(tab) > 0L, "No results with finite effects and individual q-values to plot."))

      threshold <- run$settings$q
      q <- tab$qval_individual
      q_floor <- max(.Machine$double.xmin, min(c(q[q > 0], threshold)) / 10)
      tab$log_q <- -log10(ifelse(q == 0, q_floor, q))
      tab$significance <- factor(
        ifelse(q > threshold, "Above q cutoff",
               ifelse(tab$effect >= 0, "Positive effect", "Negative effect")),
        levels = c("Positive effect", "Negative effect", "Above q cutoff")
      )
      labels <- tab[which(q <= threshold), , drop = FALSE] %>%
        dplyr::group_by(model_label) %>%
        dplyr::arrange(qval_individual, dplyr::desc(abs(effect)), .by_group = TRUE) %>%
        dplyr::slice_head(n = 5) %>%
        dplyr::ungroup()
      ggplot2::ggplot(tab, ggplot2::aes(effect, log_q, colour = significance)) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey55", linetype = "dotted") +
        ggplot2::geom_hline(yintercept = -log10(threshold), colour = "grey55", linetype = "dashed") +
        ggplot2::geom_point(alpha = 0.75, size = 2) +
        ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = feature),
                                size = 3, seed = 42, show.legend = FALSE) +
        ggplot2::facet_wrap(~ model_label, scales = "free_x", drop = FALSE) +
        ggplot2::scale_colour_manual(values = c("Positive effect" = "#B54443",
                                               "Negative effect" = "#3568A7",
                                               "Above q cutoff" = "#808080")) +
        ggplot2::labs(
          x = "Model coefficient minus the tested null value",
          y = expression(-log[10]("individual q-value")),
          colour = NULL,
          caption = paste0(
            sprintf("Dashed line marks the run's q-value cutoff (%g).", threshold),
            if (any(q == 0)) " Zero q-values are shown above the smallest positive q-value or the cutoff."
          )
        ) +
        ggplot2::theme_bw(base_size = req(input$summary_font)) +
          plot_axis_label_theme(req(input$summary_font), 10) +
        ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"), legend.position = "bottom")
    })

    output$volcano <- renderPlot({
                                   if (is.null(result())) {
                                     if (isTRUE(running())) {
                                       draw_wait_message("MaAsLin3 is running. Please wait...")
                                     } else {
                                       draw_wait_message("Click 'Run MaAsLin3' to start analysis.")
                                     }
                                     return(invisible(NULL))
                                   }
                                   volcano_plot()
                                 },
                                 width = function() volcano_dimensions()$width,
                                 height = function() if (is.null(result())) 200 else volcano_dimensions()$height)
    output$volcano_download <- downloadHandler("maaslin3-volcano.png", function(file) {
      size <- volcano_dimensions()
      ggplot2::ggsave(file, plot = volcano_plot(), device = ragg::agg_png, units = "in",
                      width = size$width / display_dpi, height = size$height / display_dpi,
                      dpi = export_dpi, bg = "white")
    }, contentType = "image/png")

    visible_results <- reactive({
      run <- current()
      tab <- run$table
      if (isTRUE(input$primary_only)) tab <- tab[tab$metadata %in% run$model$primary, , drop = FALSE]
      abundance <- tolower(tab$model) == "abundance"
      tab$effect_size <- ifelse(abundance, 2^tab$coef, exp(tab$coef))
      tab$effect_unit <- ifelse(abundance, "Relative abundance fold change", "Detection odds ratio")
      tab$reference <- ifelse(tab$metadata %in% run$model$primary, run$model$reference, "")
      tab$model_variable <- tab$metadata
      tab$test_type <- ifelse(tab$metadata %in% run$model$group, "Group test",
                              ifelse(tab$metadata %in% run$model$ordered, "Ordered contrast", ""))
      labels <- names(run$model$names)[match(tab$metadata, unname(run$model$names))]
      tab$metadata[!is.na(labels)] <- labels[!is.na(labels)]
      if (!identical(input$view, "All results")) {
        tab <- tab[tolower(tab$model) %in% tolower(input$view), , drop = FALSE]
      }
      q_col <- if (identical(input$test, "Joint")) "qval_joint" else "qval_individual"
      if (isTRUE(input$significant)) tab <- tab[which(tab[[q_col]] <= run$settings$q), , drop = FALSE]
      tab[order(tab[[q_col]], na.last = TRUE), , drop = FALSE]
    })

    output$results_width_style <- renderUI({
      width <- input$summary_width
      req(length(width) == 1L, is.finite(width), width >= 400, width <= 4000)
      tags$style(sprintf("#%s { --results-width: %dpx; }", session$ns("results_content"), as.integer(width)))
    })

    output$results <- DT::renderDT({
      tab <- visible_results()
      test_column <- if (any(nzchar(tab$test_type))) "test_type"
      columns <- c("feature", "metadata", "value", "reference", test_column, "model", "coef", "effect_size",
                   "effect_unit", "null_hypothesis", "pval_individual", "qval_individual",
                   "qval_joint", "N", "N_not_zero", "error")
      table_data <- tab[, intersect(columns, names(tab)), drop = FALSE]
      widget <- DT::datatable(
        table_data, rownames = FALSE, selection = "single",
        class = "compact stripe hover cell-border", style = "bootstrap",
        options = list(scrollX = TRUE, pageLength = 10, ordering = TRUE, order = list())
      )
      numeric_columns <- names(table_data)[vapply(table_data, is.numeric, logical(1))]
      widget <- DT::formatRound(widget, columns = numeric_columns, digits = 2)
      DT::formatStyle(widget, columns = names(table_data), fontSize = "11px", padding = "4px")
    })

    selected_association <- reactive({
      run <- current()
      rows <- input$results_rows_selected
      tab <- visible_results()
      req(length(rows) == 1L, rows %in% seq_len(nrow(tab)))
      row <- tab[rows, , drop = FALSE]
      validate(need(!grepl(":", row$name, fixed = TRUE),
                    "Association Detail is available for main effects. Interaction estimates are shown in the results table."))
      validate(need(!row$model_variable %in% run$model$special,
                    "Group tests and ordered contrasts are summarized in the results table."))
      feature <- row$feature
      display_variable <- row$metadata
      variable <- row$model_variable
      ids <- rownames(run$model$metadata)
      validate(need(
        variable %in% names(run$model$metadata) &&
          feature %in% names(run$filtered) && feature %in% names(run$transformed) &&
          all(ids %in% rownames(run$filtered)) && all(ids %in% rownames(run$transformed)),
        "Sample-level values are unavailable for this association."
      ))
      # MaAsLin3 can represent non-detections as NA in the filtered output.
      filtered <- run$filtered[ids, feature]
      observations <- data.frame(
        predictor = run$model$metadata[[variable]],
        abundance = run$transformed[ids, feature],
        detected = !is.na(filtered) & filtered > 0
      )
      matches <- with(run$table, which(
        feature == row$feature & metadata == variable &
          (value == row$value | (is.na(value) & is.na(row$value)))
      ))
      list(rows = run$table[matches, , drop = FALSE], data = observations,
           feature = feature, variable = display_variable, contrast = row$value,
           reference = row$reference)
    })

    association_fit_label <- function(detail, model) {
      fit <- detail$rows[tolower(detail$rows$model) %in% model, , drop = FALSE]
      if (nrow(fit) != 1L || !is.finite(fit$coef)) return("Model estimate unavailable")
      abundance <- identical(model, "abundance")
      effect <- if (abundance) 2^fit$coef else exp(fit$coef)
      q <- if (is.finite(fit$qval_individual)) format(fit$qval_individual, digits = 3) else "unavailable"
      warning <- if (!is.na(fit$error) && nzchar(fit$error)) paste0(" | ", fit$error) else ""
      sprintf("Adjusted %s: %s | Individual q: %s%s",
              if (abundance) "fold change" else "odds ratio",
              format(effect, digits = 3), q, warning)
    }

    volcano_dimensions <- reactive({
      width <- input$summary_width
      height <- input$summary_height
      req(isTruthy(width), isTruthy(height), is.finite(width), is.finite(height),
          width >= 400, width <= 4000, height >= 400, height <= 4000)
      list(width = as.integer(width), height = as.integer(height))
    })

    detail_theme <- function() {
      req(isTruthy(input$summary_font), is.finite(input$summary_font))
      ggplot2::theme_minimal(base_size = input$summary_font) +
        plot_axis_label_theme(input$summary_font, 10) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(face = "bold", colour = "#243746"))
    }

    output$association_detail <- renderUI({
      detail <- selected_association()
      heading <- paste(detail$feature, "—", detail$variable)
      if (!is.na(detail$contrast) && nzchar(detail$contrast) && is.factor(detail$data$predictor)) {
        reference <- if (nzchar(detail$reference)) detail$reference else levels(detail$data$predictor)[1]
        heading <- paste0(heading, " (", detail$contrast, " vs ", reference, ")")
      }
      div(
        style = "max-width: 100%; margin-top: 16px; padding: 16px; border: 1px solid #d9e2ec; border-radius: 8px;",
        h4("Association Detail"), h5(heading),
        p("Plots show observed samples. Effect sizes and individual q-values come from the adjusted model. Continuous predictor effects are per standard deviation."),
        div(
          style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 300px), 1fr)); gap: 16px;",
          div(style = "min-width: 0;", h5("Abundance"), p(association_fit_label(detail, "abundance")),
              plotOutput(session$ns("detail_abundance"), height = "300px")),
          div(style = "min-width: 0;", h5("Prevalence"), p(association_fit_label(detail, "prevalence")),
              plotOutput(session$ns("detail_prevalence"), height = "300px"))
        )
      )
    })

    output$detail_abundance <- renderPlot({
      detail <- selected_association()
      data <- detail$data[which(detail$data$detected), , drop = FALSE]
      validate(need(nrow(data) > 0L, "No detected samples."))
      plot <- ggplot2::ggplot(data, ggplot2::aes(predictor, abundance))
      if (is.factor(data$predictor)) {
        plot <- plot +
          ggplot2::geom_boxplot(ggplot2::aes(fill = predictor), width = 0.48, alpha = 0.6,
                                colour = "#264653", outlier.shape = NA, box.linewidth = 0.45, show.legend = FALSE) +
          ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.14, height = 0, seed = 42),
                              colour = "#264653", alpha = 0.55, size = 1.8) +
          ggplot2::scale_fill_viridis_d(begin = 0.25, end = 0.85) +
          ggplot2::scale_x_discrete(
            limits = function(values) levels(plot_order_factor(values, detail$variable, plot_order(), values)),
            labels = scales::label_wrap(16)
          )
      } else {
        plot <- plot + ggplot2::geom_point(colour = "#287a8b", alpha = 0.7, size = 2)
      }
      plot + ggplot2::labs(x = detail$variable, y = "Log-transformed relative abundance") +
        detail_theme()
    })

    output$detail_prevalence <- renderPlot({
      detail <- selected_association()
      data <- detail$data
      if (is.factor(data$predictor)) {
        grouped <- split(data$detected, data$predictor, drop = TRUE)
        summary <- data.frame(
          group = factor(names(grouped), levels = levels(data$predictor)),
          detected = vapply(grouped, sum, numeric(1)),
          total = lengths(grouped)
        )
        summary$fraction <- summary$detected / summary$total
        summary$label <- paste(summary$detected, summary$total, sep = " / ")
        plot <- ggplot2::ggplot(summary, ggplot2::aes(group, fraction)) +
          ggplot2::geom_col(ggplot2::aes(fill = group), width = 0.62,
                            alpha = 0.75, show.legend = FALSE) +
          ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.45, colour = "#243746",
                             fontface = "bold", size = input$summary_font / ggplot2::.pt) +
          ggplot2::scale_fill_viridis_d(begin = 0.25, end = 0.85) +
          ggplot2::scale_x_discrete(
            limits = function(values) levels(plot_order_factor(values, detail$variable, plot_order(), values)),
            labels = scales::label_wrap(16)
          ) +
          ggplot2::scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.25),
                                      limits = c(0, 1.08)) +
          ggplot2::labs(x = detail$variable, y = "Samples with detection")
      } else {
        data$presence <- factor(data$detected, levels = c(FALSE, TRUE), labels = c("Not detected", "Detected"))
        plot <- ggplot2::ggplot(data, ggplot2::aes(presence, predictor)) +
          ggplot2::geom_boxplot(ggplot2::aes(fill = presence), width = 0.48, alpha = 0.55,
                                colour = "#264653", outlier.shape = NA, box.linewidth = 0.45, show.legend = FALSE) +
          ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.14, height = 0, seed = 42),
                              colour = "#264653", alpha = 0.55, size = 1.8) +
          ggplot2::scale_fill_manual(values = c("#a7b2bc", "#64b6ac")) +
          ggplot2::labs(x = "Detection status", y = detail$variable)
      }
      plot + detail_theme()
    })
  })
}
