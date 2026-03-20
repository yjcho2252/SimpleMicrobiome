### UI 함수
mod_beta_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        width = 3,
        
        uiOutput(ns("primary_group_selector")),
        uiOutput(ns("primary_color_controls")),
        
        uiOutput(ns("secondary_group_selector")),
        hr(),
        
        h4("Plot Settings"),
        numericInput(ns("plot_width_px"), "Plot Width (pixels):", value = 600, min = 300, max = 1500, step = 50),
        numericInput(ns("plot_height_px"), "Plot Height (pixels):", value = 450, min = 300, max = 1500, step = 50),
        hr(),
        
        checkboxInput(ns("show_ellipses"), "Show Group Ellipses", value = FALSE),
        checkboxInput(ns("show_sample_names"), "Show Sample Names", value = FALSE),
        hr(),
        
        numericInput(ns("dot_size"), "Dot Size (point size):", value = 4, min = 0.5, max = 10, step = 0.5),
        
        # --- Vector/Taxa Analysis 섹션 전체 제거 ---
        # hr(),
        # h4("Vector/Taxa Analysis"),
        # ...
        # ----------------------------------------------
        
        hr(),
        
        h4("Axis Limits"),
        helpText("Leave blank to use automatic limits."),
        numericInput(ns("x_min"), "X-axis minimum:", value = NA),
        numericInput(ns("x_max"), "X-axis maximum:", value = NA),
        numericInput(ns("y_min"), "Y-axis minimum:", value = NA),
        numericInput(ns("y_max"), "Y-axis maximum:", value = NA)
      ),
      
      mainPanel(
        width = 9,
        tabsetPanel(
          id = ns("beta_tabs"),
          tabPanel("PCoA (Bray-Curtis)",
                   h4("PCoA - Bray Curtis Ordination"),
                   downloadButton(ns("download_pcoa"), "Download PCoA Plot (PNG)"),
                   tags$div(style = "margin-bottom: 15px;"),
                   uiOutput(ns("pcoa_plot_ui"))
          ),
          tabPanel("NMDS (Bray-Curtis)",
                   h4("NMDS - Bray Curtis Ordination"),
                   actionButton(ns("redraw_nmds"), "Redraw Plot"),
                   downloadButton(ns("download_nmds"), "Download NMDS Plot (PNG)"),
                   tags$div(style = "margin-bottom: 15px;"),
                   uiOutput(ns("nmds_plot_ui"))
          )
        ),
        
        hr(),
        h4("Statistical Test: PERMANOVA Results"),
        div(
          class = "alert alert-info",
          htmlOutput(ns("permanova_results_out"))
        )
      )
    )
  )
}

### Server 함수
mod_beta_server <- function(id, ps_obj, meta_cols) {
  moduleServer(id, function(input, output, session) {
    
    output$primary_group_selector <- renderUI({
      req(meta_cols())
      group_choices <- setdiff(meta_cols(), "SampleID")
      selected_col <- if (length(group_choices) > 0) group_choices[1] else NULL
      selectInput(session$ns("primary_group_var"), "Select Primary Variable (Color):",
                  choices = group_choices, selected = selected_col)
    })
    
    output$secondary_group_selector <- renderUI({
      req(meta_cols())
      group_choices <- c("None", setdiff(meta_cols(), "SampleID"))
      selected_col <- if (length(group_choices) > 1) group_choices[2] else group_choices[1]
      selectInput(session$ns("secondary_group_var"), "Select Secondary Variable (Shape):",
                  choices = group_choices, selected = "None")
    })
    
    primary_group_var <- reactive({ req(input$primary_group_var); input$primary_group_var })
    secondary_group_var <- reactive({ input$secondary_group_var })
    
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
      
      do.call(tagList, c(list(hr(), h4("Primary Group Colors")), controls, list(hr())))
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
    
    data_prepared <- reactive({
      req(ps_obj())
      ps_data <- ps_obj()
      dist_mat <- phyloseq::distance(ps_data, method = "bray")
      metadata <- as(phyloseq::sample_data(ps_data), "data.frame")
      list(ps = ps_data, dist_mat = dist_mat, metadata = metadata)
    })
    
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
    
    # ⭐ 제거됨: group_vector_data reactive
    # ⭐ 제거됨: taxa_vector_data reactive
    
    pcoa_ordination_reactive <- reactive({
      req(ps_plot_obj())
      phyloseq::ordinate(ps_plot_obj(), method = "PCoA", distance = "bray")
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
        p <- phyloseq::plot_ordination(ps_data_for_plot, ord, color = primary_group_var()) +
          ggplot2::geom_point(size = input$dot_size)
      } else {
        p <- phyloseq::plot_ordination(ps_data_for_plot, ord, color = primary_group_var(), shape = plot_shape_var()) +
          ggplot2::geom_point(size = input$dot_size)
      }
      p <- p + ggplot2::scale_color_manual(values = primary_color_map())
      
      if (input$show_ellipses) {
        p <- p + ggplot2::stat_ellipse(ggplot2::aes_string(group = primary_group_var()))
      }
      
      p <- p + ggplot2::theme_bw() +
        ggplot2::labs(title = "PCoA - Bray-Curtis",
                      color = primary_group_var(),
                      shape = if (is.null(plot_shape_var())) NULL else plot_shape_var(),
                      x = x_axis_label,
                      y = y_axis_label)
      
      plot_data <- phyloseq::plot_ordination(ps_data_for_plot, ord, justDF = TRUE)
      
      if (input$show_sample_names) {
        p <- p + ggplot2::geom_text(data = plot_data, mapping = ggplot2::aes(x = Axis.1, y = Axis.2, label = SampleID), color = "black", size = 3, hjust = 0.2, vjust = 0.5)
      }
      
      # ⭐ 제거됨: Vector/Taxa plotting logic (input$show_vectors)
      
      limits <- axis_limits()
      if (!is.null(limits$x) || !is.null(limits$y)) {
        p <- p + ggplot2::coord_cartesian(xlim = limits$x, ylim = limits$y)
      }
      
      return(p)
    })
    
    output$pcoa_plot_out <- renderPlot({ pcoa_plot_reactive() })
    
    output$download_pcoa <- downloadHandler(
      filename = function() { paste0("PCoA_BrayCurtis_", Sys.Date(), ".png") },
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
      isolate(phyloseq::ordinate(ps_plot_obj(), method = "NMDS", distance = "bray"))
    }, ignoreNULL = FALSE)
    
    nmds_plot_reactive <- reactive({
      req(ps_plot_obj(), primary_group_var(), ord <- nmds_ordination_reactive(), input$dot_size)
      ps_data_for_plot <- ps_plot_obj()
      
      if (is.null(plot_shape_var())) {
        p <- phyloseq::plot_ordination(ps_data_for_plot, ord, color = primary_group_var()) +
          ggplot2::geom_point(size = input$dot_size)
      } else {
        p <- phyloseq::plot_ordination(ps_data_for_plot, ord, color = primary_group_var(), shape = plot_shape_var()) +
          ggplot2::geom_point(size = input$dot_size)
      }
      p <- p + ggplot2::scale_color_manual(values = primary_color_map())
      
      if (input$show_ellipses) {
        p <- p + ggplot2::stat_ellipse(ggplot2::aes_string(group = primary_group_var()))
      }
      
      p <- p + ggplot2::theme_bw() +
        ggplot2::labs(title = paste0("NMDS - Bray-Curtis (Stress: ", round(ord$stress, 3), ")"),
                      color = primary_group_var(),
                      shape = if (is.null(plot_shape_var())) NULL else plot_shape_var())
      
      plot_data <- phyloseq::plot_ordination(ps_data_for_plot, ord, justDF = TRUE)
      
      if (input$show_sample_names) {
        p <- p + ggplot2::geom_text(data = plot_data, mapping = ggplot2::aes(x = NMDS1, y = NMDS2, label = SampleID), color = "black", size = 3, hjust = -0.1, vjust = 0.5)
      }
      
      # ⭐ 제거됨: Vector/Taxa plotting logic (input$show_vectors)
      
      limits <- axis_limits()
      if (!is.null(limits$x) || !is.null(limits$y)) {
        p <- p + ggplot2::coord_cartesian(xlim = limits$x, ylim = limits$y)
      }
      
      return(p)
    })
    
    output$nmds_plot_out <- renderPlot({ nmds_plot_reactive() })
    
    output$download_nmds <- downloadHandler(
      filename = function() { paste0("NMDS_BrayCurtis_", Sys.Date(), ".png") },
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
        "<strong>Distance Metric:</strong> Bray-Curtis<br>",
        "---<br>",
        "<strong>R-squared ($R^2$):</strong> ", format(permanova_result$R2, digits = 3, nsmall = 3), "<br>",
        "<strong style='color: ", ifelse(permanova_result$Pval < 0.05, "red", "black"), "; ", ifelse(permanova_result$Pval < 0.05, "font-weight: bold;", ""), "'>P-value:</strong> ", format(permanova_result$Pval, digits = 3),
        "<br><br>",
        "<em>Interpretation: $R^2$ indicates the proportion of variance explained by the group. The P-value indicates if the group centroids are significantly different.</em>"
      ))
    })
    
  })
}
