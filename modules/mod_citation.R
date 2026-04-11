## UI
mod_citation_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      h4(icon("book"), "Citation"),
      p("References are grouped by module (APA style)."),
      tags$style(HTML("
        .citation-section-title {
          margin-top: 16px;
          margin-bottom: 6px;
          font-weight: 700;
        }
        .citation-list {
          margin-top: 4px;
          padding-left: 1.25rem;
        }
        .citation-list li {
          margin-bottom: 10px;
          line-height: 1.45;
        }
      ")),
      uiOutput(ns("citation_list"))
    )
  )
}

## Server
mod_citation_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    citations_by_module <- list(
      "Data Loading / Preprocessing" = c(
        "McMurdie, P. J., & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLOS ONE, 8(4), e61217. https://doi.org/10.1371/journal.pone.0061217"
      ),
      "Taxa Barplot" = c(
        "Wickham, H. (2016). ggplot2: Elegant graphics for data analysis. Springer. https://ggplot2.tidyverse.org"
      ),
      "Taxa Comparison" = c(
        "Wickham, H. (2016). ggplot2: Elegant graphics for data analysis. Springer. https://ggplot2.tidyverse.org"
      ),
      "Alpha Diversity" = c(
        "Oksanen, J., Simpson, G. L., Blanchet, F. G., Kindt, R., Legendre, P., Minchin, P. R., O'Hara, R. B., Solymos, P., Stevens, M. H. H., Szoecs, E., & Wagner, H. (2024). vegan: Community Ecology Package (R package). https://CRAN.R-project.org/package=vegan",
        "Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., & Hornik, K. (2025). cluster: Cluster Analysis Basics and Extensions (R package). https://CRAN.R-project.org/package=cluster"
      ),
      "Beta Diversity" = c(
        "Oksanen, J., Simpson, G. L., Blanchet, F. G., Kindt, R., Legendre, P., Minchin, P. R., O'Hara, R. B., Solymos, P., Stevens, M. H. H., Szoecs, E., & Wagner, H. (2024). vegan: Community Ecology Package (R package). https://CRAN.R-project.org/package=vegan",
        "Wickham, H. (2016). ggplot2: Elegant graphics for data analysis. Springer. https://ggplot2.tidyverse.org"
      ),
      "ANCOM-BC2" = c(
        "Lin, H., & Peddada, S. D. (2020). Analysis of compositions of microbiomes with bias correction. Nature Communications, 11, Article 3514. https://doi.org/10.1038/s41467-020-17041-7"
      ),
      "MaAsLin2" = c(
        "Mallick, H., Rahnavard, A., McIver, L. J., Ma, S., Zhang, Y., Tickle, T. L., Weingart, G., Ren, B., Schwager, E. H., Thompson, K. N., Lu, Y., Waldron, L., Huttenhower, C., et al. (2021). Multivariable association discovery in population-scale meta-omics studies. PLOS Computational Biology, 17(11), e1009442. https://doi.org/10.1371/journal.pcbi.1009442"
      ),
      "Random Forest" = c(
        "Breiman, L. (2001). Random forests. Machine Learning, 45(1), 5-32. https://doi.org/10.1023/A:1010933404324",
        "Moe, Y. M., Jullum, M., & Loland, A. (2021). Explaining predictive models with dependent features using Shapley values and conditional inference. Artificial Intelligence, 298, 103502. https://doi.org/10.1016/j.artint.2021.103502",
        "Moe, Y. M., Jullum, M., & Loland, A. (2024). shapr: Explaining machine learning models with dependence-aware Shapley values (R package). https://CRAN.R-project.org/package=shapr",
        "Gu, Z., Eils, R., & Schlesner, M. (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics, 32(18), 2847-2849. https://doi.org/10.1093/bioinformatics/btw313"
      ),
      "SpiecEasi" = c(
        "Kurtz, Z. D., Muller, C. L., Miraldi, E. R., Littman, D. R., Blaser, M. J., & Bonneau, R. A. (2015). Sparse and compositionally robust inference of microbial ecological networks. PLOS Computational Biology, 11(5), e1004226. https://doi.org/10.1371/journal.pcbi.1004226",
        "Peschel, S., Muller, C. L., von Mutius, E., Boulesteix, A.-L., Depner, M., & NetCoMi Development Team. (2021). NetCoMi: Network construction and comparison for microbiome data in R. Briefings in Bioinformatics, 22(4), bbaa290. https://doi.org/10.1093/bib/bbaa290"
      ),
      "SparCC" = c(
        "Friedman, J., & Alm, E. J. (2012). Inferring correlation networks from genomic survey data. PLOS Computational Biology, 8(9), e1002687. https://doi.org/10.1371/journal.pcbi.1002687",
        "Peschel, S., Muller, C. L., von Mutius, E., Boulesteix, A.-L., Depner, M., & NetCoMi Development Team. (2021). NetCoMi: Network construction and comparison for microbiome data in R. Briefings in Bioinformatics, 22(4), bbaa290. https://doi.org/10.1093/bib/bbaa290"
      ),
      "General Visualization / Reporting" = c(
        "Wickham, H. (2016). ggplot2: Elegant graphics for data analysis. Springer. https://ggplot2.tidyverse.org"
      )
    )

    output$citation_list <- renderUI({
      tags$div(
        lapply(names(citations_by_module), function(module_name) {
          tags$div(
            tags$div(module_name, class = "citation-section-title"),
            tags$ol(
              class = "citation-list",
              lapply(citations_by_module[[module_name]], tags$li)
            )
          )
        })
      )
    })
  })
}
