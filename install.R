#!/usr/bin/env Rscript

# Install dependencies for local execution of SimpleMicrobiome.
# Usage (from project root):
#   source("install.R")

message("==> SimpleMicrobiome dependency installer")

cran_repo <- "https://cloud.r-project.org"
options(repos = c(CRAN = cran_repo))

cran_packages <- c(
  "bslib",
  "cluster",
  "dplyr",
  "DT",
  "ggpattern",
  "ggplot2",
  "ggpubr",
  "ggraph",
  "ggrepel",
  "httr2",
  "igraph",
  "randomForest",
  "RColorBrewer",
  "readr",
  "scales",
  "shapr",
  "shiny",
  "shinyWidgets",
  "tibble",
  "tidyr",
  "tidyselect",
  "tidyverse",
  "vegan"
)

bioc_packages <- c(
  "ANCOMBC",
  "ALDEx2",
  "ComplexHeatmap",
  "Maaslin2",
  "microbiome",
  "nloptr",
  "phyloseq"
)

github_packages <- c(
  "zdk123/SpiecEasi",
  "GraceYoon/SPRING",
  "stefpeschel/NetCoMi"
)

is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

install_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, is_installed, logical(1))]
  if (length(missing) == 0) {
    message("CRAN packages: all installed.")
    return(invisible(character(0)))
  }

  message("CRAN packages to install: ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE, repos = cran_repo)
  missing[!vapply(missing, is_installed, logical(1))]
}

ensure_bioc_manager <- function() {
  if (!is_installed("BiocManager")) {
    message("Installing BiocManager...")
    install.packages("BiocManager", repos = cran_repo)
  }
}

ensure_remotes <- function() {
  if (!is_installed("remotes")) {
    message("Installing remotes...")
    install.packages("remotes", repos = cran_repo)
  }
}

install_missing_bioc <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, is_installed, logical(1))]
  if (length(missing) == 0) {
    message("Bioconductor packages: all installed.")
    return(invisible(character(0)))
  }

  ensure_bioc_manager()
  message("Bioconductor packages to install: ", paste(missing, collapse = ", "))
  BiocManager::install(missing, ask = FALSE, update = FALSE)
  missing[!vapply(missing, is_installed, logical(1))]
}

install_github_if_missing <- function(spec) {
  pkg_name <- sub(".*/", "", spec)
  if (is_installed(pkg_name)) {
    message("GitHub package already installed: ", pkg_name)
    return(invisible(character(0)))
  }

  ensure_bioc_manager()
  ensure_remotes()
  message("Installing GitHub package: ", spec)

  remotes::install_github(
    spec,
    upgrade = "never",
    dependencies = TRUE,
    repos = c(c(CRAN = cran_repo), BiocManager::repositories())
  )

  if (is_installed(pkg_name)) character(0) else pkg_name
}

failed <- character(0)
failed <- c(failed, install_missing_cran(cran_packages))
failed <- c(failed, install_missing_bioc(bioc_packages))
failed <- c(failed, unlist(lapply(github_packages, install_github_if_missing), use.names = FALSE))
failed <- unique(failed)

if (length(failed) > 0) {
  message("")
  message("Some packages are still missing:")
  message("  - ", paste(failed, collapse = "\n  - "))
  message("")
  message("Please check system prerequisites (compilers, external libraries) and retry.")
} else {
  message("")
  message("All required packages are installed.")
  message("Next step: shiny::runApp(\"app.R\")")
}
