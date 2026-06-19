#!/usr/bin/env Rscript

# Install dependencies for local execution of SimpleMicrobiome.
# Keep this list aligned with app.R and modules/*.R.
# Usage (from project root):
#   source("install.R")

cran_repo <- "https://cloud.r-project.org"
options(repos = c(CRAN = cran_repo))

cran_packages <- c(
  "bslib",
  "cluster",
  "dplyr",
  "DT",
  "circlize",
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

install_pinned_cvxr <- function() {
  # 2026-06-19: ANCOMBC currently expects an older CVXR API that still exports solve().
  target_version <- "1.0-13"
  current_version <- if (is_installed("CVXR")) as.character(utils::packageVersion("CVXR")) else NA_character_
  if (identical(current_version, target_version)) {
    message("CVXR already pinned at version ", target_version, ".")
    return(invisible(character(0)))
  }

  ensure_remotes()
  if (is_installed("CVXR")) {
    message("Removing existing CVXR version: ", current_version)
    remove.packages("CVXR")
  }

  message("Installing pinned CVXR version: ", target_version)
  remotes::install_version("CVXR", version = target_version, upgrade = "never", dependencies = TRUE, repos = cran_repo)

  if (is_installed("CVXR") && identical(as.character(utils::packageVersion("CVXR")), target_version)) {
    character(0)
  } else {
    "CVXR"
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

install_ancombc <- function() {
  ensure_bioc_manager()

  if (!is_installed("CVXR")) {
    return("CVXR")
  }

  message("Reinstalling ANCOMBC against pinned CVXR...")
  BiocManager::install("ANCOMBC", ask = FALSE, update = FALSE, force = TRUE)

  if (is_installed("ANCOMBC")) character(0) else "ANCOMBC"
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
failed <- c(failed, install_pinned_cvxr())
failed <- c(failed, install_ancombc())
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
