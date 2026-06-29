#!/usr/bin/env Rscript

# Install dependencies for local execution of SimpleMicrobiome.
# Keep this list aligned with app.R and modules/*.R.
# Usage (from project root):
#   source("install.R")

cran_repo <- "https://cloud.r-project.org"
options(repos = c(CRAN = cran_repo))
runtime_dependencies <- c("Depends", "Imports", "LinkingTo")

cran_packages <- c(
  "bslib",
  "cluster",
  "CVXR",
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
  "GraceYoon/SPRING@3d641a4b939b1b3cc042c064a05000aa48266af0",
  "stefpeschel/NetCoMi"
)

package_version_or_na <- function(pkg) {
  tryCatch(
    as.character(utils::packageVersion(pkg)),
    error = function(e) NA_character_
  )
}

is_installed <- function(pkg) {
  !is.na(package_version_or_na(pkg))
}

install_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, is_installed, logical(1))]
  if (length(missing) == 0) {
    message("CRAN packages: all installed.")
    return(invisible(character(0)))
  }

  message("CRAN packages to install: ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = runtime_dependencies, repos = cran_repo)
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

install_pinned_rcpparmadillo <- function() {
  # Match the tested Docker environment used by the public app.
  install_version <- "15.0.2-2"
  check_version <- "15.0.2.2"

  current_version <- package_version_or_na("RcppArmadillo")

  if (identical(current_version, check_version)) {
    message("RcppArmadillo already pinned at version ", check_version, ".")
    return(invisible(character(0)))
  }

  ensure_remotes()
  if (is_installed("RcppArmadillo")) {
    message("Removing existing RcppArmadillo version: ", current_version)
    remove.packages("RcppArmadillo")
  }

  message("Installing pinned RcppArmadillo version: ", install_version)
  remotes::install_version(
    "RcppArmadillo",
    version = install_version,
    upgrade = "never",
    dependencies = c("Depends", "Imports", "LinkingTo"),
    repos = cran_repo
  )

  if (is_installed("RcppArmadillo") &&
      identical(package_version_or_na("RcppArmadillo"), check_version)) {
    character(0)
  } else {
    "RcppArmadillo"
  }
}

install_pinned_rbiom <- function() {
  # NetCoMi 1.2.0 expects the older rbiom API that exports unifrac().
  install_version <- "2.2.1"
  check_version <- "2.2.1"
  current_version <- package_version_or_na("rbiom")

  if (identical(current_version, check_version)) {
    message("rbiom already pinned at version ", check_version, ".")
    return(invisible(character(0)))
  }

  ensure_remotes()
  if (is_installed("rbiom")) {
    message("Removing existing rbiom version: ", current_version)
    remove.packages("rbiom")
  }

  message("Installing pinned rbiom version: ", install_version)
  remotes::install_version(
    "rbiom",
    version = install_version,
    upgrade = "never",
    dependencies = runtime_dependencies,
    repos = cran_repo
  )

  if (is_installed("rbiom") &&
      identical(package_version_or_na("rbiom"), check_version)) {
    character(0)
  } else {
    "rbiom"
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
  BiocManager::install(missing, ask = FALSE, update = FALSE, dependencies = runtime_dependencies)
  missing[!vapply(missing, is_installed, logical(1))]
}

install_ancombc <- function() {
  ensure_bioc_manager()

  if (!is_installed("CVXR")) {
    return("CVXR")
  }

  message("Installing ANCOMBC with the available CVXR version...")
  BiocManager::install("ANCOMBC", ask = FALSE, update = FALSE, force = TRUE)

  if (is_installed("ANCOMBC")) character(0) else "ANCOMBC"
}

install_github_if_missing <- function(spec) {
  repo_part <- sub("@.*$", "", spec)
  pkg_name <- sub(".*/", "", repo_part)
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
    dependencies = runtime_dependencies,
    repos = c(c(CRAN = cran_repo), BiocManager::repositories())
  )

  if (is_installed(pkg_name)) character(0) else pkg_name
}

failed <- character(0)
failed <- c(failed, install_pinned_rcpparmadillo())
failed <- c(failed, install_missing_cran(cran_packages))
failed <- c(failed, install_ancombc())
failed <- c(failed, install_missing_bioc(bioc_packages))
failed <- c(failed, install_pinned_rbiom())
if (!is_installed("bluster")) {
  ensure_bioc_manager()
  BiocManager::install("bluster", ask = FALSE, update = FALSE, type = "source")
}
if (!is_installed("bluster")) failed <- c(failed, "bluster")
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
