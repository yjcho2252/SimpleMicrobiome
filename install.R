#!/usr/bin/env Rscript

# Restore the R package environment used by SimpleMicrobiome.
#
# The package versions and remote GitHub commits are defined in renv.lock.
# Run this script from the project root with either:
#   Rscript install.R
#   source("install.R")
#
# Docker sets SIMPLEMICROBIOME_R_LIBRARY explicitly. For an interactive local
# installation, the first writable default R library is used.

options(timeout = max(1200, getOption("timeout", 60)))

find_lockfile <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_dir <- if (length(script_arg)) {
    dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE))
  } else {
    character(0)
  }

  candidates <- unique(c(
    file.path(getwd(), "renv.lock"),
    if (length(script_dir)) file.path(script_dir, "renv.lock") else character(0)
  ))
  existing <- candidates[file.exists(candidates)]

  if (!length(existing)) {
    stop(
      "renv.lock was not found. Keep install.R and renv.lock together, ",
      "or run install.R from the project root.",
      call. = FALSE
    )
  }

  normalizePath(existing[[1]], mustWork = TRUE)
}

select_install_library <- function() {
  requested <- Sys.getenv("SIMPLEMICROBIOME_R_LIBRARY", unset = "")

  if (nzchar(requested)) {
    requested <- path.expand(requested)
    dir.create(requested, recursive = TRUE, showWarnings = FALSE)
    if (file.access(requested, 2) != 0) {
      stop("R library is not writable: ", requested, call. = FALSE)
    }
    .libPaths(unique(c(requested, .libPaths())))
    return(normalizePath(requested, mustWork = TRUE))
  }

  libraries <- .libPaths()
  writable <- libraries[dir.exists(libraries) & file.access(libraries, 2) == 0]

  if (!length(writable)) {
    user_library <- path.expand(Sys.getenv("R_LIBS_USER", unset = ""))
    if (!nzchar(user_library)) {
      stop(
        "No writable R library was found. Set R_LIBS_USER or ",
        "SIMPLEMICROBIOME_R_LIBRARY to a writable directory.",
        call. = FALSE
      )
    }
    dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
    if (file.access(user_library, 2) != 0) {
      stop("R user library is not writable: ", user_library, call. = FALSE)
    }
    .libPaths(unique(c(user_library, .libPaths())))
    writable <- user_library
  }

  normalizePath(writable[[1]], mustWork = TRUE)
}

lockfile <- find_lockfile()
project_dir <- dirname(lockfile)
install_library <- select_install_library()

message("Using lockfile: ", lockfile)
message("Installing into R library: ", install_library)

# Bootstrap renv only when it is unavailable. renv then restores its recorded
# version, along with all other packages, from renv.lock.
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages(
    "renv",
    lib = install_library,
    repos = c(CRAN = "https://cloud.r-project.org"),
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

if (!requireNamespace("renv", quietly = TRUE)) {
  stop("Failed to install the renv package.", call. = FALSE)
}

lock <- renv::lockfile_read(lockfile)
locked_r_version <- lock$R$Version
current_r_version <- as.character(getRversion())
restore_repositories <- c(CRAN = "https://cloud.r-project.org")

# rocker images may define a dated Posit Package Manager snapshot in
# getOption("repos"). The lockfile was captured against CRAN, and some locked
# versions postdate that snapshot, so force the recorded CRAN repository for
# every restore operation.
options(repos = restore_repositories)

if (!identical(current_r_version, locked_r_version)) {
  stop(
    "This lockfile requires R ", locked_r_version,
    ", but the current R version is ", current_r_version, ".",
    call. = FALSE
  )
}

# ggtext is load-tested immediately after installation by renv. Restore its
# compiled gridtext dependency first so parallel restore ordering cannot leave
# ggtext temporarily unloadable.
renv::restore(
  project = project_dir,
  lockfile = lockfile,
  library = install_library,
  packages = "gridtext",
  repos = restore_repositories,
  clean = FALSE,
  prompt = FALSE
)

renv::restore(
  project = project_dir,
  lockfile = lockfile,
  library = install_library,
  repos = restore_repositories,
  exclude = c("SPRING", "NetCoMi"),
  clean = FALSE,
  prompt = FALSE
)

# NetCoMi 1.2.0 declares SPRING as an unpinned GitHub remote. Resolving that
# declaration currently selects SPRING HEAD, which requires R >= 4.6. Install
# both packages directly from the immutable commit archives recorded in
# renv.lock. Using install.packages() on the downloaded archive is intentional:
# renv::install() re-resolves the GitHub package's dependencies and can replace
# already-restored lockfile versions (for example, Rcpp) in the same batch.
pinned_github_record <- function(package) {
  record <- lock$Packages[[package]]
  required_fields <- c("RemoteUsername", "RemoteRepo", "RemoteSha")
  missing_fields <- required_fields[
    !vapply(required_fields, function(field) {
      value <- record[[field]]
      !is.null(value) && length(value) == 1 && !is.na(value) && nzchar(value)
    }, logical(1))
  ]

  if (length(missing_fields)) {
    stop(
      "The lockfile does not contain a complete pinned GitHub record for ",
      package, ": ", paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  record
}

for (package in c("SPRING", "NetCoMi")) {
  record <- pinned_github_record(package)
  specification <- paste0(
    record$RemoteUsername, "/", record$RemoteRepo, "@", record$RemoteSha
  )
  archive_url <- paste0(
    "https://api.github.com/repos/",
    record$RemoteUsername, "/", record$RemoteRepo,
    "/tarball/", record$RemoteSha
  )
  archive <- tempfile(
    pattern = paste0(package, "-", substr(record$RemoteSha, 1, 8), "-"),
    fileext = ".tar.gz"
  )

  message("Installing pinned GitHub package: ", specification)
  utils::download.file(
    url = archive_url,
    destfile = archive,
    method = "libcurl",
    mode = "wb",
    quiet = FALSE
  )
  utils::install.packages(
    pkgs = archive,
    lib = install_library,
    repos = NULL,
    type = "source",
    dependencies = FALSE,
    INSTALL_opts = "--no-multiarch"
  )
  unlink(archive)

  installed_version <- tryCatch(
    as.character(utils::packageVersion(package, lib.loc = install_library)),
    error = function(error) NA_character_
  )
  expected_version <- lock$Packages[[package]]$Version
  if (!identical(installed_version, expected_version)) {
    stop(
      package, " version verification failed: expected ", expected_version,
      ", installed ", installed_version,
      call. = FALSE
    )
  }
  tryCatch(
    loadNamespace(package),
    error = function(error) {
      stop(
        "Failed to load pinned ", package, ": ", conditionMessage(error),
        call. = FALSE
      )
    }
  )
}

required_packages <- c(
  "ANCOMBC",
  "CVXR",
  "ggpattern",
  "ggtext",
  "gridtext",
  "NetCoMi",
  "phyloseq",
  "shiny",
  "SPRING",
  "SpiecEasi"
)

for (package in required_packages) {
  message("Verifying R package: ", package)
  tryCatch(
    loadNamespace(package),
    error = function(error) {
      stop(
        "Failed to load ", package, ": ", conditionMessage(error),
        call. = FALSE
      )
    }
  )
}

message("")
message("All packages recorded in renv.lock were restored successfully.")
message("Next step: shiny::runApp(\"app.R\")")
