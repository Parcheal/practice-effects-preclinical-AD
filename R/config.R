# Shared package configuration.
# All analysis scripts in this handoff should resolve paths through this file.
#
# DATA ROOT: set the environment variable PE_DATA_ROOT to the folder that holds
# the cohort exports you obtained from A4, PREVENT-AD and ADNI (layout in
# README.md). If PE_DATA_ROOT is unset, <repo root>/data is used.
#   e.g. in R:     Sys.setenv(PE_DATA_ROOT = "/path/to/your/pe_data")
#   e.g. in shell: export PE_DATA_ROOT=/path/to/your/pe_data

`%||%` <- function(x, y) if (is.null(x)) y else x

find_package_root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(8)) {
    if (file.exists(file.path(cur, "MANIFEST.csv")) &&
        file.exists(file.path(cur, "R", "config.R"))) {
      return(normalizePath(cur, winslash = "/", mustWork = TRUE))
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("Cannot locate PE_analysis_code_package root.")
}

pkg_root <- find_package_root()

data_root <- Sys.getenv("PE_DATA_ROOT", unset = "")
if (!nzchar(data_root)) data_root <- file.path(pkg_root, "data")
data_root <- normalizePath(data_root, winslash = "/", mustWork = FALSE)

path_pkg <- function(...) file.path(pkg_root, ...)
path_data <- function(...) file.path(data_root, ...)
path_locked <- function(...) path_pkg("outputs", "locked", ...)
path_verify <- function(...) {
  p <- path_pkg("verify", ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

set.seed(42)
options(readr.show_col_types = FALSE)
