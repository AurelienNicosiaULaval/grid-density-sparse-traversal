#!/usr/bin/env Rscript

script_root <- function() {
  args <- commandArgs(FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0L) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", hit[[1L]])), "../.."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

root <- script_root()
options(grid_density_project_root = root)
source(file.path(root, "code/R/00_helpers.R"))
source_project_R(root)
ensure_project_dirs(root)

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Package 'testthat' is required to run tests. Install it with install.packages('testthat').",
       call. = FALSE)
}

test_dir <- file.path(root, "code/tests/testthat")
if (!dir.exists(test_dir)) stop("Test directory not found: ", test_dir, call. = FALSE)

testthat::test_dir(test_dir, reporter = "summary")
