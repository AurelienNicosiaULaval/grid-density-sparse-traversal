#!/usr/bin/env Rscript

# Compatibility wrapper for the estimator-family benchmark name used in the
# repository documentation. The implementation lives in
# run_estimator_family_benchmarks.R.

script_root <- function() {
  args <- commandArgs(FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0L) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", hit[[1L]])), "../.."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

root <- script_root()
target <- file.path(root, "code/scripts/run_estimator_family_benchmarks.R")
args <- commandArgs(trailingOnly = TRUE)
cmd <- file.path(R.home("bin"), "Rscript")
status <- system2(cmd, c(target, args))
if (!identical(status, 0L)) {
  stop("run_estimator_family_benchmarks.R failed with status ", status, call. = FALSE)
}
