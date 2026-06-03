#!/usr/bin/env Rscript

script_root <- function() {
  args <- commandArgs(FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0L) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", hit[[1L]])), "../.."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

parse_bool <- function(x, default = TRUE) {
  if (is.null(x) || !nzchar(x)) return(default)
  x <- tolower(trimws(x))
  if (x %in% c("1", "true", "yes", "y")) return(TRUE)
  if (x %in% c("0", "false", "no", "n")) return(FALSE)
  stop("Cannot parse boolean value: ", x, call. = FALSE)
}

root <- script_root()
set.seed(20260603)

run_full <- parse_bool(Sys.getenv("RUN_FULL", unset = "TRUE"), default = TRUE)
required_packages <- c("Rcpp", "RcppParallel", "ggplot2", "sessioninfo")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running the full reproducibility pipeline.",
    call. = FALSE
  )
}

compiler_candidates <- Sys.which(c("c++", "g++", "clang++"))
if (!any(nzchar(compiler_candidates))) {
  stop("No C++ compiler was found on PATH. A C++17 compiler is required.", call. = FALSE)
}

dirs <- file.path(root, c("code/results", "code/figures", "code/logs", "figures", "tables", "session-info"))
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

rscript <- file.path(R.home("bin"), "Rscript")
run_script <- function(script, args = character()) {
  script_path <- file.path(root, script)
  if (!file.exists(script_path)) stop("Missing script: ", script_path, call. = FALSE)
  message("\n==> Running ", script, if (length(args)) paste(" ", paste(args, collapse = " ")) else "")
  status <- system2(rscript, c(script_path, args))
  if (!identical(status, 0L)) {
    stop("Script failed with status ", status, ": ", script, call. = FALSE)
  }
  invisible(TRUE)
}

copy_figures_to_root <- function() {
  src_dir <- file.path(root, "code/figures")
  dst_dir <- file.path(root, "figures")
  dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(src_dir, all.files = FALSE, full.names = TRUE, no.. = TRUE)
  files <- files[file.info(files)$isdir == FALSE]
  if (length(files)) {
    ok <- file.copy(files, dst_dir, overwrite = TRUE)
    if (!all(ok)) stop("Failed to copy one or more figures to figures/.", call. = FALSE)
  }
  invisible(TRUE)
}

message("Project root: ", root)
message("RUN_FULL=", run_full)
if (!run_full) {
  message(
    "RUN_FULL=FALSE performs a smoke reproducibility run: correctness checks, ",
    "figure regeneration from existing CSV files, table traceability, figure sync, ",
    "and global session-info export. It does not regenerate timing benchmark CSV files."
  )
}

run_script("code/scripts/run_correctness_checks.R")

if (run_full) {
  run_script("code/scripts/run_family_benchmarks.R")
  run_script("code/scripts/run_scaling_benchmarks.R")
  run_script("code/scripts/run_occupancy_benchmarks.R")
  run_script("code/scripts/run_mode_selector_benchmarks.R")
}

run_script("code/scripts/make_estimator_family_figures.R")
run_script("code/scripts/make_scaling_figures.R")
run_script("code/scripts/check_table_traceability.R")

copy_figures_to_root()

source(file.path(root, "code/R/00_helpers.R"))
write_session_info(file.path(root, "session-info/session_info.txt"))

message("\nReproducibility pipeline completed successfully.")
message("CSV outputs: ", file.path(root, "code/results"))
message("Canonical figures: ", file.path(root, "code/figures"))
message("Reviewer figure copies: ", file.path(root, "figures"))
message("Traceability tables: ", file.path(root, "tables"))
message("Session information: ", file.path(root, "session-info/session_info.txt"))
