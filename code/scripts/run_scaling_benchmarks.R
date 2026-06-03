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
source(file.path(root, "code/R/00_helpers.R"))
source_project_R(root)
ensure_project_dirs(root)

quick <- parse_flag("--quick", default = TRUE)
out_dir <- file.path(root, "code/results")

if (quick) {
  ns <- c(500L, 1000L, 2000L, 5000L, 10000L)
  ds <- 2L:12L
  reps <- 2L
  scenario <- "normal"
  bin_width_multiplier <- 1.0
  direct_cpp_budget <- 1e8
  short_timing_inner_reps <- 10L
} else {
  ns <- c(500L, 1000L, 2000L, 5000L, 10000L, 20000L)
  ds <- 2L:15L
  reps <- 3L
  scenario <- "normal"
  bin_width_multiplier <- 1.0
  direct_cpp_budget <- 2e8
  short_timing_inner_reps <- 10L
}

message("Checking C++ backends.")
cpp_available <- compile_sparse_prefix_cpp(root, verbose = TRUE)
if (!cpp_available || !exists("lbfp_direct_cpp", mode = "function") ||
    !exists("lbfp_sparse_prefix_cpp", mode = "function")) {
  stop("Both direct and sparse-prefix C++ backends are required for scaling benchmarks.")
}

records <- list()

add_record <- function(config, fit, direct_time, direct_res, sparse_time, sparse_res) {
  metrics <- grid_occupancy_metrics(fit)
  err <- abs(direct_res$density - sparse_res$density)
  records[[length(records) + 1L]] <<- data.frame(
    scenario = config$scenario,
    n = config$n,
    d = config$d,
    rep = config$rep,
    seed = config$seed,
    bin_width_multiplier = config$bin_width_multiplier,
    G = metrics$G,
    G_occ = metrics$G_occ,
    occupancy_rate = metrics$occupancy_rate,
    theoretical_vertices = 2^fit$d,
    direct_elapsed = direct_time$elapsed,
    sparse_elapsed = sparse_time$elapsed,
    speedup_sparse_vs_direct = direct_time$elapsed / sparse_time$elapsed,
    direct_visited_mean = mean(direct_res$visited),
    direct_visited_median = stats::median(direct_res$visited),
    sparse_visited_mean = mean(sparse_res$visited),
    sparse_visited_median = stats::median(sparse_res$visited),
    sparse_visited_max = max(sparse_res$visited),
    prefix_nodes_mean = mean(sparse_res$prefix_nodes),
    prefix_nodes_median = stats::median(sparse_res$prefix_nodes),
    prefix_nodes_max = max(sparse_res$prefix_nodes),
    visited_ratio_median = stats::median(sparse_res$visited) / (2^fit$d),
    prefix_ratio_median = stats::median(sparse_res$prefix_nodes) / (2^fit$d),
    max_abs_error = max(err),
    mean_abs_error = mean(err),
    all_equal_tolerance = max(err) <= 1e-10,
    direct_timing_inner_reps = direct_time$timing_inner_reps,
    sparse_timing_inner_reps = sparse_time$timing_inner_reps,
    stringsAsFactors = FALSE
  )
}

for (n in ns) {
  for (d in ds) {
    for (rep in seq_len(reps)) {
      seed <- 20260601L + 100000L * rep + 1000L * d + n
      config <- list(
        scenario = scenario,
        n = n,
        d = d,
        rep = rep,
        seed = seed,
        bin_width_multiplier = bin_width_multiplier
      )
      message("Scaling benchmark n=", n, " d=", d, " rep=", rep)
      x <- simulate_grid_data(n = n, d = d, scenario = scenario, seed = seed)
      fit <- fit_grid_for_benchmark(x, bin_width_multiplier = bin_width_multiplier)
      order <- choose_prefix_order(fit$counts$cells, "adaptive")

      nominal_work <- nrow(x) * 2^fit$d
      if (nominal_work > direct_cpp_budget) {
        message("Skipping direct C++ because n * 2^d exceeds direct_cpp_budget.")
        next
      }

      direct_time <- elapsed_repeated(
        lbfp_eval_direct_cpp(fit, x, return_visited = TRUE),
        reps = 1L,
        inner_reps = short_timing_inner_reps
      )
      sparse_time <- elapsed_repeated(
        sparse_prefix_eval_lbfp_cpp(fit, x, order = order, return_visited = TRUE),
        reps = 1L,
        inner_reps = short_timing_inner_reps
      )

      add_record(config, fit, direct_time, direct_time$value, sparse_time, sparse_time$value)
    }
  }
}

res <- do.call(rbind, records)
write_csv_base(res, file.path(out_dir, "scaling_benchmarks.csv"))
write_session_info(file.path(out_dir, "session_info_scaling.txt"))
print(utils::head(res, 20))
message("Wrote ", file.path(out_dir, "scaling_benchmarks.csv"))
