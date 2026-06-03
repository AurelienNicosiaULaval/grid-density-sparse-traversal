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

tolerance <- as.numeric(parse_arg_value("--tolerance", default = "1e-10"))
if (!is.finite(tolerance) || tolerance <= 0) stop("--tolerance must be positive and finite.")

out_dir <- file.path(root, "code/results")
set.seed(20260531)

message("Checking C++ sparse-prefix backend availability.")
cpp_available <- compile_sparse_prefix_cpp(root, verbose = TRUE)
if (!cpp_available) {
  message("C++ backend unavailable. Correctness checks will still verify the R fallback.")
}

configs <- expand.grid(
  n = c(100L, 500L),
  d = c(2L, 3L, 5L, 8L),
  bins = c(5L, 15L),
  scenario = c("normal", "clustered"),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

records <- list()

add_record <- function(fit, config, seed, backend_used, sparse_time, sparse_density,
                       sparse_score, visited, prefix_nodes, naive_time, naive_density,
                       naive_score) {
  density_error <- abs(naive_density - sparse_density)
  loo_error <- abs(naive_score - sparse_score)
  metrics <- grid_occupancy_metrics(fit)
  records[[length(records) + 1L]] <<- data.frame(
    estimator = "LBFP",
    scenario = config$scenario,
    n = config$n,
    d = config$d,
    bins = config$bins,
    seed = seed,
    max_abs_error = max(density_error),
    mean_abs_error = mean(density_error),
    max_abs_loo_error = max(loo_error),
    mean_abs_loo_error = mean(loo_error),
    all_equal_tolerance = max(density_error) <= tolerance && max(loo_error) <= tolerance,
    tolerance = tolerance,
    backend_used = backend_used,
    elapsed_naive = naive_time$elapsed,
    elapsed_sparse_prefix = sparse_time$elapsed,
    G = metrics$G,
    G_occ = metrics$G_occ,
    occupancy_rate = metrics$occupancy_rate,
    visited_vertices_mean = mean(visited),
    visited_vertices_median = stats::median(visited),
    visited_vertices_max = max(visited),
    prefix_nodes_mean = mean(prefix_nodes),
    prefix_nodes_median = stats::median(prefix_nodes),
    prefix_nodes_max = max(prefix_nodes),
    theoretical_vertices = 2^fit$d,
    stringsAsFactors = FALSE
  )
}

for (row_id in seq_len(nrow(configs))) {
  config <- configs[row_id, ]
  seed <- 20260531L + row_id * 1009L
  message(
    "Correctness check n=", config$n,
    " d=", config$d,
    " bins=", config$bins,
    " scenario=", config$scenario
  )

  x <- simulate_grid_data(config$n, config$d, config$scenario, seed = seed)
  fit <- lbfp_fit(x, nbins = config$bins)
  order <- choose_prefix_order(fit$counts$cells, "adaptive")

  naive_time <- elapsed(lbfp_eval_naive(fit, x, return_visited = TRUE))
  naive <- naive_time$value
  self_weight <- lbfp_self_weight(fit, x)
  naive_score <- loo_from_density(naive$density, self_weight, fit$grid$volume, fit$n)

  sp_r_time <- elapsed(sparse_prefix_eval_lbfp_R(fit, x, order = order, return_visited = TRUE))
  sp_r <- sp_r_time$value
  sp_r_score <- loo_from_density(sp_r$density, self_weight, fit$grid$volume, fit$n)
  add_record(fit, config, seed, "R", sp_r_time, sp_r$density, sp_r_score,
             sp_r$visited, sp_r$prefix_nodes, naive_time, naive$density, naive_score)

  if (cpp_available && exists("lbfp_sparse_prefix_cpp", mode = "function")) {
    sp_cpp_time <- elapsed(sparse_prefix_eval_lbfp_cpp(fit, x, order = order, return_visited = TRUE))
    sp_cpp <- sp_cpp_time$value
    sp_cpp_score <- loo_from_density(sp_cpp$density, self_weight, fit$grid$volume, fit$n)
    add_record(fit, config, seed, "cpp", sp_cpp_time, sp_cpp$density, sp_cpp_score,
               sp_cpp$visited, sp_cpp$prefix_nodes, naive_time, naive$density, naive_score)
  }
}

res <- do.call(rbind, records)
write_csv_base(res, file.path(out_dir, "correctness_checks.csv"))
write_session_info(file.path(out_dir, "session_info_correctness.txt"))

print(res[, c("scenario", "n", "d", "bins", "backend_used", "max_abs_error",
              "max_abs_loo_error", "all_equal_tolerance", "elapsed_naive",
              "elapsed_sparse_prefix")])
message("Wrote ", file.path(out_dir, "correctness_checks.csv"))

if (!all(res$all_equal_tolerance)) {
  stop("Some LBFP correctness checks failed. Inspect code/results/correctness_checks.csv.")
}
