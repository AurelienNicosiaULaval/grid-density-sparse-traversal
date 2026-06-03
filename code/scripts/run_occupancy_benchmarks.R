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
source(file.path(root, "code/R/05_controlled_benchmark_utils.R"))
ensure_project_dirs(root)

out_dir <- file.path(root, "code/results")
fig_dir <- file.path(root, "code/figures")
M <- as.integer(parse_arg_value("--M", default = "1000"))
base_seed <- as.integer(parse_arg_value("--seed", default = "20260601"))
tolerance <- as.numeric(parse_arg_value("--tolerance", default = "1e-10"))

if (is.na(M) || M < 1L) stop("M must be a positive integer.", call. = FALSE)

message("Checking C++ backends.")
cpp_available <- compile_sparse_prefix_cpp(root, verbose = TRUE)
if (!cpp_available || !exists("lbfp_direct_cpp", mode = "function") ||
    !exists("lbfp_sparse_prefix_cpp", mode = "function")) {
  stop("Both direct and sparse-prefix C++ backends are required.", call. = FALSE)
}

configs <- controlled_config_grid()
records <- vector("list", nrow(configs))

for (ii in seq_len(nrow(configs))) {
  cfg <- configs[ii, ]
  d <- as.integer(cfg$d)
  side_length <- controlled_side_length(d)
  seed <- base_seed + ii * 1009L + d * 100003L
  message(
    "Occupancy benchmark ", ii, "/", nrow(configs),
    ": d=", d,
    " rho=", cfg$target_occupancy_fraction,
    " pattern=", cfg$occupancy_pattern,
    " query=", cfg$query_regime
  )

  cells <- generate_controlled_cells(
    d = d,
    side_length = side_length,
    target_occupancy_fraction = cfg$target_occupancy_fraction,
    occupancy_pattern = cfg$occupancy_pattern,
    seed = seed
  )
  fit <- make_controlled_lbfp_fit(cells, side_length = side_length)
  query <- generate_controlled_queries(
    cells = cells,
    side_length = side_length,
    M = M,
    query_regime = cfg$query_regime,
    seed = seed + 17L
  )
  order <- choose_prefix_order(fit$counts$cells, "adaptive")
  inner_reps <- controlled_timing_inner_reps(d)

  direct_time <- elapsed_repeated(
    lbfp_eval_direct_cpp(fit, query, return_visited = TRUE),
    reps = 1L,
    inner_reps = inner_reps
  )
  sparse_time <- elapsed_repeated(
    sparse_prefix_eval_lbfp_cpp(fit, query, order = order, return_visited = TRUE),
    reps = 1L,
    inner_reps = inner_reps
  )

  direct_res <- direct_time$value
  sparse_res <- sparse_time$value
  err <- abs(direct_res$density - sparse_res$density)
  metrics <- grid_occupancy_metrics(fit)
  S <- 2^d
  pbar <- mean(sparse_res$prefix_nodes)
  R_cost <- (d * metrics$G_occ + M * pbar) / (M * S)
  direct_elapsed <- direct_time$elapsed
  sparse_elapsed <- sparse_time$elapsed

  records[[ii]] <- data.frame(
    d = d,
    grid_side_length = side_length,
    G = metrics$G,
    occupancy_pattern = cfg$occupancy_pattern,
    query_regime = cfg$query_regime,
    target_occupancy_fraction = cfg$target_occupancy_fraction,
    actual_occupancy_fraction = metrics$occupancy_rate,
    G_occ = metrics$G_occ,
    M = M,
    stencil_size = S,
    pbar = pbar,
    direct_time = direct_elapsed,
    sparse_prefix_time = sparse_elapsed,
    speedup = direct_elapsed / sparse_elapsed,
    R_cost = R_cost,
    predicted_sparse = R_cost < 1,
    sparse_actually_faster = sparse_elapsed < direct_elapsed,
    prediction_correct = (R_cost < 1) == (sparse_elapsed < direct_elapsed),
    max_abs_error = max(err),
    mean_abs_error = mean(err),
    all_equal_tolerance = max(err) <= tolerance,
    prefix_nodes_median = stats::median(sparse_res$prefix_nodes),
    prefix_nodes_max = max(sparse_res$prefix_nodes),
    direct_visited_mean = mean(direct_res$visited),
    sparse_visited_mean = mean(sparse_res$visited),
    timing_inner_reps = inner_reps,
    seed = seed,
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, records)
res <- res[order(res$d, res$occupancy_pattern, res$query_regime, res$target_occupancy_fraction), ]
write_csv_base(res, file.path(out_dir, "occupancy_benchmarks.csv"))

summary_df <- stats::aggregate(
  cbind(speedup, R_cost, max_abs_error, pbar, sparse_prefix_time, direct_time) ~
    d + occupancy_pattern + query_regime,
  data = res,
  FUN = median_na
)
accuracy_df <- stats::aggregate(
  prediction_correct ~ d + occupancy_pattern + query_regime,
  data = res,
  FUN = mean
)
n_configs_df <- stats::aggregate(
  target_occupancy_fraction ~ d + occupancy_pattern + query_regime,
  data = res,
  FUN = length
)
summary_df <- merge(summary_df, accuracy_df, by = c("d", "occupancy_pattern", "query_regime"))
names(summary_df)[names(summary_df) == "prediction_correct"] <- "prediction_accuracy"
summary_df <- merge(summary_df, n_configs_df, by = c("d", "occupancy_pattern", "query_regime"))
names(summary_df)[names(summary_df) == "target_occupancy_fraction"] <- "n_configs"
summary_df$n_configs <- as.integer(summary_df$n_configs)
summary_df <- summary_df[order(summary_df$d, summary_df$occupancy_pattern, summary_df$query_regime), ]
write_csv_base(summary_df, file.path(out_dir, "occupancy_summary.csv"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  p1 <- ggplot(
    res,
    aes(x = target_occupancy_fraction, y = speedup, color = factor(d), group = interaction(d, query_regime))
  ) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey45") +
    geom_line(alpha = 0.7) +
    geom_point(size = 1.8) +
    facet_grid(occupancy_pattern ~ query_regime) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      x = "target occupancy fraction",
      y = "direct time / sparse-prefix time",
      color = "d"
    ) +
    theme_bw()
  ggsave(file.path(fig_dir, "occupancy_speedup.pdf"), p1, width = 9.5, height = 6.2)

  p2 <- ggplot(res, aes(x = R_cost, y = speedup, color = factor(d), shape = occupancy_pattern)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey45") +
    geom_hline(yintercept = 1, linetype = 2, color = "grey45") +
    geom_point(size = 2, alpha = 0.85) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      x = "predicted cost ratio R_cost",
      y = "direct time / sparse-prefix time",
      color = "d",
      shape = "occupancy pattern"
    ) +
    theme_bw()
  ggsave(file.path(fig_dir, "cost_ratio_vs_speedup.pdf"), p2, width = 8, height = 5.5)

  acc_df <- stats::aggregate(
    prediction_correct ~ d + occupancy_pattern + query_regime,
    data = res,
    FUN = mean
  )
  p3 <- ggplot(acc_df, aes(x = factor(d), y = query_regime, fill = prediction_correct)) +
    geom_tile(color = "white") +
    facet_wrap(~ occupancy_pattern) +
    scale_fill_gradient(limits = c(0, 1), low = "#fddbc7", high = "#2166ac", name = "accuracy") +
    labs(x = "dimension d", y = "query regime") +
    theme_bw()
  ggsave(file.path(fig_dir, "mode_selector_accuracy.pdf"), p3, width = 8, height = 4.8)
} else {
  pdf(file.path(fig_dir, "occupancy_speedup.pdf"), width = 9, height = 6)
  plot(res$target_occupancy_fraction, res$speedup, log = "xy",
       xlab = "target occupancy fraction", ylab = "speedup")
  abline(h = 1, lty = 2)
  dev.off()
}

write_session_info(file.path(out_dir, "session_info_occupancy.txt"))
print(utils::head(res, 20))
message("Wrote ", file.path(out_dir, "occupancy_benchmarks.csv"))
message("Wrote ", file.path(out_dir, "occupancy_summary.csv"))
