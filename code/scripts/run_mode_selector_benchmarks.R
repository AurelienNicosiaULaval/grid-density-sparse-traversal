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
pilot_size <- as.integer(parse_arg_value("--pilot-size", default = "100"))
tau <- as.numeric(parse_arg_value("--tau", default = "0.8"))
base_seed <- as.integer(parse_arg_value("--seed", default = "20260617"))
tolerance <- as.numeric(parse_arg_value("--tolerance", default = "1e-10"))

if (is.na(M) || M < 1L) stop("M must be a positive integer.", call. = FALSE)
if (is.na(pilot_size) || pilot_size < 1L) stop("pilot size must be a positive integer.", call. = FALSE)
if (!is.finite(tau) || tau <= 0) stop("tau must be positive.", call. = FALSE)

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
    "Mode-selector benchmark ", ii, "/", nrow(configs),
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
  pilot_n <- min(pilot_size, M)
  pilot_index <- sort(sample.int(M, pilot_n, replace = FALSE))
  pilot_query <- query[pilot_index, , drop = FALSE]
  S <- 2^d
  metrics <- grid_occupancy_metrics(fit)

  pilot_res <- sparse_prefix_eval_lbfp_cpp(fit, pilot_query, order = order, return_visited = TRUE)
  pbar_hat <- mean(pilot_res$prefix_nodes)
  predicted_sparse <- (d * metrics$G_occ + M * pbar_hat) < (tau * M * S)
  predicted_mode <- if (predicted_sparse) "sparse-prefix" else "direct"

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
  pbar_actual <- mean(sparse_res$prefix_nodes)
  direct_elapsed <- direct_time$elapsed
  sparse_elapsed <- sparse_time$elapsed
  actual_fastest_mode <- if (sparse_elapsed < direct_elapsed) "sparse-prefix" else "direct"
  selected_time <- if (predicted_mode == "sparse-prefix") sparse_elapsed else direct_elapsed
  oracle_time <- min(direct_elapsed, sparse_elapsed)

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
    pilot_size = pilot_n,
    tau = tau,
    stencil_size = S,
    pbar_hat = pbar_hat,
    pbar_actual = pbar_actual,
    predicted_mode = predicted_mode,
    actual_fastest_mode = actual_fastest_mode,
    selector_correct = predicted_mode == actual_fastest_mode,
    direct_time = direct_elapsed,
    sparse_prefix_time = sparse_elapsed,
    selected_time = selected_time,
    oracle_time = oracle_time,
    regret = selected_time / oracle_time,
    max_abs_error = max(err),
    mean_abs_error = mean(err),
    all_equal_tolerance = max(err) <= tolerance,
    timing_inner_reps = inner_reps,
    seed = seed,
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, records)
res <- res[order(res$d, res$occupancy_pattern, res$query_regime, res$target_occupancy_fraction), ]
write_csv_base(res, file.path(out_dir, "mode_selector_benchmarks.csv"))

summary_df <- stats::aggregate(
  cbind(regret, max_abs_error, pbar_hat, pbar_actual, selected_time, oracle_time) ~
    d + occupancy_pattern + query_regime,
  data = res,
  FUN = median_na
)
accuracy_df <- stats::aggregate(
  selector_correct ~ d + occupancy_pattern + query_regime,
  data = res,
  FUN = mean
)
n_configs_df <- stats::aggregate(
  target_occupancy_fraction ~ d + occupancy_pattern + query_regime,
  data = res,
  FUN = length
)
summary_df <- merge(summary_df, accuracy_df, by = c("d", "occupancy_pattern", "query_regime"))
names(summary_df)[names(summary_df) == "selector_correct"] <- "selector_accuracy"
summary_df <- merge(summary_df, n_configs_df, by = c("d", "occupancy_pattern", "query_regime"))
names(summary_df)[names(summary_df) == "target_occupancy_fraction"] <- "n_configs"
summary_df$n_configs <- as.integer(summary_df$n_configs)
summary_df <- summary_df[order(summary_df$d, summary_df$occupancy_pattern, summary_df$query_regime), ]
write_csv_base(summary_df, file.path(out_dir, "mode_selector_summary.csv"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  p1 <- ggplot(res, aes(x = factor(d), y = regret, color = query_regime)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey45") +
    geom_count(
      aes(size = after_stat(n)),
      position = position_dodge(width = 0.55),
      alpha = 0.75
    ) +
    scale_size_area(max_size = 8, breaks = c(1, 5, 10), name = "configurations") +
    scale_y_continuous(breaks = c(1, 1.1, 1.25, 1.5, 2.0)) +
    coord_cartesian(ylim = c(0.98, max(res$regret) * 1.08)) +
    labs(x = "dimension d", y = "selected time / oracle time", color = "query regime") +
    theme_bw()
  ggsave(file.path(fig_dir, "mode_selector_regret.pdf"), p1, width = 8.5, height = 5.2)

  mode_df <- as.data.frame(table(
    d = res$d,
    predicted = res$predicted_mode,
    fastest = res$actual_fastest_mode
  ))
  p2 <- ggplot(mode_df, aes(x = predicted, y = fastest, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 3) +
    facet_wrap(~ d) +
    scale_fill_gradient(low = "#f7f7f7", high = "#2166ac", name = "count") +
    labs(x = "selected mode", y = "fastest mode") +
    theme_bw()
  ggsave(file.path(fig_dir, "mode_selector_modes.pdf"), p2, width = 9, height = 5)

  p3 <- ggplot(res, aes(x = pbar_actual, y = pbar_hat, color = factor(d), shape = query_regime)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
    geom_point(size = 2, alpha = 0.85) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "actual mean prefix nodes", y = "pilot-estimated mean prefix nodes",
         color = "d", shape = "query regime") +
    theme_bw()
  ggsave(file.path(fig_dir, "prefix_nodes_pilot_vs_actual.pdf"), p3, width = 8, height = 5.2)
} else {
  pdf(file.path(fig_dir, "mode_selector_regret.pdf"), width = 8, height = 5)
  plot(
    jitter(as.numeric(factor(res$d)), amount = 0.08),
    res$regret,
    xaxt = "n",
    xlab = "dimension d",
    ylab = "selected time / oracle time",
    pch = 19,
    col = "grey35"
  )
  axis(1, at = seq_along(sort(unique(res$d))), labels = sort(unique(res$d)))
  abline(h = 1, lty = 2)
  dev.off()
}

write_session_info(file.path(out_dir, "session_info_mode_selector.txt"))
print(utils::head(res, 20))
message("Wrote ", file.path(out_dir, "mode_selector_benchmarks.csv"))
message("Wrote ", file.path(out_dir, "mode_selector_summary.csv"))
