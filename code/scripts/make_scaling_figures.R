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

out_dir <- file.path(root, "code/results")
fig_dir <- file.path(root, "code/figures")
csv <- file.path(out_dir, "scaling_benchmarks.csv")

if (!file.exists(csv)) {
  stop(
    "Missing code/results/scaling_benchmarks.csv. Run ",
    "`Rscript code/scripts/run_scaling_benchmarks.R --quick` first.",
    call. = FALSE
  )
}

df <- utils::read.csv(csv, stringsAsFactors = FALSE)
if (!nrow(df)) stop("scaling_benchmarks.csv exists but contains no rows.", call. = FALSE)

median_na <- function(z) {
  if (all(is.na(z))) NA_real_ else stats::median(z, na.rm = TRUE)
}

summary_df <- stats::aggregate(
  cbind(
    direct_elapsed,
    sparse_elapsed,
    speedup_sparse_vs_direct,
    occupancy_rate,
    theoretical_vertices,
    sparse_visited_median,
    prefix_nodes_median,
    visited_ratio_median,
    prefix_ratio_median,
    max_abs_error
  ) ~ n + d,
  data = df,
  FUN = median_na
)
summary_df <- summary_df[order(summary_df$n, summary_df$d), ]
write_csv_base(summary_df, file.path(out_dir, "scaling_summary.csv"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  dimension_breaks <- sort(unique(summary_df$d))
  speedup_breaks <- c(0.3, 1, 3, 10, 30)

  p1 <- ggplot(summary_df, aes(x = factor(d), y = factor(n), fill = log10(speedup_sparse_vs_direct))) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#b2182b",
      mid = "white",
      high = "#2166ac",
      midpoint = 0,
      breaks = log10(speedup_breaks),
      labels = format(speedup_breaks, trim = TRUE),
      name = "speedup",
      na.value = "grey85"
    ) +
    labs(x = "dimension d", y = "sample size n",
         title = "Sparse-prefix speedup versus direct C++ LBFP") +
    theme_bw()
  ggsave(file.path(fig_dir, "scaling_speedup_heatmap.png"), p1, width = 10, height = 6, dpi = 150)

  p2 <- ggplot(summary_df, aes(x = d, y = speedup_sparse_vs_direct, color = factor(n))) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = dimension_breaks) +
    scale_y_log10() +
    labs(x = "dimension d", y = "speedup versus direct C++", color = "n",
         title = "Speedup transition across dimension") +
    theme_bw()
  ggsave(file.path(fig_dir, "scaling_speedup_by_dimension.png"), p2, width = 10, height = 6, dpi = 150)

  p3 <- ggplot(summary_df, aes(x = d, y = prefix_ratio_median, color = factor(n))) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = dimension_breaks) +
    scale_y_log10() +
    labs(x = "dimension d", y = "median prefix nodes / 2^d", color = "n",
         title = "Observed sparse traversal cost relative to direct stencil") +
    theme_bw()
  ggsave(file.path(fig_dir, "scaling_prefix_ratio_by_dimension.png"), p3, width = 10, height = 6, dpi = 150)

  p4 <- ggplot(summary_df, aes(x = occupancy_rate, y = speedup_sparse_vs_direct, color = factor(d))) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
    geom_point(size = 2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "grid occupancy rate", y = "speedup versus direct C++", color = "d",
         title = "Speedup as a function of grid occupancy") +
    theme_bw()
  ggsave(file.path(fig_dir, "scaling_speedup_vs_occupancy.png"), p4, width = 10, height = 6, dpi = 150)
} else {
  png(file.path(fig_dir, "scaling_speedup_by_dimension.png"), width = 1000, height = 700)
  plot(summary_df$d, summary_df$speedup_sparse_vs_direct, log = "y",
       xlab = "dimension d", ylab = "speedup versus direct C++")
  abline(h = 1, lty = 2)
  dev.off()
}

write_session_info(file.path(out_dir, "session_info_scaling_figures.txt"))
message("Wrote scaling figures to ", fig_dir)
message("Wrote scaling summary to ", file.path(out_dir, "scaling_summary.csv"))
