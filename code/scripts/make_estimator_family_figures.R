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
csv <- file.path(out_dir, "estimator_family_benchmarks.csv")

if (!file.exists(csv)) {
  stop(
    "Missing code/results/estimator_family_benchmarks.csv. Run ",
    "`Rscript code/scripts/run_estimator_family_benchmarks.R --quick` first.",
    call. = FALSE
  )
}

df <- utils::read.csv(csv, stringsAsFactors = FALSE)
if (!nrow(df)) stop("estimator_family_benchmarks.csv exists but contains no rows.", call. = FALSE)

median_na <- function(z) {
  if (all(is.na(z))) NA_real_ else stats::median(z, na.rm = TRUE)
}

summary_cols <- c(
  "direct_elapsed",
  "sparse_elapsed",
  "speedup_sparse_vs_direct",
  "occupancy_rate",
  "nominal_stencil",
  "direct_visited_median",
  "sparse_visited_median",
  "prefix_nodes_median",
  "visited_ratio_median",
  "prefix_ratio_median",
  "max_abs_error"
)
groups <- unique(df[c("estimator", "n", "d")])
summary_list <- vector("list", nrow(groups))
for (i in seq_len(nrow(groups))) {
  keep <- df$estimator == groups$estimator[[i]] &
    df$n == groups$n[[i]] &
    df$d == groups$d[[i]]
  vals <- lapply(df[keep, summary_cols, drop = FALSE], median_na)
  summary_list[[i]] <- data.frame(
    groups[i, , drop = FALSE],
    as.data.frame(vals, check.names = FALSE),
    stringsAsFactors = FALSE
  )
}
summary_df <- do.call(rbind, summary_list)
summary_df <- summary_df[order(summary_df$estimator, summary_df$n, summary_df$d), ]
write_csv_base(summary_df, file.path(out_dir, "estimator_family_summary.csv"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  dimension_breaks <- sort(unique(summary_df$d))
  speedup_breaks <- c(0.3, 1, 3, 10)

  speed_df <- subset(summary_df, estimator != "histogram" &
                       is.finite(speedup_sparse_vs_direct) &
                       speedup_sparse_vs_direct > 0)
  p1 <- ggplot(speed_df, aes(x = factor(d), y = estimator, fill = log10(speedup_sparse_vs_direct))) +
    geom_tile() +
    facet_wrap(~ n) +
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
    labs(x = "dimension d", y = "estimator",
         title = "Sparse-prefix speedup across estimator families") +
    theme_bw()
  ggsave(file.path(fig_dir, "estimator_family_speedup_heatmap.png"), p1,
         width = 10, height = 6, dpi = 150)

  p2 <- ggplot(
    speed_df,
    aes(
      x = d,
      y = speedup_sparse_vs_direct,
      color = estimator,
      shape = factor(n),
      group = interaction(estimator, n)
    )
  ) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
    geom_line() +
    geom_point(size = 2) +
    scale_x_continuous(breaks = dimension_breaks) +
    scale_y_log10() +
    labs(x = "dimension d", y = "speedup versus direct evaluation", shape = "n",
         title = "Sparse-prefix speedup by estimator and sample size") +
    theme_bw()
  ggsave(file.path(fig_dir, "estimator_family_speedup_by_dimension.png"), p2,
         width = 11, height = 7, dpi = 150)

  direct_long <- data.frame(
    estimator = summary_df$estimator,
    n = summary_df$n,
    d = summary_df$d,
    mode = "direct",
    elapsed = summary_df$direct_elapsed,
    stringsAsFactors = FALSE
  )
  sparse_long <- data.frame(
    estimator = summary_df$estimator,
    n = summary_df$n,
    d = summary_df$d,
    mode = "sparse-prefix",
    elapsed = summary_df$sparse_elapsed,
    stringsAsFactors = FALSE
  )
  elapsed_df <- rbind(direct_long, sparse_long)
  elapsed_df <- subset(elapsed_df, is.finite(elapsed) & elapsed > 0)
  p2b <- ggplot(
    elapsed_df,
    aes(
      x = d,
      y = elapsed,
      color = estimator,
      linetype = mode,
      shape = factor(n),
      group = interaction(estimator, mode, n)
    )
  ) +
    geom_line() +
    geom_point(size = 2) +
    scale_x_continuous(breaks = dimension_breaks) +
    scale_y_log10() +
    labs(x = "dimension d", y = "elapsed time (s)", shape = "n",
         title = "Direct and sparse-prefix elapsed time by estimator and sample size") +
    theme_bw()
  ggsave(file.path(fig_dir, "estimator_family_elapsed_by_dimension.png"), p2b,
         width = 11, height = 7, dpi = 150)

  prefix_df <- subset(summary_df, estimator != "histogram" &
                        is.finite(prefix_ratio_median) &
                        prefix_ratio_median > 0)
  p3 <- ggplot(prefix_df, aes(x = d, y = prefix_ratio_median, color = estimator)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
    geom_line() +
    geom_point() +
    facet_wrap(~ n) +
    scale_x_continuous(breaks = dimension_breaks) +
    scale_y_log10() +
    labs(x = "dimension d", y = "median prefix nodes / nominal stencil",
         title = "Observed sparse-prefix work by estimator family") +
    theme_bw()
  ggsave(file.path(fig_dir, "estimator_family_prefix_ratio.png"), p3,
         width = 10, height = 6, dpi = 150)
} else {
  png(file.path(fig_dir, "estimator_family_speedup_heatmap.png"), width = 1000, height = 700)
  plot.new()
  title("Sparse-prefix speedup across estimator families")
  text(0.5, 0.5, "Install ggplot2 to generate estimator-family figures.")
  dev.off()
}

write_session_info(file.path(out_dir, "session_info_estimator_family_figures.txt"))
message("Wrote estimator-family figures to ", fig_dir)
message("Wrote estimator-family summary to ", file.path(out_dir, "estimator_family_summary.csv"))
