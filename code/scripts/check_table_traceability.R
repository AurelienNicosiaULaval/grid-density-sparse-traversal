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
args <- commandArgs(trailingOnly = TRUE)
strict <- !("--no-stop" %in% args)

tables_dir <- file.path(root, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(name) {
  path <- file.path(root, "code/results", name)
  if (!file.exists(path)) stop("Missing required result file: ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE)
}

checks <- list()

add_check <- function(location, quantity, source_file, computed, expected,
                      status = identical(as.character(computed), as.character(expected))) {
  checks[[length(checks) + 1L]] <<- data.frame(
    location = location,
    quantity = quantity,
    source_file = source_file,
    computed = as.character(computed),
    expected = as.character(expected),
    status = if (isTRUE(status)) "OK" else "MISMATCH",
    stringsAsFactors = FALSE
  )
}

fmt2 <- function(x) sprintf("%.2f", x)
fmt3 <- function(x) sprintf("%.3f", x)
sci <- function(x) format(x, scientific = TRUE, digits = 6)

family <- read_result("estimator_family_benchmarks.csv")
for (est in c("ASH", "LBFP", "GLBFP")) {
  sub <- subset(family, estimator == est & is.finite(speedup_sparse_vs_direct))
  max_i <- which.max(sub$speedup_sparse_vs_direct)
  expected <- switch(
    est,
    ASH = c(range = "2--6", speedup = "10.50", dmax = "6", error = "0.00"),
    LBFP = c(range = "2--6", speedup = "2.47", dmax = "6", error = "4.16e-17"),
    GLBFP = c(range = "2--6", speedup = "22.00", dmax = "6", error = "5.55e-17")
  )
  add_check("Table 2", paste(est, "tested d"), "code/results/estimator_family_benchmarks.csv",
            paste0(min(sub$d), "--", max(sub$d)), expected[["range"]])
  add_check("Table 2", paste(est, "max speedup"), "code/results/estimator_family_benchmarks.csv",
            fmt2(sub$speedup_sparse_vs_direct[[max_i]]), expected[["speedup"]])
  add_check("Table 2", paste(est, "d at max"), "code/results/estimator_family_benchmarks.csv",
            sub$d[[max_i]], expected[["dmax"]])
  err_display <- if (est == "ASH") "0.00" else sub("e-0", "e-", sprintf("%.2e", max(sub$max_abs_error, na.rm = TRUE)))
  add_check("Table 2", paste(est, "max error"), "code/results/estimator_family_benchmarks.csv",
            err_display, expected[["error"]])
}

correctness <- read_result("correctness_checks.csv")
add_check("Section 5.1", "correctness rows", "code/results/correctness_checks.csv",
          nrow(correctness), 64L)
add_check("Section 5.1", "max absolute density error", "code/results/correctness_checks.csv",
          sub("e-0", "e-", sprintf("%.2e", max(correctness$max_abs_error, na.rm = TRUE))), "8.88e-16")
add_check("Section 5.1", "family max error below 6e-17", "code/results/estimator_family_benchmarks.csv",
          sci(max(family$max_abs_error, na.rm = TRUE)), "< 6e-17",
          max(family$max_abs_error, na.rm = TRUE) < 6e-17)

scaling <- read_result("scaling_summary.csv")
expected_scaling <- data.frame(
  n = c(500L, 1000L, 2000L, 5000L, 10000L),
  d_first = c(6L, 7L, 5L, 5L, 6L),
  max_speedup = c("46.50", "33.42", "28.79", "55.87", "36.07"),
  d_at_max = c(12L, 12L, 12L, 12L, 12L),
  min_prefix_ratio = c("0.036", "0.044", "0.057", "0.078", "0.097"),
  max_error = c("1.39e-17", "1.39e-17", "3.47e-17", "1.73e-17", "2.08e-17"),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(expected_scaling))) {
  sub <- subset(scaling, n == expected_scaling$n[[i]])
  first <- min(sub$d[sub$speedup_sparse_vs_direct > 1])
  max_i <- which.max(sub$speedup_sparse_vs_direct)
  add_check("Table 3", paste(expected_scaling$n[[i]], "d_first"), "code/results/scaling_summary.csv",
            first, expected_scaling$d_first[[i]])
  add_check("Table 3", paste(expected_scaling$n[[i]], "max speedup"), "code/results/scaling_summary.csv",
            fmt2(sub$speedup_sparse_vs_direct[[max_i]]), expected_scaling$max_speedup[[i]])
  add_check("Table 3", paste(expected_scaling$n[[i]], "d at max"), "code/results/scaling_summary.csv",
            sub$d[[max_i]], expected_scaling$d_at_max[[i]])
  add_check("Table 3", paste(expected_scaling$n[[i]], "min prefix ratio"), "code/results/scaling_summary.csv",
            fmt3(min(sub$prefix_ratio_median)), expected_scaling$min_prefix_ratio[[i]])
  add_check("Table 3", paste(expected_scaling$n[[i]], "max error"), "code/results/scaling_summary.csv",
            sub("e-0", "e-", sprintf("%.2e", max(sub$max_abs_error))), expected_scaling$max_error[[i]])
}

occupancy <- read_result("occupancy_benchmarks.csv")
expected_occupancy <- data.frame(
  d = c(4L, 6L, 8L, 10L, 12L),
  median_speedup = c("1.20", "1.66", "3.61", "11.70", "10.42"),
  median_R_cost = c("0.853", "0.566", "0.286", "0.121", "0.085"),
  criterion_accuracy = c("1.00", "0.95", "0.90", "1.00", "1.00"),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(expected_occupancy))) {
  sub <- subset(occupancy, d == expected_occupancy$d[[i]])
  add_check("Table 4", paste(expected_occupancy$d[[i]], "median speedup"), "code/results/occupancy_benchmarks.csv",
            fmt2(stats::median(sub$speedup)), expected_occupancy$median_speedup[[i]])
  add_check("Table 4", paste(expected_occupancy$d[[i]], "median R_cost"), "code/results/occupancy_benchmarks.csv",
            fmt3(stats::median(sub$R_cost)), expected_occupancy$median_R_cost[[i]])
  add_check("Table 4", paste(expected_occupancy$d[[i]], "criterion accuracy"), "code/results/occupancy_benchmarks.csv",
            fmt2(mean(sub$prediction_correct)), expected_occupancy$criterion_accuracy[[i]])
}
add_check("Section 5.5", "criterion correct count", "code/results/occupancy_benchmarks.csv",
          sum(occupancy$prediction_correct), "97")
add_check("Section 5.5", "occupancy max absolute error", "code/results/occupancy_benchmarks.csv",
          sub("e-0", "e-", sprintf("%.2e", max(occupancy$max_abs_error, na.rm = TRUE))), "1.39e-17")

selector <- read_result("mode_selector_benchmarks.csv")
expected_selector <- data.frame(
  d = c(4L, 6L, 8L, 10L, 12L),
  accuracy = c("0.75", "0.95", "0.95", "1.00", "1.00"),
  median_regret = c("1.00", "1.00", "1.00", "1.00", "1.00"),
  max_regret = c("1.38", "1.00", "2.09", "1.00", "1.00"),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(expected_selector))) {
  sub <- subset(selector, d == expected_selector$d[[i]])
  add_check("Table 5", paste(expected_selector$d[[i]], "selector accuracy"), "code/results/mode_selector_benchmarks.csv",
            fmt2(mean(sub$selector_correct)), expected_selector$accuracy[[i]])
  add_check("Table 5", paste(expected_selector$d[[i]], "median regret"), "code/results/mode_selector_benchmarks.csv",
            fmt2(stats::median(sub$regret)), expected_selector$median_regret[[i]])
  add_check("Table 5", paste(expected_selector$d[[i]], "maximum regret"), "code/results/mode_selector_benchmarks.csv",
            fmt2(max(sub$regret)), expected_selector$max_regret[[i]])
}
add_check("Section 5.6", "selector correct count", "code/results/mode_selector_benchmarks.csv",
          sum(selector$selector_correct), "93")
add_check("Section 5.6", "selector median regret", "code/results/mode_selector_benchmarks.csv",
          fmt2(stats::median(selector$regret)), "1.00")
add_check("Section 5.6", "selector largest regret", "code/results/mode_selector_benchmarks.csv",
          fmt2(max(selector$regret)), "2.09")
add_check("Section 5.6", "selector max absolute error", "code/results/mode_selector_benchmarks.csv",
          sub("e-0", "e-", sprintf("%.2e", max(selector$max_abs_error, na.rm = TRUE))), "1.39e-17")

trace <- do.call(rbind, checks)
out <- file.path(tables_dir, "table_traceability.csv")
utils::write.csv(trace, out, row.names = FALSE)

bad <- trace[trace$status != "OK", , drop = FALSE]
if (nrow(bad)) {
  message("Traceability mismatches written to ", out)
  print(bad)
  if (strict) stop("Table traceability check failed.", call. = FALSE)
}

message("Wrote ", out)
message("All table traceability checks passed.")
