# 04_benchmark_utils.R
# Benchmark orchestration utilities.

simulate_grid_data <- function(n, d, scenario = c("normal", "mixture", "clustered", "uniform", "low_rank"),
                               seed = 1L, n_clusters = 10L) {
  if (length(n) != 1L || length(d) != 1L || n < 2L || d < 1L) {
    stop("n must be at least 2 and d must be at least 1.")
  }
  scenario <- match.arg(scenario)
  set.seed(seed)
  if (scenario == "normal") {
    x <- matrix(stats::rnorm(n * d), n, d)
  } else if (scenario == "mixture") {
    z <- stats::rbinom(n, 1L, 0.5)
    mu <- matrix(ifelse(z == 1L, 2, -2), n, d)
    x <- matrix(stats::rnorm(n * d), n, d) + mu
  } else if (scenario == "clustered") {
    centers <- matrix(stats::rnorm(n_clusters * d, sd = 3), n_clusters, d)
    id <- sample.int(n_clusters, n, replace = TRUE)
    x <- centers[id, , drop = FALSE] + matrix(stats::rnorm(n * d, sd = 0.15), n, d)
  } else if (scenario == "uniform") {
    x <- matrix(stats::runif(n * d, -1, 1), n, d)
  } else if (scenario == "low_rank") {
    q <- max(1L, min(3L, d))
    z <- matrix(stats::rnorm(n * q), n, q)
    A <- matrix(stats::rnorm(q * d), q, d)
    x <- z %*% A + matrix(stats::rnorm(n * d, sd = 0.05), n, d)
  }
  colnames(x) <- paste0("X", seq_len(d))
  x
}

fit_grid_for_benchmark <- function(x, bins_per_dim = NULL, bin_width_multiplier = 1) {
  x <- as_numeric_matrix(x)
  if (length(bin_width_multiplier) != 1L || !is.finite(bin_width_multiplier) || bin_width_multiplier <= 0) {
    stop("bin_width_multiplier must be a positive finite scalar.")
  }
  if (!is.null(bins_per_dim)) {
    lbfp_fit(x, nbins = bins_per_dim)
  } else {
    bw <- default_bin_width(x) * bin_width_multiplier
    lbfp_fit(x, bin_width = bw)
  }
}

empty_time <- function() {
  list(
    value = NULL,
    elapsed = NA_real_,
    user = NA_real_,
    system = NA_real_,
    elapsed_iqr = NA_real_,
    user_iqr = NA_real_,
    system_iqr = NA_real_,
    elapsed_min = NA_real_,
    elapsed_max = NA_real_,
    timing_reps = 0L,
    timing_inner_reps = 0L
  )
}

elapsed_repeated <- function(expr, reps = 1L, inner_reps = 1L) {
  reps <- as.integer(reps)
  inner_reps <- as.integer(inner_reps)
  if (length(reps) != 1L || is.na(reps) || reps < 1L) stop("reps must be a positive integer.")
  if (length(inner_reps) != 1L || is.na(inner_reps) || inner_reps < 1L) {
    stop("inner_reps must be a positive integer.")
  }

  expr_sub <- substitute(expr)
  env <- parent.frame()
  elapsed_values <- numeric(reps)
  user_values <- numeric(reps)
  system_values <- numeric(reps)
  value <- NULL

  for (rep_id in seq_len(reps)) {
    gc(FALSE)
    t <- system.time({
      for (inner_id in seq_len(inner_reps)) {
        value <- eval(expr_sub, env)
      }
    })
    elapsed_values[[rep_id]] <- unname(t[["elapsed"]]) / inner_reps
    user_values[[rep_id]] <- unname(t[["user.self"]]) / inner_reps
    system_values[[rep_id]] <- unname(t[["sys.self"]]) / inner_reps
  }

  list(
    value = value,
    elapsed = stats::median(elapsed_values),
    user = stats::median(user_values),
    system = stats::median(system_values),
    elapsed_iqr = stats::IQR(elapsed_values),
    user_iqr = stats::IQR(user_values),
    system_iqr = stats::IQR(system_values),
    elapsed_min = min(elapsed_values),
    elapsed_max = max(elapsed_values),
    timing_reps = reps,
    timing_inner_reps = inner_reps
  )
}

benchmark_lbfp_methods <- function(x, d_label = ncol(x), n_label = nrow(x),
                                   scenario = "unknown", bin_width_multiplier = 1,
                                   bins_per_dim = NULL, run_naive = TRUE,
                                   run_cpp = TRUE, run_hist = TRUE,
                                   naive_budget = 2e7, seed = 1L,
                                   tolerance = 1e-10,
                                   cpp_available = NULL,
                                   direct_cpp_budget = 5e7,
                                   timing_reps = 1L,
                                   short_timing_inner_reps = 10L) {
  x <- as_numeric_matrix(x)
  if (nrow(x) < 2L) stop("At least two observations are required.")
  if (!is.finite(tolerance) || tolerance <= 0) stop("tolerance must be positive and finite.")
  timing_reps <- as.integer(timing_reps)
  short_timing_inner_reps <- as.integer(short_timing_inner_reps)
  if (is.na(timing_reps) || timing_reps < 1L) stop("timing_reps must be a positive integer.")
  if (is.na(short_timing_inner_reps) || short_timing_inner_reps < 1L) {
    stop("short_timing_inner_reps must be a positive integer.")
  }
  if (!is.finite(direct_cpp_budget) || direct_cpp_budget <= 0) {
    stop("direct_cpp_budget must be positive and finite.")
  }

  fit_time <- elapsed({
    fit <- fit_grid_for_benchmark(x, bins_per_dim = bins_per_dim,
                                  bin_width_multiplier = bin_width_multiplier)
  })
  fit <- fit_time$value
  order <- choose_prefix_order(fit$counts$cells, "adaptive")
  theoretical_vertices <- 2^fit$d
  grid_metrics <- grid_occupancy_metrics(fit)
  bins_label <- if (is.null(bins_per_dim)) NA_integer_ else as.integer(bins_per_dim[[1L]])
  records <- list()

  add_record <- function(method, backend_used, time, visited = NULL, prefix_nodes = NULL,
                         max_abs_error = NA_real_, mean_abs_error = NA_real_,
                         notes = "") {
    all_equal_tolerance <- if (is.na(max_abs_error)) NA else isTRUE(max_abs_error <= tolerance)
    records[[length(records) + 1L]] <<- data.frame(
      scenario = scenario,
      n = n_label,
      d = d_label,
      method = method,
      backend_used = backend_used,
      elapsed = time$elapsed,
      user = time$user,
      system = time$system,
      elapsed_iqr = time$elapsed_iqr,
      user_iqr = time$user_iqr,
      system_iqr = time$system_iqr,
      elapsed_min = time$elapsed_min,
      elapsed_max = time$elapsed_max,
      timing_reps = time$timing_reps,
      timing_inner_reps = time$timing_inner_reps,
      fit_elapsed = fit_time$elapsed,
      G = grid_metrics$G,
      G_occ = grid_metrics$G_occ,
      occupancy_rate = grid_metrics$occupancy_rate,
      visited_vertices_mean = if (length(visited)) mean(visited) else NA_real_,
      visited_vertices_median = if (length(visited)) stats::median(visited) else NA_real_,
      visited_vertices_max = if (length(visited)) max(visited) else NA_real_,
      prefix_nodes_mean = if (length(prefix_nodes)) mean(prefix_nodes) else NA_real_,
      prefix_nodes_median = if (length(prefix_nodes)) stats::median(prefix_nodes) else NA_real_,
      prefix_nodes_max = if (length(prefix_nodes)) max(prefix_nodes) else NA_real_,
      theoretical_vertices = theoretical_vertices,
      nominal_stencil = theoretical_vertices,
      median_visited = if (length(visited)) stats::median(visited) else NA_real_,
      max_abs_error = max_abs_error,
      mean_abs_error = mean_abs_error,
      all_equal_tolerance = all_equal_tolerance,
      tolerance = tolerance,
      speedup = NA_real_,
      speedup_vs_naive = NA_real_,
      speedup_vs_direct_cpp = NA_real_,
      elapsed_naive = NA_real_,
      elapsed_direct_cpp = NA_real_,
      bin_width_multiplier = bin_width_multiplier,
      bins = bins_label,
      seed = seed,
      notes = notes,
      stringsAsFactors = FALSE
    )
  }

  if (run_hist) {
    message("Running sparse histogram baseline.")
    t_hist <- elapsed_repeated(hist_eval(fit, x),
                               reps = timing_reps,
                               inner_reps = short_timing_inner_reps)
    add_record("histogram_sparse_R", "R", t_hist, notes = "Histogram baseline; not an LBFP equality check.")
  }

  naive_res <- NULL
  direct_cpp_res <- NULL
  run_naive_actual <- isTRUE(run_naive) && fit$d <= 30L && nrow(x) * theoretical_vertices <= naive_budget
  if (run_naive_actual) {
    message("Running naive LBFP reference.")
    t_naive <- elapsed_repeated(lbfp_eval_naive(fit, x, return_visited = TRUE),
                                reps = timing_reps,
                                inner_reps = 1L)
    naive_res <- t_naive$value
    add_record("LBFP_naive_R", "R", t_naive,
               visited = naive_res$visited,
               max_abs_error = 0,
               mean_abs_error = 0)
  } else if (isTRUE(run_naive)) {
    message("Skipping naive LBFP reference because n * 2^d exceeds naive_budget.")
  }

  if (run_cpp) {
    if (is.null(cpp_available)) {
      cpp_available <- compile_sparse_prefix_cpp(project_root_from_script(), verbose = TRUE)
    }
    run_direct_cpp <- isTRUE(cpp_available) &&
      exists("lbfp_direct_cpp", mode = "function") &&
      fit$d <= 30L &&
      nrow(x) * theoretical_vertices <= direct_cpp_budget
    if (run_direct_cpp) {
      message("Running direct LBFP C++ backend.")
      t_direct_cpp <- elapsed_repeated(lbfp_eval_direct_cpp(fit, x, return_visited = TRUE),
                                       reps = timing_reps,
                                       inner_reps = short_timing_inner_reps)
      direct_cpp_res <- t_direct_cpp$value
      if (!is.null(naive_res)) {
        err_direct_cpp <- abs(direct_cpp_res$density - naive_res$density)
        max_err_direct_cpp <- max(err_direct_cpp)
        mean_err_direct_cpp <- mean(err_direct_cpp)
      } else {
        max_err_direct_cpp <- 0
        mean_err_direct_cpp <- 0
      }
      add_record("LBFP_direct_cpp", "cpp", t_direct_cpp,
                 visited = direct_cpp_res$visited,
                 max_abs_error = max_err_direct_cpp,
                 mean_abs_error = mean_err_direct_cpp,
                 notes = if (is.null(naive_res)) "Used as exact direct C++ reference; naive R was skipped." else "")
    } else if (isTRUE(cpp_available) && exists("lbfp_direct_cpp", mode = "function")) {
      add_record("LBFP_direct_cpp", "cpp", empty_time(),
                 notes = "Direct C++ backend skipped because n * 2^d exceeds direct_cpp_budget.")
    }
  }

  message("Running sparse-prefix LBFP R backend.")
  t_sp_r <- elapsed_repeated(sparse_prefix_eval_lbfp_R(fit, x, order = order, return_visited = TRUE),
                             reps = timing_reps,
                             inner_reps = 1L)
  sp_r <- t_sp_r$value
  reference_density <- if (!is.null(naive_res)) {
    naive_res$density
  } else if (!is.null(direct_cpp_res)) {
    direct_cpp_res$density
  } else {
    NULL
  }
  if (!is.null(reference_density)) {
    err_r <- abs(sp_r$density - reference_density)
    max_err_r <- max(err_r)
    mean_err_r <- mean(err_r)
  } else {
    max_err_r <- NA_real_
    mean_err_r <- NA_real_
  }
  add_record("LBFP_sparse_prefix_R", "R", t_sp_r,
             visited = sp_r$visited,
             prefix_nodes = sp_r$prefix_nodes,
             max_abs_error = max_err_r,
             mean_abs_error = mean_err_r)

  if (run_cpp) {
    if (is.null(cpp_available)) {
      cpp_available <- compile_sparse_prefix_cpp(project_root_from_script(), verbose = TRUE)
    }
    if (isTRUE(cpp_available) && exists("lbfp_sparse_prefix_cpp", mode = "function")) {
      message("Running sparse-prefix LBFP C++ backend.")
      t_cpp <- elapsed_repeated(sparse_prefix_eval_lbfp_cpp(fit, x, order = order, return_visited = TRUE),
                                reps = timing_reps,
                                inner_reps = short_timing_inner_reps)
      sp_cpp <- t_cpp$value
      reference <- if (!is.null(naive_res)) {
        naive_res$density
      } else if (!is.null(direct_cpp_res)) {
        direct_cpp_res$density
      } else {
        sp_r$density
      }
      err_cpp <- abs(sp_cpp$density - reference)
      add_record("LBFP_sparse_prefix_cpp", "cpp", t_cpp,
                 visited = sp_cpp$visited,
                 prefix_nodes = sp_cpp$prefix_nodes,
                 max_abs_error = max(err_cpp),
                 mean_abs_error = mean(err_cpp),
                 notes = if (is.null(naive_res)) "Compared with R sparse-prefix backend because naive was skipped." else "")
    } else {
      add_record("LBFP_sparse_prefix_cpp", "cpp", empty_time(),
                 notes = "C++ backend unavailable; R fallback was used.")
    }
  }

  add_speedups(do.call(rbind, records))
}

add_speedups <- function(df, baseline = "LBFP_naive_R") {
  if (!nrow(df)) return(df)
  if (!"speedup" %in% names(df)) df$speedup <- NA_real_
  if (!"speedup_vs_naive" %in% names(df)) df$speedup_vs_naive <- NA_real_
  if (!"speedup_vs_direct_cpp" %in% names(df)) df$speedup_vs_direct_cpp <- NA_real_
  if (!"elapsed_naive" %in% names(df)) df$elapsed_naive <- NA_real_
  if (!"elapsed_direct_cpp" %in% names(df)) df$elapsed_direct_cpp <- NA_real_
  if (!"all_equal_tolerance" %in% names(df)) df$all_equal_tolerance <- NA
  key_cols <- intersect(c("scenario", "n", "d", "bin_width_multiplier", "bins", "seed", "rep"), names(df))
  key_cols <- key_cols[vapply(df[key_cols], function(z) !all(is.na(z)), logical(1L))]
  if (!length(key_cols)) return(df)
  split_id <- interaction(df[key_cols], drop = TRUE, lex.order = TRUE)
  for (lev in levels(split_id)) {
    ii <- which(split_id == lev)
    base_idx <- ii[df$method[ii] == baseline & is.finite(df$elapsed[ii])]
    if (length(base_idx) == 1L) {
      base <- df$elapsed[base_idx]
      df$elapsed_naive[ii] <- base
      comparable <- grepl("^LBFP_sparse_prefix", df$method[ii]) |
        df$method[ii] == baseline
      passed <- is.na(df$all_equal_tolerance[ii]) | df$all_equal_tolerance[ii]
      ok <- comparable & passed & is.finite(df$elapsed[ii]) & df$elapsed[ii] > 0
      df$speedup[ii[ok]] <- base / df$elapsed[ii[ok]]
      df$speedup_vs_naive[ii[ok]] <- df$speedup[ii[ok]]
    }
    direct_idx <- ii[df$method[ii] == "LBFP_direct_cpp" & is.finite(df$elapsed[ii])]
    if (length(direct_idx) == 1L && df$elapsed[direct_idx] > 0) {
      direct_base <- df$elapsed[direct_idx]
      df$elapsed_direct_cpp[ii] <- direct_base
      comparable_cpp <- df$method[ii] %in% c("LBFP_direct_cpp", "LBFP_sparse_prefix_cpp")
      passed_cpp <- is.na(df$all_equal_tolerance[ii]) | df$all_equal_tolerance[ii]
      ok_cpp <- comparable_cpp & passed_cpp & is.finite(df$elapsed[ii]) & df$elapsed[ii] > 0
      df$speedup_vs_direct_cpp[ii[ok_cpp]] <- direct_base / df$elapsed[ii[ok_cpp]]
    }
  }
  df
}

write_csv_base <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

plot_heatmap_base <- function(df, value, file, main) {
  png(file, width = 1200, height = 800)
  on.exit(dev.off(), add = TRUE)
  op <- par(mar = c(7, 5, 4, 2))
  on.exit(par(op), add = TRUE)
  if (!nrow(df) || !value %in% names(df)) {
    plot.new(); title(main); text(0.5, 0.5, "No data")
    return(invisible(file))
  }
  df$label <- paste0(df$method, "\nn=", df$n)
  methods <- unique(df$label)
  dims <- sort(unique(df$d))
  mat <- matrix(NA_real_, nrow = length(methods), ncol = length(dims), dimnames = list(methods, dims))
  for (i in seq_len(nrow(df))) mat[df$label[i], as.character(df$d[i])] <- df[[value]][i]
  image(seq_along(dims), seq_along(methods), t(mat), axes = FALSE, xlab = "d", ylab = "method / n", main = main)
  axis(1, at = seq_along(dims), labels = dims)
  axis(2, at = seq_along(methods), labels = methods, las = 2, cex.axis = 0.75)
  invisible(file)
}
