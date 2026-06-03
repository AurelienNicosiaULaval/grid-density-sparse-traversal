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
  ns <- c(250L, 500L, 1000L, 2000L)
  ds <- 2L:6L
  reps <- 1L
  scenarios <- c("normal")
  m_value <- 2L
  timing_reps <- 1L
  short_timing_inner_reps <- 1L
  direct_budget <- 3e6
} else {
  ns <- c(500L, 1000L, 2000L)
  ds <- 2L:7L
  reps <- 3L
  scenarios <- c("normal", "mixture", "clustered")
  m_value <- 2L
  timing_reps <- as.integer(parse_arg_value("--timing-reps", default = "2"))
  short_timing_inner_reps <- 1L
  direct_budget <- 2e7
}

if (is.na(timing_reps) || timing_reps < 1L) stop("--timing-reps must be a positive integer.")

fit_estimator_family <- function(x, estimator, bin_width, m_value) {
  if (estimator == "histogram") return(hist_fit(x, bin_width = bin_width))
  if (estimator == "ASH") return(ash_fit(x, bin_width = bin_width, m = m_value))
  if (estimator == "LBFP") return(lbfp_fit(x, bin_width = bin_width))
  if (estimator == "GLBFP") return(glbfp_fit(x, bin_width = bin_width, m = m_value))
  stop("Unsupported estimator: ", estimator)
}

direct_eval_family <- function(fit, estimator, query) {
  if (estimator == "histogram") {
    dens <- hist_eval(fit, query)
    return(list(density = dens, visited = rep(1L, length(dens)), nominal_stencil = 1L))
  }
  if (estimator == "ASH") return(ash_eval_naive(fit, query, return_visited = TRUE))
  if (estimator == "LBFP") return(lbfp_eval_naive(fit, query, return_visited = TRUE))
  if (estimator == "GLBFP") return(glbfp_eval_naive(fit, query, return_visited = TRUE, max_stencil = direct_budget))
  stop("Unsupported estimator: ", estimator)
}

records <- list()

add_record <- function(config, fit, direct_time, sparse_time = NULL, sparse_res = NULL) {
  metrics <- grid_occupancy_metrics(fit)
  direct_res <- direct_time$value
  has_sparse <- !is.null(sparse_time) && !is.null(sparse_res)
  err <- if (has_sparse) abs(direct_res$density - sparse_res$density) else NA_real_
  records[[length(records) + 1L]] <<- data.frame(
    scenario = config$scenario,
    n = config$n,
    d = config$d,
    rep = config$rep,
    seed = config$seed,
    estimator = config$estimator,
    m = config$m,
    G = metrics$G,
    G_occ = metrics$G_occ,
    occupancy_rate = metrics$occupancy_rate,
    nominal_stencil = sparse_prefix_nominal_stencil(fit),
    direct_elapsed = direct_time$elapsed,
    sparse_elapsed = if (has_sparse) sparse_time$elapsed else NA_real_,
    speedup_sparse_vs_direct = if (has_sparse) direct_time$elapsed / sparse_time$elapsed else NA_real_,
    direct_visited_median = stats::median(direct_res$visited),
    sparse_visited_median = if (has_sparse) stats::median(sparse_res$visited) else NA_real_,
    prefix_nodes_median = if (has_sparse) stats::median(sparse_res$prefix_nodes) else NA_real_,
    visited_ratio_median = if (has_sparse) stats::median(sparse_res$visited) / sparse_prefix_nominal_stencil(fit) else NA_real_,
    prefix_ratio_median = if (has_sparse) stats::median(sparse_res$prefix_nodes) / sparse_prefix_nominal_stencil(fit) else NA_real_,
    max_abs_error = if (has_sparse) max(err) else NA_real_,
    mean_abs_error = if (has_sparse) mean(err) else NA_real_,
    all_equal_tolerance = if (has_sparse) max(err) <= 1e-10 else NA,
    timing_reps = direct_time$timing_reps,
    timing_inner_reps = direct_time$timing_inner_reps,
    notes = if (config$estimator == "histogram") "Histogram is direct lookup only." else "",
    stringsAsFactors = FALSE
  )
}

estimators <- c("histogram", "ASH", "LBFP", "GLBFP")

for (scenario in scenarios) {
  for (n in ns) {
    for (d in ds) {
      for (rep in seq_len(reps)) {
        seed <- 20260601L + 100000L * rep + 1000L * d + n
        message("Estimator-family benchmark scenario=", scenario, " n=", n,
                " d=", d, " rep=", rep)
        x <- simulate_grid_data(n = n, d = d, scenario = scenario, seed = seed)
        bin_width <- default_bin_width(x)

        for (estimator in estimators) {
          fit <- fit_estimator_family(x, estimator, bin_width = bin_width, m_value = m_value)
          nominal <- sparse_prefix_nominal_stencil(fit)
          if (nrow(x) * nominal > direct_budget) {
            message("Skipping ", estimator, " because n * nominal_stencil exceeds direct_budget.")
            next
          }

          config <- list(
            scenario = scenario,
            n = n,
            d = d,
            rep = rep,
            seed = seed,
            estimator = estimator,
            m = if (estimator %in% c("ASH", "GLBFP")) m_value else NA_integer_
          )

          direct_time <- elapsed_repeated(
            direct_eval_family(fit, estimator, x),
            reps = timing_reps,
            inner_reps = short_timing_inner_reps
          )

          if (estimator == "histogram") {
            add_record(config, fit, direct_time)
          } else {
            order <- choose_prefix_order(fit$counts$cells, "adaptive")
            sparse_time <- elapsed_repeated(
              sparse_prefix_eval_local_grid_R(fit, x, order = order, return_visited = TRUE),
              reps = timing_reps,
              inner_reps = short_timing_inner_reps
            )
            add_record(config, fit, direct_time, sparse_time, sparse_time$value)
          }
        }
      }
    }
  }
}

res <- do.call(rbind, records)
write_csv_base(res, file.path(out_dir, "estimator_family_benchmarks.csv"))
write_session_info(file.path(out_dir, "session_info_estimator_family.txt"))
print(utils::head(res, 20))
message("Wrote ", file.path(out_dir, "estimator_family_benchmarks.csv"))
