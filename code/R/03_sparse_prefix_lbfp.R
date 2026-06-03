# 03_sparse_prefix_lbfp.R
# Exact sparse-prefix LBFP evaluation.

choose_prefix_order <- function(cells, strategy = c("adaptive", "natural", "reverse", "random"), seed = 1L) {
  cells <- as.matrix(cells)
  d <- ncol(cells)
  if (d < 1L || nrow(cells) < 1L) stop("cells must have at least one row and one column.")
  strategy <- match.arg(strategy)
  if (strategy == "natural") return(seq_len(d))
  if (strategy == "reverse") return(rev(seq_len(d)))
  if (strategy == "random") {
    set.seed(seed)
    return(sample.int(d))
  }
  # Adaptive heuristic: put dimensions with fewer distinct occupied coordinates first.
  # These dimensions tend to reject impossible prefixes earlier.
  n_unique <- apply(cells, 2L, function(z) length(unique(z)))
  order(n_unique, decreasing = FALSE)
}

validate_prefix_order <- function(order, d) {
  if (length(order) != d) stop("order must have length d.")
  if (anyNA(order) || any(!is.finite(order))) stop("order must contain finite integer positions.")
  order <- as.integer(order)
  if (!setequal(order, seq_len(d)) || anyDuplicated(order)) {
    stop("order must be a permutation of 1:d.")
  }
  order
}

build_prefix_index <- function(cells, counts, order = NULL, sep = "|") {
  cells <- as.matrix(cells)
  storage.mode(cells) <- "integer"
  d <- ncol(cells)
  if (d < 1L || nrow(cells) < 1L) stop("cells must have at least one row and one column.")
  if (length(counts) != nrow(cells)) stop("counts length must equal number of cell rows.")
  if (anyNA(counts) || any(counts < 0L)) stop("counts must be non-negative integers.")
  if (is.null(order)) order <- choose_prefix_order(cells, "adaptive")
  order <- validate_prefix_order(order, d)
  ordered_cells <- cells[, order, drop = FALSE]
  prefix_envs <- vector("list", d)
  for (r in seq_len(d)) {
    keys <- key_from_index(ordered_cells[, seq_len(r), drop = FALSE], sep = sep)
    env <- new.env(parent = emptyenv(), hash = TRUE, size = max(29L, length(keys) * 2L))
    for (k in unique(keys)) assign(k, TRUE, envir = env)
    prefix_envs[[r]] <- env
  }
  full_keys <- key_from_index(ordered_cells, sep = sep)
  count_env <- make_count_environment(full_keys, counts)
  list(prefix_envs = prefix_envs, count_env = count_env, order = as.integer(order), sep = sep, d = d)
}

prefix_exists <- function(prefix_index, key, depth) {
  exists(key, envir = prefix_index$prefix_envs[[depth]], inherits = FALSE)
}

sparse_prefix_nominal_stencil <- function(fit) {
  estimator <- tolower(fit$estimator)
  if (estimator == "histogram") return(1)
  if (estimator == "lbfp") return(2^fit$d)
  if (estimator == "ash") return(prod(2L * fit$m - 1L))
  if (estimator == "glbfp") return(prod(2L * fit$m))
  stop("Unsupported estimator: ", fit$estimator)
}

sparse_prefix_denominator_volume <- function(fit) {
  estimator <- tolower(fit$estimator)
  if (estimator %in% c("histogram", "lbfp")) return(fit$grid$volume)
  if (estimator %in% c("ash", "glbfp")) return(fit$coarse_grid$volume)
  stop("Unsupported estimator: ", fit$estimator)
}

sparse_prefix_candidates <- function(fit, i, idx, u = NULL) {
  estimator <- tolower(fit$estimator)
  d <- fit$d

  if (estimator == "histogram") {
    return(lapply(seq_len(d), function(s) list(offset = 0L, weight = 1)))
  }

  if (estimator == "lbfp") {
    return(lapply(seq_len(d), function(s) {
      list(offset = c(0L, 1L), weight = c(1 - u[i, s], u[i, s]))
    }))
  }

  if (estimator == "ash") {
    return(lapply(seq_len(d), function(s) {
      ell <- seq.int(1L - fit$m[s], fit$m[s] - 1L)
      list(offset = ell, weight = pmax(0, 1 - abs(ell) / fit$m[s]))
    }))
  }

  if (estimator == "glbfp") {
    return(lapply(seq_len(d), function(s) {
      glbfp_weights_1d(u[i, s], fit$m[s])
    }))
  }

  stop("Unsupported estimator: ", fit$estimator)
}

sparse_prefix_eval_local_grid_R <- function(fit, query = fit$x, order = NULL,
                                            return_visited = FALSE,
                                            prefix_index = NULL) {
  query <- as_numeric_matrix(query)
  d <- fit$d
  if (ncol(query) != d) stop("query must have the same number of columns as fit$x.")
  if (is.null(prefix_index)) {
    if (is.null(order)) order <- choose_prefix_order(fit$counts$cells, "adaptive")
    prefix_index <- build_prefix_index(fit$counts$cells, fit$counts$counts, order = order)
  }
  order <- validate_prefix_order(prefix_index$order, d)
  idx <- cell_index(query, fit$grid)
  estimator <- tolower(fit$estimator)
  u <- if (estimator %in% c("lbfp", "glbfp")) relative_coord(query, fit$grid, idx) else NULL
  denom <- fit$n * sparse_prefix_denominator_volume(fit)
  nominal_stencil <- sparse_prefix_nominal_stencil(fit)
  dens <- numeric(nrow(query))
  visited <- integer(nrow(query))
  prefix_nodes <- integer(nrow(query))

  for (i in seq_len(nrow(query))) {
    candidates <- sparse_prefix_candidates(fit, i, idx, u)
    acc <- 0
    vi <- 0L
    nodes <- 0L

    rec <- function(depth, key, weight) {
      dim <- order[[depth]]
      cand <- candidates[[dim]]
      for (j in seq_along(cand$offset)) {
        w <- weight * cand$weight[[j]]
        if (w == 0) next
        coord <- idx[i, dim] + cand$offset[[j]]
        new_key <- if (depth == 1L) as.character(coord) else paste0(key, prefix_index$sep, coord)
        nodes <<- nodes + 1L
        if (!prefix_exists(prefix_index, new_key, depth)) next
        if (depth == d) {
          cc <- lookup_count(prefix_index$count_env, new_key)
          if (cc != 0L) {
            vi <<- vi + 1L
            acc <<- acc + w * cc
          }
        } else {
          rec(depth + 1L, new_key, w)
        }
      }
      invisible(NULL)
    }

    rec(1L, "", 1)
    dens[[i]] <- acc / denom
    visited[[i]] <- vi
    prefix_nodes[[i]] <- nodes
  }

  if (return_visited) {
    list(density = dens, visited = visited, prefix_nodes = prefix_nodes,
         nominal_stencil = nominal_stencil, prefix_index = prefix_index)
  } else {
    dens
  }
}

sparse_prefix_eval_lbfp_R <- function(fit, query = fit$x, order = NULL,
                                      return_visited = FALSE, prefix_index = NULL) {
  query <- as_numeric_matrix(query)
  d <- fit$d
  if (ncol(query) != d) stop("query must have the same number of columns as fit$x.")
  if (is.null(prefix_index)) {
    if (is.null(order)) order <- choose_prefix_order(fit$counts$cells, "adaptive")
    prefix_index <- build_prefix_index(fit$counts$cells, fit$counts$counts, order = order)
  }
  order <- prefix_index$order
  order <- validate_prefix_order(order, d)
  idx <- cell_index(query, fit$grid)
  u <- relative_coord(query, fit$grid, idx)
  dens <- numeric(nrow(query))
  visited <- integer(nrow(query))
  prefix_nodes <- integer(nrow(query))

  for (i in seq_len(nrow(query))) {
    acc <- 0
    vi <- 0L
    nodes <- 0L

    rec <- function(depth, key, weight) {
      dim <- order[[depth]]
      for (bit in 0:1) {
        coord <- idx[i, dim] + bit
        w <- weight * if (bit == 1L) u[i, dim] else (1 - u[i, dim])
        new_key <- if (depth == 1L) as.character(coord) else paste0(key, prefix_index$sep, coord)
        nodes <<- nodes + 1L
        if (!prefix_exists(prefix_index, new_key, depth)) next
        if (depth == d) {
          cc <- lookup_count(prefix_index$count_env, new_key)
          if (cc != 0L) {
            vi <<- vi + 1L
            acc <<- acc + w * cc
          }
        } else {
          rec(depth + 1L, new_key, w)
        }
      }
      invisible(NULL)
    }

    rec(1L, "", 1)
    dens[i] <- acc / (fit$n * fit$grid$volume)
    visited[i] <- vi
    prefix_nodes[i] <- nodes
  }

  if (return_visited) {
    list(density = dens, visited = visited, prefix_nodes = prefix_nodes,
         nominal_stencil = 2^d, prefix_index = prefix_index)
  } else {
    dens
  }
}

loo_lbfp_scores_sparse_prefix_R <- function(fit, query = fit$x, order = NULL,
                                            return_details = FALSE, prefix_index = NULL) {
  ev <- sparse_prefix_eval_lbfp_R(fit, query = query, order = order,
                                  return_visited = TRUE, prefix_index = prefix_index)
  w <- lbfp_self_weight(fit, query)
  D <- loo_from_density(ev$density, w, fit$grid$volume, fit$n)
  if (return_details) {
    ev$score <- D
    ev$self_weight <- w
    ev
  } else D
}

compile_sparse_prefix_cpp <- function(root = getwd(), verbose = TRUE) {
  cpp <- file.path(root, "code/src/sparse_prefix_lbfp.cpp")
  if (!file.exists(cpp)) {
    if (verbose) message("C++ source file not found; using R fallback.")
    return(FALSE)
  }
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    if (verbose) message("Rcpp is not installed; using R fallback.")
    return(FALSE)
  }
  if (!requireNamespace("RcppParallel", quietly = TRUE)) {
    if (verbose) message("RcppParallel is not installed; using R fallback.")
    return(FALSE)
  }
  tryCatch({
    Rcpp::sourceCpp(cpp, rebuild = FALSE, verbose = FALSE)
    TRUE
  }, error = function(e) {
    if (verbose) message("C++ backend could not be compiled: ", conditionMessage(e))
    FALSE
  })
}

sparse_prefix_eval_lbfp_cpp <- function(fit, query = fit$x, order = NULL,
                                        return_visited = FALSE) {
  if (!exists("lbfp_sparse_prefix_cpp", mode = "function")) {
    ok <- compile_sparse_prefix_cpp(project_root_from_script(), verbose = FALSE)
    if (!ok || !exists("lbfp_sparse_prefix_cpp", mode = "function")) {
      stop("C++ backend not available. Call sparse_prefix_eval_lbfp_R() instead.")
    }
  }
  if (is.null(order)) order <- choose_prefix_order(fit$counts$cells, "adaptive")
  query <- as_numeric_matrix(query)
  if (ncol(query) != fit$d) stop("query must have the same number of columns as fit$x.")
  order <- validate_prefix_order(order, fit$d)
  res <- lbfp_sparse_prefix_cpp(
    query,
    as.numeric(fit$grid$origin),
    as.numeric(fit$grid$bin_width),
    as.integer(fit$n),
    matrix(as.integer(fit$counts$cells), ncol = fit$d),
    as.integer(fit$counts$counts),
    as.integer(order)
  )
  if (return_visited) res else res$density
}

lbfp_eval_direct_cpp <- function(fit, query = fit$x, return_visited = FALSE) {
  if (!exists("lbfp_direct_cpp", mode = "function")) {
    ok <- compile_sparse_prefix_cpp(project_root_from_script(), verbose = FALSE)
    if (!ok || !exists("lbfp_direct_cpp", mode = "function")) {
      stop("Direct C++ backend not available. Call lbfp_eval_naive() instead.")
    }
  }
  query <- as_numeric_matrix(query)
  if (ncol(query) != fit$d) stop("query must have the same number of columns as fit$x.")
  res <- lbfp_direct_cpp(
    query,
    as.numeric(fit$grid$origin),
    as.numeric(fit$grid$bin_width),
    as.integer(fit$n),
    matrix(as.integer(fit$counts$cells), ncol = fit$d),
    as.integer(fit$counts$counts)
  )
  if (return_visited) res else res$density
}

loo_lbfp_scores_sparse_prefix_cpp <- function(fit, query = fit$x, order = NULL,
                                              return_details = FALSE) {
  ev <- sparse_prefix_eval_lbfp_cpp(fit, query = query, order = order, return_visited = TRUE)
  w <- lbfp_self_weight(fit, query)
  D <- loo_from_density(ev$density, w, fit$grid$volume, fit$n)
  if (return_details) {
    ev$score <- D
    ev$self_weight <- w
    ev
  } else D
}
