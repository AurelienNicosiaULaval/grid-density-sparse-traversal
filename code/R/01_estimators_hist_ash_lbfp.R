# 01_estimators_hist_ash_lbfp.R
# Reference implementations for histograms, ASH, LBFP and GLBFP-like grid estimators.

hist_fit <- function(x, bin_width = NULL, nbins = NULL, origin = NULL, grid = NULL) {
  x <- as_numeric_matrix(x)
  if (is.null(grid)) grid <- make_grid_spec(x, bin_width = bin_width, nbins = nbins, origin = origin)
  cnt <- count_grid_cells(x, grid)
  env <- make_count_environment(cnt$keys, cnt$counts)
  structure(list(
    x = x,
    grid = grid,
    counts = cnt,
    count_env = env,
    n = nrow(x),
    d = ncol(x),
    estimator = "histogram"
  ), class = "grid_density_fit")
}

hist_eval <- function(fit, query) {
  query <- as_numeric_matrix(query)
  idx <- cell_index(query, fit$grid)
  keys <- key_from_index(idx)
  counts <- lookup_count(fit$count_env, keys)
  counts / (fit$n * fit$grid$volume)
}

lbfp_fit <- function(x, bin_width = NULL, nbins = NULL, origin = NULL, grid = NULL) {
  fit <- hist_fit(x, bin_width = bin_width, nbins = nbins, origin = origin, grid = grid)
  fit$estimator <- "LBFP"
  fit
}

lbfp_eval_naive <- function(fit, query, return_visited = FALSE) {
  query <- as_numeric_matrix(query)
  if (fit$d > 30L) stop("Naive LBFP stencil is too large for d > 30.")
  idx <- cell_index(query, fit$grid)
  u <- relative_coord(query, fit$grid, idx)
  stencil <- binary_stencil(fit$d)
  out <- numeric(nrow(query))
  visited <- integer(nrow(query))
  for (i in seq_len(nrow(query))) {
    acc <- 0
    vi <- 0L
    for (r in seq_len(nrow(stencil))) {
      j <- stencil[r, ]
      w <- prod(ifelse(j == 1L, u[i, ], 1 - u[i, ]))
      key <- paste(idx[i, ] + j, collapse = "|")
      cc <- lookup_count(fit$count_env, key)
      if (cc != 0L) vi <- vi + 1L
      acc <- acc + w * cc
    }
    out[i] <- acc / (fit$n * fit$grid$volume)
    visited[i] <- vi
  }
  if (return_visited) list(density = out, visited = visited, nominal_stencil = 2^fit$d) else out
}

ash_fit <- function(x, bin_width = NULL, m = 4L, nbins = NULL, origin = NULL) {
  x <- as_numeric_matrix(x)
  d <- ncol(x)
  if (length(m) == 1L) m <- rep(as.integer(m), d)
  if (length(m) != d || any(m < 1L)) stop("m must have length 1 or d and positive entries.")
  coarse <- make_grid_spec(x, bin_width = bin_width, nbins = nbins, origin = origin)
  fine <- make_grid_spec(x, bin_width = coarse$bin_width / m, origin = coarse$origin)
  cnt <- count_grid_cells(x, fine)
  env <- make_count_environment(cnt$keys, cnt$counts)
  structure(list(
    x = x,
    coarse_grid = coarse,
    grid = fine,
    m = m,
    counts = cnt,
    count_env = env,
    n = nrow(x),
    d = d,
    estimator = "ASH"
  ), class = "grid_density_fit")
}

ash_eval_naive <- function(fit, query, return_visited = FALSE) {
  query <- as_numeric_matrix(query)
  lower <- 1L - fit$m
  upper <- fit$m - 1L
  offsets <- offset_grid(lower, upper)
  idx <- cell_index(query, fit$grid)
  out <- numeric(nrow(query))
  visited <- integer(nrow(query))
  for (i in seq_len(nrow(query))) {
    acc <- 0
    vi <- 0L
    for (r in seq_len(nrow(offsets))) {
      ell <- offsets[r, ]
      w <- prod(pmax(0, 1 - abs(ell) / fit$m))
      if (w == 0) next
      key <- paste(idx[i, ] + ell, collapse = "|")
      cc <- lookup_count(fit$count_env, key)
      if (cc != 0L) vi <- vi + 1L
      acc <- acc + w * cc
    }
    out[i] <- acc / (fit$n * fit$coarse_grid$volume)
    visited[i] <- vi
  }
  if (return_visited) list(density = out, visited = visited, nominal_stencil = nrow(offsets)) else out
}

glbfp_fit <- function(x, bin_width = NULL, m = 4L, nbins = NULL, origin = NULL) {
  # Product-form GLBFP/FP-ASH reference implementation.
  x <- as_numeric_matrix(x)
  d <- ncol(x)
  if (length(m) == 1L) m <- rep(as.integer(m), d)
  if (length(m) != d || any(m < 1L)) stop("m must have length 1 or d and positive entries.")
  coarse <- make_grid_spec(x, bin_width = bin_width, nbins = nbins, origin = origin)
  fine <- make_grid_spec(x, bin_width = coarse$bin_width / m, origin = coarse$origin)
  cnt <- count_grid_cells(x, fine)
  env <- make_count_environment(cnt$keys, cnt$counts)
  structure(list(
    x = x,
    coarse_grid = coarse,
    grid = fine,
    m = m,
    counts = cnt,
    count_env = env,
    n = nrow(x),
    d = d,
    estimator = "GLBFP"
  ), class = "grid_density_fit")
}

glbfp_weights_1d <- function(u, m) {
  ell <- seq.int(1L - m, m)
  w <- (1 - u) * pmax(0, 1 - abs(ell) / m) + u * pmax(0, 1 - abs(ell - 1L) / m)
  list(offset = ell, weight = w)
}

glbfp_eval_naive <- function(fit, query, return_visited = FALSE, max_stencil = 5e6) {
  query <- as_numeric_matrix(query)
  nominal <- prod(2L * fit$m)
  if (nominal > max_stencil) stop("GLBFP stencil too large for naive evaluation: ", nominal)
  lower <- 1L - fit$m
  upper <- fit$m
  offsets <- offset_grid(lower, upper)
  idx <- cell_index(query, fit$grid)
  u <- relative_coord(query, fit$grid, idx)
  out <- numeric(nrow(query))
  visited <- integer(nrow(query))
  for (i in seq_len(nrow(query))) {
    # Per-dimension weights indexed by offset value.
    wlist <- lapply(seq_len(fit$d), function(s) glbfp_weights_1d(u[i, s], fit$m[s]))
    acc <- 0
    vi <- 0L
    for (r in seq_len(nrow(offsets))) {
      ell <- offsets[r, ]
      w <- 1
      for (s in seq_len(fit$d)) {
        pos <- match(ell[s], wlist[[s]]$offset)
        w <- w * wlist[[s]]$weight[pos]
      }
      if (w == 0) next
      key <- paste(idx[i, ] + ell, collapse = "|")
      cc <- lookup_count(fit$count_env, key)
      if (cc != 0L) vi <- vi + 1L
      acc <- acc + w * cc
    }
    out[i] <- acc / (fit$n * fit$coarse_grid$volume)
    visited[i] <- vi
  }
  if (return_visited) list(density = out, visited = visited, nominal_stencil = nominal) else out
}
