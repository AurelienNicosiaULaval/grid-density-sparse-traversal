# Utilities for occupancy-controlled LBFP benchmarks.

controlled_side_length <- function(d) {
  lookup <- c(`4` = 10L, `6` = 6L, `8` = 4L, `10` = 3L, `12` = 3L)
  key <- as.character(as.integer(d))
  if (!key %in% names(lookup)) {
    stop("No controlled benchmark side length configured for d = ", d, call. = FALSE)
  }
  unname(lookup[[key]])
}

linear_indices_to_cells <- function(index, side_length, d) {
  index <- as.integer(index) - 1L
  out <- matrix(0L, nrow = length(index), ncol = d)
  for (s in seq_len(d)) {
    out[, s] <- index %% side_length
    index <- index %/% side_length
  }
  out
}

generate_uniform_cells <- function(d, side_length, target_occupancy_fraction, seed) {
  set.seed(seed)
  G <- side_length^d
  n_occ <- max(1L, min(as.integer(G), as.integer(round(target_occupancy_fraction * G))))
  idx <- sample.int(as.integer(G), n_occ, replace = FALSE)
  linear_indices_to_cells(idx, side_length = side_length, d = d)
}

generate_clustered_cells <- function(d, side_length, target_occupancy_fraction, seed,
                                     n_centers = NULL, cluster_probability = 0.72) {
  set.seed(seed)
  G <- side_length^d
  n_occ <- max(1L, min(as.integer(G), as.integer(round(target_occupancy_fraction * G))))
  if (is.null(n_centers)) n_centers <- max(3L, min(12L, ceiling(sqrt(n_occ))))
  centers <- matrix(sample.int(side_length, n_centers * d, replace = TRUE) - 1L,
                    nrow = n_centers, ncol = d)
  keys <- character()
  cells <- matrix(integer(), nrow = 0L, ncol = d)
  current_cluster_probability <- cluster_probability

  while (nrow(cells) < n_occ) {
    need <- n_occ - nrow(cells)
    batch <- max(1000L, min(250000L, need * 8L))
    center_id <- sample.int(n_centers, batch, replace = TRUE)
    center_cells <- centers[center_id, , drop = FALSE]
    random_cells <- matrix(sample.int(side_length, batch * d, replace = TRUE) - 1L,
                           nrow = batch, ncol = d)
    keep_center <- matrix(stats::runif(batch * d) < current_cluster_probability,
                          nrow = batch, ncol = d)
    candidate <- center_cells
    candidate[!keep_center] <- random_cells[!keep_center]
    candidate_keys <- key_from_index(candidate)
    is_new <- !candidate_keys %in% keys
    if (any(is_new)) {
      cells <- rbind(cells, candidate[is_new, , drop = FALSE])
      keys <- c(keys, candidate_keys[is_new])
    }
    if (length(keys) > n_occ) {
      cells <- cells[seq_len(n_occ), , drop = FALSE]
      keys <- keys[seq_len(n_occ)]
    }
    if (nrow(cells) < n_occ && batch == 250000L) {
      current_cluster_probability <- max(0.35, current_cluster_probability - 0.05)
    }
  }
  storage.mode(cells) <- "integer"
  cells
}

generate_controlled_cells <- function(d, side_length, target_occupancy_fraction,
                                      occupancy_pattern, seed) {
  occupancy_pattern <- match.arg(occupancy_pattern, c("uniform", "clustered"))
  if (occupancy_pattern == "uniform") {
    generate_uniform_cells(d, side_length, target_occupancy_fraction, seed)
  } else {
    generate_clustered_cells(d, side_length, target_occupancy_fraction, seed)
  }
}

make_controlled_lbfp_fit <- function(cells, side_length, counts = NULL) {
  cells <- as.matrix(cells)
  storage.mode(cells) <- "integer"
  d <- ncol(cells)
  if (is.null(counts)) counts <- rep.int(1L, nrow(cells))
  counts <- as.integer(counts)
  keys <- key_from_index(cells)
  grid <- list(
    origin = rep(0, d),
    bin_width = rep(1, d),
    d = d,
    volume = 1,
    max_index = rep(as.integer(side_length), d)
  )
  structure(list(
    x = cells + 0.5,
    grid = grid,
    counts = list(cells = cells, counts = counts, keys = keys, n_occ = length(counts)),
    count_env = make_count_environment(keys, counts),
    n = sum(counts),
    d = d,
    estimator = "LBFP"
  ), class = "grid_density_fit")
}

generate_controlled_queries <- function(cells, side_length, M, query_regime, seed) {
  query_regime <- match.arg(query_regime, c("observation_like", "random_grid"))
  set.seed(seed)
  d <- ncol(cells)
  M <- as.integer(M)
  if (query_regime == "observation_like") {
    base <- cells[sample.int(nrow(cells), M, replace = TRUE), , drop = FALSE]
  } else {
    base <- matrix(sample.int(side_length, M * d, replace = TRUE) - 1L,
                   nrow = M, ncol = d)
  }
  base + matrix(stats::runif(M * d, min = 0.05, max = 0.95), nrow = M, ncol = d)
}

controlled_timing_inner_reps <- function(d) {
  if (d <= 4L) return(50L)
  if (d <= 6L) return(25L)
  if (d <= 8L) return(10L)
  if (d <= 10L) return(5L)
  3L
}

controlled_config_grid <- function() {
  expand.grid(
    d = c(4L, 6L, 8L, 10L, 12L),
    target_occupancy_fraction = c(0.001, 0.005, 0.01, 0.05, 0.10),
    occupancy_pattern = c("uniform", "clustered"),
    query_regime = c("observation_like", "random_grid"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

median_na <- function(z) {
  if (all(is.na(z))) NA_real_ else stats::median(z, na.rm = TRUE)
}

accuracy_na <- function(z) {
  if (!length(z)) NA_real_ else mean(z, na.rm = TRUE)
}
