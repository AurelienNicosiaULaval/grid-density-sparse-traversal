# 00_helpers.R
# General utilities for grid-based density estimators.

as_numeric_matrix <- function(x, name = "x") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.vector(x)) x <- matrix(as.numeric(x), ncol = 1L)
  if (!is.matrix(x)) stop(name, " must be a numeric matrix, data frame, or vector.")
  storage.mode(x) <- "double"
  if (!is.numeric(x)) stop(name, " must be numeric.")
  if (nrow(x) < 1L || ncol(x) < 1L) stop(name, " must have at least one row and one column.")
  if (anyNA(x)) stop(name, " contains missing values.")
  if (any(!is.finite(x))) stop(name, " contains non-finite values.")
  x
}

robust_scale <- function(x) {
  x <- as_numeric_matrix(x)
  s <- apply(x, 2L, stats::sd)
  iq <- apply(x, 2L, stats::IQR) / 1.349
  out <- pmin(s, iq, na.rm = TRUE)
  out[!is.finite(out) | out <= 0] <- s[!is.finite(out) | out <= 0]
  out[!is.finite(out) | out <= 0] <- 1
  out
}

default_bin_width <- function(x, constant = 2.15, robust = TRUE) {
  x <- as_numeric_matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  scale <- if (robust) robust_scale(x) else apply(x, 2L, stats::sd)
  scale[!is.finite(scale) | scale <= 0] <- 1
  constant * scale * n^(-1 / (d + 4))
}

make_grid_spec <- function(x, bin_width = NULL, nbins = NULL, origin = NULL,
                           padding = 1e-8, robust = TRUE) {
  x <- as_numeric_matrix(x)
  d <- ncol(x)
  rng <- apply(x, 2L, range)
  if (is.null(bin_width)) {
    if (!is.null(nbins)) {
      if (length(nbins) == 1L) nbins <- rep(nbins, d)
      if (length(nbins) != d) stop("nbins must have length 1 or d.")
      if (any(!is.finite(nbins)) || any(nbins < 2L)) stop("nbins entries must be finite and at least 2.")
      span <- pmax(rng[2L, ] - rng[1L, ], .Machine$double.eps)
      bin_width <- span / pmax(as.integer(nbins) - 1L, 1L)
    } else {
      bin_width <- default_bin_width(x, robust = robust)
    }
  }
  if (length(bin_width) == 1L) bin_width <- rep(bin_width, d)
  if (length(bin_width) != d) stop("bin_width must have length 1 or d.")
  if (any(!is.finite(bin_width)) || any(bin_width <= 0)) stop("Invalid bin_width.")
  if (is.null(origin)) origin <- rng[1L, ] - padding * pmax(1, abs(rng[1L, ]))
  if (length(origin) == 1L) origin <- rep(origin, d)
  if (length(origin) != d) stop("origin must have length 1 or d.")
  if (any(!is.finite(origin))) stop("origin must contain finite values.")
  max_index <- floor((rng[2L, ] - origin) / bin_width) + 2L
  list(
    origin = as.numeric(origin),
    bin_width = as.numeric(bin_width),
    d = d,
    volume = prod(bin_width),
    max_index = as.integer(max_index)
  )
}

cell_index <- function(x, grid) {
  x <- as_numeric_matrix(x)
  if (is.null(grid$d) || is.null(grid$origin) || is.null(grid$bin_width)) {
    stop("grid must contain d, origin, and bin_width.")
  }
  if (ncol(x) != grid$d) stop("Dimension mismatch in cell_index().")
  if (length(grid$origin) != grid$d || length(grid$bin_width) != grid$d) {
    stop("grid origin and bin_width must have length d.")
  }
  idx <- floor(sweep(sweep(x, 2L, grid$origin, "-"), 2L, grid$bin_width, "/"))
  storage.mode(idx) <- "integer"
  idx
}

relative_coord <- function(x, grid, idx = NULL) {
  x <- as_numeric_matrix(x)
  if (is.null(idx)) idx <- cell_index(x, grid)
  u <- sweep(sweep(x, 2L, grid$origin + 0, "-"), 2L, grid$bin_width, "/") - idx
  u[u < 0 & u > -1e-12] <- 0
  u[u >= 1 & u < 1 + 1e-12] <- 1 - 1e-12
  pmin(pmax(u, 0), 1 - 1e-12)
}

key_from_index <- function(idx, sep = "|") {
  idx <- as.matrix(idx)
  if (ncol(idx) == 1L) return(as.character(idx[, 1L]))
  apply(idx, 1L, paste, collapse = sep)
}

split_key <- function(keys, sep = "|", d = NULL) {
  parts <- strsplit(keys, sep, fixed = TRUE)
  if (is.null(d)) d <- length(parts[[1L]])
  out <- matrix(NA_integer_, nrow = length(keys), ncol = d)
  for (i in seq_along(parts)) out[i, ] <- as.integer(parts[[i]])
  out
}

binary_stencil <- function(d) {
  if (d < 1L) stop("d must be positive.")
  if (d > 30L) stop("Refusing to create binary stencil for d > 30.")
  as.matrix(expand.grid(rep(list(0:1), d), KEEP.OUT.ATTRS = FALSE))
}

offset_grid <- function(lower, upper) {
  if (length(lower) != length(upper)) stop("lower and upper must have same length.")
  as.matrix(expand.grid(Map(seq.int, lower, upper), KEEP.OUT.ATTRS = FALSE))
}

count_grid_cells <- function(x, grid) {
  idx <- cell_index(x, grid)
  keys <- key_from_index(idx)
  tab <- table(keys)
  counts <- as.integer(tab)
  key_names <- names(tab)
  cell_mat <- split_key(key_names, d = grid$d)
  list(cells = cell_mat, counts = counts, keys = key_names, n_occ = length(counts))
}

grid_total_cells <- function(grid) {
  if (is.null(grid$max_index)) return(NA_real_)
  m <- as.numeric(grid$max_index)
  if (!length(m) || any(!is.finite(m)) || any(m <= 0)) return(NA_real_)
  out <- prod(m)
  if (!is.finite(out)) Inf else out
}

grid_occupancy_metrics <- function(fit) {
  G <- grid_total_cells(fit$grid)
  G_occ <- fit$counts$n_occ
  occupancy_rate <- if (is.finite(G) && G > 0) G_occ / G else NA_real_
  list(G = G, G_occ = G_occ, occupancy_rate = occupancy_rate)
}

ensure_project_dirs <- function(root = getwd()) {
  dirs <- file.path(root, c("code/results", "code/figures", "code/logs"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

make_count_environment <- function(keys, counts) {
  env <- new.env(parent = emptyenv(), hash = TRUE, size = max(29L, length(keys) * 2L))
  for (i in seq_along(keys)) assign(keys[[i]], counts[[i]], envir = env)
  env
}

lookup_count <- function(env, key) {
  val <- mget(key, envir = env, ifnotfound = list(0L), inherits = FALSE)
  as.integer(unlist(val, use.names = FALSE))
}

safe_density_divisor <- function(x, eps = .Machine$double.eps) {
  x[!is.finite(x) | x < eps] <- eps
  x
}

elapsed <- function(expr) {
  gc(FALSE)
  t <- system.time(value <- force(expr))
  list(value = value, elapsed = unname(t[["elapsed"]]), user = unname(t[["user.self"]]), system = unname(t[["sys.self"]]))
}

source_project_R <- function(root = getwd()) {
  files <- c(
    "code/R/00_helpers.R",
    "code/R/01_estimators_hist_ash_lbfp.R",
    "code/R/02_loo_scores.R",
    "code/R/03_sparse_prefix_lbfp.R",
    "code/R/04_benchmark_utils.R"
  )
  for (f in files) source(file.path(root, f), chdir = TRUE)
  invisible(TRUE)
}

write_session_info <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  if (requireNamespace("sessioninfo", quietly = TRUE)) {
    writeLines(capture.output(sessioninfo::session_info()), con)
  } else {
    writeLines(capture.output(sessionInfo()), con)
  }
  invisible(path)
}

project_root_from_script <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(script) == 0L) return(normalizePath(getwd()))
  normalizePath(file.path(dirname(script), "../.."))
}

parse_flag <- function(flag, default = FALSE) {
  flag %in% commandArgs(trailingOnly = TRUE) || default
}

parse_arg_value <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  pat <- paste0("^", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (!length(hit)) return(default)
  sub(pat, "", hit[[1L]])
}
