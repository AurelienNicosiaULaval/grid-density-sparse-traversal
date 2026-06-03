# 02_loo_scores.R
# Leave-one-out scores for linear grid estimators.

loo_from_density <- function(f_full, self_weight, Vb, n, eps = .Machine$double.eps) {
  f_safe <- safe_density_divisor(as.numeric(f_full), eps = eps)
  (self_weight / (Vb * f_safe) - 1) / (n - 1)
}

lbfp_self_weight <- function(fit, query) {
  query <- as_numeric_matrix(query)
  idx <- cell_index(query, fit$grid)
  u <- relative_coord(query, fit$grid, idx)
  apply(1 - u, 1L, prod)
}

glbfp_self_weight <- function(fit, query) {
  query <- as_numeric_matrix(query)
  idx <- cell_index(query, fit$grid)
  u <- relative_coord(query, fit$grid, idx)
  m <- fit$m
  apply(1 - sweep(pmin(u, 1 - u), 2L, m, "/"), 1L, prod)
}

loo_hist_scores <- function(fit, query = fit$x) {
  f <- hist_eval(fit, query)
  # Self-weight for histogram is 1 for the cell containing the observation.
  w <- rep(1, nrow(as_numeric_matrix(query)))
  loo_from_density(f, w, fit$grid$volume, fit$n)
}

loo_lbfp_scores_naive <- function(fit, query = fit$x) {
  f <- lbfp_eval_naive(fit, query)
  w <- lbfp_self_weight(fit, query)
  loo_from_density(f, w, fit$grid$volume, fit$n)
}

loo_glbfp_scores_naive <- function(fit, query = fit$x, max_stencil = 5e6) {
  f <- glbfp_eval_naive(fit, query, max_stencil = max_stencil)
  w <- glbfp_self_weight(fit, query)
  loo_from_density(f, w, fit$coarse_grid$volume, fit$n)
}

loo_lbfp_recompute_naive <- function(fit, indices = seq_len(min(20L, fit$n))) {
  # Explicitly recompute LOO density by decrementing the cell containing X_i.
  # This is only for correctness checks on small subsets.
  x <- fit$x
  f_full <- lbfp_eval_naive(fit, x[indices, , drop = FALSE])
  w <- lbfp_self_weight(fit, x[indices, , drop = FALSE])
  d_formula <- loo_from_density(f_full, w, fit$grid$volume, fit$n)
  d_recompute <- numeric(length(indices))

  for (ii in seq_along(indices)) {
    i <- indices[ii]
    xi <- x[i, , drop = FALSE]
    idx_i <- cell_index(xi, fit$grid)
    key_i <- key_from_index(idx_i)

    # Copy count environment for clarity; this is intentionally slow.
    tmp_env <- new.env(parent = emptyenv(), hash = TRUE, size = length(fit$counts$keys) * 2L)
    for (r in seq_along(fit$counts$keys)) assign(fit$counts$keys[[r]], fit$counts$counts[[r]], envir = tmp_env)
    old <- lookup_count(tmp_env, key_i)
    assign(key_i, as.integer(old - 1L), envir = tmp_env)
    tmp_fit <- fit
    tmp_fit$count_env <- tmp_env
    tmp_fit$n <- fit$n - 1L
    f_minus <- lbfp_eval_naive(tmp_fit, xi)
    d_recompute[ii] <- 1 - f_minus / f_full[ii]
  }

  data.frame(
    index = indices,
    D_formula = d_formula,
    D_recompute = d_recompute,
    abs_error = abs(d_formula - d_recompute)
  )
}

glbfp_loo_recompute_naive <- function(fit, indices = seq_len(min(10L, fit$n)), max_stencil = 1e6) {
  x <- fit$x
  f_full <- glbfp_eval_naive(fit, x[indices, , drop = FALSE], max_stencil = max_stencil)
  w <- glbfp_self_weight(fit, x[indices, , drop = FALSE])
  d_formula <- loo_from_density(f_full, w, fit$coarse_grid$volume, fit$n)
  d_recompute <- numeric(length(indices))

  for (ii in seq_along(indices)) {
    i <- indices[ii]
    xi <- x[i, , drop = FALSE]
    idx_i <- cell_index(xi, fit$grid)
    key_i <- key_from_index(idx_i)
    tmp_env <- new.env(parent = emptyenv(), hash = TRUE, size = length(fit$counts$keys) * 2L)
    for (r in seq_along(fit$counts$keys)) assign(fit$counts$keys[[r]], fit$counts$counts[[r]], envir = tmp_env)
    old <- lookup_count(tmp_env, key_i)
    assign(key_i, as.integer(old - 1L), envir = tmp_env)
    tmp_fit <- fit
    tmp_fit$count_env <- tmp_env
    tmp_fit$n <- fit$n - 1L
    f_minus <- glbfp_eval_naive(tmp_fit, xi, max_stencil = max_stencil)
    d_recompute[ii] <- 1 - f_minus / f_full[ii]
  }

  data.frame(index = indices, D_formula = d_formula, D_recompute = d_recompute,
             abs_error = abs(d_formula - d_recompute))
}
