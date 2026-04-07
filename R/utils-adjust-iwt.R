iwt_compute_row <- function(
  i,
  t0_2x,
  t_coeff_2x,
  n_perm,
  p,
  recycle,
  aggregation_strategy
) {
  js <- if (recycle) seq_len(p) else seq_len(i)
  row_vals <- numeric(if (recycle) p else i)
  for (k in seq_along(js)) {
    j <- js[k]
    inf <- j
    sup <- (p - i) + j
    t0_temp <- if (aggregation_strategy == "integral") {
      sum(t0_2x[inf:sup], na.rm = TRUE)
    } else {
      max(t0_2x[inf:sup], na.rm = TRUE)
    }
    t_temp <- if (aggregation_strategy == "integral") {
      rowSums(t_coeff_2x[, inf:sup, drop = FALSE], na.rm = TRUE)
    } else {
      apply(t_coeff_2x[, inf:sup, drop = FALSE], 1, max, na.rm = TRUE)
    }
    row_vals[k] <- sum(t_temp >= t0_temp) / n_perm
  }
  row_vals
}

iwt_compute_row_pair <- function(
  i,
  t0_2x,
  t_coeff_2x,
  n_perm,
  p,
  recycle,
  aggregation_strategy
) {
  if (i == p - i) {
    return(list(iwt_compute_row(
      i = i,
      t0_2x = t0_2x,
      t_coeff_2x = t_coeff_2x,
      n_perm = n_perm,
      p = p,
      recycle = recycle,
      aggregation_strategy = aggregation_strategy
    )))
  }
  list(
    iwt_compute_row(
      i = i,
      t0_2x = t0_2x,
      t_coeff_2x = t_coeff_2x,
      n_perm = n_perm,
      p = p,
      recycle = recycle,
      aggregation_strategy = aggregation_strategy
    ),
    iwt_compute_row(
      i = p - i,
      t0_2x = t0_2x,
      t_coeff_2x = t_coeff_2x,
      n_perm = n_perm,
      p = p,
      recycle = recycle,
      aggregation_strategy = aggregation_strategy
    )
  )
}

p_adjust_iwt <- function(
  p,
  pval,
  t0,
  t_coeff,
  n_perm,
  recycle,
  aggregation_strategy
) {
  matrice_pval_asymm <- matrix(nrow = p, ncol = p)
  matrice_pval_asymm[p, ] <- pval[seq_len(p)]
  t0_2x <- c(t0, t0)
  t_coeff_2x <- cbind(t_coeff, t_coeff)

  row_indices <- (p - 1L):1L

  perm_args <- list(
    t0_2x = t0_2x,
    t_coeff_2x = t_coeff_2x,
    n_perm = n_perm,
    p = p,
    recycle = recycle,
    aggregation_strategy = aggregation_strategy
  )

  if (mirai::daemons_set()) {
    optimized_order <- optimize_order(row_indices)
    row_tasks <- mirai::mirai_map((p - 1L):floor(p / 2), function(.i) {
      rlang::inject(iwt_compute_row_pair(.i, !!!perm_args))
    })
    row_results <- row_tasks[.progress]
    row_results <- unlist(row_results, recursive = FALSE)
    row_results <- row_results[order(optimized_order, decreasing = TRUE)]
  } else {
    row_results <- lapply(row_indices, function(.i) {
      rlang::inject(iwt_compute_row(.i, !!!perm_args))
    })
  }

  for (k in seq_along(row_indices)) {
    i <- row_indices[k]
    js <- if (recycle) seq_len(p) else seq_len(i)
    matrice_pval_asymm[i, js] <- row_results[[k]]
  }

  corrected_pval_matrix <- pval_correct_cpp(matrice_pval_asymm)
  corrected_pval <- corrected_pval_matrix[1, ]

  list(
    adjusted_pvalues = corrected_pval,
    pvalue_matrix = matrice_pval_asymm
  )
}
