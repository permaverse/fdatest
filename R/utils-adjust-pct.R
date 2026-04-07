p_adjust_pct <- function(
  partition,
  p,
  t0,
  t_coeff,
  aggregation_strategy
) {
  partition <- factor(partition)
  nintervals <- length(levels(partition))
  ntests <- 2^nintervals - 1L
  all_combs <- matrix(nrow = ntests, ncol = p)
  labels <- levels(partition)
  tt <- 1L
  for (nint in seq_len(nintervals)) {
    combinations <- utils::combn(labels, nint)
    n_comb <- dim(combinations)[2]
    for (comb in seq_len(n_comb)) {
      index <- rep(0, p)
      for (ii in seq_len(dim(combinations)[1])) {
        index <- index + as.numeric(partition == combinations[ii, comb])
      }
      all_combs[tt, ] <- index
      tt <- tt + 1L
    }
  }

  adjusted_pval <- numeric(p)
  for (test in seq_len(ntests)) {
    active <- which(all_combs[test, ] == 1)
    t0_comb <- if (aggregation_strategy == "integral") {
      sum(t0[active], na.rm = TRUE)
    } else {
      max(t0[active], na.rm = TRUE)
    }
    t_comb <- if (aggregation_strategy == "integral") {
      rowSums(t_coeff[, active, drop = FALSE], na.rm = TRUE)
    } else {
      apply(t_coeff[, active, drop = FALSE], 1, max, na.rm = TRUE)
    }
    pval_temp <- mean(t_comb >= t0_comb)
    adjusted_pval[active] <- apply(
      rbind(adjusted_pval[active], pval_temp),
      2,
      max
    )
  }

  list(
    adjusted_pvalues = adjusted_pval
  )
}
