p_adjust_global <- function(
  aggregation_strategy,
  t0,
  t_coeff,
  p
) {
  adjusted_pval <- switch(
    aggregation_strategy,
    integral = {
      t0_comb <- sum(t0)
      t_comb <- rowSums(t_coeff)
      pval_temp <- mean(t_comb >= t0_comb)
      rep(pval_temp, p)
    },
    max = {
      t0_comb <- max(t0)
      t_comb <- apply(t_coeff, 1, max)
      pval_temp <- mean(t_comb >= t0_comb)
      rep(pval_temp, p)
    }
  )

  list(
    adjusted_pvalues = adjusted_pval
  )
}
