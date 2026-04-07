p_adjust_twt <- function(pval, p, t0, t_coeff, aggregation_strategy) {
  thresholds <- c(0, sort(unique(pval)), 1)
  adjusted_pval <- pval
  pval_tmp <- rep(0, p)
  for (test in seq_along(thresholds)) {
    points_1 <- which(pval <= thresholds[test])
    t0_comb <- if (aggregation_strategy == "integral") {
      sum(t0[points_1], na.rm = TRUE)
    } else {
      max(t0[points_1], na.rm = TRUE)
    }
    t_comb <- if (aggregation_strategy == "integral") {
      rowSums(t_coeff[, points_1, drop = FALSE], na.rm = TRUE)
    } else {
      apply(t_coeff[, points_1, drop = FALSE], 1, max, na.rm = TRUE)
    }
    pval_tmp[points_1] <- mean(t_comb >= t0_comb)
    adjusted_pval <- apply(rbind(adjusted_pval, pval_tmp), 2, max)

    points_2 <- which(pval > thresholds[test])
    t0_comb <- if (aggregation_strategy == "integral") {
      sum(t0[points_2], na.rm = TRUE)
    } else {
      max(t0[points_2], na.rm = TRUE)
    }
    t_comb <- if (aggregation_strategy == "integral") {
      rowSums(t_coeff[, points_2, drop = FALSE], na.rm = TRUE)
    } else {
      apply(t_coeff[, points_2, drop = FALSE], 1, max, na.rm = TRUE)
    }
    pval_tmp[points_2] <- mean(t_comb >= t0_comb)
    adjusted_pval <- apply(rbind(adjusted_pval, pval_tmp), 2, max)
  }

  list(
    adjusted_pvalues = adjusted_pval
  )
}
