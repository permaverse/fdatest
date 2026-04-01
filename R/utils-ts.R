# Internal helper: one permutation iteration for twosample_alt_permtest.
# Returns a numeric vector of length p (the permuted test statistic row).
.twosample_one_perm <- function(coeff, n1, alternative, paired) {
  n <- nrow(coeff)
  if (paired) {
    if_perm <- stats::rbinom(n1, 1, 0.5)
    coeff_perm <- coeff
    for (couple in seq_len(n1)) {
      if (if_perm[couple] == 1L) {
        coeff_perm[c(couple, n1 + couple), ] <- coeff[c(n1 + couple, couple), ]
      }
    }
  } else {
    coeff_perm <- coeff[sample(n), ]
  }
  meandiff_perm <- colMeans(coeff_perm[seq_len(n1), , drop = FALSE]) -
    colMeans(coeff_perm[seq(n1 + 1L, n), , drop = FALSE])
  sign_diff_perm <- sign(meandiff_perm)
  sign_diff_perm[sign_diff_perm == -1L] <- 0L
  switch(
    alternative,
    two.sided = meandiff_perm^2,
    greater = (meandiff_perm * sign_diff_perm)^2,
    less = (meandiff_perm * (sign_diff_perm - 1L))^2
  )
}

# Internal helper: shared pointwise permutation test for two-sample functions
# (iwt2, twt2, fdr2, pct2). Returns list(t0, t_coeff, pval).
twosample_alt_permtest <- function(coeff, n1, n_perm, alternative, paired) {
  n <- nrow(coeff)
  p <- ncol(coeff)

  meandiff <- colMeans(coeff[seq_len(n1), , drop = FALSE], na.rm = TRUE) -
    colMeans(coeff[seq(n1 + 1L, n), , drop = FALSE], na.rm = TRUE)
  sign_diff <- sign(meandiff)
  sign_diff[sign_diff == -1L] <- 0L
  t0 <- switch(
    alternative,
    two.sided = meandiff^2,
    greater = (meandiff * sign_diff)^2,
    less = (meandiff * (sign_diff - 1L))^2
  )

  # Run permutations in parallel via mirai_map().
  # Each task returns one row of t_coeff (length p).
  perm_args <- list(
    coeff = coeff,
    n1 = n1,
    alternative = alternative,
    paired = paired
  )

  if (mirai::daemons_set()) {
    perm_tasks <- mirai::mirai_map(seq_len(n_perm), function(.x) {
      rlang::inject(.twosample_one_perm(!!!perm_args))
    })
    perm_results <- perm_tasks[.progress]
  } else {
    perm_results <- lapply(seq_len(n_perm), function(.x) {
      rlang::inject(.twosample_one_perm(!!!perm_args))
    })
  }
  t_coeff <- do.call(rbind, perm_results)

  pval <- colSums(
    t_coeff >= matrix(t0, nrow = n_perm, ncol = p, byrow = TRUE)
  ) /
    n_perm

  list(t0 = t0, t_coeff = t_coeff, pval = pval)
}
