#' Two population Threshold Wise Testing procedure
#'
#' The function implements the Threshold Wise Testing procedure for testing mean
#' differences between two functional populations. Functional data are tested
#' locally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the family-wise error rate asymptotically.
#'
#' @inherit functional_two_sample_test params return seealso
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. \emph{Biometrics} 79(2), 1119-1132.
#'
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data.
#' \emph{Journal of Nonparametric Statistics}, 29(2), 407-424.
#'
#' @export
#' @examples
#' # Performing the TWT for two populations
#' TWT_result <- TWT2(NASAtemp$paris, NASAtemp$milan)
#'
#' # Plotting the results of the TWT
#' plot(
#'   TWT_result,
#'   xrange = c(0, 12),
#'   title = 'TWT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(TWT_result$adjusted_pvalues < 0.05)
TWT2 <- function( # nolint: object_name_linter.
  data1,
  data2,
  mu = 0,
  dx = NULL,
  B = 1000L, # nolint: object_name_linter.
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  verbose = FALSE
) {
  twt2(
    data1 = data1,
    data2 = data2,
    mu = mu,
    dx = dx,
    n_perm = B,
    paired = paired,
    alternative = alternative,
    verbose = verbose
  )
}

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname TWT2
#' @export
twt2 <- function(
  data1,
  data2,
  mu = 0,
  dx = NULL,
  n_perm = 1000L,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  verbose = FALSE
) {
  alternative <- rlang::arg_match(alternative)

  inputs <- twosamples2coeffs(data1, data2, mu, dx = dx)
  coeff1 <- inputs$coeff1
  coeff2 <- inputs$coeff2
  mu_eval <- inputs$mu

  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  n <- n1 + n2
  data_eval <- rbind(coeff1, coeff2)
  p <- dim(data_eval)[2]
  coeff1 <- coeff1 - matrix(data = mu_eval, nrow = n1, ncol = p)
  coeff <- rbind(coeff1, coeff2)
  etichetta_ord <- c(rep(1, n1), rep(2, n2))

  perm_res <- twosample_alt_permtest(coeff, n1, n_perm, alternative, paired)
  t0 <- perm_res$t0
  t_coeff <- perm_res$t_coeff
  pval <- perm_res$pval

  if (verbose) {
    cli::cli_h1("Threshold-wise tests")
  }

  thresholds <- c(0, sort(unique(pval)), 1)
  adjusted_pval <- pval
  pval_tmp <- rep(0, p)
  for (test in seq_along(thresholds)) {
    points_1 <- which(pval <= thresholds[test])
    t0_comb <- sum(t0[points_1], na.rm = TRUE)
    t_comb <- rowSums(t_coeff[, points_1, drop = FALSE], na.rm = TRUE)
    pval_tmp[points_1] <- mean(t_comb >= t0_comb)
    adjusted_pval <- apply(rbind(adjusted_pval, pval_tmp), 2, max)

    points_2 <- which(pval > thresholds[test])
    t0_comb <- sum(t0[points_2])
    t_comb <- rowSums(t_coeff[, points_2, drop = FALSE], na.rm = TRUE)
    pval_tmp[points_2] <- mean(t_comb >= t0_comb)
    adjusted_pval <- apply(rbind(adjusted_pval, pval_tmp), 2, max)
  }

  out <- list(
    data = coeff,
    group_labels = etichetta_ord,
    mu = mu_eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = adjusted_pval
  )
  class(out) <- "ftwosample"
  out
}
