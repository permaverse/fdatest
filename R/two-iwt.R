#' Two population Interval Wise Testing procedure
#'
#' The function implements the Interval Wise Testing procedure for testing mean
#' differences between two functional populations. Functional data are tested
#' locally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the interval-wise error rate.
#'
#' @inherit functional_two_sample_test params return seealso
#'
#' @references
#' A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference
#' for Functional Data Controlling the Family Wise Error Rate on Intervals.
#' *Biometrics*, 73(3): 835–845.
#'
#' A. Pini and S. Vantini (2017). Interval-wise testing for functional data.
#' *Journal of Nonparametric Statistics*, 29(2), 407-424.
#'
#' @export
#' @examples
#' # Performing the IWT for two populations
#' IWT_result <- IWT2(NASAtemp$paris, NASAtemp$milan, B = 10L)
#'
#' # Plotting the results of the IWT
#' plot(
#'   IWT_result,
#'   xrange = c(0, 12),
#'   title = 'IWT results for testing mean differences'
#' )
#'
#' # Plotting the p-value heatmap
#' IWTimage(IWT_result, abscissa_range = c(0, 12))
#'
#' # Selecting the significant components at 5% level
#' which(IWT_result$adjusted_pvalues < 0.05)
IWT2 <- function( # nolint: object_name_linter.
  data1,
  data2,
  mu = 0,
  dx = NULL,
  B = 1000L, # nolint: object_name_linter.
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  verbose = FALSE,
  recycle = TRUE
) {
  iwt2(
    data1 = data1,
    data2 = data2,
    mu = mu,
    dx = dx,
    n_perm = B,
    paired = paired,
    alternative = alternative,
    verbose = verbose,
    recycle = recycle
  )
}

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname IWT2
#' @export
iwt2 <- function(
  data1,
  data2,
  mu = 0,
  dx = NULL,
  n_perm = 1000L,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  verbose = FALSE,
  recycle = TRUE
) {
  alternative <- rlang::arg_match(alternative)

  inputs <- twosamples2coeffs(data1, data2, mu, dx = dx)
  coeff1 <- inputs$coeff1
  coeff2 <- inputs$coeff2
  mu_eval <- inputs$mu

  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  p <- dim(coeff1)[2]
  n <- n1 + n2
  etichetta_ord <- c(rep(1, n1), rep(2, n2))
  coeff1 <- coeff1 - matrix(data = mu_eval, nrow = n1, ncol = p)

  coeff <- rbind(coeff1, coeff2)
  data_eval <- rbind(
    coeff1 + matrix(mu_eval, nrow = n1, ncol = p),
    coeff2
  )

  if (verbose) {
    cli::cli_h1("Point-wise tests")
  }

  perm_res <- twosample_alt_permtest(coeff, n1, n_perm, alternative, paired)
  t0 <- perm_res$t0
  t_coeff <- perm_res$t_coeff
  pval <- perm_res$pval

  if (verbose) {
    cli::cli_h1("Interval-wise tests")
  }

  matrice_pval_asymm <- matrix(nrow = p, ncol = p)
  matrice_pval_asymm[p, ] <- pval[seq_len(p)]
  t0_2x <- c(t0, t0)
  t_coeff_2x <- cbind(t_coeff, t_coeff)

  maxrow <- 1L
  if (recycle) {
    for (i in (p - 1L):maxrow) {
      for (j in seq_len(p)) {
        inf <- j
        sup <- (p - i) + j
        t0_temp <- sum(t0_2x[inf:sup])
        t_temp <- rowSums(t_coeff_2x[, inf:sup, drop = FALSE])
        matrice_pval_asymm[i, j] <- sum(t_temp >= t0_temp) / n_perm
      }
      if (verbose) {
        cli::cli_h1(
          "Creating the p-value matrix: end of row {p - i + 1} out of {p}"
        )
      }
    }
  } else {
    for (i in (p - 1L):maxrow) {
      for (j in seq_len(i)) {
        inf <- j
        sup <- (p - i) + j
        t0_temp <- sum(t0_2x[inf:sup])
        t_temp <- rowSums(t_coeff_2x[, inf:sup, drop = FALSE])
        matrice_pval_asymm[i, j] <- sum(t_temp >= t0_temp) / n_perm
      }
      if (verbose) {
        cli::cli_h1(
          "Creating the p-value matrix: end of row {p - i + 1} out of {p}"
        )
      }
    }
  }

  corrected_pval_matrix <- pval_correct(matrice_pval_asymm)
  corrected_pval <- corrected_pval_matrix[1, ]

  if (verbose) {
    cli::cli_h1("Interval-Wise Testing completed")
  }

  out <- list(
    data = data_eval,
    group_labels = etichetta_ord,
    mu = mu_eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = corrected_pval,
    pvalue_matrix = matrice_pval_asymm
  )
  class(out) <- "ftwosample"
  out
}
