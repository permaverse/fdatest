#' Two population functional Benjamini-Hochberg procedure
#'
#' The function implements the functional Benjamini Hochberg (fBH) procedure for
#' testing mean differences between two functional populations. Functional data
#' are tested locally and unadjusted and adjusted p-value functions are
#' provided. The unadjusted p-value function controls the point-wise error rate.
#' The adjusted p-value function controls the family-wise error rate
#' asymptotically.
#'
#' @inherit functional_two_sample_test params return seealso
#'
#' @references
#' Lundtorp Olsen, N., Pini, A., & Vantini, S. (2021). False discovery rate for
#' functional data \emph{TEST} 30, 784–809.
#'
#' @export
#' @examples
#' # Performing the fBH for two populations
#'
#' FDR_result <- FDR2(NASAtemp$paris, NASAtemp$milan)
#'
#' # Plotting the results of the fBH
#' plot(
#'   FDR_result,
#'   xrange = c(0, 12),
#'   title = 'FDR results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(FDR_result$adjusted_pvalues < 0.05)
FDR2 <- # nolint: object_name_linter.
  function(
    data1,
    data2,
    mu = 0,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    paired = FALSE,
    alternative = c("two.sided", "less", "greater")
  ) {
    fdr2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      dx = dx,
      n_perm = B,
      paired = paired,
      alternative = alternative
    )
  }

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname FDR2
#' @export
fdr2 <- function(
  data1,
  data2,
  mu = 0,
  dx = NULL,
  n_perm = 1000L,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater")
) {
  alternative <- rlang::arg_match(alternative)

  inputs <- twosamples2coeffs(data1, data2, mu, dx = dx)
  coeff1 <- inputs$coeff1
  coeff2 <- inputs$coeff2
  mu_eval <- inputs$mu

  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  n <- n1 + n2
  etichetta_ord <- c(rep(1, n1), rep(2, n2))
  coeff <- rbind(coeff1, coeff2)
  data_eval <- coeff

  perm_res <- twosample_alt_permtest(coeff, n1, n_perm, alternative, paired)
  pval <- perm_res$pval

  adjusted_pval <- stats::p.adjust(pval, method = "BH")

  out <- list(
    data = data_eval,
    group_labels = etichetta_ord,
    mu = mu_eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = adjusted_pval
  )
  class(out) <- "ftwosample"
  out
}
