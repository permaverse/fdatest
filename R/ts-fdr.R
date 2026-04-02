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
    alternative = c("two.sided", "less", "greater"),
    standardize = FALSE,
    verbose = FALSE
  ) {
    fdr2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      dx = dx,
      n_perm = B,
      paired = paired,
      alternative = alternative,
      standardize = standardize,
      verbose = verbose
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
  alternative = c("two.sided", "less", "greater"),
  standardize = FALSE,
  verbose = FALSE
) {
  if (verbose) {
    cli::cli_h1("Data preparation and point-wise testing")
  }

  prepped_data <- ts_prepare_data(
    data1 = data1,
    data2 = data2,
    mu = mu,
    dx = dx,
    n_perm = n_perm,
    paired = paired,
    alternative = alternative,
    standardize = standardize
  )

  data_eval <- prepped_data$data
  mu_eval <- prepped_data$mu
  group_labels <- prepped_data$group_labels
  p <- prepped_data$p

  t0 <- prepped_data$t0
  t_coeff <- prepped_data$t_coeff
  pval <- prepped_data$pval

  if (verbose) {
    cli::cli_h1("Benjamini-Hochberg FDR Testing")
  }

  adjusted_pval <- stats::p.adjust(pval, method = "BH")

  out <- list(
    data = data_eval,
    group_labels = group_labels,
    mu = mu_eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = adjusted_pval
  )
  class(out) <- "fts"
  out
}
