#' Two population functional Benjamini-Hochberg procedure
#'
#' The function implements the functional Benjamini Hochberg (fBH) procedure for
#' testing mean differences between two functional populations. Functional data
#' are tested locally and unadjusted and adjusted p-value functions are
#' provided. The unadjusted p-value function controls the point-wise error rate.
#' The adjusted p-value function controls the family-wise error rate
#' asymptotically.
#'
#' @inherit functional_two_sample_test params return
#'
#' @seealso [`global2()`], [`twt2()`], [`pct2()`], [`iwt2()`] for calling directly
#' one of the other tests, [`functional_two_sample_test()`] for calling the
#' interface test and [`plot.fts()`] for plotting the results.
#'
#' @references
#' - Lundtorp Olsen, Niels, Alessia Pini, and Simone Vantini. 2021. "False discovery
#' rate for functional data." TEST 30, 784–809.
#'
#' @export
#' @examples
#' # Performing the fBH for two populations
#'
#' FDR_result <- fdr2(NASAtemp$paris, NASAtemp$milan)
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
  functional_two_sample_test(
    data1 = data1,
    data2 = data2,
    mu = mu,
    dx = dx,
    n_perm = n_perm,
    paired = paired,
    alternative = alternative,
    standardize = standardize,
    correction = "FDR",
    verbose = verbose
  )
}

#' @param B An integer value specifying the number of permutations to use
#'   for the local testing procedure. Defaults to `1000L`.
#' @param statistic A string specifying the test statistic to use. Possible
#'   values are:
#'
#'   - `"Integral"`: Integral of the squared sample mean difference.
#'   - `"Max"`: Maximum of the squared sample mean difference.
#'   - `"Integral_std"`: Integral of the squared t-test statistic.
#'   - `"Max_std"`: Maximum of the squared t-test statistic.
#'
#'   Defaults to `"Integral"`.
#' @rdname fdr2
#' @export
FDR2 <- # nolint: object_name_linter.
  function(
    data1,
    data2,
    mu = 0,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    paired = FALSE,
    alternative = c("two.sided", "less", "greater"),
    statistic = c("Integral", "Integral_std"),
    verbose = FALSE
  ) {
    statistic <- rlang::arg_match(statistic)
    standardize <- statistic == "Integral_std"
    lifecycle::deprecate_warn(
      when = "0.2.0",
      what = "FDR2()",
      details = "Use fdr2() instead. Be mindful that the argument `statistic` has been replaced by `standardize`.",
      id = "fdatest-deprecated-fdr2"
    )
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

ts_p_adjust_fdr <- function(pval) {
  list(
    adjusted_pvalues = stats::p.adjust(pval, method = "BH")
  )
}
