#' Two population Interval Wise Testing procedure
#'
#' The function implements the Interval Wise Testing procedure for testing mean
#' differences between two functional populations. Functional data are tested
#' locally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the interval-wise error rate.
#'
#' @inherit functional_two_sample_test params return
#'
#' @seealso [`global2()`], [`twt2()`], [`pct2()`], [`fdr2()`] for calling directly
#' one of the other tests, [`functional_two_sample_test()`] for calling the
#' interface test and [`plot.fts()`] for plotting the results.
#'
#' @references
#' - Pini, Alessia, and Simone Vantini. 2016. “The interval testing procedure: a
#' general framework for inference in functional data analysis.” \emph{Biometrics} 72 (3):
#' 835–845.
#' - Pini, Alessia, and Simone Vantini. 2017. “Interval-Wise Testing for Functional
#' Data.” \emph{Journal of Nonparametric Statistics} 29 (2): 407–24.
#' - Pini, Alessia, Simone Vantini, Bianca Maria Colosimo, and Marco Grasso. 2018.
#' “Domain-Selective Functional Analysis of Variance for Supervised Statistical
#' Profile Monitoring of Signal Data.” \emph{Journal of the Royal Statistical Society
#' Series C: Applied Statistics} 67 (1): 55–81.
#' - Abramowicz, Konrad, Charlotte K Häger, Alessia Pini, Lina Schelin, Sara
#' Sjöstedt de Luna, and Simone Vantini. 2018. “Nonparametric Inference for
#' Functional-on-Scalar Linear Models Applied to Knee Kinematic Hop Data After
#' Injury of the Anterior Cruciate Ligament.” \emph{Scandinavian Journal of Statistics} 45
#' (4): 1036–61.
#'
#' @export
#' @examples
#' # Performing the IWT for two populations
#' IWT_result <- iwt2(NASAtemp$paris, NASAtemp$milan, n_perm = 10L)
#'
#' # Plotting the results of the IWT
#' plot(
#'   IWT_result,
#'   xrange = c(0, 12),
#'   title = "IWT results for testing mean differences"
#' )
#'
#' # Plotting the p-value heatmap
#' IWTimage(IWT_result, abscissa_range = c(0, 12))
#'
#' # Selecting the significant components at 5% level
#' which(IWT_result$adjusted_pvalues < 0.05)
iwt2 <- function(
  data1,
  data2,
  mu = 0,
  dx = NULL,
  n_perm = 1000L,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  standardize = FALSE,
  verbose = FALSE,
  aggregation_strategy = c("integral", "max"),
  recycle = TRUE
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
    verbose = verbose,
    correction = "IWT",
    aggregation_strategy = aggregation_strategy,
    recycle = recycle
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
#' @rdname iwt2
#' @export
IWT2 <- # nolint: object_name_linter.
  function(
    data1,
    data2,
    mu = 0,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    paired = FALSE,
    alternative = c("two.sided", "less", "greater"),
    statistic = c("Integral", "Max", "Integral_std", "Max_std"),
    verbose = FALSE,
    recycle = TRUE
  ) {
    statistic <- rlang::arg_match(statistic)
    standardize <- statistic %in% c("Integral_std", "Max_std")
    aggregation_strategy <- switch(
      statistic,
      "Integral" = "integral",
      "Max" = "max",
      "Integral_std" = "integral",
      "Max_std" = "max"
    )
    lifecycle::deprecate_warn(
      when = "0.2.0",
      what = "IWT2()",
      details = "Use iwt2() instead. Be mindful that the argument `statistic` has been replaced by `aggregation_strategy` and `standardize`.",
      id = "fdatest-deprecated-iwt2"
    )
    iwt2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      dx = dx,
      n_perm = B,
      paired = paired,
      alternative = alternative,
      standardize = standardize,
      verbose = verbose,
      aggregation_strategy = aggregation_strategy,
      recycle = recycle
    )
  }
