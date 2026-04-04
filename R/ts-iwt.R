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
#' general framework for inference in functional data analysis.” Biometrics 72 (3):
#' 835–845.
#' - Pini, Alessia, and Simone Vantini. 2017. “Interval-Wise Testing for Functional
#' Data.” Journal of Nonparametric Statistics 29 (2): 407–24.
#' - Pini, Alessia, Simone Vantini, Bianca Maria Colosimo, and Marco Grasso. 2018.
#' “Domain-Selective Functional Analysis of Variance for Supervised Statistical
#' Profile Monitoring of Signal Data.” Journal of the Royal Statistical Society
#' Series C: Applied Statistics 67 (1): 55–81.
#' - Abramowicz, Konrad, Charlotte K Häger, Alessia Pini, Lina Schelin, Sara
#' Sjöstedt de Luna, and Simone Vantini. 2018. “Nonparametric Inference for
#' Functional-on-Scalar Linear Models Applied to Knee Kinematic Hop Data After
#' Injury of the Anterior Cruciate Ligament.” Scandinavian Journal of Statistics 45
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
#'   title = 'IWT results for testing mean differences'
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

ts_p_adjust_iwt <- function(
  p,
  pval,
  t0,
  t_coeff,
  n_perm,
  recycle,
  verbose,
  aggregation_strategy
) {
  matrice_pval_asymm <- matrix(nrow = p, ncol = p)
  matrice_pval_asymm[p, ] <- pval[seq_len(p)]
  t0_2x <- c(t0, t0)
  t_coeff_2x <- cbind(t_coeff, t_coeff)

  row_indices <- (p - 1L):1L

  perm_args <- list(
    t0_2x = t0_2x,
    t_coeff_2x = t_coeff_2x,
    n_perm = n_perm,
    p = p,
    recycle = recycle,
    aggregation_strategy = aggregation_strategy
  )

  if (mirai::daemons_set()) {
    optimized_order <- optimize_order(row_indices)
    row_tasks <- mirai::mirai_map((p - 1L):floor(p / 2), function(.i) {
      rlang::inject(compute_row_pair_ts(.i, !!!perm_args))
    })
    row_results <- row_tasks[.progress]
    row_results <- unlist(row_results, recursive = FALSE)
    row_results <- row_results[order(optimized_order, decreasing = TRUE)]
  } else {
    row_results <- lapply(row_indices, function(.i) {
      rlang::inject(compute_row_ts(.i, !!!perm_args))
    })
  }

  for (k in seq_along(row_indices)) {
    i <- row_indices[k]
    js <- if (recycle) seq_len(p) else seq_len(i)
    matrice_pval_asymm[i, js] <- row_results[[k]]
    if (verbose) {
      cli::cli_h1(
        "Creating the p-value matrix: end of row {p - i + 1} out of {p}"
      )
    }
  }

  corrected_pval_matrix <- pval_correct_cpp(matrice_pval_asymm)
  corrected_pval <- corrected_pval_matrix[1, ]

  list(
    adjusted_pvalues = corrected_pval,
    pvalue_matrix = matrice_pval_asymm
  )
}
