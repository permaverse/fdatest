#' Two population Threshold Wise Testing procedure
#'
#' The function implements the Threshold Wise Testing procedure for testing mean
#' differences between two functional populations. Functional data are tested
#' locally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the family-wise error rate asymptotically.
#'
#' @inherit functional_two_sample_test params return
#'
#' @seealso [`global2()`], [`iwt2()`], [`pct2()`], [`fdr2()`] for calling directly
#' one of the other tests, [`functional_two_sample_test()`] for calling the
#' interface test and [`plot.fts()`] for plotting the results.
#'
#' @references
#' - Abramowicz, Konrad, Alessia Pini, Lina Schelin, Sara Sjöstedt de Luna,
#' Aymeric Stamm, and Simone Vantini. 2023. “Domain Selection and Familywise
#' Error Rate for Functional Data: A Unified Framework.” Biometrics 79 (2):
#' 1119–32.
#'
#' @export
#' @examples
#' # Performing the TWT for two populations
#' TWT_result <- twt2(NASAtemp$paris, NASAtemp$milan)
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
twt2 <- function(
  data1,
  data2,
  mu = 0,
  dx = NULL,
  n_perm = 1000L,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  standardize = FALSE,
  verbose = FALSE,
  aggregation_strategy = c("integral", "max")
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
    correction = "TWT",
    aggregation_strategy = aggregation_strategy
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
#' @rdname twt2
#' @export
TWT2 <- # nolint: object_name_linter.
  function(
    data1,
    data2,
    mu = 0,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    paired = FALSE,
    alternative = c("two.sided", "less", "greater"),
    statistic = c("Integral", "Max", "Integral_std", "Max_std"),
    verbose = FALSE
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
      what = "TWT2()",
      details = "Use twt2() instead. Be mindful that the argument `statistic` has been replaced by `aggregation_strategy` and `standardize`.",
      id = "fdatest-deprecated-twt2"
    )
    twt2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      dx = dx,
      n_perm = B,
      paired = paired,
      alternative = alternative,
      standardize = standardize,
      verbose = verbose,
      aggregation_strategy = aggregation_strategy
    )
  }

ts_p_adjust_twt <- function(pval, p, t0, t_coeff, aggregation_strategy) {
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
