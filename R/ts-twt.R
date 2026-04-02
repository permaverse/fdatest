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
TWT2 <- # nolint: object_name_linter.
  function(
    data1,
    data2,
    mu = 0,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    paired = FALSE,
    alternative = c("two.sided", "less", "greater"),
    standardize = FALSE,
    verbose = FALSE,
    aggregation_strategy = c("integral", "max")
  ) {
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
  standardize = FALSE,
  verbose = FALSE,
  aggregation_strategy = c("integral", "max")
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
    cli::cli_h1("Threshold-Wise Testing")
  }

  aggregation_strategy <- rlang::arg_match(aggregation_strategy)

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
