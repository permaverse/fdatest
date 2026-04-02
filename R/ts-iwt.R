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
IWT2 <- # nolint: object_name_linter.
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
    aggregation_strategy = c("integral", "max"),
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
      standardize = standardize,
      verbose = verbose,
      aggregation_strategy = aggregation_strategy,
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
  standardize = FALSE,
  verbose = FALSE,
  aggregation_strategy = c("integral", "max"),
  recycle = TRUE
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
    cli::cli_h1("Interval-Wise Testing")
  }

  aggregation_strategy <- rlang::arg_match(aggregation_strategy)

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

  if (verbose) {
    cli::cli_h1("Interval-Wise Testing completed")
  }

  out <- list(
    data = data_eval,
    group_labels = group_labels,
    mu = mu_eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = corrected_pval,
    pvalue_matrix = matrice_pval_asymm
  )
  class(out) <- "fts"
  out
}
