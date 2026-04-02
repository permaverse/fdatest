#' Two population Partition Closed Testing procedure
#'
#' The function implements the Partition Closed Testing procedure for testing
#' mean differences between two functional populations. Functional data are
#' tested locally and unadjusted and adjusted p-value functions are provided.
#' The unadjusted p-value function controls the point-wise error rate. The
#' adjusted p-value function controls the family-wise error rate asymptotically.
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
#' # Performing the PCT for two populations
#' # Choosing as partition the 4 seasons of the year
#' partition <- c(
#'   rep(1, 31 + 28 + 21),
#'   rep(2, 10 + 30 + 31 + 21),
#'   rep(3, 9 + 31 + 31 + 23),
#'   rep(4, 7 + 31 + 30 + 21),
#'   rep(1, 10)
#' )
#' partition <- factor(partition)
#'
#' PCT_result <- PCT2(NASAtemp$paris, NASAtemp$milan, partition = partition)
#'
#' # Plotting the results of the PCT
#' plot(
#'   PCT_result,
#'   xrange = c(0, 12),
#'   title = 'PCT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(PCT_result$adjusted_pvalues < 0.05)
PCT2 <- # nolint: object_name_linter.
  function(
    data1,
    data2,
    partition,
    mu = 0,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    paired = FALSE,
    alternative = c("two.sided", "less", "greater"),
    standardize = FALSE,
    verbose = FALSE,
    aggregation_strategy = c("integral", "max")
  ) {
    pct2(
      data1 = data1,
      data2 = data2,
      partition = partition,
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
#' @rdname PCT2
#' @export
pct2 <- function(
  data1,
  data2,
  partition,
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
    cli::cli_h1("Partition Closed Testing")
  }

  aggregation_strategy <- rlang::arg_match(aggregation_strategy)

  partition <- factor(partition)
  nintervals <- length(levels(partition))
  ntests <- 2^nintervals - 1L
  all_combs <- matrix(nrow = ntests, ncol = p)
  labels <- levels(partition)
  tt <- 1L
  for (nint in seq_len(nintervals)) {
    combinations <- utils::combn(labels, nint)
    n_comb <- dim(combinations)[2]
    for (comb in seq_len(n_comb)) {
      index <- rep(0, p)
      for (ii in seq_len(dim(combinations)[1])) {
        index <- index + as.numeric(partition == combinations[ii, comb])
      }
      all_combs[tt, ] <- index
      tt <- tt + 1L
    }
  }

  adjusted_pval <- numeric(p)
  for (test in seq_len(ntests)) {
    active <- which(all_combs[test, ] == 1)
    t0_comb <- if (aggregation_strategy == "integral") {
      sum(t0[active], na.rm = TRUE)
    } else {
      max(t0[active], na.rm = TRUE)
    }
    t_comb <- if (aggregation_strategy == "integral") {
      rowSums(t_coeff[, active, drop = FALSE], na.rm = TRUE)
    } else {
      apply(t_coeff[, active, drop = FALSE], 1, max, na.rm = TRUE)
    }
    pval_temp <- mean(t_comb >= t0_comb)
    adjusted_pval[active] <- apply(
      rbind(adjusted_pval[active], pval_temp),
      2,
      max
    )
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
