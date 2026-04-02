#' Two population Global Testing procedure
#'
#' The function implements the Global Testing procedure for testing mean
#' differences between two functional populations. Functional data are tested
#' locally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the interval-wise error rate.
#'
#' @inherit functional_two_sample_test params return seealso
#'
#' @references
#' A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference for
#' Functional Data Controlling the Family Wise Error Rate on Intervals.
#' Biometrics 73(3): 835–845.
#'
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data.
#' \emph{Journal of Nonparametric Statistics}, 29(2), 407-424
#'
#' @export
#' @examples
#' # Performing the Global for two populations
#' Global_result <- Global2(NASAtemp$paris, NASAtemp$milan)
#'
#' # Plotting the results of the Global
#' plot(
#'   Global_result,
#'   xrange = c(0, 12),
#'   title = 'Global results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(Global_result$adjusted_pvalues < 0.05)
Global2 <- # nolint: object_name_linter.
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
    global2(
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
#' @rdname Global2
#' @export
global2 <- function(
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

  aggregation_strategy <- rlang::arg_match(aggregation_strategy)

  adjusted_pval <- switch(
    aggregation_strategy,
    integral = {
      t0_comb <- sum(t0)
      t_comb <- rowSums(t_coeff)
      pval_temp <- mean(t_comb >= t0_comb)
      rep(pval_temp, p)
    },
    max = {
      t0_comb <- max(t0)
      t_comb <- apply(t_coeff, 1, max)
      pval_temp <- mean(t_comb >= t0_comb)
      rep(pval_temp, p)
    }
  )

  out <- list(
    data = data_eval,
    group_labels = group_labels,
    mu = mu_eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = adjusted_pval,
    global_pvalue = adjusted_pval[1]
  )
  class(out) <- "fts"
  out
}
