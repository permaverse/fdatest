#' Local testing procedures for the functional two-sample test
#'
#' @description The function implements local testing procedures for testing
#'   mean differences between two functional populations. Functional data are
#'   tested locally and unadjusted and adjusted p-value functions are provided.
#'   The unadjusted p-value function controls the point-wise error rate. The
#'   adjusted p-value function can be computed according to the following
#'   methods:
#'
#'   - global testing (controlling the FWER weakly)
#'   - interval-wise testing (controlling the interval-wise error rate)
#'   - threshold-wise testing (controlling the FWER asymptotically)
#'   - partition closed testing (controlling the FWER on a partition)
#'   - functional Benjamini Hochberg (controlling the FDR)
#'
#' @param data1 Either a numeric matrix or an object of class [`fda::fd`]
#'   specifying the data in the first sample. If the data is provided within a
#'   matrix, it should be of shape \eqn{n_1 \times J} and it should contain in
#'   each row one of the \eqn{n_1} functions in the sample and in columns the
#'   evaluation of each function on a **same** uniform grid of size \eqn{J}.
#' @param data2 Either a numeric matrix or an object of class [`fda::fd`]
#'   specifying the data in the second sample. If the data is provided within a
#'   matrix, it should be of shape \eqn{n_2 \times J} and it should contain in
#'   each row one of the \eqn{n_2} functions in the sample and in columns the
#'   evaluation of each function on a **same** uniform grid of size \eqn{J}.
#' @param correction A string specifying the correction method to perform the
#'   local functional testing procedure and adjust the p-value function. Choices
#'   are `"Global`, `"IWT"`, `"TWT"`, `"PCT"` or `"FDR"`.
#' @param mu Either a numeric value or a numeric vector or an object of class
#'   [`fda::fd`] specifying the functional mean difference under the null
#'   hypothesis. If `mu` is a constant, then a constant function is used. If
#'   `mu` is a numeric vector, it must correspond to evaluation of the mean
#'   difference function on the **same** grid that has been used to evaluate the
#'   data samples. Defaults to `0`.
#' @param dx A numeric value specifying the step of the uniform grid on which
#'   the data are evaluated. If `NULL`, the step is automatically inferred from
#'   the data. Defaults to `NULL`.
#' @param n_perm An integer value specifying the number of permutations to use
#'   for the local testing procedure. Defaults to `1000L`.
#' @param paired A boolean value specifying whether a paired test should be
#'   performed. Defaults to `FALSE`.
#' @param alternative A string specifying the type of alternative hypothesis.
#'   Choices are `"two.sided"`, `"less"` or `"greater"`. Defaults to
#'   `"two.sided"`.
#' @param standardize A boolean value specifying whether to standardize the test
#'   statistic. Defaults to `FALSE`.
#' @param verbose A boolean value specifying whether to print the progress of
#'  the computation. Defaults to `FALSE`.
#' @param aggregation_strategy A string specifying the strategy to aggregate the
#'   point-wise test statistics for the correction procedure. Possible values
#'   are `"integral"` and `"max"`. Defaults to `"integral"`.
#' @param recycle A boolean value specifying whether to recycle the test statistic
#'   values across permutations for the IWT procedure. Defaults to `TRUE`.
#' @param partition An integer vector of length \eqn{J} specifying the
#'   membership of each point of the domain to an element of the partition.
#'   Only used and **must** be set if the `correction` argument is set to
#'   `"PCT"`.
#'
#' @returns An object of class `fts` containing the following components:
#'
#'   - `data`: A numeric matrix of shape \eqn{n \times J} containing the
#'   evaluation of the \eqn{n = n_1 + n_2} functions on a **common** uniform
#'   grid of size \eqn{p}.
#'   - `group_labels`: An integer vector of size \eqn{n = n_1 + n_2} containing
#'   the group membership of each function.
#'   - `mu`: A numeric vector of shape \eqn{J} containing the evaluation of the
#'   functional mean difference under the null hypothesis on the same uniform
#'   grid used to evaluate the functional samples.
#'   - `unadjusted_pvalues`: A numeric vector of size \eqn{J} containing the
#'   evaluation of the unadjusted p-value function on the **same** uniform grid
#'   used to evaluate the functional samples.
#'   - `adjusted_pvalues`: A numeric vector of size \eqn{J} containing the
#'   evaluation of the adjusted p-value functione on the **same** uniform grid
#'   used to evaluate the functional samples.
#'
#'   Optionally, the list may contain the following components:
#'
#'   - `global_pvalue`: A numeric value containing the global p-value. Only
#'   present if the `correction` argument is set to `"Global"`.
#'   - `pvalue_matrix`: A numeric matrix of shape \eqn{p \times p} containing
#'   the p-values of the interval-wise tests. Element \eqn{i, j} contains
#'   the p-value of the test performed on the interval indexed by
#'   \eqn{j, j+1 , \dots, j+(p-i)}. Only present if the `correction` argument is
#'   set to `"IWT"`.
#'
#' @seealso See also [`plot.fts()`] for plotting the results.
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. \emph{Biometrics} 79(2), 1119-1132.
#'
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data.
#' \emph{Journal of Nonparametric Statistics}, 29(2), 407-424
#'
#' A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference for
#' Functional Data Controlling the Family Wise Error Rate on Intervals.
#' Biometrics 73(3): 835–845.
#'
#' Lundtorp Olsen, N., Pini, A., & Vantini, S. (2021). False discovery rate for
#' functional data \emph{TEST} 30, 784–809.
#'
#' @export
#' @examples
#' # Performing the TWT for two populations
#' TWT_result <- functional_two_sample_test(
#'   NASAtemp$paris, NASAtemp$milan,
#'   correction = "TWT", n_perm = 10L
#' )
#'
#' # Plotting the results of the TWT
#' plot(
#'   TWT_result,
#'   xrange = c(0, 12),
#'   title = 'TWT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(TWT_result$adjusted_pval < 0.05)
#'
#' # Performing the IWT for two populations
#' IWT_result <- functional_two_sample_test(
#'   NASAtemp$paris, NASAtemp$milan,
#'   correction = "IWT", n_perm = 10L
#' )
#'
#' # Plotting the results of the IWT
#' plot(
#'   IWT_result,
#'   xrange = c(0, 12),
#'   title = 'IWT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(IWT_result$adjusted_pval < 0.05)
functional_two_sample_test <- function(
  data1,
  data2,
  correction = c("Global", "IWT", "TWT", "PCT", "FDR"),
  mu = 0,
  dx = NULL,
  n_perm = 1000L,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  standardize = FALSE,
  verbose = FALSE,
  aggregation_strategy = c("integral", "max"),
  recycle = TRUE,
  partition = NULL
) {
  correction <- rlang::arg_match(correction)
  if (correction == "PCT" && is.null(partition)) {
    cli::cli_abort(
      'When the {.arg correction} argument is set to {.code "PCT"}, the {.arg
      partition} argument should be explicitly provided.'
    )
  }

  if (verbose) {
    cli::cli_h1("Data preparation and point-wise testing")
  }

  alternative <- rlang::arg_match(alternative)
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
    switch(
      correction,
      IWT = cli::cli_h1("P-Value Adjustment via Interval-Wise Testing"),
      TWT = cli::cli_h1("P-Value Adjustment via Threshold-Wise Testing"),
      PCT = cli::cli_h1("P-Value Adjustment via Partition Closed Testing"),
      Global = cli::cli_h1("P-Value Adjustment via Global Testing"),
      FDR = cli::cli_h1("P-Value Adjustment via Benjamini-Hochberg FDR Testing")
    )
  }

  aggregation_strategy <- rlang::arg_match(aggregation_strategy)
  adjustment_results <- switch(
    correction,
    IWT = ts_p_adjust_iwt(
      p = p,
      pval = pval,
      t0 = t0,
      t_coeff = t_coeff,
      n_perm = n_perm,
      recycle = recycle,
      verbose = verbose,
      aggregation_strategy = aggregation_strategy
    ),
    TWT = ts_p_adjust_twt(
      pval = pval,
      p = p,
      t0 = t0,
      t_coeff = t_coeff,
      aggregation_strategy = aggregation_strategy
    ),
    PCT = ts_p_adjust_pct(
      partition = partition,
      p = p,
      t0 = t0,
      t_coeff = t_coeff,
      aggregation_strategy = aggregation_strategy
    ),
    Global = ts_p_adjust_global(
      aggregation_strategy = aggregation_strategy,
      t0 = t0,
      t_coeff = t_coeff,
      p = p
    ),
    FDR = ts_p_adjust_fdr(pval = pval)
  )

  out <- list(
    data = data_eval,
    group_labels = group_labels,
    mu = mu_eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = adjustment_results$adjusted_pvalues,
    correction_method = correction
  )

  if (correction == "Global") {
    out$global_pvalue <- adjustment_results$adjusted_pvalues[1]
  }

  if (correction == "IWT") {
    out$pvalue_matrix <- adjustment_results$pvalue_matrix
  }

  class(out) <- "fts"
  out
}
