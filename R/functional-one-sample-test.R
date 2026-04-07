#' Local testing procedures for the functional one-sample test
#'
#' @description The function implements local testing procedures for testing
#'   the center of symmetry of a functional population. Functional data are
#'   tested locally and unadjusted and adjusted p-value functions are provided.
#'   The unadjusted p-value function controls the point-wise error rate. The
#'   adjusted p-value function can be computed according to the following
#'   methods:
#'
#'   - interval-wise testing (controlling the interval-wise error rate)
#'
#' @param data Either a numeric matrix or an object of class [`fda::fd`]
#'   specifying the sample data. If the data is provided within a matrix, it
#'   should be of shape \eqn{n \times J} and it should contain in each row one
#'   of the \eqn{n} functions in the sample and in columns the evaluation of
#'   each function on a **same** uniform grid of size \eqn{J}.
#' @param correction A string specifying the correction method to perform the
#'   local functional testing procedure and adjust the p-value function.
#'   Currently only `"IWT"` is available.
#' @param mu Either a numeric value or a numeric vector or an object of class
#'   [`fda::fd`] specifying the functional center of symmetry under the null
#'   hypothesis. If `mu` is a constant, then a constant function is used. If
#'   `mu` is a numeric vector, it must correspond to evaluations of the center
#'   function on the **same** grid that has been used to evaluate the data
#'   sample. Defaults to `0`.
#' @param dx A numeric value specifying the step of the uniform grid on which
#'   the data are evaluated. If `NULL`, the step is automatically inferred from
#'   the data. Defaults to `NULL`.
#' @param n_perm An integer value specifying the number of permutations to use
#'   for the local testing procedure. Defaults to `1000L`.
#' @param verbose A boolean value specifying whether to print the progress of
#'   the computation. Defaults to `FALSE`.
#' @param aggregation_strategy A string specifying the strategy to aggregate the
#'   point-wise test statistics for the correction procedure. Possible values
#'   are `"integral"` and `"max"`. Defaults to `"integral"`.
#' @param recycle A boolean value specifying whether to recycle the test
#'   statistic values across permutations for the IWT procedure. Defaults to
#'   `TRUE`.
#'
#' @returns An object of class `fos` containing the following components:
#'
#'   - `data`: A numeric matrix of shape \eqn{n \times J} containing the
#'   evaluation of the \eqn{n} functions on a uniform grid of size \eqn{J}.
#'   - `mu`: A numeric vector of shape \eqn{J} containing the evaluation of the
#'   functional center of symmetry under the null hypothesis on the same uniform
#'   grid used to evaluate the functional sample.
#'   - `unadjusted_pvalues`: A numeric vector of size \eqn{J} containing the
#'   evaluation of the unadjusted p-value function on the same uniform grid
#'   used to evaluate the functional sample.
#'   - `adjusted_pvalues`: A numeric vector of size \eqn{J} containing the
#'   evaluation of the adjusted p-value function on the same uniform grid used
#'   to evaluate the functional sample.
#'   - `correction_method`: A string containing the correction method used to
#'   compute the adjusted p-value function.
#'
#'   Optionally, the list may contain the following component:
#'
#'   - `pvalue_matrix`: A numeric matrix of shape \eqn{p \times p} containing
#'   the p-values of the interval-wise tests. Element \eqn{i, j} contains the
#'   p-value of the test performed on the interval indexed by
#'   \eqn{j, j+1, \dots, j+(p-i)}. Only present if the `correction` argument
#'   is set to `"IWT"`.
#'
#' @seealso [`iwt1()`] for calling directly the IWT test and
#'   [`plot.IWT1()`] for plotting the results.
#'
#' @references
#' For the interval-wise testing procedure:
#' - Pini, Alessia, and Simone Vantini. 2016. "The interval testing procedure: a
#' general framework for inference in functional data analysis." Biometrics 72 (3):
#' 835–845.
#' - Pini, Alessia, and Simone Vantini. 2017. "Interval-Wise Testing for Functional
#' Data." Journal of Nonparametric Statistics 29 (2): 407–24.
#'
#' @export
#' @examples
#' # Performing the IWT for one population
#' IWT_result <- functional_one_sample_test(
#'   NASAtemp$paris, mu = 4, n_perm = 10L
#' )
#'
#' # Plotting the results of the IWT
#' IWTimage(IWT_result, abscissa_range = c(0, 12))
#'
#' # Selecting the significant components at 5% level
#' which(IWT_result$adjusted_pvalues < 0.05)
functional_one_sample_test <- function(
  data,
  correction = c("IWT"),
  mu = 0,
  dx = NULL,
  n_perm = 1000L,
  verbose = FALSE,
  aggregation_strategy = c("integral", "max"),
  recycle = TRUE
) {
  correction <- rlang::arg_match(correction)

  if (verbose) {
    cli::cli_h1("Data preparation and point-wise testing")
  }

  prepped_data <- os_prepare_data(
    data = data,
    mu = mu,
    dx = dx,
    n_perm = n_perm
  )

  data_eval <- prepped_data$data
  mu_eval <- prepped_data$mu
  p <- prepped_data$p

  t0 <- prepped_data$t0
  t_coeff <- prepped_data$t_coeff
  pval <- prepped_data$pval

  if (verbose) {
    switch(
      correction,
      IWT = cli::cli_h1("P-Value Adjustment via Interval-Wise Testing")
    )
  }

  aggregation_strategy <- rlang::arg_match(aggregation_strategy)
  adjustment_results <- switch(
    correction,
    IWT = p_adjust_iwt(
      p = p,
      pval = pval,
      t0 = t0,
      t_coeff = t_coeff,
      n_perm = n_perm,
      recycle = recycle,
      aggregation_strategy = aggregation_strategy
    )
  )

  out <- list(
    data = data_eval,
    mu = mu_eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = adjustment_results$adjusted_pvalues,
    correction_method = correction
  )

  if (correction == "IWT") {
    out$pvalue_matrix <- adjustment_results$pvalue_matrix
  }

  class(out) <- "fos"
  out
}
