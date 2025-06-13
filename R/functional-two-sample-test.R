#' Two-sample local testing procedures
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
#' @inheritParams TWTaov
#' @param paired A boolean value specifying whether a paired test should be
#'   performed. Defaults to `FALSE`.
#' @param alternative A string specifying the type of alternative hypothesis.
#'   Choices are `"two.sided"`, `"less"` or `"greater"`. Defaults to
#'   `"two.sided"`.
#' @param statistic A string speicyfing the test statistic to use. Possible
#'   values are:
#'
#'   - `"Integral"`: Integral of the squared sample mean difference.
#'   - `"Max"`: Maximum of the squared sample mean difference.
#'   - `"Integral_std"`: Integral of the squared t-test statistic.
#'   - `"Max_std"`: Maximum of the squared t-test statistic.
#'
#'   Defaults to `"Integral"`.
#' @inheritParams IWTaov
#' @param partition An integer vector of length \eqn{J} specifying the
#'   membership of each point of the domain to an element of the partition.
#'   Only used and **must** be set if the `correction` argument is set to
#'   `"PCT"`.
#' @param verbose A boolean value specifying whether to print the progress of
#'  the computation. Defaults to `FALSE`.
#'
#' @returns An object of class `ftwosample` containing the following components:
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
#' @seealso See also [`plot.ftwosample()`] for plotting the results.
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
#' TWT.result <- functional_two_sample_test(
#'   NASAtemp$paris, NASAtemp$milan,
#'   correction = "TWT", B = 10L
#' )
#'
#' # Plotting the results of the TWT
#' plot(
#'   TWT.result,
#'   xrange = c(0, 12),
#'   title = 'TWT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(TWT.result$adjusted_pval < 0.05)
#'
#' # Performing the IWT for two populations
#' IWT.result <- functional_two_sample_test(
#'   NASAtemp$paris, NASAtemp$milan,
#'   correction = "IWT", B = 10L
#' )
#'
#' # Plotting the results of the IWT
#' plot(
#'   IWT.result,
#'   xrange = c(0, 12),
#'   title = 'IWT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(IWT.result$adjusted_pval < 0.05)
functional_two_sample_test <- function(
  data1,
  data2,
  correction,
  mu = 0,
  dx = NULL,
  B = 1000L,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  statistic = c("Integral", "Max", "Integral_std", "Max_std"),
  recycle = TRUE,
  partition = NULL,
  verbose = FALSE
) {
  if (correction == "PCT" && is.null(partition)) {
    cli::cli_abort(
      'When the {.arg correction} argument is set to {.code "PCT"}, the {.arg
      partition} argument should be explicitly provided.'
    )
  }

  out <- switch(
    correction,
    IWT = IWT2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      dx = dx,
      B = B,
      paired = paired,
      alternative = alternative,
      verbose = verbose,
      recycle = recycle
    ),
    TWT = TWT2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      dx = dx,
      B = B,
      paired = paired,
      alternative = alternative,
      verbose = verbose
    ),
    PCT = PCT2(
      data1 = data1,
      data2 = data2,
      partition = partition,
      mu = mu,
      dx = dx,
      B = B,
      paired = paired,
      alternative = alternative
    ),
    Global = Global2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      dx = dx,
      B = B,
      paired = paired,
      statistic = statistic
    ),
    FDR = FDR2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      dx = dx,
      B = B,
      paired = paired,
      alternative = alternative
    )
  )

  out$correction <- correction
  out
}
