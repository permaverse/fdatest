#' One population Interval Wise Testing procedure
#'
#' The function implements the Interval Wise Testing procedure for testing the
#' center of symmetry of a functional population. Functional data are tested
#' locally and  unadjusted  and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the interval-wise error rate.
#'
#' @inherit functional_one_sample_test params return
#'
#' @seealso [`functional_one_sample_test()`] for the interface function and
#'   [`plot.fos()`] and [`IWTimage()`] for plotting the results.
#'
#' @references
#' - Pini, Alessia, and Simone Vantini. 2016. "The interval testing procedure: a
#' general framework for inference in functional data analysis." Biometrics 72 (3):
#' 835–845.
#' - Pini, Alessia, and Simone Vantini. 2017. "Interval-Wise Testing for Functional
#' Data." Journal of Nonparametric Statistics 29 (2): 407–24.
#'
#' @export
#' @examples
#' # Performing the IWT for one population
#' IWT_result <- iwt1(NASAtemp$paris, mu = 4, n_perm = 10L)
#'
#' # Selecting the significant components at 5% level
#' which(IWT_result$adjusted_pvalues < 0.05)
iwt1 <- function(
  data,
  mu = 0,
  n_perm = 1000L,
  dx = NULL,
  recycle = TRUE,
  verbose = FALSE,
  aggregation_strategy = c("integral", "max")
) {
  functional_one_sample_test(
    data = data,
    mu = mu,
    dx = dx,
    n_perm = n_perm,
    verbose = verbose,
    correction = "IWT",
    aggregation_strategy = aggregation_strategy,
    recycle = recycle
  )
}

#' @param B The number of iterations of the MC algorithm to evaluate the
#'   p-values of the permutation tests. Defaults to `1000L`.
#' @rdname iwt1
#' @export
#' @examples
#' # Performing the IWT for one population
#' IWT_result <- IWT1(NASAtemp$paris, mu = 4, B = 10L)
#'
#' # Plotting the results of the IWT
#' plot(IWT_result, xrange = c(0, 12), main = "Paris temperatures")
#'
#' # Plotting the p-value heatmap
#' IWTimage(IWT_result, abscissa_range = c(0, 12))
#'
#' # Selecting the significant components at 5% level
#' which(IWT_result$adjusted_pval < 0.05)
IWT1 <- # nolint: object_name_linter.
  function(
    data,
    mu = 0,
    B = 1000L, # nolint: object_name_linter.
    dx = NULL,
    recycle = TRUE
  ) {
    result <- iwt1(data = data, mu = mu, n_perm = B, dx = dx, recycle = recycle)
    # Build backward-compatible IWT1 object so that plot.IWT1 and IWTimage
    # continue to work with the legacy class and field names.
    out <- list(
      test = "1pop",
      mu = result$mu,
      adjusted_pval = result$adjusted_pvalues,
      unadjusted_pval = result$unadjusted_pvalues,
      pval_matrix = result$pvalue_matrix,
      data_eval = result$data
    )
    class(out) <- "IWT1"
    out
  }
