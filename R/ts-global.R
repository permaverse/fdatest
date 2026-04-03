#' Two population Global Testing procedure
#'
#' The function implements the Global Testing procedure for testing mean
#' differences between two functional populations. Functional data are tested
#' locally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the interval-wise error rate.
#'
#' @inherit functional_two_sample_test params return
#'
#' @seealso [`iwt2()`], [`twt2()`], [`pct2()`], [`fdr2()`] for calling directly
#' one of the other tests, [`functional_two_sample_test()`] for calling the
#' interface test and [`plot.fts()`] for plotting the results.
#'
#' @references
#' - Hall, Peter, and Nader Tajvidi. 2002. “Permutation Tests for Equality of
#' Distributions in High-Dimensional Settings.” Biometrika 89 (2): 359–74.
#' - Pini, Alessia, Aymeric Stamm, and Simone Vantini. 2018. “Hotelling’s T2 in
#' Separable Hilbert Spaces.” Journal of Multivariate Analysis 167: 284–305.
#'
#' @export
#' @examples
#' # Performing the Global for two populations
#' Global_result <- global2(NASAtemp$paris, NASAtemp$milan)
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
#' @rdname global2
#' @export
Global2 <- # nolint: object_name_linter.
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
      what = "Global2()",
      details = "Use global2() instead. Be mindful that the argument `statistic` has been replaced by `aggregation_strategy` and `standardize`.",
      id = "fdatest-deprecated-global2"
    )
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

ts_p_adjust_global <- function(
  aggregation_strategy,
  t0,
  t_coeff,
  p
) {
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

  list(
    adjusted_pvalues = adjusted_pval
  )
}
