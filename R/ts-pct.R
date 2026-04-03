#' Two population Partition Closed Testing procedure
#'
#' The function implements the Partition Closed Testing procedure for testing
#' mean differences between two functional populations. Functional data are
#' tested locally and unadjusted and adjusted p-value functions are provided.
#' The unadjusted p-value function controls the point-wise error rate. The
#' adjusted p-value function controls the family-wise error rate asymptotically.
#'
#' @inherit functional_two_sample_test params return
#'
#' @seealso [`global2()`], [`twt2()`], [`iwt2()`], [`fdr2()`] for calling directly
#' one of the other tests, [`functional_two_sample_test()`] for calling the
#' interface test and [`plot.fts()`] for plotting the results.
#'
#' @references
#' - Vsevolozhskaya, Olga A, Mark C Greenwood, GJ Bellante, Scott L Powell, Rick
#' L Lawrence, and Kevin S Repasky. 2013. “Combining Functions and the Closure
#' Principle for Performing Follow-up Tests in Functional Analysis of Variance.”
#' Computational Statistics & Data Analysis 67: 175–84.
#' - Vsevolozhskaya, Olga, Mark Greenwood, and Dmitri Holodov. 2014. “Pairwise
#' comparison of treatment levels in functional analysis of variance with
#' application to erythrocyte hemolysis.” The Annals of Applied Statistics 8 (2):
#' 905–25. https://doi.org/10.1214/14-AOAS723.
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
    correction = "PCT",
    aggregation_strategy = aggregation_strategy,
    partition = partition
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
#' @rdname pct2
#' @export
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
      what = "PCT2()",
      details = "Use pct2() instead. Be mindful that the argument `statistic` has been replaced by `aggregation_strategy` and `standardize`.",
      id = "fdatest-deprecated-pct2"
    )
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

ts_p_adjust_pct <- function(
  partition,
  p,
  t0,
  t_coeff,
  aggregation_strategy
) {
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

  list(
    adjusted_pvalues = adjusted_pval
  )
}
