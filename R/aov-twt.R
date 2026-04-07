#' Threshold Wise Testing procedure for testing functional analysis of variance
#'
#' The function implements the Threshold Wise Testing procedure for testing mean
#' differences between several functional populations in a one-way or multi-way
#' functional analysis of variance framework. Functional data are tested locally
#' and unadjusted and adjusted p-value functions are provided. The unadjusted
#' p-value function controls the point-wise error rate. The adjusted p-value
#' function controls the threshold-wise error rate.
#'
#' @inherit functional_anova_test params return seealso
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. *Biometrics* 79(2), 1119-1132.
#'
#' D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported
#' Significance Levels. *Journal of Business & Economic Statistics* 1.4,
#' 292-298.
#'
#' B. F. J. Manly (2006). Randomization, *Bootstrap and Monte Carlo Methods
#' in Biology*. Vol. 70. CRC Press.
#'
#' @export
#' @examples
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the TWT
#' TWT_result <- TWTaov(temperature ~ groups, B = 100L)
#'
#' # Summary of the TWT results
#' summary(TWT_result)
#'
#' # Plot of the TWT results
#' layout(1)
#' plot(TWT_result)
#'
#' # All graphics on the same device
#' layout(matrix(1:4, nrow = 2, byrow = FALSE))
#' plot(
#'   TWT_result,
#'   main = "NASA data",
#'   plot_adjpval = TRUE,
#'   xlab = "Day",
#'   xrange = c(1, 365)
#' )
TWTaov <- # nolint: object_name_linter.
  function(
    formula,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    method = c("residuals", "responses")
  ) {
    twt_aov(
      formula = formula,
      dx = dx,
      n_perm = B,
      method = method
    )
  }

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname TWTaov
#' @export
twt_aov <- function(
  formula,
  dx = NULL,
  n_perm = 1000L,
  method = c("residuals", "responses")
) {
  method <- rlang::arg_match(method)
  functional_anova_test(
    formula = formula,
    correction = "TWT",
    dx = dx,
    B = n_perm,
    method = method
  )
}
