#' Global testing procedure for testing functional analysis of variance
#'
#' The function implements the Global Testing procedure for testing mean
#' differences between several functional populations in a one-way or multi-way
#' functional analysis of variance framework. Functional data are tested
#' globally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the family-wise error rate weakly. Since this is a
#' global test, the adjusted p-value function is constant.
#'
#' @inherit functional_anova_test params return seealso
#'
#' @references
#' - Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. *Biometrics* 79(2), 1119-1132.
#' - D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported
#' Significance Levels. *Journal of Business & Economic Statistics* 1.4,
#' 292-298.
#' - B. F. J. Manly (2006). Randomization, *Bootstrap and Monte Carlo
#' Methods in Biology*. Vol. 70. CRC Press.
#'
#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @export
#' @examples
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the test
#' Global_result <- global_aov(temperature ~ groups, n_perm = 1000L)
#'
#' # Summary of the test results
#' summary(Global_result)
#'
#' # Plot of the results
#' layout(1)
#' plot(Global_result)
#'
#' # All graphics on the same device
#' layout(matrix(1:4, nrow = 2, byrow = FALSE))
#' plot(
#'   Global_result,
#'   main = 'NASA data',
#'   plot.adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
global_aov <- function(
  formula,
  dx = NULL,
  n_perm = 1000L,
  method = c("residuals", "responses"),
  stat = c("Integral", "Max")
) {
  method <- rlang::arg_match(method)
  stat <- rlang::arg_match(stat)
  functional_anova_test(
    formula = formula,
    correction = "Global",
    dx = dx,
    B = n_perm,
    method = method,
    stat = stat
  )
}

#' @rdname global_aov
#' @export
Globalaov <- # nolint: object_name_linter.
  function(
    formula,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    method = c("residuals", "responses"),
    stat = c("Integral", "Max")
  ) {
    global_aov(
      formula = formula,
      dx = dx,
      n_perm = B,
      method = method,
      stat = stat
    )
  }
