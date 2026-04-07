#' Global testing procedure for testing functional-on-scalar linear models
#'
#' The function is used to fit and test functional linear models. It can be used
#' to carry out regression, and analysis of variance. It implements the global
#' testing procedure for testing the significance of the effects of scalar
#' covariates on a functional population.
#'
#' @inherit functional_lm_test params return seealso
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. \emph{Biometrics} 79(2), 1119-1132.
#'
#' D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported
#' Significance Levels. \emph{Journal of Business & Economic Statistics} 1(4),
#' 292-298.
#'
#' B. F. J. Manly (2006). Randomization, \emph{Bootstrap and Monte Carlo Methods
#' in Biology}. Vol. 70. CRC Press.
#'
#' @export
#' @examples
#' # Defining the covariates
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the IWT
#' Global_result <- Globallm(temperature ~ groups, B = 1000)
#' # Summary of the IWT results
#' summary(Global_result)
#'
#' # Plot of the IWT results
#' plot(
#'   Global_result,
#'   main = 'NASA data',
#'   plot_adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
#'
#' plot(
#'   Global_result,
#'   main = 'NASA data',
#'   plot_adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
Globallm <- # nolint: object_name_linter.
  function(
    formula,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    method = c("residuals", "responses"),
    stat = c("Integral", "Max")
  ) {
    global_lm(
      formula = formula,
      dx = dx,
      n_perm = B,
      method = method,
      stat = stat
    )
  }

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname Globallm
#' @export
global_lm <- function(
  formula,
  dx = NULL,
  n_perm = 1000L,
  method = c("residuals", "responses"),
  stat = c("Integral", "Max")
) {
  method <- rlang::arg_match(method)
  stat <- rlang::arg_match(stat)
  functional_lm_test(
    formula = formula,
    correction = "Global",
    dx = dx,
    B = n_perm,
    method = method,
    stat = stat
  )
}
