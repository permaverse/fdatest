#' Interval-wise testing procedure for testing functional-on-scalar linear
#' models
#'
#' The function is used to fit and test functional linear models. It can be used
#' to carry out regression, and analysis of variance. It implements the
#' interval-wise testing procedure (IWT) for testing the significance of the
#' effects of scalar covariates on a functional population.
#'
#' @inherit functional_lm_test params return seealso
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. \emph{Biometrics} 79(2), 1119-1132.
#' 
#' A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference for
#' Functional Data Controlling the Family Wise Error Rate on Intervals.
#' \emph{Biometrics} 73(3): 835–845.
#'
#' Pini, A., Vantini, S., Colosimo, B. M., & Grasso, M. (2018). Domain‐selective
#' functional analysis of variance for supervised statistical profile monitoring
#' of signal data. \emph{Journal of the Royal Statistical Society: Series C
#' (Applied Statistics)} 67(1), 55-81.
#'
#' Abramowicz, K., Hager, C. K., Pini, A., Schelin, L., Sjostedt de Luna, S., &
#' Vantini, S. (2018). Nonparametric inference for functional‐on‐scalar linear
#' models applied to knee kinematic hop data after injury of the anterior
#' cruciate ligament. \emph{Scandinavian Journal of Statistics} 45(4),
#' 1036-1061.
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
#' temperature <- rbind(NASAtemp$milan[, 1:100], NASAtemp$paris[, 1:100])
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the IWT
#' IWT_result <- IWTlm(temperature ~ groups, B = 2L)
#' # Summary of the IWT results
#' summary(IWT_result)
#'
#' # Plot of the IWT results
#' plot(
#'   IWT_result,
#'   main = "NASA data",
#'   plot_adjpval = TRUE,
#'   xlab = "Day",
#'   xrange = c(1, 365)
#' )
IWTlm <- # nolint: object_name_linter.
  function(
    formula,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    method = c("residuals", "responses"),
    recycle = TRUE
  ) {
    iwt_lm(
      formula = formula,
      dx = dx,
      n_perm = B,
      method = method,
      recycle = recycle
    )
  }

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname IWTlm
#' @export
iwt_lm <- function(
  formula,
  dx = NULL,
  n_perm = 1000L,
  method = c("residuals", "responses"),
  recycle = TRUE
) {
  method <- rlang::arg_match(method)
  functional_lm_test(
    formula = formula,
    correction = "IWT",
    dx = dx,
    B = n_perm,
    method = method,
    recycle = recycle
  )
}
