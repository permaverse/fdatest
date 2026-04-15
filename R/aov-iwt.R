#' Interval Wise Testing procedure for testing functional analysis of variance
#'
#' The function implements the Interval Wise Testing procedure for testing mean
#' differences between several functional populations in a one-way or multi-way
#' functional analysis of variance framework. Functional data are tested locally
#' and unadjusted and adjusted p-value functions are provided. The unadjusted
#' p-value function controls the point-wise error rate. The adjusted p-value
#' function controls the interval-wise error rate.
#'
#' @inherit functional_anova_test params return seealso
#'
#' @references
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data.
#' *Journal of Nonparametric Statistics*, 29(2), 407-424.
#'
#' Pini, A., Vantini, S., Colosimo, B. M., & Grasso, M. (2018). Domain‐selective
#' functional analysis of variance for supervised statistical profile monitoring
#' of signal data. *Journal of the Royal Statistical Society: Series C
#' (Applied Statistics)* 67(1), 55-81.
#'
#' Abramowicz, K., Hager, C. K., Pini, A., Schelin, L., Sjostedt de Luna, S., &
#' Vantini, S. (2018). Nonparametric inference for functional‐on‐scalar linear
#' models applied to knee kinematic hop data after injury of the anterior
#' cruciate ligament. *Scandinavian Journal of Statistics* 45(4),
#' 1036-1061.
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
#' # Performing the IWT
#' IWT_result <- IWTaov(temperature ~ groups, B = 10L)
#'
#' # Summary of the IWT results
#' summary(IWT_result)
#'
#' # Plot of the IWT results
#' graphics::layout(1)
#' plot(IWT_result)
#'
#' # All graphics on the same device
#' graphics::layout(matrix(1:4, nrow = 2, byrow = FALSE))
#' plot(
#'   IWT_result,
#'   main = "NASA data",
#'   plot.adjpval = TRUE,
#'   xlab = "Day",
#'   xrange = c(1, 365)
#' )
IWTaov <- # nolint: object_name_linter.
  function(
    formula,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    method = c("residuals", "responses"),
    recycle = TRUE
  ) {
    iwt_aov(
      formula = formula,
      dx = dx,
      n_perm = B,
      method = method,
      recycle = recycle
    )
  }

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname IWTaov
#' @export
iwt_aov <- function(
  formula,
  dx = NULL,
  n_perm = 1000L,
  method = c("residuals", "responses"),
  recycle = TRUE
) {
  method <- rlang::arg_match(method)
  functional_anova_test(
    formula = formula,
    correction = "IWT",
    dx = dx,
    B = n_perm,
    method = method,
    recycle = recycle
  )
}
