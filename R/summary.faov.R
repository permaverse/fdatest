#' Summarizing Functional Analysis of Variance Fits
#'
#' `summary` method for class `faov`. Returns a summary of the results of the
#' local testing procedure for functional analysis of variance: the minimum
#' adjusted p-values of the F-test on the whole model and on each factor are
#' reported.
#'
#' @param object  An object of class `faov`, usually, a result of a
#'   call to [`functional_anova_test()`].
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list of summary statistics of the fitted functional analysis of
#'   variance given in `object`, using the component `call` from its arguments,
#'   plus:
#'
#'   - `factors`: A \eqn{L \times 2} data frame with one row per factor
#'   reporting the minimum adjusted p-value of the corresponding F-test and a
#'   significance code.
#'   - `R2`: A \eqn{2 \times 1} matrix giving the range of the functional
#'   R-squared.
#'   - `ftest`: A \eqn{1 \times 2} data frame reporting the minimum adjusted
#'   p-value of the global F-test and a significance code.
#'
#' @seealso [`IWTimage()`] for the plot of p-value heatmaps and [`plot.faov()`]
#'   for the plot of analysis of variance results.
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
#' @export
#' @examples
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the IWT
#' IWT_result <- functional_anova_test(
#'   temperature ~ groups,
#'   B = 10L,
#'   correction = "IWT"
#' )
#'
#' # Summary of the IWT results
#' summary(IWT_result)
summary.faov <- function(object, ...) {
  printresult <- vector("list")
  printresult$call <- object$call
  printresult$factors <- matrix(
    data = apply(object$adjusted_pval_factors, 1, min),
    ncol = 1
  )
  var_names <- rownames(object$adjusted_pval_factors)
  rownames(printresult$factors) <- var_names
  printresult$factors <- as.data.frame(printresult$factors)
  signif <- rep("", length(var_names))
  signif[which(printresult$factors[, 1] < 0.001)] <- "***"
  signif[which(
    printresult$factors[, 1] < 0.01 &
      printresult$factors[, 1] >= 0.001
  )] <- "**"
  signif[which(
    printresult$factors[, 1] < 0.05 &
      printresult$factors[, 1] >= 0.01
  )] <- "*"
  signif[which(
    printresult$factors[, 1] < 0.1 &
      printresult$factors[, 1] >= 0.05
  )] <- "."
  printresult$factors[, 2] <- signif
  colnames(printresult$factors) <- c("Minimum p-value", "")
  printresult$R2 <- as.matrix(range(object$R2_eval))
  colnames(printresult$R2) <- "Range of functional R-squared"
  rownames(printresult$R2) <- c("Min R-squared", "Max R-squared")
  printresult$ftest <- as.matrix(min(object$adjusted_pval_F))
  printresult$ftest <- as.data.frame(printresult$ftest)
  signif_f <- ""
  signif_f[which(printresult$ftest[, 1] < 0.001)] <- "***"
  signif_f[which(
    printresult$ftest[, 1] < 0.01 &
      printresult$ftest[, 1] >= 0.001
  )] <- "**"
  signif_f[which(
    printresult$ftest[, 1] < 0.05 &
      printresult$ftest[, 1] >= 0.01
  )] <- "*"
  signif_f[which(
    printresult$ftest[, 1] < 0.1 &
      printresult$ftest[, 1] >= 0.05
  )] <- "."
  printresult$ftest[, 2] <- signif_f
  colnames(printresult$ftest) <- c("Minimum p-value", "")
  printresult
}
