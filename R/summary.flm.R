#' Summarizing Functional Linear Model Fits
#'
#' `summary` method for class `flm`. Returns a summary of the results of the
#' local testing procedure for a functional-on-scalar linear model: the minimum
#' adjusted p-values of the F-test on the whole model and the t-tests on each
#' covariate are reported.
#'
#' @param object  An object of class `flm`, usually, a result of a
#'   call to [`functional_lm_test()`].
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list of summary statistics of the fitted functional linear model
#'   given in `object`, using the component `call` from its arguments, plus:
#'
#'   - `ttest`: A \eqn{(L+1) \times 2} data frame with one row per model term
#'   (intercept plus each predictor) reporting the minimum adjusted p-value of
#'   the corresponding t-test and a significance code.
#'   - `R2`: A \eqn{2 \times 1} matrix giving the range of the functional
#'   R-squared.
#'   - `ftest`: A \eqn{1 \times 2} data frame reporting the minimum adjusted
#'   p-value of the functional F-test and a significance code.
#'
#' @seealso [`IWTimage()`] for the plot of p-value heatmaps and [`plot.flm()`]
#'   for the plot of functional linear model results.
#'
#' @references
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data.
#' \emph{Journal of Nonparametric Statistics}, 29(2), 407-424.
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
#' @export
#' @examples
#' temperature <- rbind(NASAtemp$milan[, 1:100], NASAtemp$paris[, 1:100])
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the IWT
#' IWT_result <- functional_lm_test(
#'   temperature ~ groups,
#'   B = 2L,
#'   correction = "IWT"
#' )
#'
#' # Summary of the IWT results
#' summary(IWT_result)
summary.flm <- function(object, ...) {
  printresult <- vector("list")
  printresult$call <- object$call
  printresult$ttest <- matrix(
    data = apply(object$adjusted_pval_part, 1, min),
    ncol = 1
  )
  var_names <- rownames(object$adjusted_pval_part)
  rownames(printresult$ttest) <- var_names
  printresult$ttest <- as.data.frame(printresult$ttest)
  signif <- rep("", length(var_names))
  signif[which(printresult$ttest[, 1] < 0.001)] <- "***"
  signif[which(
    printresult$ttest[, 1] < 0.01 &
      printresult$ttest[, 1] >= 0.001
  )] <- "**"
  signif[which(
    printresult$ttest[, 1] < 0.05 &
      printresult$ttest[, 1] >= 0.01
  )] <- "*"
  signif[which(
    printresult$ttest[, 1] < 0.1 &
      printresult$ttest[, 1] >= 0.05
  )] <- "."
  printresult$ttest[, 2] <- signif
  colnames(printresult$ttest) <- c("Minimum p-value", "")

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
