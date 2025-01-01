#' Summarizing Functional Analysis of Variance Fits
#' 
#' \code{summary} method for class "\code{IWTaov}". Function returning a summary
#' of the results of IWT for the test on a functional analysis of variance:
#' minimum IWT-adjusted p-values of the F-tests on the whole model and on each
#' factor are reported.
#' 
#' @param object  An object of class "\code{IWTaov}", usually, a result of a
#'   call to \code{\link{IWTaov}}.
#' @param ... Further arguments passed to or from other methods..
#' 
#' @return No value returned. The function \code{summary.IWTaov} computes and
#'   returns a list of summary statistics of the fitted functional analysis of
#'   variance given in \code{object}, using the component "\code{call}" from its
#'   arguments, plus:
#'   
#'   - `factors`: A \code{L x 1} matrix with columns for the factors of ANOVA,
#'   and corresponding (two-sided) IWT-adjusted minimum p-values of the
#'   corresponding tests of significance (i.e., the minimum p-value over all
#'   \code{p} basis components used to describe functional data).
#'   - `R2`: Range of the functional R-squared.
#'   - `ftest`: IWT-adjusted minimum p-value of functional F-test.
#' 
#' @seealso \code{\link{IWTimage}} for the plot of p-values heatmaps.
#'   \code{\link{plot.IWTaov}} for the plot of analysis of variance results. See
#'   also \code{\link{IWT1}}, \code{\link{IWT2}} to perform the ITP to test on
#'   the mean of one population and test of differences between two populations.
#'   See \code{\link{ITPaovbspline}} for functional ANOVA based on B-spline
#'   basis representation
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
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#' 
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#' 
#' # Performing the IWT
#' IWT.result <- IWTaov(temperature ~ groups, B = 10L)
#' 
#' # Summary of the IWT results
#' summary(IWT.result)
#' 
#' # Plot of the IWT results
#' layout(1)
#' plot(IWT.result)
#' 
#' # All graphics on the same device
#' layout(matrix(1:4, nrow = 2, byrow = FALSE))
#' plot(
#'   IWT.result, 
#'   main = 'NASA data', 
#'   plot_adjpval = TRUE, 
#'   xlab = 'Day', 
#'   xrange = c(1, 365)
#' )
summary.IWTaov <- function(object, ...) {
  printresult <- vector('list')
  printresult$call <- object$call
  printresult$factors <- matrix(data = apply(object$adjusted_pval_factors, 1, min),
                                ncol = 1)
  var_names <- rownames(object$adjusted_pval_factors)
  rownames(printresult$factors) <- var_names
  printresult$factors <- as.data.frame(printresult$factors)
  signif <- rep('', length(var_names))
  signif[which(printresult$factors[, 1] < 0.001)] <- '***'
  signif[which(printresult$factors[, 1] < 0.01 & printresult$factors[, 1] >= 0.001)] <- '**'
  signif[which(printresult$factors[, 1] < 0.05 & printresult$factors[, 1] >= 0.01)] <- '*'
  signif[which(printresult$factors[, 1] < 0.1 & printresult$factors[, 1] >= 0.05)] <- '.'
  printresult$factors[,2] <- signif
  colnames(printresult$factors) <- c('Minimum p-value', '')
  printresult$R2 <- as.matrix(range(object$R2.eval))
  colnames(printresult$R2) <- 'Range of functional R-squared'
  rownames(printresult$R2) <- c('Min R-squared', 'Max R-squared')
  printresult$ftest <- as.matrix(min(object$adjusted_pval_F))
  printresult$ftest <- as.data.frame(printresult$ftest)
  signif.f <- ''
  signif.f[which(printresult$ftest[, 1] < 0.001)] <- '***'
  signif.f[which(printresult$ftest[, 1] < 0.01 & printresult$ftest[, 1] >= 0.001)] <- '**'
  signif.f[which(printresult$ftest[, 1] < 0.05 & printresult$ftest[, 1] >= 0.01)] <- '*'
  signif.f[which(printresult$ftest[, 1] < 0.1 & printresult$ftest[, 1] >= 0.05)] <- '.'
  printresult$ftest[, 2] <- signif.f
  colnames(printresult$ftest) <- c('Minimum p-value', '')
  printresult
}
