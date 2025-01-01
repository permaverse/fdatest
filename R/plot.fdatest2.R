#' Plot method for fdatest results on two-population test
#' 
#' \code{plot} method for class "\code{fdatest2}". Plotting function creating a
#' graphical output of the testing procedures for the mean comparison of two
#' groups: functional data and adjusted p-values are plotted.
#' 
#' @param x  The object to be plotted. An object of class "\code{fdatest2}",
#'   usually, a result of a call to \code{\link{fdatest2}}, \code{\link{IWT2}},
#'   \code{\link{TWT2}}, \code{\link{FDR2}}, or \code{\link{PCT2}}.
#' @param xrange Range of the \code{x} axis.
#' @param alpha1 First level of significance used to select and display
#'   significant effects. Default is \code{alpha1 = 0.05}.
#' @param alpha2 Second level of significance used to select and display
#'   significant effects. Default is \code{alpha1 = 0.01}. \code{alpha1} and
#'   \code{alpha2} are s.t. \code{alpha2 < alpha1}. Otherwise the two values are
#'   switched.
#' @param ylab Label of \code{y} axis of the plot of functional data. Default is
#'   "\code{Functional Data}".
#' @param main Plot title.
#' @param lwd Line width for the plot of the adjusted p-value function. Default
#'   is \code{lwd=1}.
#' @param ylim Range of the \code{y} axis. Default is \code{NULL}, giving a plot
#'   with automatic range for functional data.
#' @param col Colors for the plot of functional data. Default is \code{col = 1}.
#' @param type line type for the plot of the adjusted p-value function. Default
#'   is type='l'.
#' @param ... Additional plotting arguments that can be used with function
#'   \code{plot}, such as \code{\link{graphical parameters}} (see
#'   \code{\link{par}}).
#' 
#' @return No value returned. The function produces a graphical output of the
#'   fdatest results:  the plot of the functional data and the one of the
#'   adjusted p-values. The portions of the domain selected as significant by
#'   the test at level \code{alpha1} and \code{alpha2} are highlighted in the
#'   plot of the adjusted p-value function and in the one of functional data by
#'   gray areas (light and dark gray, respectively).
#' 
#' @seealso \code{\link{IWTimage}} for the plot of p-values heatmaps (for IWT).
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
#' # Performing the TWT for two populations
#' TWT.result <- fdatest2(
#'   NASAtemp$paris, NASAtemp$milan, 
#'   method = "TWT", B = 10L
#' )
#'
#' # Plotting the results of the TWT
#' plot(
#'   TWT.result, 
#'   xrange = c(0, 12), 
#'   main = 'TWT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(TWT.result$adjusted_pval < 0.05)
#' 
#' # Performing the IWT for two populations
#' IWT.result <- fdatest2(
#'   NASAtemp$paris, NASAtemp$milan, 
#'   method = "IWT", B = 10L
#' )
#'
#' # Plotting the results of the IWT
#' plot(
#'   IWT.result, 
#'   xrange = c(0, 12), 
#'   main = 'IWT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(IWT.result$adjusted_pval < 0.05)
plot.fdatest2 <- function(x, xrange = c(0,1),
                      alpha1 = 0.05, alpha2 = 0.01,
                      ylab = 'Functional Data', main = NULL, 
                      lwd = 0.5, col=c(1,2), 
                      ylim = NULL, type='l', ...) {
  if (!inherits(x, "fdatest2")) stop("First argument is not a fdatest2 object.")
  if (alpha1 < alpha2) {
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  object <- x
  n <- dim(t(object$data.eval))[1]
  
  colors <- numeric(n)
  id_pop1 <- unique(object$ord_labels)[1]
  id_pop2 <- unique(object$ord_labels)[2]
  colors[which(object$ord_labels == id_pop1)] <- col[1]
  colors[which(object$ord_labels == id_pop2)] <- col[2]
  
  grDevices::devAskNewPage(ask = TRUE) 
  
  p <- length(object$unadjusted_pval)
  xmin <- xrange[1]
  xmax <- xrange[2]
  abscissa_pval = seq(xmin, xmax, len = p)
  main_data <- paste(main, ': Functional Data')
  main_data <- sub("^ : +", "", main_data)
  n_coeff <- dim(object$data.eval)[2]
  data_eval <- object$data.eval
  if (is.null(ylim)) ylim <- range(data_eval,na.rm=TRUE)
  fda::matplot(abscissa_pval, t(data_eval), type = 'l', main = main_data, 
          ylab = ylab, col = colors, lwd = lwd, ylim = ylim, ...)
  mean1 = colMeans(object$data.eval[which(object$ord_labels==id_pop1),],na.rm=TRUE)
  mean2 = colMeans(object$data.eval[which(object$ord_labels==id_pop2),],na.rm=TRUE)
  
  difference1 <- which(object$adjusted_pval < alpha1)
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min_rect <- abscissa_pval[difference1[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], col = "gray90", 
           density = -2, border = NA)
    }
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], 
         col = NULL, border = "black")
  }
  difference2 <- which(object$adjusted_pval < alpha2)
  if (length(difference2) > 0) {
    for (j in 1:length(difference2)) {
      min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], col = "gray80", 
           density = -2, border = NA)
    }
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], 
         col = NULL, border = "black")
  }
  fda::matplot(abscissa_pval, t(data_eval), type = 'l', main = main_data,
          ylab = ylab, col = colors, lwd = lwd, add = TRUE, ...)
  #graphics::matlines(abscissa_pval,cbind(mean1,mean2),col=col,lwd=2,lty=1)
  
  #  adjusted p-values
  main_p <- paste(main,': Adjusted p-values')
  main_p <- sub("^ : +", "", main_p)
  plot(abscissa_pval, object$adjusted_pval, ylim = c(0, 1),
       main = main_p, ylab = 'p-value', type=type, lwd=lwd,...)
  difference1 <- which(object$adjusted_pval < alpha1)
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min_rect <- abscissa_pval[difference1[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], col = "gray90", 
           density = -2, border = NA)
    }
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], 
         col = NULL, border = "black")
  }
  difference2 <- which(object$adjusted_pval < alpha2)
  if (length(difference2) > 0) {
    for (j in 1:length(difference2)) {
      min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], col = "gray80", 
           density = -2, border = NA)
    }
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
  }
  for (j in 0:10) {
    graphics::abline(h = j / 10, col = 'lightgray', lty = "dotted")
  }
  graphics::points(abscissa_pval, object$adjusted_pval, type=type,lwd=2)
  
  grDevices::devAskNewPage(ask = FALSE)
}
