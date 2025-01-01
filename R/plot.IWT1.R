#' Plot method for IWT results on one-population test
#' 
#' \code{plot} method for class "\code{IWT1}". Plotting function creating a
#' graphical output of the IWT for the test of the mean of one population:
#' functional data and IWT-adjusted p-values are plotted.
#' 
#' @param x  The object to be plotted. An object of class "\code{IWT1}",
#'   usually, a result of a call to \code{\link{IWT1}}.
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
#'   IWT results:  the plot of the functional data and the one of the adjusted
#'   p-values. The portions of the domain selected as significant by the test at
#'   level \code{alpha1} and \code{alpha2} are highlighted in the plot of the
#'   adjusted p-value function and in the one of functional data by gray areas
#'   (light and dark gray, respectively).
#' 
#' @seealso \code{\link{IWTimage}} for the plot of p-values heatmaps. See also
#'   \code{\link{IWT2}} to perform the ITP to test differences between two
#'   populations. See \code{\link{ITP1bspline}} for one-population test based on
#'   B-spline basis representation.
#' 
#' @references
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data.
#' \emph{Journal of Nonparametric Statistics}, 29(2), 407-424
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
#' # Performing the IWT for one population
#' IWT.result <- IWT1(NASAtemp$paris, mu = 4, B = 10L)
#'
#' # Plotting the results of the IWT
#' plot(IWT.result, xrange = c(0, 12), main = 'Paris temperatures')
#'
#' # Plotting the p-value heatmap
#' IWTimage(IWT.result,abscissa_range=c(0,12))
#'
#' # Selecting the significant components at 5% level
#' which(IWT.result$adjusted_pval < 0.05)
plot.IWT1 <- function(x, xrange = c(0,1),
                      alpha1 = 0.05, alpha2 = 0.01,
                      ylab = 'Functional Data', main = NULL, 
                      lwd = 1, col = 1, 
                      ylim = NULL,type='l', ...) {
  if (!inherits(x, "IWT1")) stop("First argument is not a IWT1 object.")
  if (alpha1 < alpha2) {
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  object <- x
  
  grDevices::devAskNewPage(ask = TRUE)
  p <- length(object$unadjusted_pval)
  n <- dim((object$data.eval))[1]
  xmin <- xrange[1]
  xmax <- xrange[2]
  abscissa_pval <- seq(xmin, xmax, len = p)
  main_data <- paste(main, ': Functional Data')
  main_data <- sub("^ : +", "", main_data)
  n_coeff <- dim(object$data.eval)[2]
  data_eval <- object$data.eval
  if (is.null(ylim)) ylim <- range(data_eval)
  fda::matplot(abscissa_pval, t(data_eval), type='l', main = main_data, 
          ylab = ylab, col = col, lwd = lwd, ylim = ylim, ...)
  difference1 <- which(object$adjusted_pval < alpha1)
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min_rect <- abscissa_pval[difference1[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
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
      min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], col = "gray80", 
           density = -2, border = NA)
    }
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], 
         col = NULL, border = "black")
  }
  fda::matplot(abscissa_pval, t(data_eval), type = 'l', main = main_data, 
          ylab = ylab, col = col, lwd = lwd, add = TRUE, ...)
  if (length(object$mu) == 1) { # mu is a constant function
    mu_eval <- rep(object$mu, p)
  } else { # mu is a functional data with no constant coefficients
    mu <- object$mu
    mu_eval <- mu
  }
  abscissa_mu <- abscissa_pval
  graphics::lines(abscissa_mu, mu_eval, col = 'gray', lwd = 2)
  # Plot adjusted p-values
  main_p <- paste(main,': Adjusted p-values')
  main_p <- sub("^ : +", "", main_p)
  plot(abscissa_pval, object$adjusted_pval, ylim = c(0, 1), 
       main = main_p, ylab = 'p-value', type=type,lwd=lwd, ...)
  difference1 <- which(object$adjusted_pval < alpha1)
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min_rect <- abscissa_pval[difference1[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
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
      min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], col = "gray80", 
           density = -2, border = NA)
    }
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], 
         col = NULL, border = "black")
  }
  for (j in 0:10) {
    graphics::abline(h = j/10, col = 'lightgray', lty = "dotted")
  }
  graphics::points(abscissa_pval, object$adjusted_pval,type=type,lwd=lwd)
  
  
  grDevices::devAskNewPage(ask = FALSE)
}