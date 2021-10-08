#' @export
###################################################################################################################################################
# Nonparametric inference for functional-on-scalar linear models applied to knee kinematic hop data after injury of the anterior cruciate ligament
# Plot IWT results in Linear Models framework
# By: A. Pini
# last modified: 09/06/2016
###################################################################################################################################################

# Inputs: 
# x:            Object of class 'IWTlm' obtained from function IWTlm
# xrange:       Range of the plotted x axis. Default is c(0,1)
# alpha1:       First level of significance used to select and display significant differences. 
#               Default is alpha1 = 0.05.
# alpha2:       Second level of significance used to select and display significant differences. 
#               Default is alpha2 = 0.01.
# plot.adjpval: A logical indicating wether the plots of adjusted p-values have to be done. 
#               Default is plot.adjpval = FALSE.
# col:          Vector of colors for the plot of functional data (first element), and functional coefficients (following elements).
#               Default is c(1,grDevices::rainbow(dim(x$adjusted_pval_part)[1]))
# ylim:         Range of the y axis. If set to NULL, the range of the data evaluations is used.
#               Default is ylim = NULL.
# ylab:         Label of y axis of the plot of functional data. Default is "Functional Data".
# main:         An overall title for the plots (it will be pasted to "Functional Data and F-test" for the first plot 
#               and "t-test" for the other plots).
# lwd:          Line width for the plot of functional data. Default is lwd=1.
# type:         Type of lines to be used for the plot of functional data. Default is type='l'.
# ...:          Additional plotting arguments that can be used with function plot, such as graphical parameters (see par).


# Value:
# No value returned. The function produces a graphical output of the IWT results: 
# the plot of the functional data, functional regression coefficients, and IWT-adjusted p-values 
# for the F-test and t-tests. The basis components selected as significant by the tests at level 
# alpha1 and alpha2 are highlighted in the plot of the corrected p-values and in the one of 
# functional data by gray areas (light and dark gray, respectively). The plot of functional data 
# reports the gray areas corresponding to a significant F-test. The plots of functional regression 
# coefficients report the gray areas corresponding to significant t-tests for the corresponding 
# covariate.

plot.IWTlm <- function(x, 
                       xrange = c(0,1), 
                       alpha1 = 0.05, 
                       alpha2 = 0.01, 
                       plot_adjpval = FALSE,
                       col = c(1,grDevices::rainbow(dim(x$adjusted_pval_part)[1])), 
                       ylim = NULL,
                       ylab='Functional Data', 
                       main=NULL, 
                       lwd = 1, 
                       type='l',
                       ...) {
  if (class(x) != "IWTlm") stop("x should be an object of the class ITPlm")
  if (alpha1 < alpha2) {
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  object <- x
  p <- length(object$unadjusted_pval_F)
  J <- p
  n <- dim(object$data.eval)[1]
  xmin <- xrange[1]
  xmax <- xrange[2]
  abscissa_pval <- seq(xmin, xmax, len = p)
  abscissa_smooth <- seq(xmin, xmax, len = J)
  grDevices::devAskNewPage(ask = TRUE)  
  main_F <- paste(main,': Functional Data and F-test')
  main_F <- sub("^ : +", "", main_F)
  fda::matplot(abscissa_smooth, t(object$data.eval), type = 'l',col = NA, main = main_F,
          ylab = ylab, ylim = ylim, lwd = lwd, ...)
  difference1 <- which(object$adjusted_pval_F < alpha1)
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
  difference2 <- which(object$adjusted_pval_F < alpha2)
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
  fda::matplot(abscissa_smooth, t(object$data.eval), type='l', col = col[1],
          add = TRUE, lwd = lwd,...)
  
  for(var in 1:(dim(object$adjusted_pval_part)[1])){
    var_name <- rownames(object$adjusted_pval_part)[var]
    main_t <- paste(main, ': t-test -', var_name, sep = ' ')
    main_t <- sub("^ : +", "", main_t)
    plot(abscissa_smooth, object$coeff.regr.eval[var,], type = 'l', 
         col = 0, ylim = range(c(0, object$coeff.regr.eval[var,])),
         lwd = 2, main = main_t, ylab = 'Regression Coefficient', ...)
    difference1 <- which(object$adjusted_pval_part[var, ] < alpha1)
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
    difference2 <- which(object$adjusted_pval_part[var,] < alpha2)
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
    graphics::lines(abscissa_smooth, object$coeff.regr.eval[var,], type = 'l', 
          col = col[var + 1], lwd = 2, ...)
    graphics::abline(h = 0, lty = 2, col = 1)
  }
  # Plot adjusted p-values
  if (plot_adjpval == TRUE) {
    main_p <- paste(main, ': Adjusted p-values - F-test')
    main_p <- sub("^ : +", "", main_p)
    plot(abscissa_pval, object$adjusted_pval_F, type=type,lwd=2, ylim = c(0,1), 
         main = main_p, ylab = 'p-value', col=0,...)
    difference1 <- which(object$adjusted_pval_F < alpha1)
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
    difference2 <- which(object$adjusted_pval_F < alpha2)
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
    for (j in 0:10) {
      graphics::abline(h = j / 10, col = 'lightgray', lty = "dotted")
    }
    graphics::lines(abscissa_pval, object$adjusted_pval_F,lwd=2, ...)
    for(var in 1:(dim(object$adjusted_pval_part)[1])){
      var_name = rownames(object$adjusted_pval_part)[var]
      main_p <- paste(main, ': Adjusted p-values - t-test -', var_name)
      main_p <- sub("^ : +", "", main_p)
      plot(abscissa_pval, object$adjusted_pval_part[var,], ylim = c(0, 1),
           main = main_p, ylab = 'p-value',col=NA, ...)
      difference1 <- which(object$adjusted_pval_part[var,] < alpha1)
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
      difference2 <- which(object$adjusted_pval_part[var,] < alpha2)
      if (length(difference2) > 0) {
        for (j in 1:length(difference2)) {
          min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
          max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
          graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], col = "gray80", density = -2, border = NA)
        }
        graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], 
             col = NULL, border = "black")
      }
      for (j in 0:10) {
        graphics::abline(h = j / 10, col = 'lightgray', lty = "dotted")
      }
      graphics::lines(abscissa_pval, object$adjusted_pval_part[var,], lwd=2, ...)
    }
  }
  grDevices::devAskNewPage(ask = FALSE)   
}
