#' @title Plot method for ITP results on two-population tests
#'
#' @description \code{plot} method for class "\code{ITP2}".
#' Plotting function creating a graphical output of the ITP for the test of comparison between two populations: 
#' functional data and ITP-adjusted p-values are plotted.
#' 
#' @param x The object to be plotted.
#' An object of class "\code{ITP2}", that is, a result of an ITP for comparison between two populations. 
#' Usually a call to \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}} or \code{\link{ITP2pafourier}}.
#' 
#' @param xrange Range of the \code{x} axis. Default is \code{xrange=c(0,1)}.
#' 
#' @param alpha1 First level of significance used to select and display significant differences. Default is \code{alpha1 = 0.05}.
#' 
#' @param alpha2 Second level of significance used to select and display significant differences. Default is \code{alpha1 = 0.01}. \code{alpha1} and \code{alpha2} are s.t. \code{alpha2 < alpha1}. Otherwise the two values are switched.
#' 
#' @param ylab Label of \code{y} axis of the plot of functional data. Default is "\code{Functional Data}".
#' 
#' @param main An overall title for the plots (it will be pasted to "\code{Functional Data}" for the first plot and "\code{adjusted p-values}" for the second plot). 
#' 
#' @param lwd Line width for the plot of functional data. 
#' 
#' @param col Color used to plot the functional data.
#' 
#' @param pch Point character for the plot of adjusted p-values.
#' 
#' @param ylim Range of the \code{y} axis.
#' 
#' @param \dots Additional plotting arguments that can be used with function \code{plot}, such as \code{\link{graphical parameters}} (see \code{\link{par}}).
#' 
#' @return No value returned. 
#' The function produces a graphical output of the ITP results:  the plot of the functional data and the one of the adjusted p-values. 
#' The basis components selected as significant by the test at level \code{alpha1} and \code{alpha2} are highlighted in the plot of the adjusted p-values and in the one of functional data (in case the test is based on a local basis, such as B-splines) by gray areas (light and dark gray, respectively). 
#' In the case of a Fourier basis with amplitude and phase decomposition, two plots of adjusted p-values are done, one for phase and one for amplitude.
#' 
#' @seealso \code{\link{ITPimage}} for the plot of p-values heatmaps. 
#' See also \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, \code{\link{ITP2pafourier}} to perform the ITP to test for differences between two populations. 
#' See \code{\link{plot.ITP1}} and \code{\link{plot.ITPlm}} for the plot method applied to the ITP results of one-population tests and a linear models, respectively.
#' 
#' @examples 
#' # Performing the ITP for two populations with the B-spline basis
#' ITP.result.bspline <- ITP2bspline(NASAtemp$milan,NASAtemp$paris,nknots=30,B=1000)
#' # Plotting the results of the ITP
#' plot(ITP.result.bspline,xlab='Day',xrange=c(1,365),main='NASA data')
#' # Selecting the significant components for the radius at 5% level
#' which(ITP.result.bspline$adjusted.pval < 0.05)
#' 
#' @references
#' #' A. Pini and S. Vantini (2017).
#' The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. Biometrics 73(3): 835–845.
#' 
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data. \emph{Journal of Nonparametric Statistics}, 29(2), 407-424
#'
#' Pini, A., Vantini, S., Colosimo, B. M., & Grasso, M. (2018). Domain‐selective functional analysis of variance for supervised statistical profile monitoring of signal data. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)} 67(1), 55-81.
#'
#' Abramowicz, K., Hager, C. K., Pini, A., Schelin, L., Sjostedt de Luna, S., & Vantini, S. (2018).
#' Nonparametric inference for functional‐on‐scalar linear models applied to knee kinematic hop data after injury of the anterior cruciate ligament. \emph{Scandinavian Journal of Statistics} 45(4), 1036-1061.
#'
#'
#' @export
#' 

plot.ITP2 <-
function(x,xrange=c(0,1),alpha1=0.05,alpha2=0.01,
                      ylab='Functional Data',main=NULL,lwd=1,col=c(1,2),pch=16,ylim=range(object$data.eval),
                      ...){
  object <- x
  if(alpha1 < alpha2){
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  graphics::par(ask=TRUE) 
  if(object$basis=='Fourier'){
    p <- length(object$pval)
    J <- dim(object$data.eval)[2]
    n <- dim(object$data.eval)[1]
    xmin <- xrange[1]
    xmax <- xrange[2]
    abscissa.pval = 1:p
    Abscissa = seq(xmin,xmax,len=J)
    main.data <- paste(main,': Functional Data')
    main.data <- sub("^ : +", "", main.data)
    colors <- numeric(n)
    colors[which(object$labels==1)] <- col[1]
    colors[which(object$labels==2)] <- col[2]
    
    fda::matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=colors,lwd=lwd,ylim=ylim,...)
    
    ################################################################
    # pval
    main.p <- paste(main,': Adjusted p-values')
    main.p <- sub("^ : +", "", main.p)
    plot(abscissa.pval,object$adjusted.pval,pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',xlab='Frequency')
    difference1 <- which(object$adjusted.pval<alpha1)
    if(length(difference1)>0){
      for(j in 1:length(difference1)){
        min.rect <- abscissa.pval[difference1[j]] - 0.5
        max.rect <- min.rect + 1
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray90',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    difference2 <- which(object$adjusted.pval<alpha2)
    if(length(difference2)>0){
      for(j in 1:length(difference2)){
        min.rect <- abscissa.pval[difference2[j]] - 0.5
        max.rect <- min.rect + 1
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray80',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    for(j in 0:10){
      graphics::abline(h=j/10,col='lightgray',lty="dotted")
    }
    graphics::points(1:p,object$adjusted.pval,pch=pch)
    
    
    
  }else if(object$basis=='B-spline'){  
    p <- length(object$pval)
    J <- dim(object$data.eval)[2]
    n <- dim(object$data.eval)[1]
    xmin <- xrange[1]
    xmax <- xrange[2]
    abscissa.pval = seq(xmin,xmax,len=p)
    Abscissa = seq(xmin,xmax,len=J)
    main.data <- paste(main,': Functional Data')
    main.data <- sub("^ : +", "", main.data)
    colors <- numeric(n)
    colors[which(object$labels==1)] <- col[1]
    colors[which(object$labels==2)] <- col[2]
    
    fda::matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=colors,lwd=lwd,ylim=ylim,...)
    difference1 <- which(object$adjusted.pval<alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray90", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    difference2 <- which(object$adjusted.pval<alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray80", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    fda::matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=colors,lwd=lwd,add=TRUE,...)
    
    
    ################################################################
    # pval
    main.p <- paste(main,': Adjusted p-values')
    main.p <- sub("^ : +", "", main.p)
    plot(abscissa.pval,object$adjusted.pval,pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',...)
    difference1 <- which(object$adjusted.pval<alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray90", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    difference2 <- which(object$adjusted.pval<alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray80", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    for(j in 0:10){
      graphics::abline(h=j/10,col='lightgray',lty="dotted")
    }
    graphics::points(abscissa.pval,object$adjusted.pval,pch=pch)
    
  }else if(object$basis=='paFourier'){
    p <- length(object$pval_phase)
    J <- dim(object$data.eval)[2]
    n <- dim(object$data.eval)[1]
    xmin <- xrange[1]
    xmax <- xrange[2]
    abscissa.pval = 1:p
    Abscissa = seq(xmin,xmax,len=J)
    main.data <- paste(main,': Functional Data')
    main.data <- sub("^ : +", "", main.data)
    colors <- numeric(n)
    colors[which(object$labels==1)] <- col[1]
    colors[which(object$labels==2)] <- col[2]
    
    fda::matplot(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=colors,lwd=lwd,ylim=ylim,...)
    
    ################################################################
    # pval phase
    main.p <- paste(main,': Adjusted p-values - phase')
    main.p <- sub("^ : +", "", main.p)
    plot(abscissa.pval,object$adjusted.pval_phase,pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',xlab='Frequency')
    difference1 <- which(object$adjusted.pval_phase<alpha1)
    if(length(difference1)>0){
      for(j in 1:length(difference1)){
        min.rect <- abscissa.pval[difference1[j]] - 0.5
        max.rect <- min.rect + 1
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray90',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    difference2 <- which(object$adjusted.pval_phase<alpha2)
    if(length(difference2)>0){
      for(j in 1:length(difference2)){
        min.rect <- abscissa.pval[difference2[j]] - 0.5
        max.rect <- min.rect + 1
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray80',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    for(j in 0:10){
      graphics::abline(h=j/10,col='lightgray',lty="dotted")
    }
    graphics::points(1:p,object$adjusted.pval_phase,pch=pch)
    
    ################################################################
    # pval amplitude
    main.p <- paste(main,': Adjusted p-values - amplitude')
    main.p <- sub("^ : +", "", main.p)
    plot(abscissa.pval,object$adjusted.pval_amplitude,pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',xlab='Frequency')
    difference1 <- which(object$adjusted.pval_amplitude<alpha1)
    if(length(difference1)>0){
      for(j in 1:length(difference1)){
        min.rect <- abscissa.pval[difference1[j]] - 0.5
        max.rect <- min.rect + 1
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray90',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    difference2 <- which(object$adjusted.pval_amplitude<alpha2)
    if(length(difference2)>0){
      for(j in 1:length(difference2)){
        min.rect <- abscissa.pval[difference2[j]] - 0.5
        max.rect <- min.rect + 1
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray80',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    for(j in 0:10){
      graphics::abline(h=j/10,col='lightgray',lty="dotted")
    }
    graphics::points(1:p,object$adjusted.pval_amplitude,pch=pch)
    
    
  }
  graphics::par(ask=FALSE)
}
