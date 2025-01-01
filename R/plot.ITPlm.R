#' Plot the results of the ITPlm function
#' 
#' @param x Results of the ITPlm function
#' @param xrange Range of the plot abscissa. The default is \code{c(0,1)}.
#' @param alpha1 Threshold for the interval-wise error rate used for the
#'   hypothesis test. The default is \code{alpha1}=0.05.
#' @param alpha2 Threshold for the interval-wise error rate used for the
#'   hypothesis test. The default is \code{alpha2}=0.01.
#' @param plot.adjpval Logical. If \code{TRUE}, the plot of the corrected
#'   p-values is displayed. The default is \code{plot.adjpval=FALSE}.
#' @param col Color of the lines in the plot. The default is
#'   \code{col=rainbow(dim(x$adjusted.pval.t)[1])}.
#' @param ylim Range of the plot ordinate. The default is `range(x$data.eval)`.
#' @param ylab Label of the plot ordinate. The default is `'Functional Data'`.
#' @param main Title of the plot. The default is \code{main=NULL}.
#' @param lwd Width of the lines in the plot. The default is \code{lwd=1}.
#' @param pch Symbol of the points in the plot. The default is \code{pch=16}.
#' @param ... Additional arguments to be passed to the plot function.
#'  
#' @return NULL
#'  
#' @export
plot.ITPlm <- function(x, 
                       xrange = c(0, 1), 
                       alpha1 = 0.05, 
                       alpha2 = 0.01, 
                       plot.adjpval = FALSE, 
                       col = c(1, grDevices::rainbow(dim(x$adjusted.pval.t)[1])), 
                       ylim = range(x$data.eval), 
                       ylab = 'Functional Data', 
                       main = NULL, 
                       lwd = 1, 
                       pch = 16, 
                       ...) {
  if(alpha1 < alpha2){
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
    
  object <- x
  p <- length(object$pval.F)
  J <- dim(object$data.eval)[2]
  n <- dim(object$data.eval)[1]
  xmin <- xrange[1]
  xmax <- xrange[2]
  abscissa.pval = seq(xmin,xmax,len=p)
  Abscissa = seq(xmin,xmax,len=J)
  graphics::par(ask=T) 
  main.f <- paste(main,': Functional Data and F-test')
  main.f <- sub("^ : +", "", main.f)
  if(length(col)< dim(object$adjusted.pval.t)[1]+1){
    col <- rep(col,(dim(object$adjusted.pval.t)[1]+1)%/%length(col)+1)
  }
  
  fda::matplot(Abscissa,t(object$data.eval),type='l',col=col[1],main=main.f,ylab=ylab,ylim=ylim,lwd=lwd,...)
  difference1 <- which(object$adjusted.pval.F < alpha1)
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
      max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
      graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray90", density = -2, border = NA)
    }
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
  }
  
  difference2 <- which(object$adjusted.pval.F < alpha2)
  if (length(difference2) > 0) {
    for (j in 1:length(difference2)) {
      min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
      max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
      graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray80", density = -2, border = NA)
    }
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
  }
  fda::matplot(Abscissa,t(object$data.eval),type='l',col=col[1],add=TRUE,lwd=lwd,...)
  
  for(var in 1:(dim(object$adjusted.pval.t)[1])){
    var.name = rownames(object$adjusted.pval.t)[var]
    main.t <- paste(main,': t-test -',var.name,sep=' ')
    main.t <- sub("^ : +", "", main.t)
    plot(Abscissa,object$coeff.regr.eval[var,],type='l',col=col[var+1],ylim=range(c(0,object$coeff.regr.eval[var,])),lwd=2,main=main.t,ylab='Regression Coefficient',...)
    difference1 <- which(object$adjusted.pval.t[var,] < alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray90", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    difference2 <- which(object$adjusted.pval.t[var,] < alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray80", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    graphics::lines(Abscissa,object$coeff.regr.eval[var,],type='l',col=col[var+1],lwd=2,...)
    #graphics::lines(ascissa,coeff.teo[1,],lty=2,add=TRUE,type='l',col=1,lwd=2)
    graphics::abline(h=0,lty=2,col=1)
  }
  #########################################################
  #plot of adjusted p-values
  if(plot.adjpval==TRUE){
    main.p <- paste(main,': Adjusted p-values - F-test')
    main.p <- sub("^ : +", "", main.p)
    Abscissa <- abscissa.pval
    plot(Abscissa,object$adjusted.pval.F,pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',...)
    difference1 <- which(object$adjusted.pval.F<alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray90", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    difference2 <- which(object$adjusted.pval.F<alpha2)
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
    graphics::points(Abscissa,object$adjusted.pval.F,pch=pch)
    
    for(var in 1:(dim(object$adjusted.pval.t)[1])){
      var.name = rownames(object$adjusted.pval.t)[var]
      main.p <- paste(main,': Adjusted p-values - t-test -',var.name)
      main.p <- sub("^ : +", "", main.p)
      plot(Abscissa,object$adjusted.pval.t[var,],pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',...)
      difference1 <- which(object$adjusted.pval.t[var,]<alpha1)
      if (length(difference1) > 0) {
        for (j in 1:length(difference1)) {
          min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
          max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
          graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray90", density = -2, border = NA)
        }
        graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
      }
      difference2 <- which(object$adjusted.pval.t[var,]<alpha2)
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
      graphics::points(Abscissa,object$adjusted.pval.t[var,],pch=pch)
    }
  }
  graphics::par(ask=F) 
}
