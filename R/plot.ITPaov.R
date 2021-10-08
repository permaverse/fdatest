#' @title Plot method for ITP results on functional ANOVA
#' 
#' @description \code{plot} method for class "\code{ITPaov}".
#' Plotting function creating a graphical output of the ITP for the test on a functional analysis of variance: 
#' functional data, and ITP-adjusted p-values of the F-tests on the whole model and on each factor are plotted.
#' 
#' @param x  The object to be plotted. An object of class "\code{ITPaov}", usually, a result of a call 
#' to \code{\link{ITPaov}}.
#' 
#' @param xrange Range of the \code{x} axis.
#' 
#' @param alpha1 First level of significance used to select and display significant effects. Default is \code{alpha1 = 0.05}.
#' 
#' @param alpha2 Second level of significance used to select and display significant effects. Default is \code{alpha1 = 0.01}. 
#' \code{alpha1} and \code{alpha2} are s.t. \code{alpha2 < alpha1}. Otherwise the two values are switched.
#' 
#' @param plot.adjpval A logical indicating wether the plots of adjusted p-values have to be done. Default is \code{plot_adjpval = FALSE}.
#' 
#' @param ylim Range of the \code{y} axis. Default is \code{NULL}, giving a plot with authomatic range for functional data.
#' 
#' @param col Colors for the plot of functional data. Default is \code{col = 1}.
#' 
#' @param ylab Label of \code{y} axis of the plot of functional data. Default is "\code{Functional Data}".
#' 
#' @param main An overall title for the plots (it will be pasted to "Functional Data and F-test" for the first plot and to factor names for the other plots).
#' 
#' @param lwd Line width for the plot of functional data. Default is \code{lwd=1}.
#' 
#' @param pch Point character for the plot of adjusted p-values. Default is \code{pch=16}.
#' 
#' @param ... Additional plotting arguments that can be used with function \code{plot}, 
#' such as \code{\link{graphical parameters}} (see \code{\link{par}}).
#' 
#' @return No value returned. 
#' The function produces a graphical output of the ITP results:  the plot of the functional data and the one of the adjusted p-values. 
#' The portions of the domain selected as significant by the test at level \code{alpha1} and \code{alpha2} are highlighted in the plot of the adjusted p-value function and in the one of functional data by gray areas (light and dark gray, respectively). 
#' The first plot reports the gray areas corresponding to a significant F-test on the whole model. The remaining plots report the gray areas corresponding to significant F-tests on each factor (with colors corresponding to the levels of the factor).
#' 
#' @seealso \code{\link{ITPimage}} for the plot of p-values heatmaps. 
#' See also \code{\link{ITP1bspline}}, \code{\link{ITP2bspline}} to perform the ITP to test on the mean of one population and test of differences between two populations. 
#' See \code{\link{IWTaov}} for functional ANOVA not based on B-spline basis representation
#' 
#' @examples 
#' temperature <- rbind(NASAtemp$milan,NASAtemp$paris)
#' groups <- c(rep(0,22),rep(1,22))
#' 
#' # Performing the ITP
#' ITP.result <- ITPaovbspline(temperature ~ groups,B=1000,nknots=20,order=3)
#' 
#' # Summary of the ITP results
#' summary(ITP.result)
#' 
#' # Plot of the ITP results
#' graphics::layout(1)
#' plot(ITP.result)
#' 
#' # All graphics on the same device
#' graphics::layout(matrix(1:4,nrow=2,byrow=FALSE))
#' plot(ITP.result,main='NASA data', plot_adjpval = TRUE,xlab='Day',xrange=c(1,365))
#' 
#' @references
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data. \emph{Journal of Nonparametric Statistics}, 29(2), 407-424
#'
#' Pini, A., Vantini, S., Colosimo, B. M., & Grasso, M. (2018). Domain‐selective functional analysis of variance for supervised statistical profile monitoring of signal data. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)} 67(1), 55-81.
#'
#' Abramowicz, K., Hager, C. K., Pini, A., Schelin, L., Sjostedt de Luna, S., & Vantini, S. (2018).
#' Nonparametric inference for functional‐on‐scalar linear models applied to knee kinematic hop data after injury of the anterior cruciate ligament. \emph{Scandinavian Journal of Statistics} 45(4), 1036-1061.
#'
#' @export

plot.ITPaov <-
function(x,xrange=c(0,1),alpha1=0.05,alpha2=0.01,plot.adjpval=FALSE,
                      ylim=range(x$data.eval),col=1,
                      ylab='Functional Data',main=NULL,lwd=1,pch=16,
                      ...){
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
  
  
  fda::matplot(Abscissa,t(object$data.eval),type='l',col=col,main=main.f,ylab=ylab,ylim=ylim,lwd=lwd,...)
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
  fda::matplot(Abscissa,t(object$data.eval),type='l',col=col,add=TRUE,lwd=lwd,...)
  
  formula <- object$call$formula
  mf <- stats::model.frame(formula)
  nvar <- dim(object$adjusted.pval.factors)[1]
  names_all <- colnames(object$design.matrix)
  interaz <- grep(':',names_all)
  
  for(var in 1:(dim(object$adjusted.pval.factors)[1])){
    var.name = rownames(object$adjusted.pval.factors)[var]
    main.t <- paste(main,': factor',var.name,sep=' ')
    main.t <- sub("^ : +", "", main.t)
    
    if(length(grep(':',var.name))>0){ # sto plottando interazione
      var12 <- strsplit(var.name,':')
      var1 <- var12[[1]][1]
      var2 <- var12[[1]][2]
      dummy.test1 <- grep(var1,names_all)
      dummy.test2 <- grep(var2,names_all)
      dummy.test <- intersect(dummy.test1,dummy.test2)
      colors <- object$design.matrix[,dummy.test]
      if(length(dim(colors))>1){
        colors <- (apply(colors,1,paste,collapse=''))
      }
      colors <- as.factor(colors)
    }else{ #sto plottando un fattore
      dummy.test <- grep(var.name,names_all)
      dummy.test <- setdiff(dummy.test,interaz)
      colors <- object$design.matrix[,dummy.test]
      if(length(dim(colors))>1){
        colors <- (apply(colors,1,paste,collapse=''))
      }
      colors <- as.factor(colors)
    }
    
    
    fda::matplot(Abscissa,t(object$data.eval),type='l',col=colors,ylim=ylim,lwd=1,main=main.t,ylab=ylab,...)
    difference1 <- which(object$adjusted.pval.factors[var,] < alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray90", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    difference2 <- which(object$adjusted.pval.factors[var,] < alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray80", density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
    }
    graphics::matlines(Abscissa,t(object$data.eval),type='l',col=colors,...)
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
    
    for(var in 1:(dim(object$adjusted.pval.factors)[1])){
      var.name = rownames(object$adjusted.pval.factors)[var]
      main.p <- paste(main,': Adjusted p-values - factor',var.name)
      main.p <- sub("^ : +", "", main.p)
      plot(Abscissa,object$adjusted.pval.factors[var,],pch=pch,ylim=c(0,1),main=main.p,ylab='p-value',...)
      difference1 <- which(object$adjusted.pval.factors[var,]<alpha1)
      if (length(difference1) > 0) {
        for (j in 1:length(difference1)) {
          min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
          max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
          graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = "gray90", density = -2, border = NA)
        }
        graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],graphics::par("usr")[4], col = NULL, border = "black")
      }
      difference2 <- which(object$adjusted.pval.factors[var,]<alpha2)
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
      graphics::points(Abscissa,object$adjusted.pval.factors[var,],pch=pch)
    }
  }
  graphics::par(ask=F) 
}
