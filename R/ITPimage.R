#' @title Heatmap plot of the Interval Testing Procedure results
#'
#' @description Plotting function creating a graphical output of the ITP: the p-value heat-map, the plot of the corrected p-values, and the plot of the functional data.
#'
#' @param ITP.result Results of the ITP, as created by \code{\link{ITP1bspline}}, \code{\link{ITP1fourier}}, \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, and \code{\link{ITP2pafourier}}.
#'
#' @param alpha Threshold for the interval-wise error rate used for the hypothesis test. The default is \code{alpha}=0.05.
#'
#' @param abscissa.range Range of the plot abscissa. The default is \code{c(0,1)}.
#'
#' @param nlevel Number of desired color levels for the p-value heatmap. The default is \code{nlevel=20}.
#'
#' @return NULL
#'
#' @seealso See \code{\link{plot.ITP1}}, \code{\link{plot.ITP2}}, \code{\link{plot.ITPlm}}, and \code{\link{plot.ITPaov}} for the plot method applied to the ITP results of one- and two-population tests, linear models, and ANOVA, respectively.
#' See also \code{\link{ITP1bspline}}, \code{\link{ITP1fourier}}, \code{\link{ITP2bspline}}, \code{\link{ITP2fourier}}, and \code{\link{ITP2pafourier}} for applying the ITP.
#'
#'
#' @examples
#' # Performing the ITP for two populations with the B-spline basis
#' ITP.result <- ITP2bspline(NASAtemp$milan,NASAtemp$paris,nknots=20,B=1000)
#'
#' # Plotting the results of the ITP
#' ITPimage(ITP.result,abscissa.range=c(0,12))
#'
#' # Selecting the significant components for the radius at 5% level
#' which(ITP.result$corrected.pval < 0.05)
#'
#' @references A. Pini and S. Vantini (2017).
#' The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. Biometrics 73(3): 835–845.
#'
#' @export



ITPimage <-
function(ITP.result,alpha=0.05,abscissa.range=c(0,1),nlevel=20){

  if(ITP.result$basis=='paFourier' & ITP.result$test=='2pop'){
    #2 plots: phase and amplitude
    graphics::par(ask=T)
    #phase:
    p <- dim(ITP.result$heatmap.matrix_phase)[1]
    min.ascissa <- 1-(p-1)/2
    max.ascissa <- p+(p-1)/2
    ascissa.grafico <- seq(min.ascissa,max.ascissa,length.out=p*4)
    ordinata.grafico <- 1:p
    colori=grDevices::rainbow(nlevel,start=0.15,end=0.67)
    colori <- colori[length(colori):1]
    ##dev.new()
    #graphics::layout(rbind(c(1,1,1,1,1,1,1,1,2),c(1,1,1,1,1,1,1,1,2),c(3,3,3,3,3,3,3,3,0),c(4,4,4,4,4,4,4,4,0)))
    graphics::layout(rbind(1:2,c(3,0),c(4,0)),widths=c(8,1),heights=c(2,1,1))
    graphics::par(mar=c(4.1, 4.1, 3, .2),cex.main=1.5,cex.lab=1.1,las=0)
    #1: heatmap
    matrice.quad <- ITP.result$heatmap.matrix_phase[,(p+1):(3*p)]
    ascissa.quad <- ascissa.grafico[(p+1):(3*p)]
    graphics::image(ascissa.quad,ordinata.grafico,t(matrice.quad[p:1,]),col=colori,ylab='Interval length',main='p-value heatmap (phase)',xlab='Abscissa',zlim=c(0,1),asp=1)
    min.plot <- graphics::par("usr")[1]
    max.plot <- graphics::par("usr")[2]
    #2: legend
    graphics::par(mar=c(4.1, 1, 3, 3),las=1)
    graphics::image(1,seq(0,1,length.out=nlevel)-0.025*seq(0,1,length.out=nlevel)+0.025*seq(1,0,length.out=nlevel),t(as.matrix(seq(0,1,length.out=nlevel))),col=colori,xaxt='n',yaxt='n',xlab='',ylab='')
    graphics::axis(4,at=seq(0,1,0.2),padj=0.4)
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    #3: corrected p values
    graphics::par(mar=c(4.1, 4.1, 3, .2),las=0)
    plot(1:p,ITP.result$corrected.pval_phase,pch=16,ylim=c(0,1),xlim=c(min.plot,max.plot),main='Corrected p-values (phase)',ylab='p-value',xlab='Component',xaxs='i')
    #gray shaded area
    difference <- which(ITP.result$corrected.pval_phase<alpha)
    abscissa.pval <- 1:p
    if(length(difference)>0){
      for(j in 1:length(difference)){
        min.rect <- abscissa.pval[difference[j]] - 0.5
        max.rect <- min.rect + 1
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray90',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    for(j in 0:10){
      graphics::abline(h=j/10,col='lightgray',lty="dotted")
    }
    graphics::points(abscissa.pval,ITP.result$corrected.pval_phase,pch=16)
    #4: functional data
    abscissa.new <- seq(abscissa.range[1],abscissa.range[2],length.out=dim(ITP.result$data.eval)[2])
    fda::matplot(abscissa.new,t(ITP.result$data.eval),col=ITP.result$labels,type='l',main='Functional data',xlab='Abscissa',ylab='Value',xaxs='i')
    #amplitude
    #dev.new()
    #graphics::layout(rbind(c(1,1,1,1,1,1,1,1,2),c(1,1,1,1,1,1,1,1,2),c(3,3,3,3,3,3,3,3,0),c(4,4,4,4,4,4,4,4,0)))
    graphics::layout(rbind(1:2,c(3,0),c(4,0)),widths=c(8,1),heights=c(2,1,1))
    graphics::par(mar=c(4.1, 4.1, 3, .2),cex.main=1.5,cex.lab=1.1,las=0)
    #1: heatmap
    matrice.quad <- ITP.result$heatmap.matrix_amplitude[,(p+1):(3*p)]
    ascissa.quad <- ascissa.grafico[(p+1):(3*p)]
    graphics::image(ascissa.quad,ordinata.grafico,t(matrice.quad[p:1,]),col=colori,ylab='Interval length',main='p-value heatmap (amplitude)',xlab='Abscissa',zlim=c(0,1),asp=1)
    min.plot <- graphics::par("usr")[1]
    max.plot <- graphics::par("usr")[2]
    #2: legend
    graphics::par(mar=c(4.1, 1, 3, 3),las=1)
    graphics::image(1,seq(0,1,length.out=nlevel)-0.025*seq(0,1,length.out=nlevel)+0.025*seq(1,0,length.out=nlevel),t(as.matrix(seq(0,1,length.out=nlevel))),col=colori,xaxt='n',yaxt='n',xlab='',ylab='')
    graphics::axis(4,at=seq(0,1,0.2),padj=0.4)
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    #3: corrected p values
    graphics::par(mar=c(4.1, 4.1, 3, .2),las=0)
    plot(1:p,ITP.result$corrected.pval_amplitude,pch=16,ylim=c(0,1),xlim=c(min.plot,max.plot),main='Corrected p-values (amplitude)',ylab='p-value',xlab='Component',xaxs='i')
    difference <- which(ITP.result$corrected.pval_amplitude<alpha)
    abscissa.pval <- 1:p
    if(length(difference)>0){
      for(j in 1:length(difference)){
        min.rect <- abscissa.pval[difference[j]] - 0.5
        max.rect <- min.rect + 1
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray90',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    for(j in 0:10){
      graphics::abline(h=j/10,col='lightgray',lty="dotted")
    }
    graphics::points(1:p,ITP.result$corrected.pval_amplitude,pch=16)
    #4: functional data
    abscissa.new <- seq(abscissa.range[1],abscissa.range[2],length.out=dim(ITP.result$data.eval)[2])
    fda::matplot(abscissa.new,t(ITP.result$data.eval),col=ITP.result$labels,type='l',main='Functional data',xlab='Abscissa',ylab='Value',xaxs='i')
    graphics::par(ask=FALSE)
  }else if(ITP.result$basis=='Fourier'){
    p <- dim(ITP.result$heatmap.matrix)[1]
    min.ascissa <- 1-(p-1)/2
    max.ascissa <- p+(p-1)/2
    ascissa.grafico <- seq(min.ascissa,max.ascissa,length.out=p*4)
    ordinata.grafico <- 1:p
    colori=grDevices::rainbow(nlevel,start=0.15,end=0.67)
    colori <- colori[length(colori):1]
    #dev.new()
    #graphics::layout(rbind(c(1,1,1,1,1,1,1,1,2),c(1,1,1,1,1,1,1,1,2),c(3,3,3,3,3,3,3,3,0),c(4,4,4,4,4,4,4,4,0)))
    graphics::layout(rbind(1:2,c(3,0),c(4,0)),widths=c(8,1),heights=c(2,1,1))
    graphics::par(mar=c(4.1, 4.1, 3, .2),cex.main=1.5,cex.lab=1.1,las=0)
    #1: heatmap
    matrice.quad <- ITP.result$heatmap.matrix[,(p+1):(3*p)]
    ascissa.quad <- ascissa.grafico[(p+1):(3*p)]
    graphics::image(ascissa.quad,ordinata.grafico,t(matrice.quad[p:1,]),col=colori,ylab='Interval length',main='p-value heatmap',xlab='Abscissa',zlim=c(0,1),asp=1)
    min.plot <- graphics::par("usr")[1]
    max.plot <- graphics::par("usr")[2]
    #2: legend
    graphics::par(mar=c(4.1, 1, 3, 3),las=1)
    graphics::image(1,seq(0,1,length.out=nlevel)-0.025*seq(0,1,length.out=nlevel)+0.025*seq(1,0,length.out=nlevel),t(as.matrix(seq(0,1,length.out=nlevel))),col=colori,xaxt='n',yaxt='n',xlab='',ylab='')
    graphics::axis(4,at=seq(0,1,0.2),padj=0.4)
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    #3: corrected p values
    graphics::par(mar=c(4.1, 4.1, 3, .2),las=0)
    plot(1:p,ITP.result$corrected.pval,pch=16,ylim=c(0,1),xlim=c(min.plot,max.plot),main='Corrected p-values',ylab='p-value',xlab='Component',xaxs='i')
    difference <- which(ITP.result$corrected.pval<alpha)
    abscissa.pval <- 1:p
    if(length(difference)>0){
      for(j in 1:length(difference)){
        min.rect <- abscissa.pval[difference[j]] - 0.5
        max.rect <- min.rect + 1
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray90',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    for(j in 0:10){
      graphics::abline(h=j/10,col='lightgray',lty="dotted")
    }
    graphics::points(1:p,ITP.result$corrected.pval,pch=16)
    #4: functional data
    abscissa.new <- seq(abscissa.range[1],abscissa.range[2],length.out=dim(ITP.result$data.eval)[2])
    fda::matplot(abscissa.new,t(ITP.result$data.eval),col=ITP.result$labels,type='l',main='Functional data',xlab='Abscissa',ylab='Value',xaxs='i')
    if(ITP.result$test=='1pop'){
      if(length(ITP.result$mu)==1){
        abscissa.mu <- abscissa.new
        mu <- rep(ITP.result$mu,1000)
      }else{
        abscissa.mu <- seq(abscissa.range[1],abscissa.range[2],length.out=length(ITP.result$mu))
        mu <- ITP.result$mu
      }
      graphics::lines(abscissa.mu,mu,col='gray')
    }
  }else if(ITP.result$basis=='B-spline'){
    min.ascissa <- abscissa.range[1]-(abscissa.range[2]-abscissa.range[1])/2
    max.ascissa <- abscissa.range[2]+(abscissa.range[2]-abscissa.range[1])/2
    p <- dim(ITP.result$heatmap.matrix)[1]
    ordinata.grafico <- seq(abscissa.range[1],abscissa.range[2],length.out=p) - abscissa.range[1]
    colori=grDevices::rainbow(nlevel,start=0.15,end=0.67)
    colori <- colori[length(colori):1]
    #dev.new()
    #graphics::layout(rbind(c(1,1,1,1,1,1,1,1,2),c(1,1,1,1,1,1,1,1,2),c(3,3,3,3,3,3,3,3,0),c(4,4,4,4,4,4,4,4,0)))
    graphics::layout(rbind(1:2,c(3,0),c(4,0)),widths=c(8,1),heights=c(2,1,1))
    #1: heatmap
    graphics::par(mar=c(4.1, 4.1, 3, .2),cex.main=1.5,cex.lab=1.1,las=0)
    matrice.quad <- ITP.result$heatmap.matrix[,(p+1):(3*p)]
    ascissa.quad <- seq(abscissa.range[1],abscissa.range[2],length.out=p*2)
    graphics::image(ascissa.quad,ordinata.grafico,t(matrice.quad[p:1,]),col=colori,ylab='Interval length',main='p-value heatmap',xlab='Abscissa',zlim=c(0,1),asp=1)
    min.plot <- graphics::par("usr")[1]
    max.plot <- graphics::par("usr")[2]
    #2: legend
    graphics::par(mar=c(4.1, 1, 3, 3),las=1)
    graphics::image(1,seq(0,1,length.out=nlevel)-0.025*seq(0,1,length.out=nlevel)+0.025*seq(1,0,length.out=nlevel),t(as.matrix(seq(0,1,length.out=nlevel))),col=colori,xaxt='n',yaxt='n',xlab='',ylab='')
    graphics::axis(4,at=seq(0,1,0.2),padj=0.4)
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    #3: corrected p values
    graphics::par(mar=c(4.1, 4.1, 3, .2),las=0)
    abscissa.pval <- seq(abscissa.range[1],abscissa.range[2],length.out=p)
    plot(abscissa.pval,ITP.result$corrected.pval,pch=16,ylim=c(0,1),xlim=c(min.plot,max.plot),main='Corrected p-values',ylab='p-value',xlab='Component',xaxs='i')
    difference <- which(ITP.result$corrected.pval<alpha)
    if(length(difference) >0){
      for(j in 1:length(difference)){
        min.rect <- abscissa.pval[difference[j]] - (abscissa.pval[2]-abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2]-abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray90',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    for(j in 0:10){
      graphics::abline(h=j/10,col='lightgray',lty="dotted")
    }
    graphics::points(abscissa.pval,ITP.result$corrected.pval,pch=16)
    #3: functional data
    abscissa.new <- seq(abscissa.range[1],abscissa.range[2],length.out=dim(ITP.result$data.eval)[2])
    fda::matplot(abscissa.new,t(ITP.result$data.eval),col=ITP.result$labels,type='l',xlim=c(min.plot,max.plot),main='Functional data',xlab='Abscissa',ylab='Value',xaxs='i')

    if(length(difference) >0){
      for(j in 1:length(difference)){
        min.rect <- abscissa.pval[difference[j]] - (abscissa.pval[2]-abscissa.pval[1])/2
        max.rect <- min.rect + (abscissa.pval[2]-abscissa.pval[1])
        graphics::rect(min.rect, graphics::par("usr")[3], max.rect, graphics::par("usr")[4], col = 'gray90',density=-2,border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = NULL,border='black')
    }
    fda::matplot(abscissa.new,t(ITP.result$data.eval),col=ITP.result$labels,type='l',add=TRUE)

    if(ITP.result$test=='1pop'){
      if(length(ITP.result$mu)==1){
        abscissa.mu <- abscissa.new
        mu <- rep(ITP.result$mu,1000)
      }else{
        abscissa.mu <- seq(abscissa.range[1],abscissa.range[2],length.out=length(ITP.result$mu))
        mu <- ITP.result$mu
      }
      graphics::lines(abscissa.mu,mu,col='blue')
    }
  }
}
