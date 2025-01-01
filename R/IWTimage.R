#' Heatmap plot of the Interval Wise Testing Procedure results
#'
#' Plotting function creating a graphical output of the IWT: the p-value
#' heat-map, the plot of the corrected p-values, and the plot of the functional
#' data.
#' 
#' @param IWT_result Results of the ITP, as created by \code{\link{IWT1}},
#'   \code{\link{IWT2}}, \code{\link{IWTaov}}, and \code{\link{IWTlm}}.
#' @param alpha Threshold for the interval-wise error rate used for the
#'   hypothesis test. The default is \code{alpha=0.05}.
#' @param abscissa_range Range of the plot abscissa. The default is \code{c(0,1)}.
#' @param nlevel Number of desired color levels for the p-value heatmap. The
#'   default is \code{nlevel=20}.
#' @param plot_unadjusted Flag indicating if the unadjusted p-value function has
#'   to be added to the plots. The default is \code{FALSE}.
#' 
#' @return No value returned.
#' 
#' @seealso See \code{\link{plot.IWT1}}, \code{\link{plot.fdatest2}},
#'   \code{\link{plot.IWTlm}}, and \code{\link{plot.IWTaov}} for the plot method
#'   applied to the IWT results of one- and two-population tests, linear models,
#'   and ANOVA, respectively.
#' 
#' @references
#' Pini, A., & Vantini, S. (2018). Interval-wise testing for functional data.
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
#' # Performing the IWT for two populations 
#' IWT.result <- IWT2(NASAtemp$milan, NASAtemp$paris, B = 10L)
#' 
#' # Plotting the results of the IWT
#' IWTimage(IWT.result, abscissa_range = c(0, 12))
#' 
#' # Selecting the significant components for the radius at 5\% level
#' which(IWT.result$corrected.pval < 0.05)
IWTimage <- function(IWT_result, 
                     alpha = 0.05,  
                     abscissa_range = c(0, 1), 
                     nlevel = 20, 
                     plot_unadjusted = FALSE) {
  min_ascissa <- abscissa_range[1] - (abscissa_range[2] - abscissa_range[1]) / 2
  max_ascissa <- abscissa_range[2] + (abscissa_range[2] - abscissa_range[1]) / 2
  
  if(inherits(IWT_result, "IWT1") | inherits(IWT_result, "IWT2")){
    p <- length(IWT_result$adjusted_pval)
    
    # heatmap matrix
    heatmap_matrix <- matrix(nrow=p,ncol=4*p)
    for(i in 0:(p-1)){
      for(j in 1:(2*p)){
        heatmap_matrix[p-i,j+i+p] <- IWT_result$pval_matrix[p-i,(j+1)%/%2]
        if(j+i>2*p-i){
          heatmap_matrix[p-i,j+i-p] <- IWT_result$pval_matrix[p-i,(j+1)%/%2]
        }
      }
    }
    
    ordinata_grafico <- seq(abscissa_range[1], abscissa_range[2], length.out = p) - abscissa_range[1]
    colori <- grDevices::rainbow(nlevel, start = 0.15, end = 0.67)
    colori <- colori[length(colori):1]
    graphics::layout(rbind(1:2, c(3, 0), c(4, 0)), widths = c(8, 1), heights = c(2, 1, 1))
    # 1: heatmap
    graphics::par(mar = c(4.1, 4.1, 3, .2), cex.main = 1.5, cex.lab = 1.1, las = 0)
    matrice_quad <- heatmap_matrix[,(p + 1):(3 * p)]
    ascissa_quad <- seq(abscissa_range[1], abscissa_range[2], length.out = p * 2)
    graphics::image(ascissa_quad, ordinata_grafico,t(matrice_quad[p:1, ]), 
          col = colori, ylab = 'Interval length', main = 'p-value heatmap', 
          xlab = 'Abscissa', zlim = c(0, 1), asp = 1)
    min_plot <- graphics::par("usr")[1]
    max_plot <- graphics::par("usr")[2]
    # 2: legend
    # graphics::par(mar = c(4.1, 1, 3, 3), las = 1) # margins too large!
    graphics::par(mar = c(4.1, 1, 3, 2.4), las = 1)
    graphics::image(1, seq(0, 1, length.out = nlevel) - 0.025 * seq(0, 1, length.out=nlevel) 
          + 0.025 * seq(1, 0, length.out = nlevel), 
          t(as.matrix(seq(0, 1, length.out = nlevel))), col = colori,
          xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    graphics::axis(4, at = seq(0, 1, 0.2), padj = 0.4)
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], 
         col = NULL, border = 'black')
    # 3: adjusted p-values
    # graphics::par(mar = c(4.1, 4.1, 3, .2), las = 0) # margins too large!
    graphics::par(mar = c(2.1, 4.1, 3, .2), las = 0)
    abscissa_pval <- seq(abscissa_range[1], abscissa_range[2], length.out = p)
    plot(abscissa_pval, IWT_result$adjusted_pval, ylim = c(0, 1), 
         xlim = c(min_plot, max_plot), main='Adjusted p-value function', 
         ylab = 'p-value', xlab = 'Abscissa', xaxs = 'i',type='l',lwd=2,lty=1)
    if(plot_unadjusted==TRUE){
      graphics::lines(abscissa_pval, IWT_result$unadjusted_pval,lwd=2,lty=2)
    }
    difference <- which(IWT_result$adjusted_pval < alpha)
    if (length(difference) > 0) {
      for (j in 1:length(difference)) {
        min_rect <- abscissa_pval[difference[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], 
             col = 'gray90', density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], 
           col = NULL, border = 'black')      
    }
    for (j in 0:10) {
      graphics::abline(h = j / 10, col = 'lightgray', lty = "dotted")
    }
    graphics::lines(abscissa_pval, IWT_result$adjusted_pval, pch = 16,lwd=2,lty=1)
    if(plot_unadjusted==TRUE){
      graphics::lines(abscissa_pval, IWT_result$unadjusted_pval,lwd=2,lty=2)
    }
    # 4: functional data
    if (IWT_result$test == '2pop') {
      col_pop <- IWT_result$ord_labels
    } else { # 1pop
      col_pop <- 1 # black
    }
    fda::matplot(abscissa_pval, t(IWT_result$data.eval), col = col_pop, 
            type = 'l', 
            main = 'Functional data', xlab = 'Abscissa', 
            ylab = 'Value', xaxs = 'i',xlim = c(min_plot, max_plot))
    if (length(difference) > 0) {
      for (j in 1:length(difference)) {
        min_rect <- abscissa_pval[difference[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], col = 'gray90',
             density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], 
           col = NULL, border = 'black')      
    }
    fda::matplot(abscissa_pval, t(IWT_result$data.eval), col = col_pop, type = 'l', add = TRUE)
    if (IWT_result$test == '1pop'){
      if (length(IWT_result$mu) == 1){
        mu_eval <- rep(IWT_result$mu, p)
      }
      graphics::lines(abscissa_pval, mu_eval, col = 'blue')
    }
  }else if(inherits(IWT_result, "IWTaov")){
    p <- length(IWT_result$adjusted_pval_F)
    
    # heatmap matrix
    nvar <- dim(IWT_result$adjusted_pval_factor)[1]
    heatmap_matrix_F <- matrix(nrow=p,ncol=4*p)
    heatmap_matrix_factor <- array(dim = c(nvar, p, 4 * p))
    for (i in 0:(p-1)){
      for (j in 1:(2*p)){
        heatmap_matrix_F[p-i, j+i+p] <- IWT_result$pval_matrix_F[p-i,(j+1) %/% 2]
        if (j+i > 2*p-i){
          heatmap_matrix_F[p-i, j+i-p] <- IWT_result$pval_matrix_F[p - i,(j + 1) %/% 2]
        }
        for (ii in 1:nvar){
          heatmap_matrix_factor[ii, p - i, j + i + p] <- IWT_result$pval_matrix_factor[ii,p - i,(j + 1) %/% 2]
          if(j + i > 2 * p - i){
            heatmap_matrix_factor[ii, p - i, j + i - p] <- IWT_result$pval_matrix_factor[ii,p - i,(j + 1) %/% 2]
          }
        }
      }
    }
    
    
    ordinata_grafico <- seq(abscissa_range[1], abscissa_range[2], length.out = p) - abscissa_range[1]
    colori <- grDevices::rainbow(nlevel, start = 0.15, end = 0.67)
    colori <- colori[length(colori):1]
    matrix_plot <- rbind(1:2, c(3, 0), c(4, 0))
    for(ii in 1:nvar){
      start <- max(matrix_plot) 
      matrix_temp <- rbind((1:2)+start, c(3+start, 0), c(4+start, 0))
      matrix_plot <- cbind(matrix_plot,matrix_temp)
    }
    graphics::layout(matrix_plot, widths = rep(c(8/(nvar+1), 1),nvar+1), heights = c(2, 1, 1))
    # 1: heatmap
    graphics::par(mar = c(4.1, 4.1, 3, .2), cex.main = 1.5, cex.lab = 1.1, las = 0)
    
    # F test
    matrice_quad_F <- heatmap_matrix_F[,(p + 1):(3 * p)]
    ascissa_quad <- seq(abscissa_range[1], abscissa_range[2], length.out = p * 2)
    graphics::image(ascissa_quad, ordinata_grafico,t(matrice_quad_F[p:1, ]), 
          col = colori, ylab = 'Interval length', main = 'p-value heatmap, F test', 
          xlab = 'Abscissa', zlim = c(0, 1), asp = 1)
    min_plot <- graphics::par("usr")[1]
    max_plot <- graphics::par("usr")[2]
    # 2: legend
    # graphics::par(mar = c(4.1, 1, 3, 3), las = 1) # margins too large!
    graphics::par(mar = c(4.1, 1, 3, 2.4), las = 1)
    graphics::image(1, seq(0, 1, length.out = nlevel) - 0.025 * seq(0, 1, length.out=nlevel) 
          + 0.025 * seq(1, 0, length.out = nlevel), 
          t(as.matrix(seq(0, 1, length.out = nlevel))), col = colori,
          xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    graphics::axis(4, at = seq(0, 1, 0.2), padj = 0.4)
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], 
         col = NULL, border = 'black')
    # 3: adjusted p-values
    # graphics::par(mar = c(4.1, 4.1, 3, .2), las = 0) # margins too large!
    graphics::par(mar = c(2.1, 4.1, 3, .2), las = 0)
    abscissa_pval <- seq(abscissa_range[1], abscissa_range[2], length.out = p)
    plot(abscissa_pval, IWT_result$adjusted_pval_F, ylim = c(0, 1), 
         xlim = c(min_plot, max_plot), main='Adjusted p-value function', 
         ylab = 'p-value', xlab = 'Abscissa', xaxs = 'i',type='l',lwd=2,lty=1)
    if(plot_unadjusted==TRUE){
      graphics::lines(abscissa_pval, IWT_result$unadjusted_pval_F,lwd=2,lty=2)
    }
    difference <- which(IWT_result$adjusted_pval_F < alpha)
    if (length(difference) > 0) {
      for (j in 1:length(difference)) {
        min_rect <- abscissa_pval[difference[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], 
             col = 'gray90', density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], 
           col = NULL, border = 'black')      
    }
    for (j in 0:10) {
      graphics::abline(h = j / 10, col = 'lightgray', lty = "dotted")
    }
    graphics::lines(abscissa_pval, IWT_result$adjusted_pval_F, pch = 16,lwd=2,lty=1)
    if(plot_unadjusted==TRUE){
      graphics::lines(abscissa_pval, IWT_result$unadjusted_pval_F,lwd=2,lty=2)
    }
    # 4: functional data
    col_pop <- factor(apply(IWT_result$design_matrix, 1, paste, collapse = ''))
    fda::matplot(abscissa_pval, t(IWT_result$data.eval), col = col_pop, 
            type = 'l', 
            main = 'Functional data', xlab = 'Abscissa', 
            ylab = 'Value', xaxs = 'i',xlim = c(min_plot, max_plot))
    if (length(difference) > 0) {
      for (j in 1:length(difference)) {
        min_rect <- abscissa_pval[difference[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], col = 'gray90',
             density = -2, border = NA)
      }
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], 
           col = NULL, border = 'black')      
    }
    fda::matplot(abscissa_pval, t(IWT_result$data.eval), col = col_pop, type = 'l', add = TRUE)
    
    factor.names <- rownames(IWT_result$unadjusted_pval_factor)
    all.names <- colnames(IWT_result$design_matrix)
    
    interaz <- grep(':', all.names)
    for(var in 1:nvar){
      graphics::par(mar = c(4.1, 4.1, 3, .2), cex.main = 1.5, cex.lab = 1.1, las = 0)
      
      # factors
      matrice_quad_factor <- heatmap_matrix_factor[var,,(p + 1):(3 * p)]
      graphics::image(ascissa_quad, ordinata_grafico,t(matrice_quad_factor[p:1, ]), 
            col = colori, ylab = 'Interval length', main = paste('p-value heatmap, factor:',factor.names[var]), 
            xlab = 'Abscissa', zlim = c(0, 1), asp = 1)
      min_plot <- graphics::par("usr")[1]
      max_plot <- graphics::par("usr")[2]
      # 2: legend
      # graphics::par(mar = c(4.1, 1, 3, 3), las = 1) # margins too large!
      graphics::par(mar = c(4.1, 1, 3, 2.4), las = 1)
      graphics::image(1, seq(0, 1, length.out = nlevel) - 0.025 * seq(0, 1, length.out=nlevel) 
            + 0.025 * seq(1, 0, length.out = nlevel), 
            t(as.matrix(seq(0, 1, length.out = nlevel))), col = colori,
            xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
      graphics::axis(4, at = seq(0, 1, 0.2), padj = 0.4)
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], 
           col = NULL, border = 'black')
      # 3: adjusted p-values
      # graphics::par(mar = c(4.1, 4.1, 3, .2), las = 0) # margins too large!
      graphics::par(mar = c(2.1, 4.1, 3, .2), las = 0)
      abscissa_pval <- seq(abscissa_range[1], abscissa_range[2], length.out = p)
      plot(abscissa_pval, IWT_result$adjusted_pval_factor[var,], ylim = c(0, 1), 
           xlim = c(min_plot, max_plot), main='Adjusted p-value function', 
           ylab = 'p-value', xlab = 'Abscissa', xaxs = 'i',type='l',lwd=2,lty=1)
      if(plot_unadjusted==TRUE){
        graphics::lines(abscissa_pval, IWT_result$unadjusted_pval_factor[var,],lwd=2,lty=2)
      }
      difference <- which(IWT_result$adjusted_pval_factor[var,] < alpha)
      if (length(difference) > 0) {
        for (j in 1:length(difference)) {
          min_rect <- abscissa_pval[difference[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
          max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
          graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], 
               col = 'gray90', density = -2, border = NA)
        }
        graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], 
             col = NULL, border = 'black')      
      }
      for (j in 0:10) {
        graphics::abline(h = j / 10, col = 'lightgray', lty = "dotted")
      }
      graphics::lines(abscissa_pval, IWT_result$adjusted_pval_factor[var,], pch = 16,lwd=2,lty=1)
      if(plot_unadjusted==TRUE){
        graphics::lines(abscissa_pval, IWT_result$unadjusted_pval_factor[var,],lwd=2,lty=2)
      }
      # 4: functional data
      var_name <- factor.names[var]
      if (length(grep(':', var_name)) > 0) { # Plot interaction
        var12 <- strsplit(var_name, ':')
        var1 <- var12[[1]][1]
        var2 <- var12[[1]][2]
        dummy_test1 <- grep(var1, all.names)
        dummy_test2 <- grep(var2, all.names)
        dummy_test <- intersect(dummy_test1, dummy_test2)
        colors <- IWT_result$design_matrix[, dummy_test]
        if (length(dim(colors)) > 1) {
          colors <- (apply(colors, 1, paste, collapse = ''))
        }
        colors <- as.factor(colors)
      } else { # Plot of a factor
        dummy_test <- grep(var_name, all.names)
        dummy_test <- setdiff(dummy_test, interaz)
        colors <- IWT_result$design_matrix[, dummy_test]
        if (length(dim(colors)) > 1) {
          colors <- (apply(colors, 1, paste, collapse = ''))
        }
        colors <- as.factor(colors)
      }
      
      fda::matplot(abscissa_pval, t(IWT_result$data.eval), col = colors, 
              type = 'l', 
              main = 'Functional data', xlab = 'Abscissa', 
              ylab = 'Value', xaxs = 'i',xlim = c(min_plot, max_plot))
      if (length(difference) > 0) {
        for (j in 1:length(difference)) {
          min_rect <- abscissa_pval[difference[j]] - (abscissa_pval[2] - abscissa_pval[1]) / 2
          max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
          graphics::rect(min_rect, graphics::par("usr")[3], max_rect, graphics::par("usr")[4], col = 'gray90',
               density = -2, border = NA)
        }
        graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], 
             col = NULL, border = 'black')      
      }
      fda::matplot(abscissa_pval, t(IWT_result$data.eval), col = colors,type = 'l',add=TRUE)
    }
  }
}
