#' @export
plot.IWT1 <- function(x, xrange = c(0,1),
                      alpha1 = 0.05, alpha2 = 0.01,
                      ylab = 'Functional Data', main = NULL, 
                      lwd = 1, col = 1, 
                      ylim = NULL,type='l', ...) {
  if (class(x) != "IWT1") stop("First argument is not a IWT1 object.")
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