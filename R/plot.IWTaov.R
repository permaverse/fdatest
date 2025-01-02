#' Plot method for IWT results on functional ANOVA
#' 
#' \code{plot} method for class "\code{IWTaov}". Plotting function creating a
#' graphical output of the IWT for the test on a functional analysis of
#' variance: functional data, and IWT-adjusted p-values of the F-tests on the
#' whole model and on each factor are plotted.
#' 
#' @param x  The object to be plotted. An object of class "\code{IWTaov}",
#'   usually, a result of a call to \code{\link{IWTaov}}.
#' @param xrange Range of the \code{x} axis.
#' @param alpha1 First level of significance used to select and display
#'   significant effects. Default is \code{alpha1 = 0.05}.
#' @param alpha2 Second level of significance used to select and display
#'   significant effects. Default is \code{alpha1 = 0.01}. \code{alpha1} and
#'   \code{alpha2} are s.t. \code{alpha2 < alpha1}. Otherwise the two values are
#'   switched.
#' @param plot_adjpval A logical indicating wether the plots of adjusted
#'   p-values have to be done. Default is \code{plot_adjpval = FALSE}.
#' @param ylim Range of the \code{y} axis. Default is \code{NULL}, giving a plot
#'   with authomatic range for functional data.
#' @param col Colors for the plot of functional data. Default is \code{col = 1}.
#' @param ylab Label of \code{y} axis of the plot of functional data. Default is
#'   "\code{Functional Data}".
#' @param main An overall title for the plots (it will be pasted to "Functional
#'   Data and F-test" for the first plot and to factor names for the other
#'   plots).
#' @param lwd Line width for the plot of the adjusted p-value function. Default
#'   is \code{lwd=1}.
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
#'   (light and dark gray, respectively). The first plot reports the gray areas
#'   corresponding to a significant F-test on the whole model. The remaining
#'   plots report the gray areas corresponding to significant F-tests on each
#'   factor (with colors corresponding to the levels of the factor).
#' 
#' @seealso \code{\link{IWTimage}} for the plot of p-values heatmaps. See also
#'   \code{\link{IWT1}}, \code{\link{IWT2}} to perform the ITP to test on the
#'   mean of one population and test of differences between two populations. See
#'   \code{\link{ITPaovbspline}} for functional ANOVA based on B-spline basis
#'   representation
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
#' graphics::layout(1)
#' plot(IWT.result)
#' 
#' # All graphics on the same device
#' graphics::layout(matrix(1:4, nrow = 2, byrow = FALSE))
#' plot(
#'   IWT.result, 
#'   main = 'NASA data', 
#'   plot_adjpval = TRUE, 
#'   xlab = 'Day', 
#'   xrange = c(1, 365)
#' )
plot.IWTaov <- function(x, 
                        xrange = c(0, 1), 
                        alpha1 = 0.05, 
                        alpha2 = 0.01, 
                        plot_adjpval = FALSE,
                        ylim = NULL, 
                        col = 1,
                        ylab = "Functional Data", 
                        main = NULL, 
                        lwd = 0.5, 
                        type = "l", 
                        ...) {
  if (!inherits(x, "IWTaov"))
    stop("x should be an object of the class IWTaov")
  
  if (alpha1 < alpha2) {
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  
  object <- x
  nvar <- dim(object$adjusted_pval_factors)[1]
  p <- length(object$unadjusted_pval_F)
  n <- dim(t(object$data.eval))[1]
  xmin <- xrange[1]
  xmax <- xrange[2]
  abscissa_pval <- seq(xmin, xmax, len = p)
  grDevices::devAskNewPage(ask = TRUE)
  main_f <- paste(main, ": Functional Data and F-test")
  main_f <- sub("^ : +", "", main_f)
  
  if (is.null(ylim))
    ylim <- range(object$data.eval)
  
  if (nvar > 1) {
    fda::matplot(
      abscissa_pval,
      t(object$data.eval),
      type = 'l',
      col = col,
      main = main_f,
      ylab = ylab,
      ylim = ylim,
      lwd = lwd,
      ...
    )
    
    difference1 <- which(object$adjusted_pval_F < alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min_rect <- abscissa_pval[difference1[j]] - 
          (abscissa_pval[2] - abscissa_pval[1]) / 2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        graphics::rect(
          min_rect,
          graphics::par("usr")[3],
          max_rect,
          graphics::par("usr")[4],
          col = "gray90",
          density = -2,
          border = NA
        )
      }
      graphics::rect(
        graphics::par("usr")[1],
        graphics::par("usr")[3],
        graphics::par("usr")[2],
        graphics::par("usr")[4],
        col = NULL,
        border = "black"
      )
    }
    
    difference2 <- which(object$adjusted_pval_F < alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min_rect <- abscissa_pval[difference2[j]] - 
          (abscissa_pval[2] - abscissa_pval[1]) / 2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        graphics::rect(
          min_rect,
          graphics::par("usr")[3],
          max_rect,
          graphics::par("usr")[4],
          col = "gray80",
          density = -2,
          border = NA
        )
      }
      graphics::rect(
        graphics::par("usr")[1],
        graphics::par("usr")[3],
        graphics::par("usr")[2],
        graphics::par("usr")[4],
        col = NULL,
        border = "black"
      )
    }
    
    fda::matplot(
      abscissa_pval,
      t(object$data.eval),
      type = 'l',
      col = col,
      add = TRUE,
      lwd = lwd,
      ...
    )
  }
  
  names_all <- colnames(object$design_matrix)
  interaz <- grep(':', names_all)
  for (var in 1:(dim(object$adjusted_pval_factors)[1])) {
    var_name <- rownames(object$adjusted_pval_factors)[var]
    main_t <- paste(main, ': factor', var_name, sep = ' ')
    main_t <- sub("^ : +", "", main_t)
    if (length(grep(':', var_name)) > 0) { # Plot interaction
      var12 <- strsplit(var_name, ':')
      var1 <- var12[[1]][1]
      var2 <- var12[[1]][2]
      dummy_test1 <- grep(var1, names_all)
      dummy_test2 <- grep(var2, names_all)
      dummy_test <- intersect(dummy_test1, dummy_test2)
      colors <- object$design_matrix[, dummy_test]
      if (length(dim(colors)) > 1) {
        colors <- apply(colors, 1, paste, collapse = '')
      }
      colors <- as.factor(colors)
    } else { # Plot of a factor
      dummy_test <- grep(var_name, names_all)
      dummy_test <- setdiff(dummy_test, interaz)
      colors <- object$design_matrix[, dummy_test]
      if (length(dim(colors)) > 1) {
        colors <- apply(colors, 1, paste, collapse = '')
      }
      colors <- as.factor(colors)
    }
    
    fda::matplot(
      abscissa_pval,
      t(object$data.eval),
      type = 'l',
      col = colors,
      ylim = ylim,
      lwd = 1,
      main = main_t,
      ylab = ylab,
      ...
    )
    
    difference1 <- which(object$adjusted_pval_factors[var, ] < alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min_rect <- abscissa_pval[difference1[j]] - 
          (abscissa_pval[2] - abscissa_pval[1]) / 2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        graphics::rect(
          min_rect,
          graphics::par("usr")[3],
          max_rect,
          graphics::par("usr")[4],
          col = "gray90",
          density = -2,
          border = NA
        )
      }
      graphics::rect(
        graphics::par("usr")[1],
        graphics::par("usr")[3],
        graphics::par("usr")[2],
        graphics::par("usr")[4],
        col = NULL,
        border = "black"
      )
    }
    
    difference2 <- which(object$adjusted_pval_factors[var, ] < alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min_rect <- abscissa_pval[difference2[j]] - 
          (abscissa_pval[2] - abscissa_pval[1]) / 2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        graphics::rect(
          min_rect,
          graphics::par("usr")[3],
          max_rect,
          graphics::par("usr")[4],
          col = "gray80",
          density = -2,
          border = NA
        )
      }
      graphics::rect(
        graphics::par("usr")[1],
        graphics::par("usr")[3],
        graphics::par("usr")[2],
        graphics::par("usr")[4],
        col = NULL,
        border = "black"
      )
    }
    graphics::matlines(
      abscissa_pval,
      t(object$data.eval),
      type = 'l',
      col = colors,
      ...
    )
    graphics::abline(h = 0, lty = 2, col = 1)
  }
  
  # Plot adjusted p-values
  if (plot_adjpval) {
    main_p <- paste(main, ': Adjusted p-values - F-test')
    main_p <- sub("^ : +", "", main_p)
    plot(
      abscissa_pval,
      object$adjusted_pval_F,
      ylim = c(0, 1),
      main = main_p,
      ylab = 'p-value',
      type = type,
      lwd = 2,
      ...
    )
    
    difference1 <- which(object$adjusted_pval_F < alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min_rect <- abscissa_pval[difference1[j]] - 
          (abscissa_pval[2] - abscissa_pval[1]) / 2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        graphics::rect(
          min_rect,
          graphics::par("usr")[3],
          max_rect,
          graphics::par("usr")[4],
          col = "gray90",
          density = -2,
          border = NA
        )
      }
      graphics::rect(
        graphics::par("usr")[1],
        graphics::par("usr")[3],
        graphics::par("usr")[2],
        graphics::par("usr")[4],
        col = NULL,
        border = "black"
      )
    }
    
    difference2 <- which(object$adjusted_pval_F < alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min_rect <- abscissa_pval[difference2[j]] - 
          (abscissa_pval[2] - abscissa_pval[1]) / 2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        graphics::rect(
          min_rect,
          graphics::par("usr")[3],
          max_rect,
          graphics::par("usr")[4],
          col = "gray80",
          density = -2,
          border = NA
        )
      }
      graphics::rect(
        graphics::par("usr")[1],
        graphics::par("usr")[3],
        graphics::par("usr")[2],
        graphics::par("usr")[4],
        col = NULL,
        border = "black"
      )
    }
    
    for (j in 0:10) {
      graphics::abline(h = j / 10, col = 'lightgray', lty = "dotted")
    }
    
    graphics::lines(
      abscissa_pval,
      object$adjusted_pval_F,
      lwd = 2,
      type = type,
      ...
    )
    
    for (var in 1:(dim(object$adjusted_pval_factors)[1])) {
      var_name <- rownames(object$adjusted_pval_factors)[var]
      main_p <- paste(main, ': Adjusted p-values - factor', var_name)
      main_p <- sub("^ : +", "", main_p)
      plot(
        abscissa_pval,
        object$adjusted_pval_factors[var, ],
        ylim = c(0, 1),
        main = main_p,
        ylab = 'p-value',
        lwd = 2,
        type = type,
        ...
      )
      
      difference1 <- which(object$adjusted_pval_factors[var, ] < alpha1)
      if (length(difference1) > 0) {
        for (j in 1:length(difference1)) {
          min_rect <- abscissa_pval[difference1[j]] - 
            (abscissa_pval[2] - abscissa_pval[1]) / 2
          max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
          graphics::rect(
            min_rect,
            graphics::par("usr")[3],
            max_rect,
            graphics::par("usr")[4],
            col = "gray90",
            density = -2,
            border = NA
          )
        }
        graphics::rect(
          graphics::par("usr")[1],
          graphics::par("usr")[3],
          graphics::par("usr")[2],
          graphics::par("usr")[4],
          col = NULL,
          border = "black"
        )
      }
      
      difference2 <- which(object$adjusted_pval_factors[var, ] < alpha2)
      if (length(difference2) > 0) {
        for (j in 1:length(difference2)) {
          min_rect <- abscissa_pval[difference2[j]] - 
            (abscissa_pval[2] - abscissa_pval[1]) / 2
          max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
          graphics::rect(
            min_rect,
            graphics::par("usr")[3],
            max_rect,
            graphics::par("usr")[4],
            col = "gray80",
            density = -2,
            border = NA
          )
        }
        graphics::rect(
          graphics::par("usr")[1],
          graphics::par("usr")[3],
          graphics::par("usr")[2],
          graphics::par("usr")[4],
          col = NULL,
          border = "black"
        )
      }
      
      for (j in 0:10) {
        graphics::abline(h = j / 10, col = 'lightgray', lty = "dotted")
      }
      
      graphics::lines(
        abscissa_pval,
        object$adjusted_pval_factors[var, ],
        type = type,
        lwd = 2,
        ...
      )
    }
  }
  
  grDevices::devAskNewPage(ask = FALSE)  
}
