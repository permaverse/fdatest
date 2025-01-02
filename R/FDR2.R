#' Two population functional Benjamini-Hochberg procedure
#'
#' The function implements the functional Benjamini Hochberg (fBH) procedure for
#' testing mean differences between two functional populations. Functional data
#' are tested locally and unadjusted and adjusted p-value functions are
#' provided. The unadjusted p-value function controls the point-wise error rate.
#' The adjusted p-value function controls the family-wise error rate
#' asymptotically.
#'
#' @param data1 First population's data. Either pointwise evaluations of the
#'   functional data set on a uniform grid, or a \code{fd} object from the
#'   package \code{fda}. If pointwise evaluations are provided, \code{data2} is
#'   a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns
#'   and \code{n1} units on rows.
#' @param data2 Second population's data. Either pointwise evaluations of the
#'   functional data set on a uniform grid, or a \code{fd} object from the
#'   package \code{fda}. If pointwise evaluations are provided, \code{data2} is
#'   a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns
#'   and \code{n2} units on rows.
#' @param mu Functional mean difference under the null hypothesis. Three
#'   possibilities are available for \code{mu}: a constant (in this case, a
#'   constant function is used); a \code{J}-dimensional vector containing the
#'   evaluations on the same grid which \code{data} are evaluated; a \code{fd}
#'   object from the package \code{fda} containing one function. The default is
#'   \code{mu=0}.
#' @param B The number of iterations of the MC algorithm to evaluate the
#'   p-values of the permutation tests. The defualt is \code{B=1000}.
#' @param paired A logical indicating whether a paired test has to be performed.
#'   Default is \code{FALSE}.
#' @param dx Used only if a \code{fd} object is provided. In this case,
#'   \code{dx} is the size of the discretization step of the grid  used to
#'   evaluate functional data. If set to \code{NULL}, a grid of size 100 is
#'   used. Default is \code{NULL}.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of "\code{two.sided}" (default), "\code{greater}" or
#'   "\code{less}".
#'
#' @return An object of class `fdatest2` containing the following components:
#'   
#'   - `test`: String vector indicating the type of test performed. In this case
#'   equal to \code{"2pop"}.
#'   - `mu`: Evaluation on a grid of the functional mean difference under the
#'   null hypothesis (as entered by the user).
#'   - `unadjusted_pval`: Evaluation on a grid of the unadjusted p-value
#'   function.
#'   - `adjusted_pval`: Evaluation on a grid of the adjusted p-value function.
#'   - `data.eval`: Evaluation on a grid of the functional data.
#'   - `ord_labels`: Vector of labels indicating the group membership of
#'   `data.eval`.
#'
#' @seealso See also \code{\link{plot.fdatest2}} for plotting the results.
#'
#' @references
#' Lundtorp Olsen, N., Pini, A., & Vantini, S. (2021). False discovery rate for
#' functional data \emph{TEST} 30, 784–809.
#'
#' @export
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#'
#' # Performing the fBH for two populations
#' 
#' FDR.result <- FDR2(NASAtemp$paris, NASAtemp$milan)
#'
#' # Plotting the results of the fBH
#' plot(
#'   FDR.result, 
#'   xrange = c(0, 12), 
#'   main = 'FDR results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(FDR.result$adjusted_pval < 0.05)
FDR2 <- function(data1, data2, 
                 mu = 0, 
                 B = 1000L, 
                 paired = FALSE, 
                 dx = NULL, 
                 alternative = "two.sided") {
  alternative <- rlang::arg_match(alternative, values = AVAILABLE_ALTERNATIVES())
  
  inputs <- twosamples2coeffs(data1, data2, mu, dx = dx)
  coeff1 <- inputs$coeff1
  coeff2 <- inputs$coeff2
  mu.eval <- inputs$mu
  
  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  J <- dim(coeff1)[2]
  n <- n1 + n2
  etichetta_ord <- c(rep(1, n1), rep(2, n2))
  
  eval <- coeff <- rbind(coeff1, coeff2)
  p <- dim(coeff)[2]
  
  data.eval <- eval
  
  #univariate permutations
  
  meandiff <- colMeans(coeff[1:n1, , drop = FALSE], na.rm = TRUE) - 
    colMeans(coeff[(n1 + 1):n, , drop = FALSE], na.rm = TRUE)
  sign.diff <- sign(meandiff)
  sign.diff[which(sign.diff == -1)] <- 0
  T0 <- switch(
    alternative,
    two.sided =  (meandiff)^2, 
    greater   =  (meandiff * sign.diff)^2, 
    less      =  (meandiff * (sign.diff - 1))^2
  )
  
  T_coeff <- matrix(ncol = p, nrow = B)
  for (perm in 1:B) {
    if (paired) {
      if.perm <- stats::rbinom(n1, 1, 0.5)
      coeff_perm <- coeff
      for (couple in 1:n1) {
        if (if.perm[couple] == 1) {
          coeff_perm[c(couple, n1 + couple), ] <- coeff[c(n1 + couple, couple), ]
        }
      }
    } else {
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni, ]
    }
    
    meandiff <- colMeans(coeff_perm[1:n1, , drop = FALSE], na.rm = TRUE) - 
      colMeans(coeff_perm[(n1 + 1):n, , drop = FALSE], na.rm = TRUE)
    sign.diff <- sign(meandiff)
    sign.diff[which(sign.diff == -1)] <- 0
    T_coeff[perm, ] <- switch(
      alternative,
      two.sided =  (meandiff)^2, 
      greater   =  (meandiff * sign.diff)^2, 
      less      =  (meandiff * (sign.diff - 1))^2
    )
  }
  
  pval <- numeric(p)
  for (i in 1:p) {
    pval[i] <- sum(T_coeff[, i] >= T0[i]) / B
  }
  
  #combination
  adjusted.pval <- stats::p.adjust(pval, method = 'BH')
  
  out <- list(
    test = '2pop',
    mu = mu.eval,
    adjusted_pval = adjusted.pval,
    unadjusted_pval = pval,
    data.eval = data.eval,
    ord_labels = etichetta_ord
  )
  class(out) <- 'fdatest2'
  out
}

