#' Two population Partition Closed Testing procedure
#'
#' The function implements the Partition Closed Testing procedure for testing
#' mean differences between two functional populations. Functional data are
#' tested locally and unadjusted and adjusted p-value functions are provided.
#' The unadjusted p-value function controls the point-wise error rate. The
#' adjusted p-value function controls the family-wise error rate asymptotically.
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
#' @param paired Flag indicating whether a paired test has to be performed.
#'   Default is \code{FALSE}.
#' @param dx Used only if a \code{fd} object is provided. In this case,
#'   \code{dx} is the size of the discretization step of the grid  used to
#'   evaluate functional data. If set to \code{NULL}, a grid of size 100 is
#'   used. Default is \code{NULL}.
#' @param partition Vector of length \code{J} containing the labels assigning
#'   each point of the domain to an element of the partition.
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
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. \emph{Biometrics} 79(2), 1119-1132.
#'
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data.
#' \emph{Journal of Nonparametric Statistics}, 29(2), 407-424.
#'
#' @export
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#'
#' # Performing the PCT for two populations
#' # Choosing as partition the 4 seasons of the year
#' partition <- c(
#'   rep(1, 31 + 28 + 21), 
#'   rep(2, 10 + 30 + 31 + 21), 
#'   rep(3, 9 + 31 + 31 + 23), 
#'   rep(4, 7 + 31 + 30 + 21), 
#'   rep(1, 10)
#' )
#' partition <- factor(partition)
#' 
#' PCT.result <- PCT2(NASAtemp$paris, NASAtemp$milan, partition = partition)
#'
#' # Plotting the results of the PCT
#' plot(
#'   PCT.result, 
#'   xrange = c(0, 12), 
#'   main = 'PCT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(PCT.result$adjusted_pval < 0.05)
PCT2 <- function(data1, data2, partition, 
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
  
  #splines coefficients:
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
    two.sided = (meandiff)^2,
    greater   = (meandiff*sign.diff)^2,
    less      = (meandiff*(sign.diff-1))^2
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
      two.sided = (meandiff)^2,
      greater   = (meandiff*sign.diff)^2,
      less      = (meandiff*(sign.diff-1))^2
    )
  }
  
  pval <- numeric(p)
  for (i in 1:p) {
    pval[i] <- sum(T_coeff[, i] >= T0[i]) / B
  }
  
  #combination
  partition <- factor(partition)
  nintervals <- length(levels(partition))
  ntests <- 2^nintervals - 1
  all_combs <- matrix(nrow = ntests, ncol = p)
  labels <- levels(partition)
  tt <- 1
  for (nint in 1:nintervals) {
    combinations <- utils::combn(labels, nint)
    n.comb <- dim(combinations)[2]
    for (comb in 1:n.comb) {
      index <- rep(0, p)
      for (ii in 1:dim(combinations)[1]) {
        index <- index + as.numeric(partition == combinations[ii, comb])
      }
      all_combs[tt, ] <- index
      tt <- tt + 1
    }
  }
  
  #interval-wise tests
  adjusted.pval <- numeric(p)
  responsible.test <- matrix(nrow = p, ncol = p)
  for (test in 1:ntests) {
    T0_comb <- sum(T0[which(all_combs[test, ] == 1)])
    T_comb <- (rowSums(T_coeff[, which(all_combs[test, ] == 1), drop = FALSE]))
    pval.temp <- mean(T_comb >= T0_comb)
    indexes <- which(all_combs[test, ] == 1)
    max <- apply(rbind(adjusted.pval[indexes], pval.temp), 2, which.max)
    adjusted.pval[indexes] <- apply(rbind(adjusted.pval[indexes], pval.temp), 2, max)
    if (2 %in% max) {
      responsible.test[indexes[which(max == 2)] , ] <- matrix(
        data = all_combs[test, ],
        nrow = sum((max == 2)),
        ncol = p,
        byrow = TRUE
      )
    }
  }
  
  out <- list(
    test = '2pop',
    mu = mu.eval,
    adjusted_pval = adjusted.pval,
    unadjusted_pval = pval,
    data.eval = data.eval,
    ord_labels = etichetta_ord
  )
  class(out) <- "fdatest2"
  out
}
