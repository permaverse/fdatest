#' Two population Global Testing procedure
#'
#' The function implements the Global Testing procedure for testing mean
#' differences between two functional populations. Functional data are tested
#' locally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the interval-wise error rate.
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
#' @param stat Test statistic used for the global test. Possible values are:
#'   \code{"Integral"}: integral of the squared sample mean difference;
#'   \code{"Max"}: maximum of the squared sample mean difference;
#'   \code{"Integral_std"}: integral of the squared t-test statistic;
#'   \code{"Max_std"}: maximum of the squared t-test statistic. Default is
#'   \code{"Integral"}.
#'
#' @return An object of class `fdatest2`, containing the following components:
#'   
#'   - `test`: String vector indicating the type of test performed. In this case
#'   equal to \code{"2pop"}.
#'   - `mu`: Evaluation on a grid of the functional mean difference under the
#'   null hypothesis (as entered by the user).
#'   - `unadjusted_pval`: Evaluation on a grid of the unadjusted p-value
#'   function (it is a constant function according to the global testing
#'   procedure).
#'   - `adjusted_pval`: Evaluation on a grid of the adjusted p-value function.
#'   - `data.eval`: Evaluation on a grid of the functional data.
#'   - `ord_labels`: Vector of labels indicating the group membership of
#'   `data.eval`.
#'
#' @seealso See also \code{\link{IWT2}} for local inference. See
#'   \code{\link{plot.fdatest2}} for plotting the results.
#'
#' @references
#' A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference for
#' Functional Data Controlling the Family Wise Error Rate on Intervals.
#' Biometrics 73(3): 835–845.
#'
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data.
#' \emph{Journal of Nonparametric Statistics}, 29(2), 407-424
#'
#' @export
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#'
#' # Performing the Global for two populations
#' Global.result <- Global2(NASAtemp$paris, NASAtemp$milan)
#'
#' # Plotting the results of the Global
#' plot(
#'   Global.result, 
#'   xrange = c(0, 12), 
#'   main = 'Global results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(Global.result$adjusted_pval < 0.05)
Global2 <- function(data1, data2, 
                    mu = 0, 
                    B = 1000L, 
                    paired = FALSE, 
                    dx = NULL, 
                    stat = 'Integral') {
  inputs <- twosamples2coeffs(data1, data2, mu, dx = dx)
  coeff1 <- inputs$coeff1
  coeff2 <- inputs$coeff2
  mu.eval <- inputs$mu
  
  possible_statistics <- c("Integral", "Integral_std", "Max", "Max_std")
  if (!(stat %in% possible_statistics)) {
    stop(paste0(
      'Possible statistics are ',
      paste0(possible_statistics, collapse = ', ')
    ))
  }
  
  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  J <- dim(coeff1)[2]
  n <- n1+n2
  etichetta_ord <- c(rep(1, n1), rep(2, n2))
  
  #splines coefficients:
  eval <- coeff <- rbind(coeff1,coeff2)
  p <- dim(coeff)[2]
  
  data.eval <- eval
  
  #univariate permutations
  meandiff2 <- (colMeans(coeff[1:n1, ]) - colMeans(coeff[(n1 + 1):n, ]))^2
  S1 <- stats::cov(coeff[1:n1, ])
  S2 <- stats::cov(coeff[(n1 + 1):n, ])
  Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  T0 <- switch(
    stat,
    Integral     = meandiff2,
    Max          = meandiff2,
    Integral_std = meandiff2 / diag(Sp),
    Max_std      = meandiff2 / diag(Sp)
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
    
    meandiff2_perm <- (colMeans(coeff_perm[1:n1, ]) - colMeans(coeff_perm[(n1 + 1):n, ]))^2
    S1_perm <- stats::cov(coeff_perm[1:n1, ])
    S2_perm <- stats::cov(coeff_perm[(n1 + 1):n, ])
    Sp_perm <- ((n1 - 1) * S1_perm + (n2 - 1) * S2_perm) / (n1 + n2 - 2)
    T_coeff[perm, ] <- switch(
      stat,
      Integral     = meandiff2_perm,
      Max          = meandiff2_perm,
      Integral_std = meandiff2_perm / diag(Sp_perm), 
      Max_std      = meandiff2_perm / diag(Sp_perm)
    )
    
  }
  pval <- numeric(p)
  for (i in 1:p) {
    pval[i] <- sum(T_coeff[, i] >= T0[i]) / B
  }
  
  #combination
  all_combs <- rbind(rep(1, p))
  ntests <- 1
  adjusted.pval <- numeric(p)
  
  if (stat =='Integral' || stat == 'Integral_std') {
    T0_comb <- sum(T0[which(all_combs[1, ] == 1)])
    T_comb <- (rowSums(T_coeff[, which(all_combs[1, ] == 1), drop = FALSE]))
    pval.temp <- mean(T_comb >= T0_comb)
    indexes <- which(all_combs[1, ] == 1)
    adjusted.pval[indexes] <- pval.temp
  } else if (stat == 'Max' || stat == 'Max_std') {
    T0_comb <- max(T0[which(all_combs[1, ] == 1)])
    T_comb <- (apply(T_coeff[, which(all_combs[1, ] == 1)], 1, max))
    pval.temp <- mean(T_comb >= T0_comb)
    indexes <- which(all_combs[1, ] == 1)
    adjusted.pval[indexes] <- pval.temp
  }
  
  out <- list(
    test = '2pop', 
    mu = mu.eval,
    adjusted_pval = adjusted.pval,
    unadjusted_pval = pval,
    data.eval=data.eval,
    ord_labels = etichetta_ord,
    global_pval = adjusted.pval[1]
  )
  class(out) <- 'fdatest2'
  out
}
