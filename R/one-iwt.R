#' One population Interval Wise Testing procedure
#'
#' The function implements the Interval Wise Testing procedure for testing the
#' center of symmetry of a functional population. Functional data are tested
#' locally and  unadjusted  and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the interval-wise error rate.
#'
#' @param data Either pointwise evaluations of the functional data set on a
#'   uniform grid, or an \code{\link[fda]{fd}}. If pointwise evaluations are
#'   provided, `data` is a matrix of dimensions `c(n, J)`, with `J` evaluations
#'   on columns and `n` units on rows.
#' @param mu The center of symmetry under the null hypothesis. Three
#'   possibilities are available for `mu`:
#'   
#' - a constant (in this case, a constant function is used);
#' - a `J`-dimensional vector containing the evaluations on the same grid which
#' `data` are evaluated;
#' - a \code{\link[fda]{fd}} object containing one function.
#' Defaults to `0`.
#' @param B The number of iterations of the MC algorithm to evaluate the
#'   p-values of the permutation tests. Defaults to `1000L`.
#' @param dx Used only if an \code{\link[fda]{fd}} object is provided. In this
#'   case, `dx` is the size of the discretization step of the gridused to
#'   evaluate functional data. If set to `NULL`, a grid of size `100L` is used.
#'   Defaults to `NULL`.
#' @param recycle Flag used to decide whether the recycled version of the IWT
#'   should be used (see Pini and Vantini, 2017 for details). Defaults to
#'   `TRUE`.
#'
#' @return An object of class \code{\link{IWT1}}, which is a list containing at
#'   least the following components:
#' 
#' - `test`: String vector indicating the type of test performed. In this case
#' equal to `"1pop"`.
#' - `mu`: Evaluation on a grid of the center of symmetry under the null
#' hypothesis (as entered by the user).
#' - `unadjusted_pval`: Evaluation on a grid of the unadjusted p-value function.
#' - `pval_matrix`: Matrix of dimensions `c(p, p)` of the p-values of the
#' multivariate tests. The element `(i,j)` of matrix `pval.matrix` contains the
#' p-value of the joint NPC test of the components `(j,j+1,...,j+(p-i))`.
#' - `adjusted_pval`: Evaluation on a grid of the adjusted p-value function.
#' - `data.eval`: Evaluation on a grid of the functional data.
#'
#' @seealso See also \code{\link{plot.IWT1}} and \code{\link{IWTimage}} for
#'   plotting the results.
#'
#' @references 
#' A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference
#' for Functional Data Controlling the Family Wise Error Rate on Intervals.
#' *Biometrics*, 73(3): 835–845.
#' 
#' A. Pini and S. Vantini (2017). Interval-wise testing for functional data.
#' *Journal of Nonparametric Statistics*, 29(2), 407-424.
#'
#' @export
#' @examples
#' # Performing the IWT for one population
#' IWT.result <- IWT1(NASAtemp$paris, mu = 4, B = 10L)
#'
#' # Plotting the results of the IWT
#' plot(IWT.result, xrange = c(0, 12), main = 'Paris temperatures')
#'
#' # Plotting the p-value heatmap
#' IWTimage(IWT.result, abscissa_range = c(0, 12))
#'
#' # Selecting the significant components at 5% level
#' which(IWT.result$adjusted_pval < 0.05)
IWT1 <- function(data, mu = 0, B = 1000, dx = NULL, recycle = TRUE) {
  # data preprocessing
  inputs <- onesample2coeffs(data, mu, dx = dx)
  coeff <- inputs$coeff
  mu.eval <- inputs$mu

  n <- dim(coeff)[1]
  p <- dim(coeff)[2]
  data.eval <- coeff <- coeff - matrix(
    data = mu.eval,
    nrow = n,
    ncol = p,
    byrow = TRUE
  )

  #univariate permutations
  cli::cli_h1("Point-wise tests")
  
  T0 <- abs(colMeans(coeff))^2  #sample mean
  T_coeff <- matrix(ncol = p, nrow = B)
  for (perm in 1:B) {
    signs <- stats::rbinom(n, 1, 0.5) * 2 - 1
    coeff_perm <- coeff * signs
    T_coeff[perm, ] <- abs(colMeans(coeff_perm))^2
  }
  pval <- numeric(p)
  for (i in 1:p) {
    pval[i] <- sum(T_coeff[, i] >= T0[i]) / B
  }

  #combination
  cli::cli_h1("Interval-wise tests")

  #asymmetric combination matrix:
  matrice_pval_asymm <- matrix(nrow = p, ncol = p)
  matrice_pval_asymm[p, ] <- pval[1:p]
  T0_2x <- c(T0, T0)
  T_coeff_2x <- cbind(T_coeff, T_coeff)

  maxrow <- 1

  if (recycle) {
    for (i in (p - 1):maxrow) { # rows
      for (j in 1:p) { # columns
        inf <- j
        sup <- (p - i) + j
        T0_temp <- sum(T0_2x[inf:sup])
        T_temp <- rowSums(T_coeff_2x[, inf:sup])
        pval_temp <- sum(T_temp >= T0_temp) / B
        matrice_pval_asymm[i, j] <- pval_temp
      }
      cli::cli_h1("creating the p-value matrix: end of row {p - i + 1} out of {p}")
    }
  } else { # without recycling
    for (i in (p - 1):maxrow) { # rows
      for (j in 1:i) { # columns
        inf <- j
        sup <- (p - i) + j
        T0_temp <- sum(T0_2x[inf:sup])
        T_temp <- rowSums(T_coeff_2x[, inf:sup])
        pval_temp <- sum(T_temp >= T0_temp) / B
        matrice_pval_asymm[i, j] <- pval_temp
      }
      cli::cli_h1("Creating the p-value matrix: end of row {p - i + 1} out of {p}")
    }
  }

  corrected.pval.matrix <- pval_correct(matrice_pval_asymm)
  corrected.pval <- corrected.pval.matrix[1, ]
  
  cli::cli_h1("Interval-Wise Testing completed")
  
  out <- list(
    test = '1pop',
    mu = mu.eval,
    adjusted_pval = corrected.pval,
    unadjusted_pval = pval,
    pval_matrix = matrice_pval_asymm,
    data.eval = data.eval
  )
  class(out) <- 'IWT1'
  out
}
