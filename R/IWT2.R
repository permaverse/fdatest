#' Two population Interval Wise Testing procedure
#'
#' The function implements the Interval Wise Testing procedure for testing mean
#' differences between two functional populations. Functional data are tested
#' locally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the interval-wise error rate.
#'
#' @inherit functional_two_sample_test params return seealso
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
#' # Performing the IWT for two populations
#' IWT.result <- IWT2(NASAtemp$paris, NASAtemp$milan, B = 10L)
#'
#' # Plotting the results of the IWT
#' plot(
#'   IWT.result, 
#'   xrange = c(0, 12), 
#'   main = 'IWT results for testing mean differences'
#' )
#'
#' # Plotting the p-value heatmap
#' IWTimage(IWT.result, abscissa_range = c(0, 12))
#'
#' # Selecting the significant components at 5% level
#' which(IWT.result$adjusted_pvalues < 0.05)
IWT2 <- function(data1, data2, 
                 mu = 0, 
                 dx = NULL, 
                 B = 1000L, 
                 paired = FALSE, 
                 alternative = c("two.sided", "less", "greater"), 
                 verbose = FALSE,
                 recycle = TRUE) {
  alternative <- rlang::arg_match(alternative)

  # data preprocessing
  inputs <- twosamples2coeffs(data1, data2, mu, dx = dx)
  coeff1 <- inputs$coeff1
  coeff2 <- inputs$coeff2
  mu.eval <- inputs$mu

  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  p <- dim(coeff1)[2]
  n <- n1 + n2
  etichetta_ord <- c(rep(1, n1), rep(2, n2))
  coeff1 <- coeff1 - matrix(data = mu.eval, nrow = n1, ncol = p)

  #splines coefficients:
  eval <- coeff <- rbind(coeff1, coeff2)

  data.eval <- eval
  data.eval[1:n1, ] <- data.eval[1:n1, ] + matrix(
    data = mu.eval, nrow = n1, ncol = p
  )

  if (verbose)
    cli::cli_h1("Point-wise tests")
  
  #univariate permutations
  meandiff <- colMeans(coeff[1:n1, , drop = FALSE], na.rm = TRUE) - 
    colMeans(coeff[(n1 + 1):n, , drop = FALSE], na.rm = TRUE)
  sign.diff <- sign(meandiff)
  sign.diff[which(sign.diff == -1)] <- 0
  T0 <- switch(
    alternative,
    two.sided = (meandiff)^2,
    greater   = (meandiff * sign.diff)^2,
    less      = (meandiff * (sign.diff - 1))^2
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
      greater   = (meandiff * sign.diff)^2,
      less      = (meandiff * (sign.diff - 1))^2
    )
  }
  
  pval <- numeric(p)
  for (i in 1:p) {
    pval[i] <- sum(T_coeff[, i] >= T0[i]) / B
  }

  #combination
  if (verbose)
    cli::cli_h1("Interval-wise tests")

  #asymmetric combination matrix:
  matrice_pval_asymm <- matrix(nrow = p, ncol = p)
  matrice_pval_asymm[p,] <- pval[1:p]
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
      
      if (verbose)
        cli::cli_h1("Creating the p-value matrix: end of row {p - i + 1} out of {p}")
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
      
      if (verbose)
        cli::cli_h1("Creating the p-value matrix: end of row {p - i + 1} out of {p}")
    }
  }

  corrected.pval.matrix <- pval_correct(matrice_pval_asymm)
  corrected.pval <- corrected.pval.matrix[1, ]

  if (verbose)
    cli::cli_h1("Interval-Wise Testing completed")
  
  out <- list(
    data = data.eval,
    group_labels = etichetta_ord,
    mu = mu.eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = corrected.pval,
    pvalue_matrix = matrice_pval_asymm
  )
  class(out) <- 'ftwosample'
  out
}
