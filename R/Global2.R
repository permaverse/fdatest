#' Two population Global Testing procedure
#'
#' The function implements the Global Testing procedure for testing mean
#' differences between two functional populations. Functional data are tested
#' locally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the interval-wise error rate.
#'
#' @inherit functional_two_sample_test params return seealso
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
#' which(Global.result$adjusted_pvalues < 0.05)
Global2 <- function(data1, data2, 
                    mu = 0, 
                    dx = NULL, 
                    B = 1000L, 
                    paired = FALSE, 
                    statistic = c("Integral", "Max", "Integral_std", "Max_std")) {
  statistic <- rlang::arg_match(statistic)
  inputs <- twosamples2coeffs(data1, data2, mu, dx = dx)
  coeff1 <- inputs$coeff1
  coeff2 <- inputs$coeff2
  mu.eval <- inputs$mu
  
  # Check the statistic
  
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
    statistic,
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
      statistic,
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
  
  if (statistic =='Integral' || statistic == 'Integral_std') {
    T0_comb <- sum(T0[which(all_combs[1, ] == 1)])
    T_comb <- (rowSums(T_coeff[, which(all_combs[1, ] == 1), drop = FALSE]))
    pval.temp <- mean(T_comb >= T0_comb)
    indexes <- which(all_combs[1, ] == 1)
    adjusted.pval[indexes] <- pval.temp
  } else if (statistic == 'Max' || statistic == 'Max_std') {
    T0_comb <- max(T0[which(all_combs[1, ] == 1)])
    T_comb <- (apply(T_coeff[, which(all_combs[1, ] == 1)], 1, max))
    pval.temp <- mean(T_comb >= T0_comb)
    indexes <- which(all_combs[1, ] == 1)
    adjusted.pval[indexes] <- pval.temp
  }
  
  out <- list(
    data = data.eval,
    group_labels = etichetta_ord,
    mu = mu.eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = adjusted.pval,
    global_pvalue = adjusted.pval[1]
  )
  class(out) <- 'ftwosample'
  out
}
