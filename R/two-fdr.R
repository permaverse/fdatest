#' Two population functional Benjamini-Hochberg procedure
#'
#' The function implements the functional Benjamini Hochberg (fBH) procedure for
#' testing mean differences between two functional populations. Functional data
#' are tested locally and unadjusted and adjusted p-value functions are
#' provided. The unadjusted p-value function controls the point-wise error rate.
#' The adjusted p-value function controls the family-wise error rate
#' asymptotically.
#'
#' @inherit functional_two_sample_test params return seealso
#'
#' @references
#' Lundtorp Olsen, N., Pini, A., & Vantini, S. (2021). False discovery rate for
#' functional data \emph{TEST} 30, 784–809.
#'
#' @export
#' @examples
#' # Performing the fBH for two populations
#'
#' FDR.result <- FDR2(NASAtemp$paris, NASAtemp$milan)
#'
#' # Plotting the results of the fBH
#' plot(
#'   FDR.result,
#'   xrange = c(0, 12),
#'   title = 'FDR results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(FDR.result$adjusted_pvalues < 0.05)
FDR2 <- function(
  data1,
  data2,
  mu = 0,
  dx = NULL,
  B = 1000L,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater")
) {
  alternative <- rlang::arg_match(alternative)

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
    two.sided = (meandiff)^2,
    greater = (meandiff * sign.diff)^2,
    less = (meandiff * (sign.diff - 1))^2
  )

  T_coeff <- matrix(ncol = p, nrow = B)
  for (perm in 1:B) {
    if (paired) {
      if.perm <- stats::rbinom(n1, 1, 0.5)
      coeff_perm <- coeff
      for (couple in 1:n1) {
        if (if.perm[couple] == 1) {
          coeff_perm[c(couple, n1 + couple), ] <- coeff[
            c(n1 + couple, couple),
          ]
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
      greater = (meandiff * sign.diff)^2,
      less = (meandiff * (sign.diff - 1))^2
    )
  }

  pval <- numeric(p)
  for (i in 1:p) {
    pval[i] <- sum(T_coeff[, i] >= T0[i]) / B
  }

  #combination
  adjusted.pval <- stats::p.adjust(pval, method = 'BH')

  out <- list(
    data = data.eval,
    group_labels = etichetta_ord,
    mu = mu.eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = adjusted.pval
  )
  class(out) <- 'ftwosample'
  out
}
