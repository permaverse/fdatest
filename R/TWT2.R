#' Two population Threshold Wise Testing procedure
#'
#' The function implements the Threshold Wise Testing procedure for testing mean
#' differences between two functional populations. Functional data are tested
#' locally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the family-wise error rate asymptotically.
#'
#' @inherit functional_two_sample_test params return seealso
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
#' # Performing the TWT for two populations
#' TWT.result <- TWT2(NASAtemp$paris, NASAtemp$milan)
#'
#' # Plotting the results of the TWT
#' plot(
#'   TWT.result,
#'   xrange = c(0, 12),
#'   title = 'TWT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(TWT.result$adjusted_pvalues < 0.05)
TWT2 <- function(
  data1,
  data2,
  mu = 0,
  dx = NULL,
  B = 1000L,
  paired = FALSE,
  alternative = c("two.sided", "less", "greater"),
  verbose = FALSE
) {
  alternative <- rlang::arg_match(alternative)

  inputs <- twosamples2coeffs(data1, data2, mu, dx = dx)
  coeff1 <- inputs$coeff1
  coeff2 <- inputs$coeff2
  mu.eval <- inputs$mu

  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  n <- n1 + n2
  data.eval <- rbind(coeff1, coeff2)
  p <- dim(data.eval)[2]
  coeff1 <- coeff1 - matrix(data = mu.eval, nrow = n1, ncol = p)

  coeff <- rbind(coeff1, coeff2)
  etichetta_ord <- c(rep(1, n1), rep(2, n2))

  # First part:
  # univariate permutation test for each point
  # this is the computation that needs to be done for each voxel
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
    # loop on random permutations
    if (paired) {
      # paired test (for brain data we will not need it)
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
      # unpaired test
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

  # p-value computation
  pval <- numeric(p)
  for (i in 1:p) {
    pval[i] <- sum(T_coeff[, i] >= T0[i]) / B
  }

  # Second part:
  # combination into subsets
  if (verbose) cli::cli_h1("Threshold-wise tests")

  thresholds <- c(0, sort(unique(pval)), 1)
  adjusted.pval <- pval # we initialize the adjusted p-value as unadjusted one
  pval.tmp <- rep(0, p) # inizialize p-value vector resulting from combined test
  for (test in 1:length(thresholds)) {
    # test below threshold
    points.1 <- which(pval <= thresholds[test])
    T0_comb <- sum(T0[points.1], na.rm = TRUE) # combined test statistic
    T_comb <- (rowSums(T_coeff[, points.1, drop = FALSE], na.rm = TRUE))
    pval.test <- mean(T_comb >= T0_comb)
    pval.tmp[points.1] <- pval.test
    # compute maximum
    adjusted.pval <- apply(rbind(adjusted.pval, pval.tmp), 2, max)

    # test above threshold
    points.2 <- which(pval > thresholds[test])
    T0_comb <- sum(T0[points.2]) # combined test statistic
    T_comb <- (rowSums(T_coeff[, points.2, drop = FALSE], na.rm = TRUE))
    pval.test <- mean(T_comb >= T0_comb)
    pval.tmp[points.2] <- pval.test
    # compute maximum
    adjusted.pval <- apply(rbind(adjusted.pval, pval.tmp), 2, max)
  }

  out <- list(
    data = coeff,
    group_labels = etichetta_ord,
    mu = mu.eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = adjusted.pval
  )
  class(out) <- "ftwosample"
  out
}
