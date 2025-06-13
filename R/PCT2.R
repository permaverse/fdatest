#' Two population Partition Closed Testing procedure
#'
#' The function implements the Partition Closed Testing procedure for testing
#' mean differences between two functional populations. Functional data are
#' tested locally and unadjusted and adjusted p-value functions are provided.
#' The unadjusted p-value function controls the point-wise error rate. The
#' adjusted p-value function controls the family-wise error rate asymptotically.
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
#'   title = 'PCT results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(PCT.result$adjusted_pvalues < 0.05)
PCT2 <- function(
  data1,
  data2,
  partition,
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
    adjusted.pval[indexes] <- apply(
      rbind(adjusted.pval[indexes], pval.temp),
      2,
      max
    )
    if (2 %in% max) {
      responsible.test[indexes[which(max == 2)], ] <- matrix(
        data = all_combs[test, ],
        nrow = sum((max == 2)),
        ncol = p,
        byrow = TRUE
      )
    }
  }

  out <- list(
    data = data.eval,
    group_labels = etichetta_ord,
    mu = mu.eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = adjusted.pval
  )
  class(out) <- "ftwosample"
  out
}
