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
#' Global_result <- Global2(NASAtemp$paris, NASAtemp$milan)
#'
#' # Plotting the results of the Global
#' plot(
#'   Global_result,
#'   xrange = c(0, 12),
#'   title = 'Global results for testing mean differences'
#' )
#'
#' # Selecting the significant components at 5% level
#' which(Global_result$adjusted_pvalues < 0.05)
Global2 <- # nolint: object_name_linter.
  function(
    data1,
    data2,
    mu = 0,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    paired = FALSE,
    statistic = c("Integral", "Max", "Integral_std", "Max_std")
  ) {
    global2(
      data1 = data1,
      data2 = data2,
      mu = mu,
      dx = dx,
      n_perm = B,
      paired = paired,
      statistic = statistic
    )
  }

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname Global2
#' @export
global2 <- function(
  data1,
  data2,
  mu = 0,
  dx = NULL,
  n_perm = 1000L,
  paired = FALSE,
  statistic = c("Integral", "Max", "Integral_std", "Max_std")
) {
  statistic <- rlang::arg_match(statistic)
  inputs <- twosamples2coeffs(data1, data2, mu, dx = dx)
  coeff1 <- inputs$coeff1
  coeff2 <- inputs$coeff2
  mu_eval <- inputs$mu

  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  n <- n1 + n2
  etichetta_ord <- c(rep(1, n1), rep(2, n2))
  coeff <- rbind(coeff1, coeff2)
  p <- dim(coeff)[2]
  data_eval <- coeff

  meandiff2 <- (colMeans(coeff[seq_len(n1), ]) -
    colMeans(coeff[(n1 + 1):n, ]))^2 # nolint: indentation_linter.
  s1 <- stats::cov(coeff[seq_len(n1), ])
  s2 <- stats::cov(coeff[(n1 + 1):n, ])
  sp <- ((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2)
  t0 <- switch(
    statistic,
    Integral = meandiff2,
    Max = meandiff2,
    Integral_std = meandiff2 / diag(sp),
    Max_std = meandiff2 / diag(sp)
  )

  t_coeff <- matrix(ncol = p, nrow = n_perm)
  for (perm in seq_len(n_perm)) {
    if (paired) {
      if_perm <- stats::rbinom(n1, 1, 0.5)
      coeff_perm <- coeff
      for (couple in seq_len(n1)) {
        if (if_perm[couple] == 1) {
          coeff_perm[c(couple, n1 + couple), ] <- coeff[
            c(n1 + couple, couple),
          ]
        }
      }
    } else {
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni, ]
    }

    meandiff2_perm <- (colMeans(coeff_perm[seq_len(n1), ]) -
      colMeans(coeff_perm[(n1 + 1):n, ]))^2 # nolint: indentation_linter.
    s1_perm <- stats::cov(coeff_perm[seq_len(n1), ])
    s2_perm <- stats::cov(coeff_perm[(n1 + 1):n, ])
    sp_perm <- ((n1 - 1) * s1_perm + (n2 - 1) * s2_perm) / (n1 + n2 - 2)
    t_coeff[perm, ] <- switch(
      statistic,
      Integral = meandiff2_perm,
      Max = meandiff2_perm,
      Integral_std = meandiff2_perm / diag(sp_perm),
      Max_std = meandiff2_perm / diag(sp_perm)
    )
  }

  pval <- numeric(p)
  for (i in seq_len(p)) {
    pval[i] <- sum(t_coeff[, i] >= t0[i]) / n_perm
  }

  adjusted_pval <- if (statistic %in% c("Integral", "Integral_std")) {
    t0_comb <- sum(t0)
    t_comb <- rowSums(t_coeff)
    pval_temp <- mean(t_comb >= t0_comb)
    rep(pval_temp, p)
  } else {
    t0_comb <- max(t0)
    t_comb <- apply(t_coeff, 1, max)
    pval_temp <- mean(t_comb >= t0_comb)
    rep(pval_temp, p)
  }

  out <- list(
    data = data_eval,
    group_labels = etichetta_ord,
    mu = mu_eval,
    unadjusted_pvalues = pval,
    adjusted_pvalues = adjusted_pval,
    global_pvalue = adjusted_pval[1]
  )
  class(out) <- "ftwosample"
  out
}
