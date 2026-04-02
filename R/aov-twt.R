#' Threshold Wise Testing procedure for testing functional analysis of variance
#'
#' The function implements the Threshold Wise Testing procedure for testing mean
#' differences between several functional populations in a one-way or multi-way
#' functional analysis of variance framework. Functional data are tested locally
#' and unadjusted and adjusted p-value functions are provided. The unadjusted
#' p-value function controls the point-wise error rate. The adjusted p-value
#' function controls the threshold-wise error rate.
#'
#' @inherit functional_anova_test params return seealso
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. *Biometrics* 79(2), 1119-1132.
#'
#' D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported
#' Significance Levels. *Journal of Business & Economic Statistics* 1.4,
#' 292-298.
#'
#' B. F. J. Manly (2006). Randomization, *Bootstrap and Monte Carlo Methods
#' in Biology*. Vol. 70. CRC Press.
#'
#' @export
#' @examples
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the TWT
#' TWT_result <- TWTaov(temperature ~ groups, B = 100L)
#'
#' # Summary of the TWT results
#' summary(TWT_result)
#'
#' # Plot of the TWT results
#' layout(1)
#' plot(TWT_result)
#'
#' # All graphics on the same device
#' layout(matrix(1:4, nrow = 2, byrow = FALSE))
#' plot(
#'   TWT_result,
#'   main = 'NASA data',
#'   plot_adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
TWTaov <- # nolint: object_name_linter.
  function(
    formula,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    method = c("residuals", "responses")
  ) {
    twt_aov(
      formula = formula,
      dx = dx,
      n_perm = B,
      method = method
    )
  }

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname TWTaov
#' @export
twt_aov <- function(
  formula,
  dx = NULL,
  n_perm = 1000L,
  method = c("residuals", "responses")
) {
  method <- rlang::arg_match(method)
  cl <- match.call()

  cli::cli_h1("Point-wise tests")
  res <- aov_permtest(formula, dx, n_perm, method)
  coeff <- res$coeff
  n <- res$n
  p <- res$p
  nvar <- res$nvar
  var_names <- res$var_names
  design_matrix <- res$design_matrix
  regr0 <- res$regr0
  t0_part <- res$t0_part
  t0_glob <- res$t0_glob
  t_glob <- res$t_glob
  t_part <- res$t_part
  pval_glob <- res$pval_glob
  pval_part <- res$pval_part

  cli::cli_h1("Threshold-wise tests")

  # F-test
  thresholds <- c(0, sort(unique(pval_glob)), 1)
  adjusted_pval_glob <- pval_glob
  pval_tmp <- rep(0, p)
  for (test in seq_along(thresholds)) {
    points_1 <- which(pval_glob <= thresholds[test])
    t0_comb <- sum(t0_glob[points_1], na.rm = TRUE)
    t_comb <- rowSums(t_glob[, points_1, drop = FALSE], na.rm = TRUE)
    pval_tmp[points_1] <- mean(t_comb >= t0_comb)
    adjusted_pval_glob <- apply(rbind(adjusted_pval_glob, pval_tmp), 2, max)

    points_2 <- which(pval_glob > thresholds[test])
    t0_comb <- sum(t0_glob[points_2])
    t_comb <- rowSums(t_glob[, points_2, drop = FALSE], na.rm = TRUE)
    pval_tmp[points_2] <- mean(t_comb >= t0_comb)
    adjusted_pval_glob <- apply(rbind(adjusted_pval_glob, pval_tmp), 2, max)
  }

  # F-tests on single factors
  thresholds <- c(0, sort(unique(as.numeric(pval_part))), 1)
  adjusted_pval_part <- pval_part

  for (ii in seq_len(nvar)) {
    pval_tmp <- rep(0, p)
    for (test in seq_along(thresholds)) {
      points_1 <- which(pval_part[ii, ] <= thresholds[test])
      t0_comb <- sum(t0_part[ii, points_1], na.rm = TRUE)
      t_comb <- rowSums(t_part[, ii, points_1, drop = FALSE], na.rm = TRUE)
      pval_tmp[points_1] <- mean(t_comb >= t0_comb)
      adjusted_pval_part[ii, ] <- apply(
        rbind(adjusted_pval_part[ii, ], pval_tmp),
        2,
        max
      )

      points_2 <- which(pval_part[ii, ] > thresholds[test])
      t0_comb <- sum(t0_part[ii, points_2])
      t_comb <- rowSums(t_part[, ii, points_2, drop = FALSE], na.rm = TRUE)
      pval_tmp[points_2] <- mean(t_comb >= t0_comb)
      adjusted_pval_part[ii, ] <- apply(
        rbind(adjusted_pval_part[ii, ], pval_tmp),
        2,
        max
      )
    }
  }

  coeff_t <- regr0$coeff
  fitted_t <- regr0$fitted.values

  rownames(adjusted_pval_part) <- var_names
  rownames(coeff_t) <- colnames(design_matrix)
  rownames(pval_part) <- var_names

  residuals_t <- coeff - fitted_t
  ybar_t <- colMeans(coeff)
  r2_t <- colSums(
    (fitted_t - matrix(ybar_t, nrow = n, ncol = p, byrow = TRUE))^2
  ) /
    colSums(
      (coeff - matrix(ybar_t, nrow = n, ncol = p, byrow = TRUE))^2
    )

  cli::cli_h1("Threshold-Wise Testing completed")

  out <- list(
    call = cl,
    design_matrix = design_matrix,
    unadjusted_pval_F = pval_glob,
    adjusted_pval_F = adjusted_pval_glob,
    unadjusted_pval_factors = pval_part,
    adjusted_pval_factors = adjusted_pval_part,
    data_eval = coeff,
    coeff_regr_eval = coeff_t,
    fitted_eval = fitted_t,
    residuals_eval = residuals_t,
    R2_eval = r2_t
  )
  class(out) <- "faov"
  out
}
