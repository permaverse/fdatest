#' Global testing procedure for testing functional analysis of variance
#'
#' The function implements the Global Testing procedure for testing mean
#' differences between several functional populations in a one-way or multi-way
#' functional analysis of variance framework. Functional data are tested
#' globally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the family-wise error rate weakly. Since this is a
#' global test, the adjusted p-value function is constant.
#'
#' @inherit functional_anova_test params return seealso
#'
#' @references
#' - Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. *Biometrics* 79(2), 1119-1132.
#' - D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported
#' Significance Levels. *Journal of Business & Economic Statistics* 1.4,
#' 292-298.
#' - B. F. J. Manly (2006). Randomization, *Bootstrap and Monte Carlo
#' Methods in Biology*. Vol. 70. CRC Press.
#'
#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @export
#' @examples
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the test
#' Global_result <- global_aov(temperature ~ groups, n_perm = 1000L)
#'
#' # Summary of the test results
#' summary(Global_result)
#'
#' # Plot of the results
#' layout(1)
#' plot(Global_result)
#'
#' # All graphics on the same device
#' layout(matrix(1:4, nrow = 2, byrow = FALSE))
#' plot(
#'   Global_result,
#'   main = 'NASA data',
#'   plot.adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
global_aov <- function(
  formula,
  dx = NULL,
  n_perm = 1000L,
  method = c("residuals", "responses"),
  stat = c("Integral", "Max")
) {
  method <- rlang::arg_match(method)
  stat <- rlang::arg_match(stat)
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

  cli::cli_h1("Global test")
  if (stat == "Integral") {
    t0_comb <- sum(t0_glob)
    t_comb <- rowSums(t_glob)
    global_pval_f <- sum(t_comb >= t0_comb) / n_perm

    global_pval_factors <- numeric(nvar)
    for (ii in seq_len(nvar)) {
      t0_comb <- sum(t0_part[ii, ])
      t_comb <- rowSums(t_part[, ii, , drop = FALSE])
      global_pval_factors[ii] <- sum(t_comb >= t0_comb) / n_perm
    }
  } else {
    t0_comb <- max(t0_glob)
    t_comb <- apply(t_glob, 1, max)
    global_pval_f <- sum(t_comb >= t0_comb) / n_perm

    global_pval_factors <- numeric(nvar)
    for (ii in seq_len(nvar)) {
      t0_comb <- max(t0_part[ii, ])
      t_comb <- apply(t_part[, ii, , drop = FALSE], 1, max)
      global_pval_factors[ii] <- sum(t_comb >= t0_comb) / n_perm
    }
  }

  corrected_pval_glob <- rep(global_pval_f, p)
  corrected_pval_part <- matrix(nrow = nvar, ncol = p)
  for (ii in seq_len(nvar)) {
    corrected_pval_part[ii, ] <- rep(global_pval_factors[ii], p)
  }

  coeff_t <- regr0$coeff
  fitted_t <- regr0$fitted.values

  rownames(corrected_pval_part) <- var_names
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

  cli::cli_h1("Global Testing completed")

  out <- list(
    call = cl,
    design_matrix = design_matrix,
    unadjusted_pval_F = pval_glob,
    adjusted_pval_F = corrected_pval_glob,
    unadjusted_pval_factors = pval_part,
    adjusted_pval_factors = corrected_pval_part,
    Global_pval_F = global_pval_f,
    Global_pval_factors = global_pval_factors,
    data_eval = coeff,
    coeff_regr_eval = coeff_t,
    fitted_eval = fitted_t,
    residuals_eval = residuals_t,
    R2_eval = r2_t
  )
  class(out) <- "faov"
  out
}

#' @rdname global_aov
#' @export
Globalaov <- # nolint: object_name_linter.
  function(
    formula,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    method = c("residuals", "responses"),
    stat = c("Integral", "Max")
  ) {
    global_aov(
      formula = formula,
      dx = dx,
      n_perm = B,
      method = method,
      stat = stat
    )
  }
