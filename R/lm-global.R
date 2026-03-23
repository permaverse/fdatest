#' Global testing procedure for testing functional-on-scalar linear models
#'
#' The function is used to fit and test functional linear models. It can be used
#' to carry out regression, and analysis of variance. It implements the global
#' testing procedure for testing the significance of the effects of scalar
#' covariates on a functional population.
#'
#' @inherit functional_lm_test params return seealso
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified
#' framework. \emph{Biometrics} 79(2), 1119-1132.
#'
#' D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported
#' Significance Levels. \emph{Journal of Business & Economic Statistics} 1(4),
#' 292-298.
#'
#' B. F. J. Manly (2006). Randomization, \emph{Bootstrap and Monte Carlo Methods
#' in Biology}. Vol. 70. CRC Press.
#'
#' @export
#' @examples
#' # Defining the covariates
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the IWT
#' Global_result <- Globallm(temperature ~ groups, B = 1000)
#' # Summary of the IWT results
#' summary(Global_result)
#'
#' # Plot of the IWT results
#' plot(
#'   Global_result,
#'   main = 'NASA data',
#'   plot_adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
#'
#' plot(
#'   Global_result,
#'   main = 'NASA data',
#'   plot_adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
Globallm <- function( # nolint: object_name_linter.
  formula,
  dx = NULL,
  B = 1000L, # nolint: object_name_linter.
  method = c("residuals", "responses"),
  stat = c("Integral", "Max")
) {
  global_lm(
    formula = formula,
    dx = dx,
    n_perm = B,
    method = method,
    stat = stat
  )
}

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname Globallm
#' @export
global_lm <- function(
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
  res <- lm_permtest(formula, dx, n_perm, method)
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
    t0_temp <- sum(t0_glob[seq_len(p)])
    t_temp <- rowSums(t_glob[, seq_len(p), drop = FALSE])
    global_pval_f <- sum(t_temp >= t0_temp) / n_perm

    global_pval_part <- numeric(nvar + 1L)
    for (ii in seq_len(nvar + 1L)) {
      t0_temp <- sum(t0_part[ii, seq_len(p)])
      t_temp <- rowSums(t_part[, ii, seq_len(p), drop = FALSE])
      global_pval_part[ii] <- sum(t_temp >= t0_temp) / n_perm
    }

    corrected_pval_glob <- rep(global_pval_f, p)
    corrected_pval_part <- matrix(nrow = nvar + 1L, ncol = p)
    for (ii in seq_len(nvar + 1L)) {
      corrected_pval_part[ii, ] <- rep(global_pval_part[ii], p)
    }
  } else {
    t0_temp <- max(t0_glob)
    t_temp <- apply(t_glob, 1, max)
    global_pval_f <- sum(t_temp >= t0_temp) / n_perm

    global_pval_part <- numeric(nvar + 1L)
    for (ii in seq_len(nvar + 1L)) {
      t0_temp <- max(t0_part[ii, ])
      t_temp <- apply(t_part[, ii, ], 1, max)
      global_pval_part[ii] <- sum(t_temp >= t0_temp) / n_perm
    }

    corrected_pval_glob <- rep(global_pval_f, p)
    corrected_pval_part <- matrix(nrow = nvar + 1L, ncol = p)
    for (ii in seq_len(nvar + 1L)) {
      corrected_pval_part[ii, ] <- rep(global_pval_part[ii], p)
    }
  }

  coeff_t <- regr0$coeff
  fitted_t <- regr0$fitted

  rownames(corrected_pval_part) <- var_names
  rownames(coeff_t) <- var_names
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
    unadjusted_pval_part = pval_part,
    adjusted_pval_part = corrected_pval_part,
    Global_pval_F = global_pval_f,
    Global_pval_part = global_pval_part,
    data_eval = coeff,
    coeff_regr_eval = coeff_t,
    fitted_eval = fitted_t,
    residuals_eval = residuals_t,
    R2_eval = r2_t
  )
  class(out) <- "flm"
  out
}
