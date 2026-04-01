#' Interval Wise Testing procedure for testing functional analysis of variance
#'
#' The function implements the Interval Wise Testing procedure for testing mean
#' differences between several functional populations in a one-way or multi-way
#' functional analysis of variance framework. Functional data are tested locally
#' and unadjusted and adjusted p-value functions are provided. The unadjusted
#' p-value function controls the point-wise error rate. The adjusted p-value
#' function controls the interval-wise error rate.
#'
#' @inherit functional_anova_test params return seealso
#'
#' @references
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data.
#' *Journal of Nonparametric Statistics*, 29(2), 407-424.
#'
#' Pini, A., Vantini, S., Colosimo, B. M., & Grasso, M. (2018). Domain‐selective
#' functional analysis of variance for supervised statistical profile monitoring
#' of signal data. *Journal of the Royal Statistical Society: Series C
#' (Applied Statistics)* 67(1), 55-81.
#'
#' Abramowicz, K., Hager, C. K., Pini, A., Schelin, L., Sjostedt de Luna, S., &
#' Vantini, S. (2018). Nonparametric inference for functional‐on‐scalar linear
#' models applied to knee kinematic hop data after injury of the anterior
#' cruciate ligament. *Scandinavian Journal of Statistics* 45(4),
#' 1036-1061.
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
#' # Performing the IWT
#' IWT_result <- IWTaov(temperature ~ groups, B = 10L)
#'
#' # Summary of the IWT results
#' summary(IWT_result)
#'
#' # Plot of the IWT results
#' graphics::layout(1)
#' plot(IWT_result)
#'
#' # All graphics on the same device
#' graphics::layout(matrix(1:4, nrow = 2, byrow = FALSE))
#' plot(
#'   IWT_result,
#'   main = 'NASA data',
#'   plot.adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
IWTaov <- # nolint: object_name_linter.
  function(
    formula,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    method = c("residuals", "responses"),
    recycle = TRUE
  ) {
    iwt_aov(
      formula = formula,
      dx = dx,
      n_perm = B,
      method = method,
      recycle = recycle
    )
  }

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname IWTaov
#' @export
iwt_aov <- function(
  formula,
  dx = NULL,
  n_perm = 1000L,
  method = c("residuals", "responses"),
  recycle = TRUE
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

  cli::cli_h1("Interval-wise tests")

  matrice_pval_asymm_glob <- matrix(nrow = p, ncol = p)
  matrice_pval_asymm_glob[p, ] <- pval_glob[seq_len(p)]
  t0_2x_glob <- c(t0_glob, t0_glob)
  t_2x_glob <- cbind(t_glob, t_glob)

  matrice_pval_asymm_part <- array(dim = c(nvar, p, p))
  t0_2x_part <- cbind(t0_part, t0_part)
  t_2x_part <- array(dim = c(n_perm, nvar, p * 2L))
  for (ii in seq_len(nvar)) {
    matrice_pval_asymm_part[ii, p, ] <- pval_part[ii, seq_len(p)]
    t_2x_part[, ii, ] <- cbind(t_part[, ii, ], t_part[, ii, ])
  }

  row_indices <- (p - 1L):1L

  perm_args <- list(
    t0_2x_glob = t0_2x_glob,
    t_2x_glob = t_2x_glob,
    t0_2x_part = t0_2x_part,
    t_2x_part = t_2x_part,
    n_perm = n_perm,
    p = p,
    nvar = nvar,
    recycle = recycle
  )

  if (mirai::daemons_set()) {
    optimized_order <- optimize_order(row_indices)
    row_tasks <- mirai::mirai_map(
      (p - 1L):floor((p - 1L) / 2),
      function(.i) {
        rlang::inject(compute_row_pair_aov(.i, !!!perm_args))
      }
    )
    row_results <- row_tasks[.progress]
    row_results <- unlist(row_results, recursive = FALSE)
    row_results <- row_results[order(
      optimized_order,
      decreasing = TRUE
    )]
  } else {
    row_results <- lapply(row_indices, function(.i) {
      rlang::inject(compute_row_aov(.i, !!!perm_args))
    })
  }

  for (k in seq_along(row_indices)) {
    i <- row_indices[k]
    js <- if (recycle) seq_len(p) else seq_len(i)
    matrice_pval_asymm_glob[i, js] <- row_results[[k]]$glob
    for (ii in seq_len(nvar)) {
      matrice_pval_asymm_part[ii, i, js] <- row_results[[k]]$part[ii, ]
    }
    cli::cli_h1(
      "Creating the p-value matrix: end of row {p - i + 1} out of {p}"
    )
  }

  corrected_pval_matrix_glob <- pval_correct_cpp(matrice_pval_asymm_glob)
  corrected_pval_glob <- corrected_pval_matrix_glob[1, ]

  corrected_pval_part <- matrix(nrow = nvar, ncol = p)
  corrected_pval_matrix_part <- array(dim = c(nvar, p, p))
  for (ii in seq_len(nvar)) {
    corrected_pval_matrix_part[ii, , ] <- pval_correct_cpp(
      matrice_pval_asymm_part[ii, , ]
    )
    corrected_pval_part[ii, ] <- corrected_pval_matrix_part[ii, 1, ]
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

  cli::cli_h1("Interval-Wise Testing completed")

  out <- list(
    call = cl,
    design_matrix = design_matrix,
    unadjusted_pval_F = pval_glob,
    pval_matrix_F = matrice_pval_asymm_glob,
    adjusted_pval_F = corrected_pval_glob,
    unadjusted_pval_factors = pval_part,
    pval_matrix_factors = matrice_pval_asymm_part,
    adjusted_pval_factors = corrected_pval_part,
    data_eval = coeff,
    coeff_regr_eval = coeff_t,
    fitted_eval = fitted_t,
    residuals_eval = residuals_t,
    R2_eval = r2_t
  )
  class(out) <- "fanova"
  out
}
