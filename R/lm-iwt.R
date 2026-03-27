#' Interval-wise testing procedure for testing functional-on-scalar linear
#' models
#'
#' The function is used to fit and test functional linear models. It can be used
#' to carry out regression, and analysis of variance. It implements the
#' interval-wise testing procedure (IWT) for testing the significance of the
#' effects of scalar covariates on a functional population.
#'
#' @inherit functional_lm_test params return seealso
#'
#' @references
#' A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference for
#' Functional Data Controlling the Family Wise Error Rate on Intervals.
#' Biometrics 73(3): 835–845.
#'
#' Pini, A., Vantini, S., Colosimo, B. M., & Grasso, M. (2018). Domain‐selective
#' functional analysis of variance for supervised statistical profile monitoring
#' of signal data. \emph{Journal of the Royal Statistical Society: Series C
#' (Applied Statistics)} 67(1), 55-81.
#'
#' Abramowicz, K., Hager, C. K., Pini, A., Schelin, L., Sjostedt de Luna, S., &
#' Vantini, S. (2018). Nonparametric inference for functional‐on‐scalar linear
#' models applied to knee kinematic hop data after injury of the anterior
#' cruciate ligament. \emph{Scandinavian Journal of Statistics} 45(4),
#' 1036-1061.
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
#' temperature <- rbind(NASAtemp$milan[, 1:100], NASAtemp$paris[, 1:100])
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the IWT
#' IWT_result <- IWTlm(temperature ~ groups, B = 2L)
#' # Summary of the IWT results
#' summary(IWT_result)
#'
#' # Plot of the IWT results
#' plot(
#'   IWT_result,
#'   main = 'NASA data',
#'   plot_adjpval = TRUE,
#'   xlab = 'Day',
#'   xrange = c(1, 365)
#' )
IWTlm <- # nolint: object_name_linter.
  function(
    formula,
    dx = NULL,
    B = 1000L, # nolint: object_name_linter.
    method = c("residuals", "responses"),
    recycle = TRUE
  ) {
    iwt_lm(
      formula = formula,
      dx = dx,
      n_perm = B,
      method = method,
      recycle = recycle
    )
  }

#' @param n_perm An integer value specifying the number of permutations for the
#'   permutation tests. Defaults to `1000L`.
#' @rdname IWTlm
#' @export
iwt_lm <- function(
  formula,
  dx = NULL,
  n_perm = 1000L,
  method = c("residuals", "responses"),
  recycle = TRUE
) {
  method <- rlang::arg_match(method)
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

  cli::cli_h1("Interval-wise tests")

  matrice_pval_asymm_glob <- matrix(nrow = p, ncol = p)
  matrice_pval_asymm_glob[p, ] <- pval_glob[seq_len(p)]
  t0_2x_glob <- c(t0_glob, t0_glob)
  t_2x_glob <- cbind(t_glob, t_glob)

  matrice_pval_asymm_part <- array(dim = c(nvar + 1L, p, p))
  t0_2x_part <- cbind(t0_part, t0_part)
  t_2x_part <- array(dim = c(n_perm, nvar + 1L, p * 2L))
  for (ii in seq_len(nvar + 1L)) {
    matrice_pval_asymm_part[ii, p, ] <- pval_part[ii, seq_len(p)]
    t_2x_part[, ii, ] <- cbind(t_part[, ii, ], t_part[, ii, ])
  }

  row_indices <- (p - 1L):1L

  compute_row_lm <- function(i) {
    js <- if (recycle) seq_len(p) else seq_len(i)
    glob_vals <- numeric(length(js))
    part_vals <- matrix(nrow = nvar + 1L, ncol = length(js))
    for (k in seq_along(js)) {
      j <- js[k]
      inf <- j
      sup <- (p - i) + j
      t0_temp <- sum(t0_2x_glob[inf:sup])
      t_temp <- rowSums(t_2x_glob[, inf:sup, drop = FALSE])
      glob_vals[k] <- sum(t_temp >= t0_temp) / n_perm
      for (ii in seq_len(nvar + 1L)) {
        t0_temp <- sum(t0_2x_part[ii, inf:sup])
        t_temp <- rowSums(t_2x_part[, ii, inf:sup, drop = FALSE])
        part_vals[ii, k] <- sum(t_temp >= t0_temp) / n_perm
      }
    }
    list(glob = glob_vals, part = part_vals)
  }

  if (mirai::daemons_set()) {
    perm_args <- list(
      compute_row_lm = compute_row_lm,
      t0_2x_glob = t0_2x_glob,
      t_2x_glob = t_2x_glob,
      t0_2x_part = t0_2x_part,
      t_2x_part = t_2x_part,
      n_perm = n_perm,
      p = p,
      nvar = nvar,
      recycle = recycle
    )
    row_tasks <- mirai::mirai_map(
      row_indices,
      \(.i) compute_row_lm(.i),
      .args = perm_args
    )
    row_results <- row_tasks[]
  } else {
    row_results <- lapply(row_indices, compute_row_lm)
  }

  for (k in seq_along(row_indices)) {
    i <- row_indices[k]
    js <- if (recycle) seq_len(p) else seq_len(i)
    matrice_pval_asymm_glob[i, js] <- row_results[[k]]$glob
    for (ii in seq_len(nvar + 1L)) {
      matrice_pval_asymm_part[ii, i, js] <- row_results[[k]]$part[ii, ]
    }
    cli::cli_h1(
      "Creating the p-value matrix: end of row {p - i + 1} out of {p}"
    )
  }

  corrected_pval_matrix_glob <- pval_correct(matrice_pval_asymm_glob)
  corrected_pval_glob <- corrected_pval_matrix_glob[1, ]

  corrected_pval_part <- matrix(nrow = nvar + 1L, ncol = p)
  corrected_pval_matrix_part <- array(dim = c(nvar + 1L, p, p))
  for (ii in seq_len(nvar + 1L)) {
    corrected_pval_matrix_part[ii, , ] <- pval_correct(
      matrice_pval_asymm_part[ii, , ]
    )
    corrected_pval_part[ii, ] <- corrected_pval_matrix_part[ii, 1, ]
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

  cli::cli_h1("Interval-Wise Testing completed")

  out <- list(
    call = cl,
    design_matrix = design_matrix,
    unadjusted_pval_F = pval_glob,
    pval_matrix_F = matrice_pval_asymm_glob,
    adjusted_pval_F = corrected_pval_glob,
    unadjusted_pval_part = pval_part,
    pval_matrix_part = matrice_pval_asymm_part,
    adjusted_pval_part = corrected_pval_part,
    data_eval = coeff,
    coeff_regr_eval = coeff_t,
    fitted_eval = fitted_t,
    residuals_eval = residuals_t,
    R2_eval = r2_t
  )
  class(out) <- "flm"
  out
}
