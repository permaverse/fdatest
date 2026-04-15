#' Local testing procedures for functional-on-scalar linear models
#'
#' @description Implements local testing procedures for testing the
#'   significance of the effects of scalar covariates on a functional response
#'   in a functional-on-scalar linear model. Functional data are tested locally
#'   and unadjusted and adjusted p-value functions are provided. The unadjusted
#'   p-value function controls the point-wise error rate. The adjusted p-value
#'   function can be computed according to the following methods:
#'
#'   - global testing (controlling the FWER weakly)
#'   - interval-wise testing (controlling the interval-wise error rate)
#'   - threshold-wise testing (controlling the FWER asymptotically)
#'
#' @inheritParams functional_anova_test
#'
#' @returns An object of class `flm` containing the following components:
#'
#'   - `call`: The matched call.
#'   - `design_matrix`: The design matrix of the functional-on-scalar linear
#'   model.
#'   - `unadjusted_pval_F`: A numeric vector of length \eqn{J} containing the
#'   unadjusted p-value function of the global F-test evaluated on the grid.
#'   - `adjusted_pval_F`: A numeric vector of length \eqn{J} containing the
#'   adjusted p-value function of the global F-test evaluated on the grid.
#'   - `unadjusted_pval_part`: A numeric matrix with one row per model term
#'   containing the unadjusted p-value functions of the per-term t-tests.
#'   - `adjusted_pval_part`: A numeric matrix with one row per model term
#'   containing the adjusted p-value functions of the per-term t-tests.
#'   - `data_eval`: A numeric matrix containing the functional response
#'   evaluated on the grid.
#'   - `coeff_regr_eval`: A numeric matrix containing the functional regression
#'   coefficients evaluated on the grid.
#'   - `fitted_eval`: A numeric matrix containing the fitted values of the
#'   functional regression evaluated on the grid.
#'   - `residuals_eval`: A numeric matrix containing the residuals of the
#'   functional regression evaluated on the grid.
#'   - `R2_eval`: A numeric vector containing the functional R-squared evaluated
#'   on the grid.
#'
#'   Optionally, the list may contain the following components:
#'
#'   - `pval_matrix_F`: A matrix of dimensions \eqn{p \times p} of p-values of
#'   the interval-wise F-tests. Element \eqn{(i,j)} contains the p-value of the
#'   test on the interval \eqn{(j, j+1, \ldots, j+(p-i))}. Present only if
#'   `correction` is `"IWT"`.
#'   - `pval_matrix_part`: An array of dimensions \eqn{(L+1) \times p \times p}
#'   of p-values of the per-term interval-wise t-tests. Element \eqn{(l,i,j)}
#'   contains the p-value of the joint test on term \eqn{l} and interval
#'   \eqn{(j, j+1, \ldots, j+(p-i))}. Present only if `correction` is `"IWT"`.
#'   - `global_pval_F`: Global p-value of the overall F-test. Present only if
#'   `correction` is `"Global"`.
#'   - `global_pval_part`: A numeric vector of global p-values of the per-term
#'   t-tests. Present only if `correction` is `"Global"`.
#'
#' @seealso [`iwt_lm()`], [`twt_lm()`] and [`global_lm()`] for calling a
#'   specific correction directly. [`plot.flm()`] for plotting the results and
#'   [`summary.flm()`] for summarizing the results.
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
#' @export
#' @examples
#' # Defining the covariates
#' temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
#' groups <- c(rep(0, 22), rep(1, 22))
#'
#' # Performing the TWT
#' TWT_result <- functional_lm_test(
#'   temperature ~ groups,
#'   correction = "TWT",
#'   B = 10L
#' )
#' # Summary of the TWT results
#' summary(TWT_result)
functional_lm_test <- function(
  formula,
  correction = c("IWT", "TWT", "Global"),
  dx = NULL,
  B = 1000L, # nolint: object_name_linter.
  method = c("residuals", "responses"),
  recycle = TRUE,
  stat = c("Integral", "Max")
) {
  correction <- rlang::arg_match(correction)
  method <- rlang::arg_match(method)
  stat <- rlang::arg_match(stat)
  cl <- match.call()

  aggregation_strategy <- switch(stat, Integral = "integral", Max = "max")

  cli::cli_h1("Point-wise tests")
  res <- lm_prepare_data(formula, dx, B, method)

  coeff <- res$coeff
  n <- res$n
  p <- res$p
  nvar <- res$nvar
  var_names <- res$var_names
  design_matrix <- res$design_matrix
  regr0 <- res$regr0
  t0_glob <- res$t0_glob
  t_glob <- res$t_glob
  pval_glob <- res$pval_glob
  t0_part <- res$t0_part
  t_part <- res$t_part
  pval_part <- res$pval_part

  cli::cli_h1(switch(
    correction,
    IWT = "Interval-wise tests",
    TWT = "Threshold-wise tests",
    Global = "Global test"
  ))

  # Adjust global F-test p-values
  adj_glob <- switch(
    correction,
    IWT = p_adjust_iwt(
      p = p,
      pval = pval_glob,
      t0 = t0_glob,
      t_coeff = t_glob,
      n_perm = B,
      recycle = recycle,
      aggregation_strategy = aggregation_strategy
    ),
    TWT = p_adjust_twt(
      pval = pval_glob,
      p = p,
      t0 = t0_glob,
      t_coeff = t_glob,
      aggregation_strategy = aggregation_strategy
    ),
    Global = p_adjust_global(
      aggregation_strategy = aggregation_strategy,
      t0 = t0_glob,
      t_coeff = t_glob,
      p = p
    )
  )

  # Adjust per-variable p-values — one call to p_adjust_xx() per variable
  # (nvar + 1: intercept plus each predictor)
  adj_part_list <- lapply(seq_len(nvar + 1L), function(ii) {
    switch(
      correction,
      IWT = p_adjust_iwt(
        p = p,
        pval = pval_part[ii, ],
        t0 = t0_part[ii, ],
        t_coeff = t_part[, ii, ],
        n_perm = B,
        recycle = recycle,
        aggregation_strategy = aggregation_strategy
      ),
      TWT = p_adjust_twt(
        pval = pval_part[ii, ],
        p = p,
        t0 = t0_part[ii, ],
        t_coeff = t_part[, ii, ],
        aggregation_strategy = aggregation_strategy
      ),
      Global = p_adjust_global(
        aggregation_strategy = aggregation_strategy,
        t0 = t0_part[ii, ],
        t_coeff = t_part[, ii, ],
        p = p
      )
    )
  })

  adjusted_pval_part <- matrix(nrow = nvar + 1L, ncol = p)
  for (ii in seq_len(nvar + 1L)) {
    adjusted_pval_part[ii, ] <- adj_part_list[[ii]]$adjusted_pvalues
  }
  rownames(adjusted_pval_part) <- var_names
  rownames(pval_part) <- var_names

  coeff_t <- regr0$coeff
  fitted_t <- regr0$fitted
  rownames(coeff_t) <- var_names

  residuals_t <- coeff - fitted_t
  ybar_t <- colMeans(coeff)
  r2_t <- colSums(
    (fitted_t - matrix(ybar_t, nrow = n, ncol = p, byrow = TRUE))^2
  ) /
    colSums(
      (coeff - matrix(ybar_t, nrow = n, ncol = p, byrow = TRUE))^2
    )

  cli::cli_h1(switch(
    correction,
    IWT = "Interval-Wise Testing completed",
    TWT = "Threshold-Wise Testing completed",
    Global = "Global Testing completed"
  ))

  out <- list(
    call = cl,
    design_matrix = design_matrix,
    unadjusted_pval_F = pval_glob,
    adjusted_pval_F = adj_glob$adjusted_pvalues,
    unadjusted_pval_part = pval_part,
    adjusted_pval_part = adjusted_pval_part,
    data_eval = coeff,
    coeff_regr_eval = coeff_t,
    fitted_eval = fitted_t,
    residuals_eval = residuals_t,
    R2_eval = r2_t,
    correction = correction
  )

  if (correction == "IWT") {
    out$pval_matrix_F <- adj_glob$pvalue_matrix
    pval_matrix_part <- array(dim = c(nvar + 1L, p, p))
    for (ii in seq_len(nvar + 1L)) {
      pval_matrix_part[ii, , ] <- adj_part_list[[ii]]$pvalue_matrix
    }
    out$pval_matrix_part <- pval_matrix_part
  }

  if (correction == "Global") {
    out$global_pval_F <- adj_glob$adjusted_pvalues[1]
    out$global_pval_part <- vapply(
      adj_part_list,
      function(x) x$adjusted_pvalues[1],
      numeric(1L)
    )
  }

  class(out) <- "flm"
  out
}
