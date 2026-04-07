#' Local testing procedures for functional-on-scalar linear models
#'
#' @description The function implements local testing procedures for testing
#'   mean differences between multiple functional populations. Functional data are
#'   tested locally and unadjusted and adjusted p-value functions are provided.
#'   The unadjusted p-value function controls the point-wise error rate. The
#'   adjusted p-value function can be computed according to the following
#'   methods:
#'
#'   - global testing (controlling the FWER weakly)
#'   - interval-wise testing (controlling the interval-wise error rate)
#'   - threshold-wise testing (controlling the FWER asymptotically)
#'   - partition closed testing (controlling the FWER on a partition)
#'   - functional Benjamini Hochberg (controlling the FDR)
#'
#' @inheritParams functional_anova_test
#'
#' @returns An object of class `flm` containing the following components:
#'
#'   - `call`: The matched call.
#'   - `design_matrix`: The design matrix of the functional-on-scalar linear
#'   model.
#'   - `unadjusted_pval_F`: Evaluation on a grid of the unadjusted p-value
#'   function of the functional F-test.
#'   - `adjusted_pval_F`: Evaluation on a grid of the adjusted p-value function
#'   of the functional F-test.
#'   - `unadjusted_pval_part`: Evaluation on a grid of the unadjusted p-value
#'   function of the functional F-tests on each factor of the analysis of
#'   variance (rows).
#'   - `adjusted_pval_part`: Adjusted p-values of the functional F-tests on
#'   each factor of the analysis of variance (rows) and each basis coefficient
#'   (columns).
#'   - `data_eval`: Evaluation on a fine uniform grid of the functional data
#'   obtained through the basis expansion.
#'   - `coeff_regr_eval`: Evaluation on a fine uniform grid of the functional
#'   regression coefficients.
#'   - `fitted_eval`: Evaluation on a fine uniform grid of the fitted values of
#'   the functional regression.
#'   - `residuals_eval`: Evaluation on a fine uniform grid of the residuals of
#'   the functional regression.
#'   - `R2_eval`: Evaluation on a fine uniform grid of the functional R-squared
#'   of the regression.
#'
#'   Optionally, the list may contain the following components:
#'
#'   - `pval_matrix_F`: Matrix of dimensions \code{c(p,p)} of the p-values of
#'   the intervalwise F-tests. The element \eqn{(i,j)} of matrix `pval_matrix`
#'   contains the p-value of the test of interval indexed by
#'   \eqn{(j,j+1,...,j+(p-i))}; this component is present only if `correction`
#'   is set to `"IWT"`.
#'   - `pval_matrix_part`: Array of dimensions `c(L+1,p,p)` of the p-values
#'   of the multivariate F-tests on factors. The element \eqn{(l,i,j)} of array
#'   `pval_matrix_part` contains the p-value of the joint NPC test on factor `l` of
#'   the components \eqn{(j,j+1,...,j+(p-i))}; this component is present only if
#'   `correction` is set to `"IWT"`.
#'   - `Global_pval_F`: Global p-value of the overall test F; this component is
#'   present only if `correction` is set to `"Global"`.
#'   - `Global_pval_part`: Global p-value of test F involving each factor
#'   separately; this component is present only if `correction` is set to `"Global"`.
#'
#' @seealso See also [`plot.flm()`] for plotting the results and
#'  [`summary.flm()`] for summarizing the results of the functional analysis
#'  of variance.
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
#' TWT_result <- TWTlm(temperature ~ groups, B = 100L)
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
    out$Global_pval_F <- adj_glob$adjusted_pvalues[1]
    out$Global_pval_part <- vapply(
      adj_part_list,
      function(x) x$adjusted_pvalues[1],
      numeric(1L)
    )
  }

  class(out) <- "flm"
  out
}
