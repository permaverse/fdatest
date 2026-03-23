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
#'   - `data.eval`: Evaluation on a fine uniform grid of the functional data
#'   obtained through the basis expansion.
#'   - `coeff.regr.eval`: Evaluation on a fine uniform grid of the functional
#'   regression coefficients.
#'   - `fitted.eval`: Evaluation on a fine uniform grid of the fitted values of
#'   the functional regression.
#'   - `residuals.eval`: Evaluation on a fine uniform grid of the residuals of
#'   the functional regression.
#'   - `R2.eval`: Evaluation on a fine uniform grid of the functional R-squared
#'   of the regression.
#'
#'   Optionally, the list may contain the following components:
#'
#'   - `pval_matrix_F`: Matrix of dimensions \code{c(p,p)} of the p-values of
#'   the intervalwise F-tests. The element \eqn{(i,j)} of matrix `pval.matrix`
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
#' TWT.result <- TWTlm(temperature ~ groups, B = 100L)
#' # Summary of the TWT results
#' summary(TWT.result)
functional_lm_test <- function(
  formula,
  correction,
  dx = NULL,
  B = 1000L,
  method = c("residuals", "responses"),
  recycle = TRUE,
  stat = c("Integral", "Max")
) {
  out <- switch(
    correction,
    IWT = IWTaov(
      formula = formula,
      dx = dx,
      B = B,
      method = method,
      recycle = recycle
    ),
    TWT = TWTaov(
      formula = formula,
      dx = dx,
      B = B,
      method = method
    ),
    Global = Globalaov(
      formula = formula,
      dx = dx,
      B = B,
      method = method,
      stat = stat
    )
  )

  out$correction <- correction
  out
}
