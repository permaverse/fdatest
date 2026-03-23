# Local testing procedures for functional-on-scalar linear models

The function implements local testing procedures for testing mean
differences between multiple functional populations. Functional data are
tested locally and unadjusted and adjusted p-value functions are
provided. The unadjusted p-value function controls the point-wise error
rate. The adjusted p-value function can be computed according to the
following methods:

- global testing (controlling the FWER weakly)

- interval-wise testing (controlling the interval-wise error rate)

- threshold-wise testing (controlling the FWER asymptotically)

- partition closed testing (controlling the FWER on a partition)

- functional Benjamini Hochberg (controlling the FDR)

## Usage

``` r
functional_lm_test(
  formula,
  correction,
  dx = NULL,
  B = 1000L,
  method = c("residuals", "responses"),
  recycle = TRUE,
  stat = c("Integral", "Max")
)
```

## Arguments

- formula:

  An object of class
  [`stats::formula`](https://rdrr.io/r/stats/formula.html) (or one that
  can be coerced to that class) specifying the model to be fitted in a
  symbolic fashion. The output variable (left-hand side) of the formula
  can be either a matrix of dimension \\n \times J\\ containing the
  pointwise evaluations of \\n\\ functions on the **same** grid of \\J\\
  points, or an object of class
  [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html).

- correction:

  A string specifying the method used to calculate the adjusted p-value
  function. Choices are either `"Global"` for global testing, `"IWT"`
  for interval-wise testing or `"TWT"` for threshold-wise testing.

- dx:

  A numeric value specifying the discretization step of the grid used to
  evaluate functional data when it is provided as objects of class
  [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html). Defaults to `NULL`,
  in which case a default value of `0.01` is used which corresponds to a
  grid of size `100L`. Unused if functional data is provided in the form
  of matrices.

- B:

  An integer value specifying the number of iterations of the MC
  algorithm to evaluate the p-value of the permutation tests. Defaults
  to `1000L`.

- method:

  A string specifying the method used to calculate the p-value of
  permutation tests. Choices are either `"residuals"` which performs
  permutation of residuals under the reduced model according to the
  Freedman and Lane scheme or `"responses"`, which performs permutation
  of the responses, according to the Manly scheme. Defaults to
  `"residuals"`.

- recycle:

  A boolean value specifying whether the recycled version of the
  interval-wise testing procedure should be used. See Pini and
  Vantini (2017) for details. Defaults to `TRUE`.

- stat:

  A string specifying the test statistic used for the global test.
  Choices are either `"Integral"`, in which case the statistic is
  defined as the integral of the F-test statistic over the domain, or
  `"Max"`, in which case the statistic is defined as the maximum of the
  F-test statistic over the domain. Defaults to `"Integral"`.

## Value

An object of class `flm` containing the following components:

- `call`: The matched call.

- `design_matrix`: The design matrix of the functional-on-scalar linear
  model.

- `unadjusted_pval_F`: Evaluation on a grid of the unadjusted p-value
  function of the functional F-test.

- `adjusted_pval_F`: Evaluation on a grid of the adjusted p-value
  function of the functional F-test.

- `unadjusted_pval_part`: Evaluation on a grid of the unadjusted p-value
  function of the functional F-tests on each factor of the analysis of
  variance (rows).

- `adjusted_pval_part`: Adjusted p-values of the functional F-tests on
  each factor of the analysis of variance (rows) and each basis
  coefficient (columns).

- `data_eval`: Evaluation on a fine uniform grid of the functional data
  obtained through the basis expansion.

- `coeff_regr_eval`: Evaluation on a fine uniform grid of the functional
  regression coefficients.

- `fitted_eval`: Evaluation on a fine uniform grid of the fitted values
  of the functional regression.

- `residuals_eval`: Evaluation on a fine uniform grid of the residuals
  of the functional regression.

- `R2_eval`: Evaluation on a fine uniform grid of the functional
  R-squared of the regression.

Optionally, the list may contain the following components:

- `pval_matrix_F`: Matrix of dimensions `c(p,p)` of the p-values of the
  intervalwise F-tests. The element \\(i,j)\\ of matrix `pval_matrix`
  contains the p-value of the test of interval indexed by
  \\(j,j+1,...,j+(p-i))\\; this component is present only if
  `correction` is set to `"IWT"`.

- `pval_matrix_part`: Array of dimensions `c(L+1,p,p)` of the p-values
  of the multivariate F-tests on factors. The element \\(l,i,j)\\ of
  array `pval_matrix_part` contains the p-value of the joint NPC test on
  factor `l` of the components \\(j,j+1,...,j+(p-i))\\; this component
  is present only if `correction` is set to `"IWT"`.

- `Global_pval_F`: Global p-value of the overall test F; this component
  is present only if `correction` is set to `"Global"`.

- `Global_pval_part`: Global p-value of test F involving each factor
  separately; this component is present only if `correction` is set to
  `"Global"`.

## References

Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
“Domain selection and familywise error rate for functional data: A
unified framework. *Biometrics* 79(2), 1119-1132.

D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of
Reported Significance Levels. *Journal of Business & Economic
Statistics* 1.4, 292-298.

B. F. J. Manly (2006). Randomization, *Bootstrap and Monte Carlo Methods
in Biology*. Vol. 70. CRC Press.

Pini, A., & Vantini, S. (2017). Interval-wise testing for functional
data. *Journal of Nonparametric Statistics*, 29(2), 407-424.

Pini, A., Vantini, S., Colosimo, B. M., & Grasso, M. (2018).
Domain‐selective functional analysis of variance for supervised
statistical profile monitoring of signal data. *Journal of the Royal
Statistical Society: Series C (Applied Statistics)* 67(1), 55-81.

Abramowicz, K., Hager, C. K., Pini, A., Schelin, L., Sjostedt de Luna,
S., & Vantini, S. (2018). Nonparametric inference for
functional‐on‐scalar linear models applied to knee kinematic hop data
after injury of the anterior cruciate ligament. *Scandinavian Journal of
Statistics* 45(4), 1036-1061.

## See also

See also
[`plot.flm()`](https://permaverse.github.io/fdatest/reference/plot.flm.md)
for plotting the results and
[`summary.flm()`](https://permaverse.github.io/fdatest/reference/summary.flm.md)
for summarizing the results of the functional analysis of variance.

## Examples

``` r
# Defining the covariates
temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
groups <- c(rep(0, 22), rep(1, 22))

# Performing the TWT
TWT_result <- TWTlm(temperature ~ groups, B = 100L)
#> 
#> ── Point-wise tests ────────────────────────────────────────────────────────────
#> 
#> ── Threshold-wise tests ────────────────────────────────────────────────────────
#> 
#> ── Threshold-Wise Testing completed ────────────────────────────────────────────
# Summary of the TWT results
summary(TWT_result)
#> $call
#> twt_lm(formula = formula, dx = dx, n_perm = B, method = method)
#> 
#> $ttest
#>             Minimum p-value    
#> (Intercept)               0 ***
#> groups                    0 ***
#> 
#> $R2
#>               Range of functional R-squared
#> Min R-squared                  3.390203e-05
#> Max R-squared                  5.399620e-01
#> 
#> $ftest
#>   Minimum p-value    
#> 1               0 ***
#> 
```
