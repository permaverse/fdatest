# Interval Wise Testing procedure for testing functional analysis of variance

The function implements the Interval Wise Testing procedure for testing
mean differences between several functional populations in a one-way or
multi-way functional analysis of variance framework. Functional data are
tested locally and unadjusted and adjusted p-value functions are
provided. The unadjusted p-value function controls the point-wise error
rate. The adjusted p-value function controls the interval-wise error
rate.

## Usage

``` r
IWTaov(
  formula,
  dx = NULL,
  B = 1000L,
  method = c("residuals", "responses"),
  recycle = TRUE
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

## Value

An object of class `fanova` containing the following components:

- `call`: The matched call.

- `design_matrix`: The design matrix of the functional-on-scalar linear
  model.

- `unadjusted_pval_F`: Evaluation on a grid of the unadjusted p-value
  function of the functional F-test.

- `adjusted_pval_F`: Evaluation on a grid of the adjusted p-value
  function of the functional F-test.

- `unadjusted_pval_factors`: Evaluation on a grid of the unadjusted
  p-value function of the functional F-tests on each factor of the
  analysis of variance (rows).

- `adjusted_pval_factors`: Adjusted p-values of the functional F-tests
  on each factor of the analysis of variance (rows) and each basis
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

- `pval_matrix_factors`: Array of dimensions `c(L+1,p,p)` of the
  p-values of the multivariate F-tests on factors. The element
  \\(l,i,j)\\ of array `pval_matrix` contains the p-value of the joint
  NPC test on factor `l` of the components \\(j,j+1,...,j+(p-i))\\; this
  component is present only if `correction` is set to `"IWT"`.

- `heatmap_matrix_F`: Heatmap matrix of p-values of functional F-test
  (used only for plots); this component is present only if `correction`
  is set to `"IWT"`.

- `heatmap_matrix_factors`: Heatmap matrix of p-values of functional
  F-tests on each factor of the analysis of variance (used only for
  plots); this component is present only if `correction` is set to
  `"IWT"`.

- `Global_pval_F`: Global p-value of the overall test F; this component
  is present only if `correction` is set to `"Global"`.

- `Global_pval_factors`: Global p-value of test F involving each factor
  separately; this component is present only if `correction` is set to
  `"Global"`.

## References

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

D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of
Reported Significance Levels. *Journal of Business & Economic
Statistics* 1.4, 292-298.

B. F. J. Manly (2006). Randomization, *Bootstrap and Monte Carlo Methods
in Biology*. Vol. 70. CRC Press.

## See also

See also
[`plot.fanova()`](https://permaverse.github.io/fdatest/reference/plot.fanova.md)
for plotting the results and
[`summary.fanova()`](https://permaverse.github.io/fdatest/reference/summary.fanova.md)
for summarizing the results of the functional analysis of variance.

## Examples

``` r
temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
groups <- c(rep(0, 22), rep(1, 22))

# Performing the IWT
IWT_result <- IWTaov(temperature ~ groups, B = 10L)
#> Error in eval(predvars, data, env): object 'groups' not found

# Summary of the IWT results
summary(IWT_result)
#> Error: object 'IWT_result' not found

# Plot of the IWT results
graphics::layout(1)
plot(IWT_result)
#> Error: object 'IWT_result' not found

# All graphics on the same device
graphics::layout(matrix(1:4, nrow = 2, byrow = FALSE))
plot(
  IWT_result,
  main = 'NASA data',
  plot.adjpval = TRUE,
  xlab = 'Day',
  xrange = c(1, 365)
)
#> Error: object 'IWT_result' not found
```
