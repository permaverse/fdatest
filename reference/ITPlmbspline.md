# Interval Testing Procedure for testing Functional-on-Scalar Linear Models with B-spline basis

The function is used to fit and test functional linear models. It can be
used to carry out regression, and analysis of variance. It implements
the Interval Testing Procedure for testing the significance of the
effects of scalar covariates on a functional population evaluated on a
uniform grid. Data are represented by means of the B-spline basis and
the significance of each basis coefficient is tested with an
interval-wise control of the Family Wise Error Rate. The default
parameters of the basis expansion lead to the piece-wise interpolating
function.

## Usage

``` r
ITPlmbspline(
  formula,
  order = 2,
  nknots = dim(stats::model.response(stats::model.frame(formula)))[2],
  B = 10000,
  method = "residuals"
)
```

## Arguments

- formula:

  An object of class "[`formula`](https://rdrr.io/r/stats/formula.html)"
  (or one that can be coerced to that class): a symbolic description of
  the model to be fitted.

- order:

  Order of the B-spline basis expansion. The default is `order=2`.

- nknots:

  Number of knots of the B-spline basis expansion. The default is
  `nknots=dim(data1)[2]`.

- B:

  The number of iterations of the MC algorithm to evaluate the p-values
  of the permutation tests. The defualt is `B=1000`.

- method:

  Permutation method used to calculate the p-value of permutation tests.
  Choose "`residuals`" for the permutations of residuals under the
  reduced model, according to the Freedman and Lane scheme, and
  "`responses`" for the permutation of the responses, according to the
  Manly scheme.

## Value

`ITPlmbspline` returns an object of
[`class`](https://rdrr.io/r/base/class.html) "`ITPlm`". The function
`summary` is used to obtain and print a summary of the results. An
object of class "`ITPlm`" is a list containing at least the following
components:

- call:

  The matched call.

- design_matrix:

  The design matrix of the functional-on-scalar linear model.

- basis:

  String vector indicating the basis used for the first phase of the
  algorithm. In this case equal to `"B-spline"`.

- coeff:

  Matrix of dimensions `c(n,p)` of the `p` coefficients of the B-spline
  basis expansion. Rows are associated to units and columns to the basis
  index.

- coeff_regr:

  Matrix of dimensions `c(L+1,p)` of the `p` coefficients of the
  B-spline basis expansion of the intercept (first row) and the `L`
  effects of the covariates specified in `formula`. Columns are
  associated to the basis index.

- pval_F:

  Unadjusted p-values of the functional F-test for each basis
  coefficient.

- pval_matrix_F:

  Matrix of dimensions `c(p,p)` of the p-values of the multivariate
  F-tests. The element `(i,j)` of matrix `pval_matrix` contains the
  p-value of the joint NPC test of the components `(j,j+1,...,j+(p-i))`.

- adjusted_pval_F:

  Adjusted p-values of the functional F-test for each basis coefficient.

- pval_t:

  Unadjusted p-values of the functional t-tests for each partial
  regression coefficient including the intercept (rows) and each basis
  coefficient (columns).

- pval_matrix_t:

  Array of dimensions `c(L+1,p,p)` of the p-values of the multivariate
  t-tests. The element `(l,i,j)` of array `pval_matrix` contains the
  p-value of the joint NPC test on covariate `l` of the components
  `(j,j+1,...,j+(p-i))`.

- adjusted_pval_t:

  adjusted p-values of the functional t-tests for each partial
  regression coefficient including the intercept (rows) and each basis
  coefficient (columns).

- data_eval:

  Evaluation on a fine uniform grid of the functional data obtained
  through the basis expansion.

- coeff_regr_eval:

  Evaluation on a fine uniform grid of the functional regression
  coefficients.

- fitted_eval:

  Evaluation on a fine uniform grid of the fitted values of the
  functional regression.

- residuals_eval:

  Evaluation on a fine uniform grid of the residuals of the functional
  regression.

- R2_eval:

  Evaluation on a fine uniform grid of the functional R-squared of the
  regression.

- heatmap_matrix_F:

  Heatmap matrix of p-values of functional F-test (used only for plots).

- heatmap_matrix_t:

  Heatmap matrix of p-values of functional t-tests (used only for
  plots).

## References

A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference
for Functional Data Controlling the Family Wise Error Rate on Intervals.
Biometrics 73(3): 835–845.

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
Statistics* 1(4), 292-298.

B. F. J. Manly (2006). Randomization, *Bootstrap and Monte Carlo Methods
in Biology*. Vol. 70. CRC Press.

## Examples

``` r
# Defining the covariates
temperature <- rbind(NASAtemp$milan,NASAtemp$paris)
groups <- c(rep(0,22),rep(1,22))

# Performing the ITP
#ITP_result <- ITPlmbspline(temperature ~ groups,B=100,nknots=20)
# Summary of the ITP results
#summary(ITP_result)

# Plot of the ITP results
#plot(ITP_result,main='NASA data', plot_adjpval = TRUE,xlab='Day',xrange=c(1,365))

#plot(ITP_result,main='NASA data', plot_adjpval = TRUE,xlab='Day',xrange=c(1,365))
```
