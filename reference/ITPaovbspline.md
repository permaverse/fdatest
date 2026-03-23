# Interval Testing Procedure for testing unctional analysis of variance with B-spline basis

The function implements the Interval Testing Procedure for testing for
significant differences between several functional population evaluated
on a uniform grid, in a functional analysis of variance setting. Data
are represented by means of the B-spline basis and the significance of
each basis coefficient is tested with an interval-wise control of the
Family Wise Error Rate. The default parameters of the basis expansion
lead to the piece-wise interpolating function.

## Usage

``` r
ITPaovbspline(
  formula,
  order = 2,
  nknots = dim(stats::model.response(stats::model.frame(formula)))[2],
  B = 1000,
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

`ITPaovbspline` returns an object of
[`class`](https://rdrr.io/r/base/class.html) "`ITPaov`". The function
`summary` is used to obtain and print a summary of the results. An
object of class "`ITPaov`" is a list containing at least the following
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

- pval_factors:

  Unadjusted p-values of the functional F-tests on each factor of the
  analysis of variance, separately (rows) and each basis coefficient
  (columns).

- pval_matrix_factors:

  Array of dimensions `c(L+1,p,p)` of the p-values of the multivariate
  F-tests on factors. The element `(l,i,j)` of array `pval_matrix`
  contains the p-value of the joint NPC test on factor `l` of the
  components `(j,j+1,...,j+(p-i))`.

- adjusted_pval_factors:

  adjusted p-values of the functional F-tests on each factor of the
  analysis of variance (rows) and each basis coefficient (columns).

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

- heatmap_matrix_factors:

  Heatmap matrix of p-values of functional F-tests on each factor of the
  analysis of variance (used only for plots).

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
Statistics* 1.4, 292-298.

B. F. J. Manly (2006). Randomization, *Bootstrap and Monte Carlo Methods
in Biology*. Vol. 70. CRC Press.

## Examples

``` r
temperature <- rbind(NASAtemp$milan,NASAtemp$paris)
groups <- c(rep(0,22),rep(1,22))

# Performing the ITP
ITP_result <- ITPaovbspline(temperature ~ groups,B=5,nknots=20,order=3)
#> Warning: `ITPaovbspline()` was deprecated in fdatest 2.2.0.
#> ℹ Please use `IWTaov()` instead.
#> 
#> ── Point-wise tests ────────────────────────────────────────────────────────────
#> 
#> ── Interval-wise tests ─────────────────────────────────────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 2 out of 365 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 3 out of 365 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 4 out of 365 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 5 out of 365 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 6 out of 365 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 7 out of 365 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 8 out of 365 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 9 out of 365 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 10 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 11 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 12 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 13 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 14 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 15 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 16 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 17 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 18 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 19 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 20 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 21 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 22 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 23 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 24 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 25 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 26 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 27 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 28 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 29 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 30 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 31 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 32 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 33 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 34 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 35 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 36 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 37 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 38 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 39 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 40 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 41 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 42 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 43 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 44 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 45 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 46 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 47 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 48 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 49 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 50 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 51 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 52 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 53 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 54 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 55 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 56 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 57 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 58 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 59 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 60 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 61 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 62 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 63 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 64 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 65 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 66 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 67 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 68 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 69 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 70 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 71 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 72 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 73 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 74 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 75 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 76 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 77 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 78 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 79 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 80 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 81 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 82 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 83 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 84 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 85 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 86 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 87 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 88 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 89 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 90 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 91 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 92 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 93 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 94 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 95 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 96 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 97 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 98 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 99 out of 365 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 100 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 101 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 102 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 103 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 104 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 105 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 106 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 107 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 108 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 109 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 110 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 111 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 112 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 113 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 114 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 115 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 116 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 117 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 118 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 119 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 120 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 121 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 122 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 123 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 124 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 125 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 126 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 127 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 128 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 129 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 130 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 131 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 132 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 133 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 134 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 135 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 136 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 137 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 138 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 139 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 140 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 141 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 142 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 143 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 144 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 145 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 146 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 147 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 148 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 149 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 150 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 151 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 152 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 153 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 154 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 155 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 156 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 157 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 158 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 159 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 160 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 161 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 162 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 163 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 164 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 165 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 166 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 167 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 168 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 169 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 170 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 171 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 172 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 173 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 174 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 175 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 176 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 177 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 178 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 179 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 180 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 181 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 182 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 183 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 184 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 185 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 186 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 187 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 188 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 189 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 190 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 191 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 192 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 193 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 194 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 195 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 196 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 197 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 198 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 199 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 200 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 201 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 202 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 203 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 204 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 205 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 206 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 207 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 208 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 209 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 210 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 211 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 212 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 213 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 214 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 215 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 216 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 217 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 218 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 219 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 220 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 221 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 222 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 223 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 224 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 225 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 226 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 227 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 228 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 229 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 230 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 231 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 232 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 233 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 234 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 235 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 236 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 237 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 238 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 239 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 240 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 241 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 242 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 243 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 244 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 245 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 246 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 247 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 248 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 249 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 250 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 251 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 252 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 253 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 254 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 255 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 256 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 257 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 258 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 259 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 260 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 261 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 262 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 263 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 264 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 265 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 266 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 267 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 268 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 269 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 270 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 271 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 272 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 273 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 274 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 275 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 276 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 277 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 278 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 279 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 280 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 281 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 282 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 283 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 284 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 285 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 286 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 287 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 288 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 289 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 290 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 291 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 292 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 293 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 294 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 295 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 296 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 297 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 298 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 299 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 300 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 301 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 302 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 303 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 304 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 305 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 306 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 307 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 308 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 309 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 310 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 311 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 312 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 313 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 314 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 315 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 316 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 317 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 318 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 319 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 320 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 321 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 322 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 323 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 324 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 325 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 326 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 327 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 328 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 329 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 330 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 331 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 332 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 333 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 334 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 335 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 336 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 337 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 338 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 339 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 340 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 341 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 342 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 343 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 344 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 345 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 346 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 347 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 348 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 349 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 350 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 351 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 352 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 353 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 354 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 355 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 356 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 357 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 358 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 359 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 360 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 361 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 362 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 363 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 364 out of 365 ──────────────────────
#> 
#> ── Creating the p-value matrix: end of row 365 out of 365 ──────────────────────
#> 
#> ── Interval-Wise Testing completed ─────────────────────────────────────────────

# Summary of the ITP results
summary(ITP_result)
#> $call
#> iwt_aov(formula = formula, dx = dx, n_perm = B, method = method, 
#>     recycle = recycle)
#> 
#> $factors
#>        Minimum p-value    
#> groups               0 ***
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

# Plot of the ITP results
plot(
  ITP_result,
  main = "NASA data",
  plot_adjpval = TRUE,
  xlab = 'Day',
  xrange = c(1, 365)
)

```
