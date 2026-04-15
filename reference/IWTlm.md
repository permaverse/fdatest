# Interval-wise testing procedure for testing functional-on-scalar linear models

The function is used to fit and test functional linear models. It can be
used to carry out regression, and analysis of variance. It implements
the interval-wise testing procedure (IWT) for testing the significance
of the effects of scalar covariates on a functional population.

## Usage

``` r
IWTlm(
  formula,
  dx = NULL,
  B = 1000L,
  method = c("residuals", "responses"),
  recycle = TRUE
)

iwt_lm(
  formula,
  dx = NULL,
  n_perm = 1000L,
  method = c("residuals", "responses"),
  recycle = TRUE
)
```

## Arguments

- formula:

  An object of class
  [`stats::formula`](https://rdrr.io/r/stats/formula.html) (or one that
  can be coerced to that class) specifying the model to be fitted in a
  symbolic fashion. The response (left-hand side) can be either a matrix
  of dimension \\n \times J\\ containing the pointwise evaluations of
  \\n\\ functions on the **same** grid of \\J\\ points, or an object of
  class [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html).

- dx:

  A numeric value specifying the discretization step of the grid used to
  evaluate functional data when it is provided as objects of class
  [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html). Defaults to `NULL`,
  in which case a default value of `0.01` is used which corresponds to a
  grid of size `100L`. Unused if functional data is provided in the form
  of matrices.

- B:

  An integer value specifying the number of permutations used to
  evaluate the p-values of the permutation tests. Defaults to `1000L`.
  Passed as `n_perm` in
  [`iwt_aov()`](https://permaverse.github.io/fdatest/reference/IWTaov.md),
  [`twt_aov()`](https://permaverse.github.io/fdatest/reference/TWTaov.md)
  and
  [`global_aov()`](https://permaverse.github.io/fdatest/reference/global_aov.md).

- method:

  A string specifying the permutation scheme. `"residuals"` permutes
  residuals under the reduced model (Freedman-Lane scheme);
  `"responses"` permutes the responses (Manly scheme). Defaults to
  `"residuals"`.

- recycle:

  A boolean value specifying whether the recycled version of the
  interval-wise testing procedure should be used. See Pini and
  Vantini (2017) for details. Defaults to `TRUE`.

- n_perm:

  An integer value specifying the number of permutations for the
  permutation tests. Defaults to `1000L`.

## Value

An object of class `flm` containing the following components:

- `call`: The matched call.

- `design_matrix`: The design matrix of the functional-on-scalar linear
  model.

- `unadjusted_pval_F`: A numeric vector of length \\J\\ containing the
  unadjusted p-value function of the global F-test evaluated on the
  grid.

- `adjusted_pval_F`: A numeric vector of length \\J\\ containing the
  adjusted p-value function of the global F-test evaluated on the grid.

- `unadjusted_pval_part`: A numeric matrix with one row per model term
  containing the unadjusted p-value functions of the per-term t-tests.

- `adjusted_pval_part`: A numeric matrix with one row per model term
  containing the adjusted p-value functions of the per-term t-tests.

- `data_eval`: A numeric matrix containing the functional response
  evaluated on the grid.

- `coeff_regr_eval`: A numeric matrix containing the functional
  regression coefficients evaluated on the grid.

- `fitted_eval`: A numeric matrix containing the fitted values of the
  functional regression evaluated on the grid.

- `residuals_eval`: A numeric matrix containing the residuals of the
  functional regression evaluated on the grid.

- `R2_eval`: A numeric vector containing the functional R-squared
  evaluated on the grid.

Optionally, the list may contain the following components:

- `pval_matrix_F`: A matrix of dimensions \\p \times p\\ of p-values of
  the interval-wise F-tests. Element \\(i,j)\\ contains the p-value of
  the test on the interval \\(j, j+1, \ldots, j+(p-i))\\. Present only
  if `correction` is `"IWT"`.

- `pval_matrix_part`: An array of dimensions \\(L+1) \times p \times p\\
  of p-values of the per-term interval-wise t-tests. Element \\(l,i,j)\\
  contains the p-value of the joint test on term \\l\\ and interval
  \\(j, j+1, \ldots, j+(p-i))\\. Present only if `correction` is
  `"IWT"`.

- `global_pval_F`: Global p-value of the overall F-test. Present only if
  `correction` is `"Global"`.

- `global_pval_part`: A numeric vector of global p-values of the
  per-term t-tests. Present only if `correction` is `"Global"`.

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

## See also

`iwt_lm()`,
[`twt_lm()`](https://permaverse.github.io/fdatest/reference/TWTlm.md)
and
[`global_lm()`](https://permaverse.github.io/fdatest/reference/Globallm.md)
for calling a specific correction directly.
[`plot.flm()`](https://permaverse.github.io/fdatest/reference/plot.flm.md)
for plotting the results and
[`summary.flm()`](https://permaverse.github.io/fdatest/reference/summary.flm.md)
for summarizing the results.

## Examples

``` r
# Defining the covariates
temperature <- rbind(NASAtemp$milan[, 1:100], NASAtemp$paris[, 1:100])
groups <- c(rep(0, 22), rep(1, 22))

# Performing the IWT
IWT_result <- IWTlm(temperature ~ groups, B = 2L)
#> 
#> ── Point-wise tests ────────────────────────────────────────────────────────────
#> 
#> ── Interval-wise tests ─────────────────────────────────────────────────────────
#> 
#> ── Interval-Wise Testing completed ─────────────────────────────────────────────
# Summary of the IWT results
summary(IWT_result)
#> $call
#> functional_lm_test(formula = formula, correction = "IWT", dx = dx, 
#>     B = n_perm, method = method, recycle = recycle)
#> 
#> $ttest
#>             Minimum p-value    
#> (Intercept)               0 ***
#> groups                    0 ***
#> 
#> $R2
#>               Range of functional R-squared
#> Min R-squared                  0.0003189364
#> Max R-squared                  0.2476892354
#> 
#> $ftest
#>   Minimum p-value    
#> 1               0 ***
#> 

# Plot of the IWT results
plot(
  IWT_result,
  main = "NASA data",
  plot_adjpval = TRUE,
  xlab = "Day",
  xrange = c(1, 365)
)
```
