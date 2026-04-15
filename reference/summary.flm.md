# Summarizing Functional Linear Model Fits

`summary` method for class `flm`. Returns a summary of the results of
the local testing procedure for a functional-on-scalar linear model: the
minimum adjusted p-values of the F-test on the whole model and the
t-tests on each covariate are reported.

## Usage

``` r
# S3 method for class 'flm'
summary(object, ...)
```

## Arguments

- object:

  An object of class `flm`, usually, a result of a call to
  [`functional_lm_test()`](https://permaverse.github.io/fdatest/reference/functional_lm_test.md).

- ...:

  Further arguments passed to or from other methods.

## Value

A list of summary statistics of the fitted functional linear model given
in `object`, using the component `call` from its arguments, plus:

- `ttest`: A \\(L+1) \times 2\\ data frame with one row per model term
  (intercept plus each predictor) reporting the minimum adjusted p-value
  of the corresponding t-test and a significance code.

- `R2`: A \\2 \times 1\\ matrix giving the range of the functional
  R-squared.

- `ftest`: A \\1 \times 2\\ data frame reporting the minimum adjusted
  p-value of the functional F-test and a significance code.

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

## See also

[`IWTimage()`](https://permaverse.github.io/fdatest/reference/IWTimage.md)
for the plot of p-value heatmaps and
[`plot.flm()`](https://permaverse.github.io/fdatest/reference/plot.flm.md)
for the plot of functional linear model results.

## Examples

``` r
temperature <- rbind(NASAtemp$milan[, 1:100], NASAtemp$paris[, 1:100])
groups <- c(rep(0, 22), rep(1, 22))

# Performing the IWT
IWT_result <- functional_lm_test(
  temperature ~ groups,
  B = 2L,
  correction = "IWT"
)
#> 
#> ── Point-wise tests ────────────────────────────────────────────────────────────
#> 
#> ── Interval-wise tests ─────────────────────────────────────────────────────────
#> 
#> ── Interval-Wise Testing completed ─────────────────────────────────────────────

# Summary of the IWT results
summary(IWT_result)
#> $call
#> functional_lm_test(formula = temperature ~ groups, correction = "IWT", 
#>     B = 2L)
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
```
