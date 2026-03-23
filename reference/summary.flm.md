# Summarizing Functional Linear Model Fits

`summary` method for class `flm`. Function returning a summary of the
results of IWT for the test on a functional analysis of variance:
minimum IWT-adjusted p-values of the F-tests on the whole model and on
each factor are reported.

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

No value returned. The function
[`summary.fanova()`](https://permaverse.github.io/fdatest/reference/summary.fanova.md)
computes and returns a list of summary statistics of the fitted
functional linear model given in `object`, using the component `call`
from its arguments, plus:

- `ttest`: A \\(L+1) \times 1\\ matrix with columns for the functional
  regression coefficients, and corresponding (two-sided) IWT-adjusted
  minimum p-values of t-tests (i.e., the minimum p-value over all \\p\\
  basis components used to describe functional data).

- `R2`: Range of the functional R-squared.

- `ftest`: IWT-adjusted minimum p-value of functional F-test.

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
for the plot of p-values heatmaps and
[`plot.fanova()`](https://permaverse.github.io/fdatest/reference/plot.fanova.md)
for the plot of analysis of variance results.

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
#> ── Creating the p-value matrix: end of row 2 out of 100 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 3 out of 100 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 4 out of 100 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 5 out of 100 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 6 out of 100 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 7 out of 100 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 8 out of 100 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 9 out of 100 ────────────────────────
#> 
#> ── Creating the p-value matrix: end of row 10 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 11 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 12 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 13 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 14 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 15 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 16 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 17 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 18 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 19 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 20 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 21 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 22 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 23 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 24 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 25 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 26 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 27 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 28 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 29 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 30 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 31 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 32 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 33 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 34 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 35 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 36 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 37 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 38 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 39 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 40 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 41 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 42 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 43 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 44 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 45 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 46 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 47 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 48 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 49 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 50 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 51 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 52 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 53 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 54 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 55 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 56 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 57 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 58 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 59 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 60 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 61 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 62 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 63 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 64 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 65 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 66 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 67 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 68 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 69 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 70 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 71 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 72 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 73 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 74 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 75 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 76 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 77 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 78 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 79 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 80 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 81 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 82 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 83 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 84 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 85 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 86 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 87 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 88 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 89 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 90 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 91 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 92 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 93 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 94 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 95 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 96 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 97 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 98 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 99 out of 100 ───────────────────────
#> 
#> ── Creating the p-value matrix: end of row 100 out of 100 ──────────────────────
#> 
#> ── Interval-Wise Testing completed ─────────────────────────────────────────────

# Summary of the IWT results
summary(IWT_result)
#> $call
#> iwt_lm(formula = formula, dx = dx, n_perm = B, method = method, 
#>     recycle = recycle)
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
