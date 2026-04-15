# Summarizing Functional Analysis of Variance Fits

`summary` method for class `faov`. Returns a summary of the results of
the local testing procedure for functional analysis of variance: the
minimum adjusted p-values of the F-test on the whole model and on each
factor are reported.

## Usage

``` r
# S3 method for class 'faov'
summary(object, ...)
```

## Arguments

- object:

  An object of class `faov`, usually, a result of a call to
  [`functional_anova_test()`](https://permaverse.github.io/fdatest/reference/functional_anova_test.md).

- ...:

  Further arguments passed to or from other methods.

## Value

A list of summary statistics of the fitted functional analysis of
variance given in `object`, using the component `call` from its
arguments, plus:

- `factors`: A \\L \times 2\\ data frame with one row per factor
  reporting the minimum adjusted p-value of the corresponding F-test and
  a significance code.

- `R2`: A \\2 \times 1\\ matrix giving the range of the functional
  R-squared.

- `ftest`: A \\1 \times 2\\ data frame reporting the minimum adjusted
  p-value of the global F-test and a significance code.

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
[`plot.faov()`](https://permaverse.github.io/fdatest/reference/plot.faov.md)
for the plot of analysis of variance results.

## Examples

``` r
temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
groups <- c(rep(0, 22), rep(1, 22))

# Performing the IWT
IWT_result <- functional_anova_test(
  temperature ~ groups,
  B = 10L,
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
#> functional_anova_test(formula = temperature ~ groups, correction = "IWT", 
#>     B = 10L)
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
```
