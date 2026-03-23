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
temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
groups <- c(rep(0, 22), rep(1, 22))

# Performing the IWT
IWT_result <- functional_lm_test(
  temperature ~ groups,
  B = 2L,
  correction = "IWT"
)
#> Error in eval(predvars, data, env): object 'groups' not found

# Summary of the IWT results
summary(IWT_result)
#> Error: object 'IWT_result' not found
```
