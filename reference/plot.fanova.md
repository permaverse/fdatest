# Plot for Functional ANOVA Test Results

The `S3` methods `autoplot.fanova()` and `plot.fanova()` are methods for
plotting results of functional analysis of variance tests. They
visualize the functional data and the adjusted p-values obtained from
the testing procedures for mean comparison of multiple groups. The plots
highlight significant effects at two levels of significance, `alpha1`
and `alpha2`, using shaded areas.

## Usage

``` r
# S3 method for class 'fanova'
autoplot(
  object,
  xrange = c(0, 1),
  alpha1 = 0.05,
  alpha2 = 0.01,
  plot_adjpval = FALSE,
  ylim = NULL,
  col = 1,
  ylabel = "Functional Data",
  title = NULL,
  linewidth = 0.5,
  type = "l",
  ...
)

# S3 method for class 'fanova'
plot(
  x,
  xrange = c(0, 1),
  alpha1 = 0.05,
  alpha2 = 0.01,
  plot_adjpval = FALSE,
  ylim = NULL,
  col = 1,
  ylab = "Functional Data",
  main = NULL,
  lwd = 0.5,
  type = "l",
  ...
)
```

## Arguments

- object, x:

  An object of class `fanova`, usually a result of a call to
  [`functional_anova_test()`](https://permaverse.github.io/fdatest/reference/functional_anova_test.md),
  [`iwt_aov()`](https://permaverse.github.io/fdatest/reference/IWTaov.md),
  [`twt_aov()`](https://permaverse.github.io/fdatest/reference/TWTaov.md)
  or
  [`global_aov()`](https://permaverse.github.io/fdatest/reference/global_aov.md).

- xrange:

  A length-2 numeric vector specifying the range of the x-axis for the
  plots. Defaults to `c(0, 1)`. This should match the domain of the
  functional data.

- alpha1:

  A numeric value specifying the first level of significance used to
  select and display significant effects. Defaults to `alpha1 = 0.05`.

- alpha2:

  A numeric value specifying the second level of significance used to
  select and display significant effects. Defaults to `alpha2 = 0.01`.

- plot_adjpval:

  A boolean value specifying whether the plots of adjusted p-values
  should be displayed. Defaults to `FALSE`.

- ylim:

  A 2-length numeric vector specifying the range of the y-axis. Defaults
  to `NULL`, which determines automatically the range from functional
  data.

- col:

  An integer specifying the color for the plot of functional data.
  Defaults to `1L`.

- ylabel, ylab:

  A string specifying the label of the y-axis of the functional data
  plot. Defaults to `"Functional Data"`.

- title, main:

  A string specifying the title of the functional data plot. Defaults to
  `NULL` in which case no title is displayed.

- linewidth, lwd:

  A numeric value specifying the width of the line for the functional
  data plot. Note that the line width for the adjusted p-value plot will
  be twice this value. Defaults to `0.5`.

- type:

  A string specifying the type of plot for the functional data. Defaults
  to `"l"` for lines.

- ...:

  Other arguments passed to specific methods. Not used in this function.

## Value

The `autoplot.fanova()` function creates a ggplot object that displays
the functional data and the adjusted p-values. The significant intervals
at levels `alpha1` and `alpha2` are highlighted in the plots. The
`plot.fanova()` function is a wrapper around `autoplot.fanova()` that
prints the plot directly.

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
for the plot of p-values heatmaps (for IWT).

## Examples

``` r
temperature <- rbind(NASAtemp$milan, NASAtemp$paris)
groups <- c(rep(0, 22), rep(1, 22))

# Performing the TWT
TWT_result <- functional_anova_test(
  temperature ~ groups,
  correction = "TWT",
  B = 5L
)
#> 
#> ── Point-wise tests ────────────────────────────────────────────────────────────
#> 
#> ── Threshold-wise tests ────────────────────────────────────────────────────────
#> 
#> ── Threshold-Wise Testing completed ────────────────────────────────────────────

# Plotting the results of the TWT
plot(
  TWT_result,
  xrange = c(0, 12),
  main = 'TWT results for testing mean differences'
)
```
