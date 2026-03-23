# Heatmap plot of the Interval Wise Testing Procedure results

Plotting function creating a graphical output of the IWT: the p-value
heat-map, the plot of the corrected p-values, and the plot of the
functional data.

## Usage

``` r
IWTimage(
  IWT_result,
  alpha = 0.05,
  abscissa_range = c(0, 1),
  nlevel = 20L,
  plot_unadjusted = FALSE
)
```

## Arguments

- IWT_result:

  Results of the IWT, as created by
  [`IWT1`](https://permaverse.github.io/fdatest/reference/IWT1.md),
  [`IWT2`](https://permaverse.github.io/fdatest/reference/IWT2.md),
  [`IWTaov`](https://permaverse.github.io/fdatest/reference/IWTaov.md),
  and
  [`IWTlm`](https://permaverse.github.io/fdatest/reference/IWTlm.md).

- alpha:

  Threshold for the interval-wise error rate used for the hypothesis
  test. The default is `alpha=0.05`.

- abscissa_range:

  Range of the plot abscissa. The default is `c(0,1)`.

- nlevel:

  Number of desired color levels for the p-value heatmap. The default is
  `nlevel=20`.

- plot_unadjusted:

  Flag indicating if the unadjusted p-value function has to be added to
  the plots. The default is `FALSE`.

## Value

No value returned.

## References

Pini, A., & Vantini, S. (2018). Interval-wise testing for functional
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

See
[`plot.IWT1()`](https://permaverse.github.io/fdatest/reference/plot.IWT1.md),
[`plot.ftwosample()`](https://permaverse.github.io/fdatest/reference/plot.ftwosample.md),
[`plot.flm()`](https://permaverse.github.io/fdatest/reference/plot.flm.md)
and
[`plot.fanova()`](https://permaverse.github.io/fdatest/reference/plot.fanova.md)
for the plot method applied to the IWT results of one- and
two-population tests, linear models, and ANOVA, respectively.

## Examples

``` r
# Performing the IWT for two populations
IWT_result <- IWT2(NASAtemp$milan, NASAtemp$paris, B = 10L)

# Plotting the results of the IWT
IWTimage(IWT_result, abscissa_range = c(0, 12))

# Selecting the significant components for the radius at 5\% level
which(IWT_result$corrected_pval < 0.05)
#> integer(0)
```
