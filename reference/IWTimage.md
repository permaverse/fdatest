# Heatmap plot of the Interval Wise Testing Procedure results

Plotting function creating a ggplot2 graphical output of the IWT: the
p-value heat-map, the adjusted p-value function, and the functional
data, assembled via patchwork.

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
  [`functional_one_sample_test()`](https://permaverse.github.io/fdatest/reference/functional_one_sample_test.md),
  [`iwt1()`](https://permaverse.github.io/fdatest/reference/IWT1.md),
  [`functional_two_sample_test()`](https://permaverse.github.io/fdatest/reference/functional_two_sample_test.md),
  [`iwt2()`](https://permaverse.github.io/fdatest/reference/IWT2.md), or
  the legacy functions
  [`IWT1()`](https://permaverse.github.io/fdatest/reference/IWT1.md) and
  [`IWTaov()`](https://permaverse.github.io/fdatest/reference/IWTaov.md).
  When using
  [`functional_two_sample_test()`](https://permaverse.github.io/fdatest/reference/functional_two_sample_test.md)
  or [`iwt2()`](https://permaverse.github.io/fdatest/reference/IWT2.md),
  `correction` must be `"IWT"`.

- alpha:

  Threshold for the interval-wise error rate used for the hypothesis
  test. Regions where the adjusted p-value is below `alpha` are
  highlighted. The default is `alpha = 0.05`.

- abscissa_range:

  Range of the plot abscissa. The default is `c(0, 1)`.

- nlevel:

  Number of desired color levels for the p-value heatmap. The default is
  `nlevel = 20`.

- plot_unadjusted:

  Flag indicating if the unadjusted p-value function has to be overlaid
  (dashed line) on the adjusted p-value panel. The default is `FALSE`.

## Value

An object of class `patchwork` containing the assembled ggplot2 panels,
returned invisibly. The plot is also printed as a side effect.

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
[`plot.fos()`](https://permaverse.github.io/fdatest/reference/plot.fos.md),
[`plot.fts()`](https://permaverse.github.io/fdatest/reference/plot.fts.md),
[`plot.flm()`](https://permaverse.github.io/fdatest/reference/plot.flm.md)
and
[`plot.faov()`](https://permaverse.github.io/fdatest/reference/plot.faov.md)
for the plot method applied to the IWT results of one- and
two-population tests, linear models, and ANOVA, respectively.

## Examples

``` r
# Performing the IWT for one population
IWT_result <- functional_one_sample_test(
  NASAtemp$paris, mu = 4, n_perm = 10L
)

# Plotting the results of the IWT
IWTimage(IWT_result, abscissa_range = c(0, 12))


# Selecting the significant components at 5% level
which(IWT_result$adjusted_pvalues < 0.05)
#>   [1]   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
#>  [19]  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36
#>  [37]  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54
#>  [55]  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72
#>  [73]  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90
#>  [91]  91  92  93  94  95  96  97  98  99 100 101 102 103 104 105 106 107 108
#> [109] 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126
#> [127] 127 128 129 130 131 132 133 134 135 258 259 260 261 262 263 264 265 266
#> [145] 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284
#> [163] 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302
#> [181] 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320
#> [199] 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338
#> [217] 339 340 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 356
#> [235] 357 358 359 360 361 362 363 364 365
```
